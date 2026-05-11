nextflow.enable.dsl = 2

/*
  Stages (ordered):
  1) data_download.sh        — parallel per (sample, SRA)
  2) quality_control.sh      — parallel per sample
  3) read_alignment.sh       — STAR genome index (once) + STARsolo per sample (parallel)
  4) per_sample_qc_and_cell_filtering.py — parallel per sample
  5) integration_and_batch_correction.py — single merge + Harmony
  6) differential_analysis.R — uses integrated h5ad
  7) visualize_sc_rna_data.R — figures from integrated h5ad (via zellkonverter if no RDS)
*/

workflow {
    ch_samples = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample.toString(), row.sra.toString()) }

    DOWNLOAD_SAMPLES(ch_samples)

    TRIM_GALORE(DOWNLOAD_SAMPLES.out)

    STAR_INDEX(
        channel.of(
            tuple(
                file(params.ref_fasta, checkIfExists: false),
                file(params.ref_gtf, checkIfExists: false)
            )
        )
    )

    ch_star_index = STAR_INDEX.out.map { it }.first()

    ch_align_in = TRIM_GALORE.out
        .combine(ch_star_index)
        .map { sample, trim_dir, index_dir ->
            tuple(
                sample,
                trim_dir,
                index_dir,
                file(params.star_whitelist, checkIfExists: false)
            )
        }

    STAR_SOLO(ch_align_in)

    PER_SAMPLE_QC(STAR_SOLO.out)

    ch_t1 = PER_SAMPLE_QC.out.filter { sample, _h5 -> sample == params.tumor1_sample }.map { _s, h -> h }.first()
    ch_t2 = PER_SAMPLE_QC.out.filter { sample, _h5 -> sample == params.tumor2_sample }.map { _s, h -> h }.first()
    ch_n  = PER_SAMPLE_QC.out.filter { sample, _h5 -> sample == params.normal_sample }.map { _s, h -> h }.first()

    INTEGRATION(ch_t1, ch_t2, ch_n)

    DIFFERENTIAL(INTEGRATION.out)

    // Run visualization only after differential analysis completes (same integrated h5ad + DE outputs)
    VISUALIZE(INTEGRATION.out.combine(DIFFERENTIAL.out))
}

process DOWNLOAD_SAMPLES {
    tag "$sample"
    label "io"
    cpus 2
    memory '8 GB'
    time '24.h'

    input:
    tuple val(sample), val(sra)

    output:
    tuple val(sample), path("${sample}_R1.fastq.gz"), path("${sample}_R2.fastq.gz")

    script:
    """
    set -euo pipefail
    prefetch "${sra}"
    fasterq-dump --split-files "${sra}"
    if command -v pigz &>/dev/null; then
      pigz -p ${task.cpus} *.fastq 2>/dev/null || true
    else
      for f in *.fastq; do
        [ -f "\$f" ] || continue
        gzip -f "\$f"
      done
    fi
    shopt -s nullglob
    r1=( *_1.fastq.gz )
    r2=( *_2.fastq.gz )
    if [ \${#r1[@]} -lt 1 ] || [ \${#r2[@]} -lt 1 ]; then
      echo "Could not find *_1.fastq.gz / *_2.fastq.gz after download." >&2
      ls -la >&2 || true
      exit 1
    fi
    mv "\${r1[0]}" "${sample}_R1.fastq.gz"
    mv "\${r2[0]}" "${sample}_R2.fastq.gz"
    """
}

process TRIM_GALORE {
    tag "$sample"
    label "cpu_medium"
    cpus { params.trim_cores as int }
    memory '16 GB'
    time '24.h'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("trimmed")

    script:
    """
    set -euo pipefail
    mkdir -p trimmed
    bash "${params.scripts_dir}/quality_control.sh" "${r1}" "${r2}" trimmed ${task.cpus}
    """
}

process STAR_INDEX {
    label "star_index"
    cpus { params.star_index_threads as int }
    memory '64 GB'
    time '48.h'

    input:
    tuple path(ref_fasta), path(ref_gtf)

    output:
    path("${params.star_genome_dir}")

    script:
    """
    set -euo pipefail
    bash "${params.scripts_dir}/read_alignment.sh" index \\
      "${params.star_genome_dir}" \\
      "${ref_fasta}" \\
      "${ref_gtf}" \\
      ${task.cpus}
    """
}

process STAR_SOLO {
    tag "$sample"
    label "star_align"
    cpus { params.star_align_threads as int }
    memory '64 GB'
    time '48.h'

    input:
    tuple val(sample), path(trim_dir), path(genome_index), path(whitelist)

    output:
    tuple val(sample), path("${sample}")

    script:
    """
    set -euo pipefail
    R1=\$(ls "${trim_dir}"/*_val_1.fq.gz 2>/dev/null | head -1)
    R2=\$(ls "${trim_dir}"/*_val_2.fq.gz 2>/dev/null | head -1)
    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
      echo "Trimmed reads not found (expected *_val_1.fq.gz / *_val_2.fq.gz in ${trim_dir})." >&2
      ls -la "${trim_dir}" >&2 || true
      exit 1
    fi
    mkdir -p "${sample}"
    bash "${params.scripts_dir}/read_alignment.sh" align \\
      "${genome_index}" \\
      "${whitelist}" \\
      "\$R2" \\
      "\$R1" \\
      "${sample}/" \\
      ${task.cpus}
    """
}

process PER_SAMPLE_QC {
    tag "$sample"
    label "scanpy"
    cpus 4
    memory '32 GB'
    time '12.h'

    input:
    tuple val(sample), path(star_out)

    output:
    tuple val(sample), path("${sample}_filtered.h5ad")

    script:
    def mtx_dir = "${star_out}/Solo.out/GeneFull/raw"
    """
    set -euo pipefail
    python3 "${params.scripts_dir}/per_sample_qc_and_cell_filtering.py" \\
      --mtx-dir "${mtx_dir}" \\
      --sample "${sample}" \\
      --out "${sample}_filtered.h5ad"
    """
}

process INTEGRATION {
    label "integration"
    cpus 8
    memory '64 GB'
    time '24.h'

    input:
    path tumor1_h5
    path tumor2_h5
    path normal_h5

    output:
    path("integrated_adata.h5ad")

    script:
    """
    set -euo pipefail
    python3 "${params.scripts_dir}/integration_and_batch_correction.py" \\
      --tumor1 "${tumor1_h5}" \\
      --tumor2 "${tumor2_h5}" \\
      --normal "${normal_h5}" \\
      --out integrated_adata.h5ad
    """
}

process DIFFERENTIAL {
    label "deseq2"
    cpus 4
    memory '32 GB'
    time '12.h'

    input:
    path(integrated_h5ad)

    output:
    path("deseq2_results")

    script:
    """
    set -euo pipefail
    mkdir -p deseq2_results
    Rscript "${params.scripts_dir}/differential_analysis.R" \\
      "${integrated_h5ad}" \\
      deseq2_results
    """
}

process VISUALIZE {
    label "r_vis"
    cpus 4
    memory '32 GB'
    time '12.h'

    input:
    tuple path(integrated_h5ad), path(deseq2_results)

    output:
    path("figures")

    script:
    """
    set -euo pipefail
    export ANNDATA_H5AD="${integrated_h5ad}"
    export CLUSTER_COL="leiden"
    export SAMPLE_COL="sample"
    export CONDITION_COL="condition"
    export FORCE_H5AD=1
    mkdir -p figures
    Rscript "${params.scripts_dir}/visualize_sc_rna_data.R" \\
      "${integrated_h5ad}" \\
      figures
    """
}
