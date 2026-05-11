#!/usr/bin/env bash
# Launch the single-cell Nextflow pipeline on SLURM (Nextflow submits tasks to SLURM).
# Usage: sbatch run_slurm.sh   OR   bash run_slurm.sh
#
# Edit #SBATCH lines and NEXTFLOW_* variables for your site.

#SBATCH --job-name=sc_scrna_nf
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

# Repo root (directory containing Single_cell_tumour_data/ and single_cell_scrna/)
PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${PIPELINE_ROOT}"

mkdir -p logs work

module purge 2>/dev/null || true
# module load nextflow  # uncomment on your cluster
# module load singularity  # if using -with-singularity

export NXF_OPTS="${NXF_OPTS:--Xms500M -Xmx4G}"

# Optional: central Nextflow work directory on scratch
# export NXF_WORK="/scratch/${USER}/nxf-work-sc"

nextflow run "${PIPELINE_ROOT}/single_cell_scrna/main.nf" \
  -profile slurm \
  -work-dir "${PIPELINE_ROOT}/work" \
  --ref_fasta "${PIPELINE_ROOT}/refs/Sscrofa11.1.dna.fa" \
  --ref_gtf "${PIPELINE_ROOT}/refs/Sus_scrofa.Sscrofa11.1.gtf" \
  --star_whitelist "${PIPELINE_ROOT}/refs/3M-february-2018.txt" \
  --samplesheet "${PIPELINE_ROOT}/single_cell_scrna/assets/samplesheet.csv" \
  --slurm_partition "${SLURM_PARTITION:-compute}" \
  ${SLURM_ACCOUNT:+--slurm_account "${SLURM_ACCOUNT}"} \
  "$@"
