#!/bin/bash
#SBATCH --job-name=BioInspyredAI_PeptideGA
#SBATCH --output=../log/ga_%j.log
#SBATCH --error=../log/ga_%j.log
#SBATCH --partition=edu-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 0-24:00:00

echo "Job started on $(hostname) at $(date)"

# 1. Carica l'ambiente (adatta il path al tuo miniconda)
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate bioai

RECEPTOR=${RECEPTOR:-"../resources/pdbqt/2P3D.pdbqt"}
TMP_BASE="../resources/tmp"

echo "Lancio Job $SLURM_JOB_ID con Recettore: $RECEPTOR"
echo "Temp Base: $TMP_BASE"

python -u main.py                                               \
    2P3D                                                        \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        10                                  \
    --population_size       5                                   \
    --generations           3                                   \
    --initial_mutation_rate 0.30                                \
    --final_mutation_rate   0.05                                \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result_$SLURM_JOB_ID.txt \
    --center_x=8.084                                            \
    --center_y=-13.829                                          \
    --center_z=-0.140                                           \
    --size_x                32                                  \
    --size_y                32                                  \
    --size_z                32                                  \
    --exhaustiveness        8                                   \
    --vina_exe_path         vina                                \
    --no_delete

echo "Job finished at $(date)"
