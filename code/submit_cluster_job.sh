#!/bin/bash
#SBATCH --job-name=BioInspyredAI_PeptideGA_v7
#SBATCH --output=../logs/ga_%j.log
#SBATCH --error=../logs/ga_%j.log
#SBATCH --partition=edu-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 0-24:00:00

echo "Job started on $(hostname) at $(date)"

# Load the environment
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate bioai

cd $HOME/bioAI_project/code

RECEPTOR=${RECEPTOR:-"../resources/pdbqt/2P3D_no_ligand.pdbqt"}
TMP_BASE="../resources/tmp"

echo "Launching Job $SLURM_JOB_ID with Receptor: $RECEPTOR"
echo "Base Temp: $TMP_BASE"

python3 -u main.py                                              \
    2P3D_no_ligand                                              \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        6                                   \
    --population_size       50                                  \
    --generations           50                                  \
    --initial_mutation_rate 0.60                                \
    --final_mutation_rate   0.15                                \
    --hydrophobicity_weight 0.015                               \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result                   \
    --deadline              23:56:00                            \
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
