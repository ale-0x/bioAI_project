#!/bin/bash
#SBATCH --job-name=clem_bioai
#SBATCH --output=bioai_%j.log
#SBATCH --error=bioai_%j.log
#SBATCH --partition=edu-long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -t 0-24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user='calabrese.clemente@studenti.unitn.it'

echo "Job started on $(hostname) at $(date)"

# 1. Carica l'ambiente (adatta il path al tuo miniconda)
source $HOME/anaconda3/etc/profile.d/conda.sh
conda activate bioai

cd $HOME/bioAI_project/code

RECEPTOR=${RECEPTOR:-"../resources/pdbqt/2P3D.pdbqt"}
TMP_BASE="../resources/tmp"

echo "Lancio Job $SLURM_JOB_ID con Recettore: $RECEPTOR"
echo "Temp Base: $TMP_BASE"

python3 -u main.py                                              \
    2P3D_no_ligand                                                        \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        6                                   \
    --population_size       5                                   \
    --generations           8                                   \
    --initial_mutation_rate 0.30                                \
    --final_mutation_rate   0.05                                \
    --hydrophobicity_weight 0.10                                \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result                   \
    --center_x=8.084                                            \
    --center_y=-13.829                                          \
    --center_z=-0.140                                           \
    --size_x                32                                  \
    --size_y                32                                  \
    --size_z                32                                  \
    --exhaustiveness        4                                   \
    --vina_exe_path         /home/clemente.calabrese/.conda/envs/bioai/bin/vina                                \
    --no_delete

echo "Job finished at $(date)"
