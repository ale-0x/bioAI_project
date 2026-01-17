# Find the last used number in the "bioai_ga_XXXXXXX" folders
TMP_BASE="../resources/tmp"

# Extract ID and remove leading zeros using sed
LAST_ID=$(ls -1 "$TMP_BASE" 2>/dev/null | grep "bioai_ga_" | sed 's/bioai_ga_//' | sort -n | tail -1)

# ID Calculation and Cleanup (Existing logic kept to calculate NEW_ID)
if [ -z "$LAST_ID" ]; then
    NEW_ID=1
    CLEAN_ID=0
    echo "No previous runs found. Starting from ID: 0000001"
else
    CLEAN_ID=$(echo $LAST_ID | sed 's/^0*//')
    [ -z "$CLEAN_ID" ] && CLEAN_ID=0
    NEW_ID=$((CLEAN_ID + 1))
fi

CHECKPOINT=""
RESUME_ARG=""

if [ "$1" = "--checkpoint" ]; then
    # If the second argument ($2) exists, use that as the path
    if [ -n "$2" ]; then
        CHECKPOINT="$2"
    else
        # Otherwise use automatic logic (based on CLEAN_ID calculated above)
        CHECKPOINT_JOB_ID=$(printf "%07d" $CLEAN_ID)
        CHECKPOINT="../results/result_${CHECKPOINT_JOB_ID}/checkpoint_${CHECKPOINT_JOB_ID}.pkl"
    fi
    
    echo "Previous run found: $CLEAN_ID"
    echo "I set up a resume from: $CHECKPOINT"
    
    # Create the complete string for the Python command
    RESUME_ARG="--resume $CHECKPOINT"
fi


# Format back to 7 digits
SLURM_JOB_ID=$(printf "%07d" $NEW_ID)

RECEPTOR=${RECEPTOR:-"../resources/pdbqt/2P3D.pdbqt"}
SLURM_CPUS_PER_TASK="6"

if [ -n "$CHECKPOINT" ]; then
    echo "---------------------------------------"
    echo "CURRENT JOB ID       : $SLURM_JOB_ID"
    echo "RESUME POSSIBLE FROM : $CHECKPOINT"
    echo "---------------------------------------"
    echo ""
fi

echo "Launching Job $SLURM_JOB_ID with Receptor: $RECEPTOR"
echo "Base Temp: $TMP_BASE"

python3 -u main.py                                              \
    2P3D                                                        \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        5                                   \
    --population_size       4                                   \
    --generations           10                                  \
    --initial_mutation_rate 0.60                                \
    --final_mutation_rate   0.05                                \
    --hydrophobicity_weight 0.02                                \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result                   \
    --deadline              00:04:00                            \
    --center_x=8.084                                            \
    --center_y=-13.829                                          \
    --center_z=-0.140                                           \
    --size_x                32                                  \
    --size_y                32                                  \
    --size_z                32                                  \
    --exhaustiveness        1                                   \
    --vina_exe_path         vina                                \
    --no_delete                                                 \
    $RESUME_ARG

echo "Job finished at $(date)"
