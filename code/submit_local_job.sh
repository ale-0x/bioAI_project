# 1. Trova l'ultimo numero utilizzato nelle cartelle "bioai_ga_XXXXXXX"
# ls -1: elenca i file/cartelle
# grep: filtra solo quelle con il prefisso giusto
# sed: estrae solo la parte numerica finale
# sort -n | tail -1: prende il numero più grande
TMP_BASE="../resources/tmp"

# 1. Estrae l'ID e rimuove gli zeri iniziali usando sed
LAST_ID=$(ls -1 "$TMP_BASE" 2>/dev/null | grep "bioai_ga_" | sed 's/bioai_ga_//' | sort -n | tail -1)

if [ -z "$LAST_ID" ]; then
    NEW_ID=1
else
    # Pulizia: rimuove gli zeri iniziali per evitare l'errore ottale
    # 's/^0*//' dice a sed di cancellare tutti gli zeri all'inizio della stringa
    CLEAN_ID=$(echo $LAST_ID | sed 's/^0*//')
    
    # Se dopo la pulizia la stringa è vuota (era "0000000"), impostiamo a 0
    [ -z "$CLEAN_ID" ] && CLEAN_ID=0
    
    NEW_ID=$((CLEAN_ID + 1))
fi

# 3. Formatta di nuovo a 7 cifre
SLURM_JOB_ID=$(printf "%07d" $NEW_ID)

RECEPTOR=${RECEPTOR:-"../resources/pdbqt/2P3D.pdbqt"}
SLURM_CPUS_PER_TASK="6"

echo "Lancio Job $SLURM_JOB_ID con Recettore: $RECEPTOR"
echo "Temp Base: $TMP_BASE"

python3 -u main.py                                              \
    2P3D                                                        \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        5                                   \
    --population_size       4                                   \
    --generations           2                                   \
    --initial_mutation_rate 0.30                                \
    --final_mutation_rate   0.05                                \
    --hydrophobicity_weight 0.10                                \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result                   \
    --deadline              24:00:00                            \
    --center_x=8.084                                            \
    --center_y=-13.829                                          \
    --center_z=-0.140                                           \
    --size_x                32                                  \
    --size_y                32                                  \
    --size_z                32                                  \
    --exhaustiveness        1                                   \
    --vina_exe_path         vina                                \
    --no_delete

echo "Job finished at $(date)"
