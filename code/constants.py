# constants.py
from Bio.Align import substitution_matrices

# I 20 aminoacidi canonici (codici a una lettera)
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

# Carica la matrice BLOSUM62
try:
    BLOSUM62 = substitution_matrices.load("BLOSUM62")
except ImportError:
    print("ATTENZIONE: Biopython non installato o BLOSUM62 non trovata. La mutazione userà probabilità uniformi.")
    BLOSUM62 = None

# --- Preparazione dell'originale .pdb ---
PROTEIN_NAME      = "2P3D"
INPUT_PDB_FILE    = f"../resources/pdb/{PROTEIN_NAME}.pdb"
OUTPUT_PDBQT_FILE = f"../resources/pdbqt/{PROTEIN_NAME}.pdbqt"

# --- Parametri del GA ---
PEPTIDE_LENGTH  = 10                                # Lunghezza fissa del peptide da evolvere
POPULATION_SIZE = 1                                # Dimensione della popolazione per generazione (50?)
MAX_GENERATIONS = 0                               # Numero massimo di generazioni

# --- Parametri Adattivi (Strategia Esplorazione -> Sfruttamento) ---
INITIAL_MUTATION_RATE = 0.30                        # 30% all'inizio (Alto caos/esplorazione)
FINAL_MUTATION_RATE   = 0.05                        # 5% alla fine (Raffinamento/sfruttamento)

# --- Configurazione Vina Docking ---
# ATTENZIONE: Questi valori devono essere specifici per la proteina (es. 7CAM)
RECEPTOR_FILE                = OUTPUT_PDBQT_FILE        # Il recettore preparato
CENTER_X, CENTER_Y, CENTER_Z = 8.084, -13.829, -0.140   # Coordinate centro tasca
SIZE_X, SIZE_Y, SIZE_Z       = 20, 20, 20               # Dimensioni box (Angstrom)
EXHAUSTIVENESS               = 7                        # Precisione di ricerca (8 è default, 32 è lento ma preciso)
VINA_EXE_PATH                = "vina"                   # Assumi che 'vina' sia nel PATH
TEMP_DOCKING                 = "../resources/tmp/temp_docking"
