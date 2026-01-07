# constants.py

__author__  = "Alex Callegaro"
__version__ = '1.0.0'

from Bio.Align import substitution_matrices

HYDROPHOBICITY_WEIGHT = 0.1  # Peso della somma di idrofobicità nella funzione di fitness

# Hydrophobicity index at pH 7 
# (from Monera, O. D., et al. Relationship of Sidechain Hydrophobicity and Alpha-Helical Propensity on the Stability of the Single-Stranded Amphipathic Alpha-Helix. J Pept Sci. (1995).)
amino_hydrophobicity_ph7 = {
    'L': 97,   # Leucine
    'I': 99,   # Isoleucine
    'F': 100,  # Phenylalanine
    'W': 97,   # Tryptophan
    'V': 76,   # Valine
    'M': 74,   # Methionine
    'C': 49,   # Cysteine
    'Y': 63,   # Tyrosine
    'A': 41,   # Alanine
    'T': 13,   # Threonine
    'E': -31,  # Glutamate
    'H': 8,    # Histidine
    'G': 0,    # Glycine
    'S': -5,   # Serine
    'Q': -10,  # Glutamine
    'D': -55,  # Aspartate
    'R': -14,  # Arginine
    'K': -23,  # Lysine
    'N': -28,  # Asparagine
    'P': 5     # Proline (using the 6.5 value provided)
}

# Calcola la MEDIA dei valori di idrofobicità per la sequenza
def get_hydrophobicity(peptide: str) -> float:
    """Calcola l'idrofobicità media di una sequenza peptidica."""
    average_hydrophobicity = sum(amino_hydrophobicity_ph7.get(aa, 0) for aa in peptide) / len(peptide)
    return average_hydrophobicity

SEED = 77

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
PEPTIDE_LENGTH  = 10                               # Lunghezza fissa del peptide da evolvere
POPULATION_SIZE = 40                                # Dimensione della popolazione per generazione (50?)
MAX_GENERATIONS = 30                               # Numero massimo di generazioni

# --- Parametri Adattivi (Strategia Esplorazione -> Sfruttamento) ---
CROSSOVER_PROBABILITY = 1.0                         # Probabilità che il crossover avvenga: 100% fa sempre il crossover
INITIAL_MUTATION_RATE = 0.30                        # 30% all'inizio (Alto caos/esplorazione)
FINAL_MUTATION_RATE   = 0.05                        # 5% alla fine (Raffinamento/sfruttamento)

# --- Configurazione Vina Docking ---
# ATTENZIONE: Questi valori devono essere specifici per la proteina (es. 7CAM)
RECEPTOR_FILE                = OUTPUT_PDBQT_FILE        # Il recettore preparato
CENTER_X, CENTER_Y, CENTER_Z = 8.084, -13.829, -0.140   # Coordinate centro tasca
SIZE_X, SIZE_Y, SIZE_Z       = 32, 32, 32               # Dimensioni box (Angstrom)
CPUS                         = 8
EXHAUSTIVENESS               = 8                        # Precisione di ricerca (8 è default, 32 è lento ma preciso)
VINA_EXE_PATH                = "vina"                   # Assumi che 'vina' sia nel PATH
