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

# --- Parametri del GA ---
PEPTIDE_LENGTH  = 10        # Lunghezza fissa del peptide da evolvere
POPULATION_SIZE = 10        # Dimensione della popolazione per generazione (50?)
MAX_GENERATIONS = 100       # Numero massimo di generazioni
MUTATION_RATE   = 0.05      # Probabilità di mutazione per singolo aminoacido