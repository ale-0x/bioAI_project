# constants.py

__author__  = "Alex Callegaro and Clemente Calabrese"
__version__ = '1.2.1'

from Bio.Align import substitution_matrices
from pathlib   import Path


SEED = 77

# The 20 canonical amino acids (one-letter codes)
AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

# Load the BLOSUM62 matrix
try:
    BLOSUM62 = substitution_matrices.load("BLOSUM62")
except ImportError:
    from utils import print_warning
    print_warning("Biopython not installed or BLOSUM62 not found. Mutation will use uniform probabilities.")
    BLOSUM62 = None

# --- Preparing the original .pdb ---
PROTEIN_NAME      = "2P3D"
INPUT_PDB_FILE    = Path(f"../resources/pdb/{PROTEIN_NAME}.pdb")
OUTPUT_PDBQT_FILE = Path(f"../resources/pdbqt/{PROTEIN_NAME}.pdbqt")
SAFETY_MARGIN     = 30

# --- Genetic Algorithm Parameters ---
PEPTIDE_LENGTH  = 10        # Fixed length of the peptide to evolve
POPULATION_SIZE = 50        # Population size per generation
MAX_GENERATIONS = 60        # Maximum number of generations

# --- Adaptive Parameters ---
CROSSOVER_PROBABILITY = 1.00                        # Probability of crossover occurring: 100% fa sempre il crossover
INITIAL_MUTATION_RATE = 0.70                        # 70% at the beginning (High Chaos/Exploration)
FINAL_MUTATION_RATE   = 0.20                        # 20% at the end (Refinement/Exploitation)

# --- Vina Docking Setup ---
# CAUTION: These values ​​must be protein specific (e.g. 2P3D)
RECEPTOR_FILE                = OUTPUT_PDBQT_FILE        # The prepared receptor
CENTER_X, CENTER_Y, CENTER_Z = 8.084, -13.829, -0.140   # Pocket center coordinates
SIZE_X, SIZE_Y, SIZE_Z       = 32, 32, 32               # Box dimensions (Angstroms)
CPUS                         = 8
EXHAUSTIVENESS               = 8                        # Search precision (8 is the default, 32 is slow but precise)
VINA_EXE_PATH                = Path("vina")             # Assume 'vina' is in the PATH

# --- Hydrophobicity parameters ---
HYDROPHOBICITY_WEIGHT = 0.02        # Weight of the sum of hydrophobicity in the fitness function
AMINO_HYDROPHOBICITY_PH7 = {        # Hydrophobicity index at pH 7 (from Monera, O. D., et al. Relationship of Sidechain Hydrophobicity and Alpha-Helical Propensity on the Stability of the Single-Stranded Amphipathic Alpha-Helix. J Pept Sci. (1995).)
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

# --- Parameters for producing IEEE-style plots ---

IEEE_FIGURE_WIDTH  = 3.5            # Figure width in inches
IEEE_FIGURE_HEIGHT = 2.5            # Figure height in inches

# Complete IEEE style configuration
IEEE_PLOT_PARAMS = {
    "font.family"    : "serif",
    "font.serif"     : ["Times New Roman", "Liberation Serif", "DejaVu Serif", "serif"],
    "font.size"      : 10,
    "axes.labelsize" : 10,
    "axes.titlesize" : 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.5,
    "figure.dpi"     : 300,
    "savefig.bbox"   : "tight",
}

