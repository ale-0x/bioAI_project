# ga_problem.py
import random
from typing import List, Dict, Any, Union
import inspyred.ec
from constants import AMINO_ACIDS, PEPTIDE_LENGTH

# --- Funzione di Generazione della Popolazione Iniziale ---

def peptide_generator(random: random.Random, args: Dict[str, Any]) -> str:
    """
    Generatore della sequenza peptidica iniziale per la popolazione.

    Crea una sequenza casuale di aminoacidi di lunghezza fissa 
    (definita da `PEPTIDE_LENGTH` o dagli argomenti).

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale fornita dal framework `inspyred`.
    args : `Dict[str, Any]`
        Dizionario di argomenti opzionali. Può contenere la chiave 
        `'peptide_length'` (`int`) per specificare la lunghezza.

    Returns
    -------
    `str`
        Una sequenza peptidica casuale (stringa di codici a una lettera).

    Examples
    --------
    >>> generator = peptide_generator(random_instance, {'peptide_length': 5})
    >>> print(generator)
    'ALTSV' # Esempio di output casuale
    """
    length: int = args.get('peptide_length', PEPTIDE_LENGTH)

    return "".join(random.choice(AMINO_ACIDS) for _ in range(length))



# --- Funzione di Valutazione (Evaluation/Fitness) ---

def evaluate_peptide_binding(candidates: List[str], args: Dict[str, Any]) -> List[float]:
    """
    Valutazione della Fitness: Stima l'energia di legame (Binding Affinity).

    Questa funzione calcola la fitness (energia di legame) per ogni sequenza candidata. 
    L'obiettivo è massimizzare questo valore (che è negativo per l'energia di legame di Vina).

    Parameters
    ----------
    candidates : `List[str]`
        Una lista di sequenze peptidiche (stringhe) da valutare.
    args : `Dict[str, Any]`
        Dizionario di argomenti opzionali (es. posizione della tasca di docking).

    Returns
    -------
    `List[float]`
        Una lista di valori di fitness (energia di legame Vina stimata, in kcal/mol)
        per ogni sequenza candidata.

    Notes
    -----
    ATTUALMENTE È UN PLACEHOLDER. L'implementazione finale dovrà integrare:
    1. Sequence-to-Structure Prediction (es. PepFold) per generare il PDB.
    2. Preparazione dei file PDBQT.
    3. Esecuzione del Docking con AutoDock Vina o un'alternativa.
    4. Estrazione del punteggio di energia.
    """
    
    # --- LOGICA ATTUALE (PLACEHOLDER CASUALE) ---
    rand_placeholder = random.Random(42) 
    
    # Valori casuali tra -10.0 (buona affinità) e -2.0 (scarsa affinità)
    fitnesses: List[float] = [rand_placeholder.uniform(-10.0, -2.0) for _ in candidates]
    
    return fitnesses
    # ---------------------------------------------