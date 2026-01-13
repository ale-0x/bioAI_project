#  peptide_operators.py
import random

from inspyred.ec import Individual
from typing import List, Tuple, Dict, Any, TYPE_CHECKING

import constants as C

from utils import print_verbose, print_warning

# Solo per l'hinting, evitando dipendenze cicliche o problemi di runtime
if TYPE_CHECKING:
    from inspyred.ec.archivers import Archiver
    from inspyred.ec.ec        import EvolutionaryComputation



# --- Operatore di Crossover ---

def single_point_crossover(
    random: random.Random, 
    parent1: Any, 
    parent2: Any
) -> List[str]:
    """
    Esegue il Crossover a Punto Singolo tra due sequenze peptidiche.

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale.
    parent1 : `Individual` o `str`
        Il primo genitore.
    parent2 : `Individual` o `str`
        Il secondo genitore.

    Returns
    -------
    `List[str]`
        Una lista contenente le due nuove sequenze peptidiche (stringhe).
    """
    
    # Estraggo la stringa candidata in modo sicuro (come dalla correzione precedente)
    parent1_str: str = parent1.candidate if isinstance(parent1, Individual) else str(parent1)
    parent2_str: str = parent2.candidate if isinstance(parent2, Individual) else str(parent2)

    length: int = len(parent1_str)

    if length <= 1:
        return [parent1_str, parent2_str]

    crossover_point: int = random.randint(1, length - 1)

    child1: str = parent1_str[:crossover_point] + parent2_str[crossover_point:]
    child2: str = parent2_str[:crossover_point] + parent1_str[crossover_point:]

    return [child1, child2]



# --- Logica per Mutazione BLOSUM ---

def get_blosum_weights(
        original_aa: str,
        args       : Dict[str, Any]
) -> Tuple[List[str], List[float]]:
    """
    Genera i pesi di probabilità per la sostituzione di un aminoacido basandosi 
    sulla matrice BLOSUM62.

    I punteggi BLOSUM62 sono trasformati in pesi positivi per essere utilizzati 
    nella selezione pesata casuale, favorendo mutazioni biologicamente plausibili.

    Parameters
    ----------
    original_aa : `str`
        L'aminoacido di partenza (codice a una lettera) che deve essere mutato.

    Returns
    -------
    `Tuple[List[str], List[float]]`
        Una tupla contenente:
        - Una lista di stringhe con tutti gli aminoacidi candidati alla sostituzione.
        - Una lista di float che rappresenta il peso (probabilità relativa) di ciascun 
          aminoacido target.
    
    Notes
    -----
    Viene applicata una trasformazione `weight = max(1, (score + 5)^2)` per convertire i 
    punteggi negativi BLOSUM in probabilità positive e amplificare le mutazioni 
    biologicamente favorite.
    """
    # Ritorna alla selezione uniforme se la matrice non è disponibile
    if C.BLOSUM62 is None:
        return list(C.AMINO_ACIDS), [1.0] * len(C.AMINO_ACIDS)
        
    targets: List[str]   = []
    weights: List[float] = []
    
    for target_aa in C.AMINO_ACIDS:
        try:
            # Recupera il punteggio
            # Utilizza .get per gestire le coppie (target, original) se non sono simmetriche o presenti
            score: int = C.BLOSUM62.get((original_aa, target_aa), 0)
        except (KeyError, ValueError):
            continue

        # Trasformazione del punteggio in un peso positivo (stessa logica di prima)
        OFFSET  : int   = 5
        EXPONENT: int   = 2
        weight  : float = float(max(1, (score + OFFSET) ** EXPONENT)) 
        
        targets.append(target_aa)
        weights.append(weight)

    return targets, weights



# --- Funzione Helper per il Rate Adattivo ---

def get_adaptive_mutation_rate(args: Dict[str, Any]) -> float:
    """
    Calcola il tasso di mutazione adattivo basandosi sul progresso del GA.

    Il rate decade linearmente da `INITIAL_MUTATION_RATE` a `FINAL_MUTATION_RATE` 
    man mano che le generazioni avanzano. Questo implementa una strategia 
    "Esplorazione" (alto rate) iniziale e "Sfruttamento" (basso rate) finale.

    Parameters
    ----------
    args : `Dict[str, Any]`
        Dizionario di argomenti, deve contenere la chiave '_ec' (EvolutionaryComputation) 
        e 'max_generations' per il calcolo.

    Returns
    -------
    `float`
        Il tasso di mutazione corrente.
    
    Examples
    --------
    # Se il GA è a 50 su 100 generazioni (a metà)
    # Rate = 0.30 - (0.30 - 0.05) * 0.5 = 0.175
    """
    # print_verbose("get_adaptive_mutation_rate:", args)
    # inspyred passa l'oggetto EvolutionaryComputation in args['_ec']
    ec = args.get('_ec', None)
    
    if ec is None:
        return C.FINAL_MUTATION_RATE # Fallback sicuro
        
    current_gen = ec.num_generations
    max_gen = args.get('max_generations', C.MAX_GENERATIONS) # Deve corrispondere a MAX_GENERATIONS
    
    if max_gen == 0: return C.FINAL_MUTATION_RATE

    # Progresso da 0.0 a 1.0
    progress = min(1.0, current_gen / max_gen)
    
    # Interpolazione lineare (Lerp)
    # Rate = Start - (Start - End) * Progress
    current_rate = C.INITIAL_MUTATION_RATE - (C.INITIAL_MUTATION_RATE - C.FINAL_MUTATION_RATE) * progress
    
    return current_rate


# --- Operatore di Mutazione ---

def blosum_peptide_mutator(
        random   : random.Random, 
        candidate: str,
        args     : Dict[str, Any]
) -> str:
    """
    Mutatore che applica mutazioni agli aminoacidi di una sequenza peptidica.

    La probabilità di mutare un aminoacido in un altro è ponderata dalla 
    matrice di sostituzione BLOSUM62, favorendo le sostituzioni conservative.
    Il tasso di mutazione viene adattato dinamicamente con l'avanzare delle generazioni.

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale.
    candidate : `str`
        La sequenza peptidica da mutare.
    args : `Dict[str, Any]`
        Dizionario di argomenti, usato per recuperare il rate adattivo.

    Returns
    -------
    `str`
        La sequenza peptidica mutata.

    Examples
    --------
    >>> seq = 'LVTA'
    >>> mutated_seq = blosum_peptide_mutator(random_instance, seq, {'mutation_rate': 0.1})
    >>> print(mutated_seq) 
    'IVTA' # Esempio di mutazione da L a I, favorita da BLOSUM.
    """
    # print_verbose("blosum_peptide_mutator:", args)
    mutation_probability: float     = get_adaptive_mutation_rate(args)
    mutated_sequence    : List[str] = list(candidate)
    
    for i in range(len(mutated_sequence)):
        if random.random() <= mutation_probability:
            original_aa        : str = mutated_sequence[i]
            targets, weights         = get_blosum_weights(original_aa, args) # Passiamo args anche qui
            new_amino_acid     : str = random.choices(targets, weights = weights, k = 1)[0]
            mutated_sequence[i]      = new_amino_acid

    return "".join(mutated_sequence)


# --- Operatore a Catena ---

def peptide_chain_variator(
    random: random.Random, 
    candidates: List[Any], # Lista di Individual o str, a seconda di inspyred
    args: Dict[str, Any]
) -> List[str]:
    """
    Combina Crossover e Mutazione in un'unica catena di variatori.

    Questo operatore è usato come singolo `ea.variator` per garantire 
    che l'oggetto passato al core di `inspyred` sia una funzione e non una tupla,
    evitando l'AttributeError riscontrato.
    
    Il 100% dei candidati selezionati per la variazione (genitori) viene 
    sottoposto a crossover, e il 100% dei figli risultanti viene sottoposto a mutazione.

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale.
    candidates : `List[Individual]`
        Lista dei candidati (genitori) selezionati per la variazione.
    args : `Dict[str, Any]`
        Dizionario di argomenti, deve contenere il rate di mutazione.

    Returns
    -------
    `List[str]`
        Lista di sequenze peptidiche (figli) variate.
    """
    # print_verbose("peptide_chain_variator:", args)

    # --- DEBUG START ---
    print_verbose(f"DEBUG: Input Candidates: {len(candidates)}")
    # --- DEBUG END ---
    # Inizializza la lista per i nuovi figli
    new_population: List[str] = []
    
    # Itera sulla lista dei candidati (genitori) a passi di 2
    for i in range(0, len(candidates), 2):
        
        # Gestisce il caso in cui ci sia un solo candidato rimanente (lista dispari)
        if i + 1 >= len(candidates):
            # Se c'è un solo genitore, lo mutiamo e lo aggiungiamo alla nuova popolazione
            lonely_parent = candidates[i]
            # Assumiamo che il genitore solitario sia una stringa (o lo estraiamo)
            lonely_seq = lonely_parent.candidate if isinstance(lonely_parent, Individual) else str(lonely_parent)
            
            mutated_child = blosum_peptide_mutator(random, lonely_seq, args)
            new_population.append(mutated_child)
            break 
            
        parent1 = candidates[i]
        parent2 = candidates[i+1]
        
        # 1. Esegui il Crossover in base a una probabilità
        if random.random() <= C.CROSSOVER_PROBABILITY:
            children_sequences: List[str] = single_point_crossover(random, parent1, parent2)
        else:
            # Se non avviene il crossover, i figli sono copie dei genitori
            parent1_str = parent1.candidate if isinstance(parent1, Individual) else str(parent1)
            parent2_str = parent2.candidate if isinstance(parent2, Individual) else str(parent2)
            children_sequences = [parent1_str, parent2_str]
        
        # 2. Esegui la Mutazione su ciascun figlio
        for child_seq in children_sequences:
            mutated_child: str = blosum_peptide_mutator(random, child_seq, args)
            new_population.append(mutated_child)
        
    # --- DEBUG START ---
    print_verbose(f"DEBUG: Output Children: {len(new_population)}")
    if len(new_population) < len(candidates):
        print_warning("!!! ALLARME: Sto perdendo individui per strada! !!!")
    # --- DEBUG END ---
            
    return new_population

# --- Hydrophobucity ---
def get_hydrophobicity(peptide: str) -> float:
    """Calcola l'idrofobicità media di una sequenza peptidica."""
    average_hydrophobicity = sum(C.AMINO_HYDROPHOBICITY_PH7.get(aa, 0) for aa in peptide) / len(peptide)
    return average_hydrophobicity
