#  peptide_operators.py
import random

from inspyred.ec import Individual
from typing      import List, Tuple, Dict, Any


import constants as C

from utils import print_verbose, print_warning

# --- Crossover Operator ---

def single_point_crossover(
    random: random.Random, 
    parent1: Any, 
    parent2: Any
) -> List[str]:
    """
    Performs Single-Point Crossover between two peptide sequences.

    Selects a random point in the sequence and swaps the segments 
    between the two parents to create two new children.

    Parameters
    ----------
    random : `random.Random`
        The random number generator instance.
    parent1 : `Individual` or `str`
        The first parent peptide sequence.
    parent2 : `Individual` or `str`
        The second parent peptide sequence.

    Returns
    -------
    `List[str]`
        A list containing the two generated child sequences.

    Examples
    --------
    >>> childs = single_point_crossover(rnd, "AAAAA", "BBBBB")
    >>> print(childs)
    ['AABBB', 'BBAAA']
    """
    parent1_str: str = parent1.candidate if isinstance(parent1, Individual) else str(parent1)
    parent2_str: str = parent2.candidate if isinstance(parent2, Individual) else str(parent2)

    length: int = len(parent1_str)

    if length <= 1:
        return [parent1_str, parent2_str]

    crossover_point: int = random.randint(1, length - 1)

    child1: str = parent1_str[:crossover_point] + parent2_str[crossover_point:]
    child2: str = parent2_str[:crossover_point] + parent1_str[crossover_point:]

    return [child1, child2]



# --- Logic for BLOSUM Mutation ---

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
    # Return to uniform selection if matrix is ​​not available
    if C.BLOSUM62 is None:
        return list(C.AMINO_ACIDS), [1.0] * len(C.AMINO_ACIDS)
        
    targets: List[str]   = []
    weights: List[float] = []
    
    for target_aa in C.AMINO_ACIDS:
        try:
            # Retrieve the score
            # Use .get to handle pairs (target, original) if they are not symmetric or present
            score: int = C.BLOSUM62.get((original_aa, target_aa), 0)
        except (KeyError, ValueError):
            continue

        # Transforming the score into a positive weight (same logic as before)
        OFFSET  : int   = 5
        EXPONENT: int   = 2
        weight  : float = float(max(1, (score + OFFSET) ** EXPONENT)) 
        
        targets.append(target_aa)
        weights.append(weight)

    return targets, weights



# --- Adaptive Rate Helper Function ---

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
        e 'real_max_gen' per il calcolo.

    Returns
    -------
    `float`
        Il tasso di mutazione corrente.
    
    Examples
    --------
    # Se il GA è a 50 su 100 generazioni (a metà)
    # Rate = 0.30 - (0.30 - 0.05) * 0.5 = 0.175
    """
    # inspyred passes the EvolutionaryComputation object into args['_ec']
    ec                    = args.get('_ec', None)
    max_gen               = args.get('real_max_gen', C.MAX_GENERATIONS) # Must match MAX_GENERATIONS
    offset                = args.get('generation_offset', 0)
    initial_mutation_rate = args.get('initial_mutation_rate', C.INITIAL_MUTATION_RATE)
    final_mutation_rate   = args.get('final_mutation_rate', C.FINAL_MUTATION_RATE)
    
    current_gen     = ec.num_generations if ec else 0
    real_generation = current_gen + offset
    
    if ec is None:   return final_mutation_rate # Secure fallback  
    if max_gen == 0: return final_mutation_rate

    # Progresso da 0.0 a 1.0
    progress = min(1.0, real_generation / max_gen)
    
    # Linear Interpolation (Lerp)
    # Rate = Start - (Start - End) * Progress
    current_rate = initial_mutation_rate - (initial_mutation_rate - final_mutation_rate) * progress

    print_verbose(f"[get_adaptive_mutation_rate] current_rate: {current_rate}", to_print = args.get('verbose', False), lock = args.get('print_lock', None))
    return max(final_mutation_rate, current_rate)


# --- Mutation Operator ---

def blosum_peptide_mutator(
        random   : random.Random, 
        candidate: str,
        args     : Dict[str, Any]
) -> str:
    """
    Mutates a peptide sequence using probabilities derived from the BLOSUM62 matrix.

    Iterates through the amino acids of the sequence and, based on the mutation rate,
    replaces residues with chemically likely substitutes defined by the substitution matrix.

    Parameters
    ----------
    random : `random.Random`
        The random number generator instance.
    candidate : `str`
        The peptide sequence to mutate.
    args : `Dict[str, Any]`
        Arguments containing 'mutation_rate' and optionally the 'blosum_matrix'.

    Returns
    -------
    `str`
        The mutated peptide sequence.

    Examples
    --------
    >>> mutated = blosum_peptide_mutator(rnd, "ACDEF", {'mutation_rate': 0.1})
    >>> print(mutated)
    'ACDWF'
    """
    mutation_probability: float     = get_adaptive_mutation_rate(args)
    mutated_sequence    : List[str] = list(candidate)
    
    for i in range(len(mutated_sequence)):
        if random.random() <= mutation_probability:
            original_aa        : str = mutated_sequence[i]
            targets, weights         = get_blosum_weights(original_aa, args)
            new_amino_acid     : str = random.choices(targets, weights = weights, k = 1)[0]
            mutated_sequence[i]      = new_amino_acid

    return "".join(mutated_sequence)


# --- Chain Operator ---

def peptide_chain_variator(
    random: random.Random, 
    candidates: List[Any], # List of Individual or str, depending on inspyred
    args: Dict[str, Any]
) -> List[str]:
    """
    The main genetic variator operator for the evolutionary algorithm.

    Applies both Crossover (recombination) and Mutation to the population.
    1. Pairs parents sequentially.
    2. Applies crossover with probability `CROSSOVER_PROBABILITY`.
    3. Applies mutation to the resulting offspring using adaptive mutation rates.

    Parameters
    ----------
    random : `random.Random`
        The random number generator instance.
    candidates : `List[Individual]`
        The list of parent individuals selected for reproduction.
    args : `Dict[str, Any]`
        Arguments containing global parameters (probabilities, mutation rates).

    Returns
    -------
    `List[str]`
        The new list of offspring sequences to form the next generation.

    Examples
    --------
    >>> offspring = peptide_chain_variator(rnd, parents, args)
    """
    verbose = args.get('verbose', False)
    lock    = args.get('print_lock', None)

    print_verbose(f"[peptide_chain_variator] Input Candidates: {len(candidates)}", to_print = verbose, lock = lock)
    # Initialize the list for new children
    new_population: List[str] = []
    
    # Iterates over the list of candidates (parents) in steps of 2
    for i in range(0, len(candidates), 2):
        
        # Handles the case where there is only one remaining candidate (odd list)
        if i + 1 >= len(candidates):
            # If there is only one parent, we mutate it and add it to the new population
            lonely_parent = candidates[i]
            # We assume the lone parent is a string (or we extract it)
            lonely_seq = lonely_parent.candidate if isinstance(lonely_parent, Individual) else str(lonely_parent)
            
            mutated_child = blosum_peptide_mutator(random, lonely_seq, args)
            new_population.append(mutated_child)
            break 
            
        parent1 = candidates[i]
        parent2 = candidates[i+1]
        
        # Perform Crossover based on a probability
        if random.random() <= C.CROSSOVER_PROBABILITY:
            children_sequences: List[str] = single_point_crossover(random, parent1, parent2)
        else:
            # If crossover does not occur, the children are copies of their parents.
            parent1_str = parent1.candidate if isinstance(parent1, Individual) else str(parent1)
            parent2_str = parent2.candidate if isinstance(parent2, Individual) else str(parent2)
            children_sequences = [parent1_str, parent2_str]
        
        # Perform Mutation on each child
        for child_seq in children_sequences:
            mutated_child: str = blosum_peptide_mutator(random, child_seq, args)
            new_population.append(mutated_child)
        
    print_verbose(f"[peptide_chain_variator] Output Children: {len(new_population)}", to_print = verbose, lock = lock)
    if len(new_population) < len(candidates):
        print_warning(f"I'm losing people along the way! {len(new_population)} candidates  < {len(new_population)} children", lock = lock)
            
    return new_population

# --- Hydrophobucity ---
def get_hydrophobicity(peptide: str) -> float:
    """
    Calculates the average hydrophobicity index of a peptide sequence.

    Uses a specific scale (e.g., Monera et al. at pH 7) to sum the hydrophobicity
    values of all residues and computes the average.

    Parameters
    ----------
    peptide : `str`
        The input peptide sequence (e.g., "WWPYWW").

    Returns
    -------
    `float`
        The average hydrophobicity value. Positive values indicate hydrophobic, 
        negative values indicate hydrophilic.

    Examples
    --------
    >>> val = get_hydrophobicity("WWWW")
    >>> print(val)
    97.0
    """
    average_hydrophobicity = sum(C.AMINO_HYDROPHOBICITY_PH7.get(aa, 0) for aa in peptide) / len(peptide)
    return average_hydrophobicity
