# peptide_operators.py
import random

from constants   import AMINO_ACIDS, BLOSUM62

from inspyred.ec import Individual                          # Importa Individual per l'hinting corretto
from typing      import List, Tuple, Dict, Any, TYPE_CHECKING

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
    if BLOSUM62 is None:
        return list(AMINO_ACIDS), [1.0] * len(AMINO_ACIDS)
        
    targets: List[str]   = []
    weights: List[float] = []
    
    for target_aa in AMINO_ACIDS:
        try:
            # Recupera il punteggio
            # Utilizza .get per gestire le coppie (target, original) se non sono simmetriche o presenti
            score: int = BLOSUM62.get((original_aa, target_aa), 0)
        except (KeyError, ValueError):
            continue

        # Trasformazione del punteggio in un peso positivo (stessa logica di prima)
        OFFSET  : int   = 5
        EXPONENT: int   = 2
        weight  : float = float(max(1, (score + OFFSET) ** EXPONENT)) 
        
        targets.append(target_aa)
        weights.append(weight)

    return targets, weights



# --- Operatore di Mutazione ---

def blosum_peptide_mutator(
        random   : random.Random, 
        candidate: str,
        args     : Dict[str, Any]
) -> str:
    """
    Esegue la Mutazione a Sostituzione con probabilità pesate da BLOSUM62.

    Per ogni posizione nella sequenza, viene applicata una mutazione con 
    probabilità `mutation_rate`. Se la mutazione avviene, l'aminoacido 
    di sostituzione è scelto in base alla distribuzione pesata da BLOSUM62.

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale fornita dal framework `inspyred`.
    candidate : `str`
        La sequenza peptidica (stringa) da mutare.
    args : `Dict[str, Any]`
        Dizionario di argomenti, deve contenere la chiave `'mutation_rate'` (`float`).

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
    mutation_probability: float = args.get('mutation_rate') 
    
    if mutation_probability is None:
        from constants import MUTATION_RATE
        mutation_probability = MUTATION_RATE
        
    mutated_sequence: List[str] = list(candidate)
    
    for i in range(len(mutated_sequence)):
        if random.random() < mutation_probability:
            original_aa: str = mutated_sequence[i]

            targets, weights = get_blosum_weights(original_aa)
            
            # random.choices seleziona un elemento basato sui pesi
            new_amino_acid: str = random.choices(targets, weights=weights, k=1)[0]
            
            mutated_sequence[i] = new_amino_acid

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
        
        # 1. Esegui il Crossover (Restituisce una lista di 2 stringhe)
        children_sequences: List[str] = single_point_crossover(random, parent1, parent2)
        
        # 2. Esegui la Mutazione su ciascun figlio
        for child_seq in children_sequences:
            mutated_child: str = blosum_peptide_mutator(random, child_seq, args)
            new_population.append(mutated_child)
            
    return new_population
