# main.py
import random
import inspyred.ec
from inspyred.ec import terminators, observers
from typing import Optional, Any

# Importa i moduli del progetto
from constants import POPULATION_SIZE, MAX_GENERATIONS, PEPTIDE_LENGTH
from ga_problem import peptide_generator, evaluate_peptide_binding
from peptide_operators import peptide_chain_variator, single_point_crossover, blosum_peptide_mutator
from inspyred.ec import Individual # Necessario per l'hinting di ritorno

def run_peptide_ga() -> Individual:
    """
    Esegue l'Algoritmo Genetico (GA) per l'ottimizzazione della sequenza peptidica.

    Configura e lancia il framework `inspyred`, utilizzando i generatori, 
    i variatori (crossover e mutazione) e gli evaluatori specifici del problema.

    Parameters
    ----------
    Nessuno. I parametri del GA sono definiti in `constants.py`.

    Returns
    -------
    `inspyred.ec.Individual`
        L'individuo (sequenza peptidica) migliore trovato alla fine dell'evoluzione,
        con la massima fitness (energia di legame minima).

    Examples
    --------
    >>> best_peptide = run_peptide_ga()
    >>> print(f"Best sequence: {best_peptide.candidate}")
    Best sequence: LWRTAVI # Esempio di miglior sequenza
    """
    
    rand = random.Random()
    rand.seed(77) # Seed per la riproducibilità

    print(f"--- Avvio Algoritmo Genetico (Peptide Length: {PEPTIDE_LENGTH}) ---")
    print(f"Popolazione: {POPULATION_SIZE}, Generazioni: {MAX_GENERATIONS}")

    # 1. Setup dell'Algoritmo Evolutivo (EA)
    ea = inspyred.ec.EvolutionaryComputation(rand)
    
    # Imposta i parametri e le funzioni
    ea.selector = inspyred.ec.selectors.tournament_selection
    ea.replacer = inspyred.ec.replacers.generational_replacement
    
    # 2. Assegna le tue funzioni di Crossover e Mutazione
    ea.variator = peptide_chain_variator
    
    # 3. Assegna le funzioni di Generazione e Valutazione
    ea.generator = peptide_generator
    ea.evaluator = evaluate_peptide_binding
    
    # 4. Imposta il terminatore e l'osservatore
    ea.terminator = terminators.generation_termination
    
    # Esegui l'EA
    final_pop = ea.evolve(
        generator        = peptide_generator,
        evaluator        = evaluate_peptide_binding,
        num_elites       = 1,                                       # Conserva il miglior individuo dalla generazione precedente
        
        generator_params = {'peptide_length': PEPTIDE_LENGTH},      # Parametri passati in args ai variatori
        variator_params  = {'max_generations': MAX_GENERATIONS},    # Parametri passati a peptide_chain_variator (max_generations qui per il calcolo adattivo)
        num_generations  = MAX_GENERATIONS,
        pop_size         = POPULATION_SIZE,
        maximize         = False,                                   # Vina: Più negativo è meglio. Quindi vogliamo minimizzare
    )

    best_individual: Individual = max(final_pop)
    
    print("\n--- Risultati Finali ---")
    print(f"Miglior sequenza trovata: {best_individual.candidate}")
    print(f"Miglior fitness (Energia di legame Vina stimata): {best_individual.fitness:.3f} kcal/mol")

    return best_individual

if __name__ == '__main__':
    best = run_peptide_ga()
