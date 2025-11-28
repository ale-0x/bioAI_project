# main.py
import random
import inspyred.ec
from inspyred.ec import terminators, observers
from typing import Optional, Any

# Importa i moduli del progetto
from constants import POPULATION_SIZE, MAX_GENERATIONS, MUTATION_RATE, PEPTIDE_LENGTH
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
    print(f"Popolazione: {POPULATION_SIZE}, Generazioni: {MAX_GENERATIONS}, Mutazione/AA: {MUTATION_RATE}")

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
    # ea.observer   = observers.best_observer
    
    # Esegui l'EA
    final_pop = ea.evolve(
        generator=peptide_generator,
        evaluator=evaluate_peptide_binding,
        
        generator_params = {'peptide_length': PEPTIDE_LENGTH},
        variator_params  = {'mutation_rate' : MUTATION_RATE},   # Parametri passati a peptide_chain_variator
        num_generations  = MAX_GENERATIONS,
        pop_size         = POPULATION_SIZE,
        maximize         = True 
    )

    best_individual: Individual = max(final_pop)
    
    print("\n--- Risultati Finali ---")
    print(f"Miglior sequenza trovata: {best_individual.candidate}")
    print(f"Miglior fitness (Energia di legame Vina stimata): {best_individual.fitness:.3f} kcal/mol")

    return best_individual

if __name__ == '__main__':
    best = run_peptide_ga()
    