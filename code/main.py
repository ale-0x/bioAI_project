# main.py
import argparse
import inspyred.ec
import multiprocessing
import os
import random
import shutil
import matplotlib.pyplot as plt

from inspyred.ec            import Individual, terminators
from inspyred.ec.evaluators import parallel_evaluation_mp
from sys                    import argv
from typing                 import Any, Optional

import constants as C

from ga_problem             import evaluate_peptide_binding, peptide_generator
from utils                  import COL_LIGHT_BLUE, color_text, print_info, print_error
from peptide_operators      import blosum_peptide_mutator, peptide_chain_variator, single_point_crossover

def parse_arguments() -> argparse.Namespace:
    """
    Handles all input and initial settings.
    
    Return
    ----------
    argparse.Namespace
        Simple object for storing attributes.
    """

    description = "Evolutionary Algorithm for Peptide-Protein Docking Optimization using AutoDock Vina."
    parser = argparse.ArgumentParser(description = description, prog = argv[0])

    parser.add_argument('receptor_name'                , action = 'store'                                                      , help = f'Il nome identificativo o codice PDB della proteina bersaglio (recettore).')
    parser.add_argument('receptor_file'                , action = 'store'                                                      , help = f'Il percorso relativo al file della proteina già preparato in formato PDBQT (necessario per Vina, include cariche e idrogeni).')
    
    parser.add_argument('-j', '--job_id'               , action = 'store'     , type = str  , default = "local_test"           , help = f'Identificativo univoco del Job (es. SLURM_JOB_ID). Usato per creare cartelle temporanee univoche.')
    parser.add_argument('-c', '--cpus'                 , action = 'store'     , type = int  , default = C.CPUS                 , help = f'Numero totale di CPU da utilizzare. [Default: {C.CPUS}]')
    parser.add_argument('-l', '--peptide_length'       , action = 'store'     , type = int  , default = C.PEPTIDE_LENGTH       , help = f'Lunghezza della sequenza peptidica in aminoacidi. [Default: {C.PEPTIDE_LENGTH}]')
    parser.add_argument('-n', '--population_size'      , action = 'store'     , type = int  , default = C.POPULATION_SIZE      , help = f'Dimensione della popolazione per generazione. [Default: {C.POPULATION_SIZE}]')
    parser.add_argument('-g', '--generations'          , action = 'store'     , type = int  , default = C.MAX_GENERATIONS      , help = f'Numero massimo di generazioni. [Default: {C.MAX_GENERATIONS}]')
    parser.add_argument(      '--initial_mutation_rate', action = 'store'     , type = float, default = C.INITIAL_MUTATION_RATE, help = f'Probabilità che un aminoacido muti all\'inizio della simulazione. [Default: {C.INITIAL_MUTATION_RATE}]')
    parser.add_argument(      '--final_mutation_rate'  , action = 'store'     , type = float, default = C.FINAL_MUTATION_RATE  , help = f'Probabilità di mutazione alla fine del processo. [Default: {C.FINAL_MUTATION_RATE}]')
    parser.add_argument(      '--temp_dir_base'        , action = 'store'     , type = str  , default = "../resources/tmp"     , help = f'La cartella principale destinata a contenere tutti i dati temporanei generati durante l\'esecuzione. [Default: "../resources/tmp"]')
    parser.add_argument('-o', '--output'               , action = 'store'     , type = str  , default = "result.txt"           , help = f'Il nome del file finale in cui verranno scritti i risultati migliori (sequenza e fitness) al termine dell\'evoluzione. [Default: "result.txt"]')

    # Parametri Docking di Vina (Sovrascrivono constants.py)
    parser.add_argument('-x', '--center_x'             , action = 'store'     , type = float, default = C.CENTER_X             , help = f'Coordinata cartesiana X del centro esatto della Grid Box (la zona della proteina dove si cercherà il legame con il peptide). [Default: {C.CENTER_X}]')
    parser.add_argument('-y', '--center_y'             , action = 'store'     , type = float, default = C.CENTER_Y             , help = f'Coordinata cartesiana Z del centro esatto della Grid Box (la zona della proteina dove si cercherà il legame con il peptide). [Default: {C.CENTER_Y}]')
    parser.add_argument('-z', '--center_z'             , action = 'store'     , type = float, default = C.CENTER_Z             , help = f'Coordinata cartesiana Y del centro esatto della Grid Box (la zona della proteina dove si cercherà il legame con il peptide). [Default: {C.CENTER_Z}]')
    parser.add_argument('-X', '--size_x'               , action = 'store'     , type = int  , default = C.SIZE_X               , help = f'Dimensione (in Ångström) della scatola di ricerca l\'asse X. [Default: {C.SIZE_X}]')
    parser.add_argument('-Y', '--size_y'               , action = 'store'     , type = int  , default = C.SIZE_Y               , help = f'Dimensione (in Ångström) della scatola di ricerca l\'asse Y. [Default: {C.SIZE_Y}]')
    parser.add_argument('-Z', '--size_z'               , action = 'store'     , type = int  , default = C.SIZE_Z               , help = f'Dimensione (in Ångström) della scatola di ricerca l\'asse Z. [Default: {C.SIZE_Z}]')
    parser.add_argument('-e', '--exhaustiveness'       , action = 'store'     , type = int  , default = C.EXHAUSTIVENESS       , help = f'Definisce quanto deve essere accurata la ricerca globale di Vina. [Default: {C.EXHAUSTIVENESS}]')
    parser.add_argument(      '--vina_exe_path'        , action = 'store'     , type = str  , default = C.VINA_EXE_PATH        , help = f'Il comando o il percorso assoluto per invocare l\'eseguibile di AutoDock Vina nel sistema. [Default: {C.VINA_EXE_PATH}]')

    parser.add_argument(      '--no_delete'            , action = 'store_true'              , default = 'False'                , help = 'Se impostato, lo script non cancellerà i file temporanei creati durante e dopo i calcoli.')
    parser.add_argument('-v', '--verbose'              , action = 'store_true'              , default = 'False'                , help = 'Mostrerà molte più informazioni nel terminale.')
    
    parser.add_argument(      '--version'              , action = 'version'   , version = '%(prog)s - Version {} created by {}'.format(C.__version__, C.__author__), help = "Output version information and exit.")
    
    return parser.parse_args()

def print_arguments(arguments: argparse.Namespace, name = "", to_print: bool = True) -> str | None:
    output = list()
    output.append(color_text(COL_LIGHT_BLUE, "\n" + "="*40))
    output.append(color_text(COL_LIGHT_BLUE, f" CONFIGURAZIONE JOB: {name}"))
    output.append(color_text(COL_LIGHT_BLUE, "="*40))
        
    for key, value in vars(arguments).items():
        output.append(color_text(COL_LIGHT_BLUE, f"{key.ljust(22)}: {value}"))
    output.append(color_text(COL_LIGHT_BLUE, "="*40 + "\n"))

    if to_print:
        print("\n".join(output))
    else:
        return "\n".join(output)


def run_peptide_ga() -> None:
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

    options = parse_arguments()
    RECEPTOR_NAME        : str   = options.receptor_name
    RECEPTOR_FILE        : str   = options.receptor_file
    
    JOB_ID               : str   = options.job_id
    CPUS                 : int   = options.cpus
    PEPTIDE_LENGTH       : int   = options.peptide_length
    POPULATION_SIZE      : int   = options.population_size
    MAX_GENERATIONS      : int   = options.generations
    INITIAL_MUTATION_RATE: float = options.initial_mutation_rate
    FINAL_MUTATION_RATE  : float = options.final_mutation_rate
    TEMP_DIR_BASE        : str   = options.temp_dir_base
    OUTPUT               : str   = options.output
    
    CENTER_X             : float = options.center_x
    CENTER_Y             : float = options.center_y
    CENTER_Z             : float = options.center_z
    SIZE_X               : int   = options.size_x
    SIZE_Y               : int   = options.size_y
    SIZE_Z               : int   = options.size_z
    EXHAUSTIVENESS       : int   = options.exhaustiveness
    VINA_EXE_PATH        : str   = options.vina_exe_path

    NO_DELETE            : bool  = options.no_delete
    VERBOSE              : bool  = options.verbose

    print_arguments(options, JOB_ID)

    rand = random.Random()
    rand.seed(C.SEED)                # Seed per la riproducibilità

    job_temp_dir = os.path.join(TEMP_DIR_BASE, f"bioai_ga_{JOB_ID}")
    os.makedirs(job_temp_dir, exist_ok = True)

    print(f"--- Avvio Algoritmo Genetico (Peptide Length: {PEPTIDE_LENGTH}) ---")
    print(f"Popolazione: {POPULATION_SIZE}, Generazioni: {MAX_GENERATIONS}")

    # 1. Setup dell'Algoritmo Evolutivo (EA)
    ea = inspyred.ec.EvolutionaryComputation(rand)

    #1.1 Aggiungi un Observer
    ea.observer = inspyred.ec.observers.file_observer # add required arguments
    
    # Imposta i parametri e le funzioni
    ea.selector = inspyred.ec.selectors.tournament_selection
    ea.replacer = inspyred.ec.replacers.plus_replacement
    
    # 2. Assegna le tue funzioni di Crossover e Mutazione
    ea.variator = peptide_chain_variator
    
    # 3. Assegna le funzioni di Generazione e Valutazione
    ea.generator = peptide_generator
    # ea.evaluator = evaluate_peptide_binding
    cpus = multiprocessing.cpu_count()
    print(f"Utilizzo {cpus} CPU in parallelo per il docking.")
    
    ea.evaluator = parallel_evaluation_mp
    
    # 4. Imposta il terminatore e l'osservatore
    ea.terminator = terminators.generation_termination

    ga_config = {
        'receptor_file'        : RECEPTOR_FILE,
        'initial_mutation_rate': INITIAL_MUTATION_RATE,
        'final_mutation_rate'  : FINAL_MUTATION_RATE,
        'center_x'             : CENTER_X,
        'center_y'             : CENTER_Y,
        'center_z'             : CENTER_Z,
        'size_x'               : SIZE_X,                # Se vuoi puoi parametrizzare anche questi
        'size_y'               : SIZE_Y,
        'size_z'               : SIZE_Z,
        'temp_dir'             : job_temp_dir,          # Passiamo la cartella creata sopra
        'exhaustiveness'       : EXHAUSTIVENESS,
        'vina_exe_path'        : VINA_EXE_PATH,
        'no_delete'            : NO_DELETE,
        'verbose'              : VERBOSE,
        'max_generations'      : MAX_GENERATIONS        # <--- Fondamentale per terminators e variators
    }

    statistics_file = open(f"ga_observer_{JOB_ID}.csv", "w")
    individuals_file = open(f"ga_individuals_{JOB_ID}.csv", "w")
    
    # Esegui l'EA
    try:
        final_pop = ea.evolve(
            generator        = peptide_generator,
            evaluator        = parallel_evaluation_mp,
            observer         = ea.observer,
            mp_evaluator     = evaluate_peptide_binding,
            mp_num_cpus      = cpus,
            num_elites       = 1,                                       # Conserva il miglior individuo dalla generazione precedente
            
            statistics_file  = statistics_file,
            individuals_file = individuals_file,                        # Parametri passati all'observer
            generator_params = {'peptide_length': PEPTIDE_LENGTH},      # Parametri passati in args ai variatori
            variator_params  = {},                                      # Parametri passati a peptide_chain_variator (max_generations qui per il calcolo adattivo)
            num_generations  = MAX_GENERATIONS,
            pop_size         = POPULATION_SIZE,
            maximize         = False,                                   # Vina: Più negativo è meglio. Quindi vogliamo minimizzare

            **ga_config
        )

        print(color_text(COL_LIGHT_BLUE, "\n" + "=" * 40))
        print(color_text(COL_LIGHT_BLUE, f"           POPOLAZIONE FINALE           "))
        print(color_text(COL_LIGHT_BLUE, "=" * 40))
        for i, individual in enumerate(final_pop):
            assert isinstance(individual, Individual)
            print(color_text(COL_LIGHT_BLUE, f" Individual {i}\n\t- Candidate : {individual.candidate}"))
            print(color_text(COL_LIGHT_BLUE, f"\t- Fitness   : {individual.fitness}"))
            print(color_text(COL_LIGHT_BLUE, f"\t- Birthdate : {individual.birthdate}\n"))
            if i < len(final_pop) - 1:
                print(color_text(COL_LIGHT_BLUE, "-" * 40))
        print(color_text(COL_LIGHT_BLUE, "=" * 40))
        
        best_individual: Individual = min(final_pop) # choose best fitness (minimize energy)
        
        print("\n--- Risultati Finali ---")
        print(f"Miglior sequenza trovata: {best_individual.candidate}")
        print(f"Miglior fitness (Energia di legame Vina stimata): {best_individual.fitness:.3f} kcal/mol")

        # Salva output finale
        with open(OUTPUT, "w") as f:
            f.write(f"ID       : {JOB_ID}\n")
            f.write(f"Receptor : {RECEPTOR_FILE}\n")
            f.write(f"Sequence : {best_individual.candidate}\n")
            f.write(f"Fitness  : {best_individual.fitness}\n")

            f.write("")
            f.write(print_arguments(options, JOB_ID, False))

        statistics_file.close()
        individuals_file.close()

        plot_file = f"plots/generation_plot_{JOB_ID}.png"
        os.makedirs("plots", exist_ok=True)  # Ensure the directory exists
        inspyred.ec.analysis.generation_plot(open(f"ga_observer_{JOB_ID}.csv", 'r'))
        plt.savefig(plot_file)
        plt.close()
        print(f"Generation plot saved to {plot_file}")

    # except Exception as e:
        # print(e)

    finally:
        # Rimuove l'intera cartella temporanea del job alla fine
        if os.path.exists(job_temp_dir) and not NO_DELETE:
            print(f"Pulizia cartella temporanea: {job_temp_dir} ...", end = "")
            shutil.rmtree(job_temp_dir)
            print("Done")
    # Close the files opened for statistics and individuals

if __name__ == '__main__':
    run_peptide_ga()
