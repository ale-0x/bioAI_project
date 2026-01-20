# main.py
import argparse
import inspyred.ec
import multiprocessing
import pickle
import random
import shutil
import time

from datetime               import datetime, timedelta
from inspyred.ec            import Individual, terminators
from inspyred.ec.evaluators import parallel_evaluation_mp
from pathlib                import Path
from sys                    import argv

import constants as C

from ga_problem             import evaluate_peptide_binding, peptide_generator
from utils                  import custom_file_observer, generation_tracker_observer, checkpoint_observer, print_arguments, print_error, print_info, print_verbose, print_warning, restore_context, time_based_termination
from peptide_operators      import peptide_chain_variator, get_hydrophobicity
from plots                  import plot_energy_vs_hydrophobicity, plot_energy_vs_hydrophobicity_2, plot_observer_statistics, plot_observer_statistics_2

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

    parser.add_argument('receptor_name'                , action = 'store'                                                      , help = f'The identifying name or PDB code of the target protein (receptor).')
    parser.add_argument('receptor_file'                , action = 'store'                                                      , help = f'The relative path to the protein file already prepared in PDBQT format (required for Vina, includes charges and hydrogens).')
    
    parser.add_argument('-j', '--job_id'               , action = 'store'     , type = str  , default = "local_test"           , help = f'Unique job identifier (e.g., SLURM_JOB_ID). Used to create unique temporary folders.')
    parser.add_argument('-c', '--cpus'                 , action = 'store'     , type = int  , default = C.CPUS                 , help = f'Total number of CPUs to use. [Default: {C.CPUS}]')
    parser.add_argument('-l', '--peptide_length'       , action = 'store'     , type = int  , default = C.PEPTIDE_LENGTH       , help = f'Length of the peptide sequence in amino acids. [Default: {C.PEPTIDE_LENGTH}]')
    parser.add_argument('-n', '--population_size'      , action = 'store'     , type = int  , default = C.POPULATION_SIZE      , help = f'Population size by generation. [Default: {C.POPULATION_SIZE}]')
    parser.add_argument('-g', '--generations'          , action = 'store'     , type = int  , default = C.MAX_GENERATIONS      , help = f'Maximum number of generations. [Default: {C.MAX_GENERATIONS}]')
    parser.add_argument(      '--initial_mutation_rate', action = 'store'     , type = float, default = C.INITIAL_MUTATION_RATE, help = f'Probability that an amino acid will mutate at the start of the simulation. [Default: {C.INITIAL_MUTATION_RATE}]')
    parser.add_argument(      '--final_mutation_rate'  , action = 'store'     , type = float, default = C.FINAL_MUTATION_RATE  , help = f'Probability of mutation at the end of the process. [Default: {C.FINAL_MUTATION_RATE}]')
    parser.add_argument(      '--hydrophobicity_weight', action = 'store'     , type = float, default = C.HYDROPHOBICITY_WEIGHT, help = f'Weight of hydrophobicity. [Default: {C.HYDROPHOBICITY_WEIGHT}]')
    parser.add_argument(      '--temp_dir_base'        , action = 'store'     , type = str  , default = "../resources/tmp"     , help = f'The main folder intended to contain all temporary data generated during execution. [Default: "../resources/tmp"]')
    parser.add_argument('-o', '--output'               , action = 'store'     , type = str  , default = "result.txt"           , help = f'The name of the final folder where the best results (sequence and fitness) will be written at the end of the evolution. [Default: "result_job_id"]')
    parser.add_argument('-d', '--deadline'             , action = 'store'     , type = str  , default = "23:55:00"             , help = f'The maximum time limit for the job to run (format HH:MM:SS). The script will attempt to finish execution before this time. [Default: "23:55:00"]')
    parser.add_argument('-r', '--resume'               , action = 'store'     , type = str  , default = None                   , help = f'Path to the .pkl file from which to resume the simulation')

    # Parametri Docking di Vina (Sovrascrivono constants.py)
    parser.add_argument('-x', '--center_x'             , action = 'store'     , type = float, default = C.CENTER_X             , help = f'Cartesian X coordinate of the exact center of the Grid Box (the region of the protein where the peptide will be bound). [Default: {C.CENTER_X}]')
    parser.add_argument('-y', '--center_y'             , action = 'store'     , type = float, default = C.CENTER_Y             , help = f'Cartesian Z coordinate of the exact center of the Grid Box (the region of the protein where the peptide will be bound). [Default: {C.CENTER_Y}]')
    parser.add_argument('-z', '--center_z'             , action = 'store'     , type = float, default = C.CENTER_Z             , help = f'Cartesian Y coordinate of the exact center of the Grid Box (the region of the protein where the peptide will be bound). [Default: {C.CENTER_Z}]')
    parser.add_argument('-X', '--size_x'               , action = 'store'     , type = int  , default = C.SIZE_X               , help = f'Size (in Ångström) of the search box on the X axis. [Default: {C.SIZE_X}]')
    parser.add_argument('-Y', '--size_y'               , action = 'store'     , type = int  , default = C.SIZE_Y               , help = f'Size (in Ångström) of the Y-axis search box. [Default: {C.SIZE_Y}]')
    parser.add_argument('-Z', '--size_z'               , action = 'store'     , type = int  , default = C.SIZE_Z               , help = f'Size (in Ångström) of the search box on the Z axis. [Default: {C.SIZE_Z}]')
    parser.add_argument('-e', '--exhaustiveness'       , action = 'store'     , type = int  , default = C.EXHAUSTIVENESS       , help = f'Defines how accurate Vina\'s global search should be. [Default: {C.EXHAUSTIVENESS}]')
    parser.add_argument(      '--vina_exe_path'        , action = 'store'     , type = str  , default = C.VINA_EXE_PATH        , help = f'The command or absolute path to invoke the AutoDock Vina executable on the system. [Default: {C.VINA_EXE_PATH}]')

    parser.add_argument(      '--no_delete'            , action = 'store_true'              , default = 'False'                , help = 'If set, the script will not delete temporary files created during and after calculations.')
    parser.add_argument('-v', '--verbose'              , action = 'store_true'              , default = 'False'                , help = 'It will show much more information in the terminal.')
    
    parser.add_argument(      '--version'              , action = 'version'   , version = '%(prog)s - Version {} created by {}'.format(C.__version__, C.__author__), help = "Output version information and exit.")
    
    return parser.parse_args()

def run_peptide_ga() -> None:
    """
    Main execution function for the Peptide Evolutionary Algorithm.

    Orchestrates the entire workflow:
    1. Sets up the output directories and logging.
    2. Initializes the Random Number Generator.
    3. Configures the `inspyred` Evolutionary Computation (EC) engine.
    4. Defines the pipeline: Generator, Evaluator, Variator, Replacer, Terminator.
    5. Runs the evolution for `MAX_GENERATIONS`.
    6. Saves the results and generates analysis plots.
    """
    START_TIME = time.time()

    options = parse_arguments()

    RECEPTOR_NAME        : str         = options.receptor_name
    RECEPTOR_FILE        : Path        = Path(options.receptor_file)
    if not RECEPTOR_FILE.exists():
        print_error(f"The receptor file does not exist: {RECEPTOR_FILE.resolve()}", code = -1)
        
    JOB_ID               : str         = options.job_id
    CPUS                 : int         = options.cpus
    PEPTIDE_LENGTH       : int         = options.peptide_length
    POPULATION_SIZE      : int         = options.population_size
    MAX_GENERATIONS      : int         = options.generations
    INITIAL_MUTATION_RATE: float       = options.initial_mutation_rate
    FINAL_MUTATION_RATE  : float       = options.final_mutation_rate
    HYDROPHOBICITY_WEIGHT: float       = options.hydrophobicity_weight
    TEMP_DIR_BASE        : Path        = Path(options.temp_dir_base)
    OUTPUT               : Path        = Path(f"{options.output}_{JOB_ID}")
    RESUME               : Path | None = Path(options.resume) if options.resume else None

    try:
        h, m, s = map(int, options.deadline.split(':')) # HH:MM:SS format
        DEADLINE : timedelta = timedelta(hours = h, minutes = m, seconds = s)
        
        # Let's calculate the absolute final timestamp
        END_TIME : float = START_TIME + DEADLINE.total_seconds()
        
        BUFFER_SECONDS = 0 
        GLOBAL_VINA_DEADLINE = END_TIME - BUFFER_SECONDS
        
    except ValueError:
        print_error("Deadline format error. Use HH:MM:SS (e.g., 23:55:00).", code = -1)

    CENTER_X             : float       = options.center_x
    CENTER_Y             : float       = options.center_y
    CENTER_Z             : float       = options.center_z
    SIZE_X               : int         = options.size_x
    SIZE_Y               : int         = options.size_y
    SIZE_Z               : int         = options.size_z
    EXHAUSTIVENESS       : int         = options.exhaustiveness
    VINA_EXE_PATH        : Path        = Path(options.vina_exe_path)

    NO_DELETE            : bool        = options.no_delete
    VERBOSE              : bool        = options.verbose

    print_arguments(options, JOB_ID)

    rand = random.Random()
    # rand.seed(C.SEED)                # Seed for reproducibility
    

    job_temp_dir = TEMP_DIR_BASE / f"bioai_ga_{JOB_ID}"
    if not job_temp_dir.exists():
        print_verbose(f"[run_peptide_ga] Creating a temporary directory '{job_temp_dir}' ...", to_print = VERBOSE)
        job_temp_dir.mkdir(parents = True, exist_ok = True)
        if job_temp_dir.exists():
            print_verbose(f"[run_peptide_ga] Done.", to_print = VERBOSE)
        else:
            print_error(f"[run_peptide_ga] Temporary directory '{job_temp_dir}' not created!", code = -1)
    
    if not OUTPUT.exists():
        print_verbose(f"[run_peptide_ga] Creating output directory '{OUTPUT}' ...", to_print = VERBOSE)
        OUTPUT.mkdir(parents = True, exist_ok = True)
        if OUTPUT.exists():
            print_verbose(f"[run_peptide_ga] Done.", to_print = VERBOSE)
        else:
            print_error(f"[run_peptide_ga] Output directory '{OUTPUT}' not created!", code = -1)
    

    print(f"--- Start Genetic Algorithm (Peptide Length: {PEPTIDE_LENGTH}) ---")
    print(f"Start Time             : {datetime.fromtimestamp(START_TIME).strftime("%d/%m/%Y - %H:%M:%S")}")
    print(f"Max Evolution Duration : {datetime.fromtimestamp(END_TIME).strftime("%d/%m/%Y - %H:%M:%S")}")
    print(f"Population Dimension   : {POPULATION_SIZE}")
    print(f"Max Generations        : {MAX_GENERATIONS}")

    # --- CHECKPOINT / RESUME LOGIC ---
    initial_population = []
    generations_done   = 0
    initial_cache_data = {}
    checkpoint_path    = OUTPUT / f"checkpoint_{JOB_ID}.pkl"

    if RESUME:
        if RESUME.exists():
            print(f"RESUME mode detected by: {RESUME}")
            try:
                loaded_gen, loaded_pop, loaded_cache, loaded_rand = restore_context(
                    checkpoint_path = RESUME,
                    new_temp_dir    = job_temp_dir,
                    new_job_id      = JOB_ID,
                    new_results_dir = OUTPUT
                )

                print_info(f"[RESUME] Restore completed successfully.")
                generations_done   = loaded_gen
                initial_population = loaded_pop
                initial_cache_data = loaded_cache
                if loaded_rand:
                    rand.setstate(loaded_rand)
                    print_info(f"[RESUME] Previous random component restored.")
                
                print_info(f"[RESUME] Restore {len(initial_population)} solutions from Gen {generations_done}")
            except Exception as e:
                print_error(f"Error reading checkpoint: {e}", code = -1)
        else:
            print_error(f"Resume file '{RESUME}' not found. I'm starting from scratch.")

    # Calculate how much is missing
    gens_to_run = MAX_GENERATIONS - generations_done

    # Genetic Algorithm (GA) Setup
    ea = inspyred.ec.GA(rand)

    ea.observer = [
        # inspyred.ec.observers.file_observer,
        custom_file_observer,
        generation_tracker_observer,
        checkpoint_observer,
    ]
    
    ea.selector = inspyred.ec.selectors.tournament_selection
    ea.replacer = inspyred.ec.replacers.plus_replacement
    
    # Assign your Crossover and Mutation functions
    ea.variator = peptide_chain_variator
    
    # Assign Generation and Evaluation functions
    ea.generator = peptide_generator
    # ea.evaluator = evaluate_peptide_binding
    cpus = multiprocessing.cpu_count()
    print(f"Using {cpus} CPU in parallel for docking.")
    
    ea.evaluator = parallel_evaluation_mp
    
    # Set terminator and observer
    ea.terminator = [
        time_based_termination,
        terminators.generation_termination,
    ]

    with multiprocessing.Manager() as manager:

        multiprocessing_cache = manager.dict(initial_cache_data)
        generation_num        = manager.Value('i', -1)  # Shared generation counter
        print_lock            = manager.Lock()
        
        ga_config = {
            'receptor_file'        : RECEPTOR_FILE,
            'initial_mutation_rate': INITIAL_MUTATION_RATE,
            'final_mutation_rate'  : FINAL_MUTATION_RATE,
            'hydrophobicity_weight': HYDROPHOBICITY_WEIGHT,
            'center_x'             : CENTER_X,
            'center_y'             : CENTER_Y,
            'center_z'             : CENTER_Z,
            'size_x'               : SIZE_X,
            'size_y'               : SIZE_Y,
            'size_z'               : SIZE_Z,
            'temp_dir'             : job_temp_dir,
            'exhaustiveness'       : EXHAUSTIVENESS,
            'vina_exe_path'        : VINA_EXE_PATH,
            'no_delete'            : NO_DELETE,
            'verbose'              : VERBOSE,
            'real_max_gen'         : MAX_GENERATIONS,
            'max_generations'      : gens_to_run,        # <--- Essential for terminators and variators
            'generation_offset'    : generations_done,
            'job_id'               : JOB_ID,
            'peptide_length'       : PEPTIDE_LENGTH,
            'multiprocessing_cache': multiprocessing_cache,
            'global_deadline'      : GLOBAL_VINA_DEADLINE,
            'tournament_size'      : 2,
            'generation_num'       : generation_num,
            'print_lock'           : print_lock,
            'results_dir'          : OUTPUT,
            'checkpoint_path'      : checkpoint_path,
            'rand'                 : rand,
        }


        statistics_filepath  = OUTPUT / f"ga_observer_{JOB_ID}.csv"
        individuals_filepath = OUTPUT / f"ga_individuals_{JOB_ID}.csv"

        statistics_file  = statistics_filepath .open('a' if statistics_filepath .exists() and RESUME else 'w')
        individuals_file = individuals_filepath.open('a' if individuals_filepath.exists() and RESUME else 'w')
        
        # Run the GA
        try:
            final_pop = ea.evolve(
                generator        = peptide_generator,                             # Initial solution generation function
                evaluator        = parallel_evaluation_mp,                        # Solution evaluation function
                observer         = ea.observer,                                   # Observation function at each generation
                mp_evaluator     = evaluate_peptide_binding,                      # Custom parallel evaluation function
                mp_num_cpus      = cpus,                                          # Number of CPUs to use for multiprocessing
                num_elites       = 0,                                             # Best n individuals from the previous generation to retain
                      
                statistics_file  = statistics_file,                               
                individuals_file = individuals_file,                              # Parameters passed to the observer
                generator_params = {'peptide_length': PEPTIDE_LENGTH},            # Parameters passed in args to the variators
                variator_params  = {},                                            # Parameters passed to peptide_chain_variator (max_generations here for adaptive computation)
                num_generations  = gens_to_run,                                   # Number of generations
                pop_size         = POPULATION_SIZE,                               # Population size
                num_offspring    = POPULATION_SIZE,                               # Number of offspring per generation
                maximize         = False,                                         # Vina: The more negative the better. So we want to minimize
                seeds            = [ind.candidate for ind in initial_population], # Initial solutions in resume mode

                **ga_config
            )

            print("\n" + "=" * 40)
            print(f"            FINAL POPULATION            ")
            print("=" * 40)
            for i, individual in enumerate(final_pop):
                assert isinstance(individual, Individual)
                print(f" Individual {i}\n\t- Candidate        : {individual.candidate}")
                print(f"\t- Fitness          : {individual.fitness}")
                print(f"\t- Binding Energy   : {individual.fitness - (get_hydrophobicity(individual.candidate) * HYDROPHOBICITY_WEIGHT)}")
                print(f"\t- Hydrophobicity   : {get_hydrophobicity(individual.candidate) * HYDROPHOBICITY_WEIGHT}")
                print(f"\t- Birthdate        : {datetime.fromtimestamp(individual.birthdate).strftime("%d/%m/%Y - %H:%M:%S")}\n")
                if i < len(final_pop) - 1:
                    print("-" * 40)
            print("=" * 40)
            
            best_individual: Individual = min(final_pop, key = lambda x: x.fitness - (get_hydrophobicity(x.candidate) * HYDROPHOBICITY_WEIGHT)) # choose best fitness (minimize energy)
            
            print("\n --- Final Results ---")

            raw_cand = str(best_individual.candidate)
            clean_seq = raw_cand.split(" : ")[0].strip()

            # Calculating correct values ​​using the clean sequence
            hydro_val = get_hydrophobicity(clean_seq)
            binding_en = best_individual.fitness - (hydro_val * HYDROPHOBICITY_WEIGHT)
            weighted_hydro = hydro_val * HYDROPHOBICITY_WEIGHT
            birth_date_str = datetime.fromtimestamp(best_individual.birthdate).strftime("%d/%m/%Y - %H:%M:%S")

            print(f"Sequence       : {clean_seq}")
            print(f"Fitness        : {best_individual.fitness}")
            print(f"Binding Energy : {binding_en}")
            print(f"Hydrophobicity : {weighted_hydro}")
            print(f"Birthdate      : {birth_date_str}")

            with (OUTPUT / f"result_{JOB_ID}.txt").open("w") as f:
                f.write(f"ID             : {JOB_ID}\n")
                f.write(f"Receptor Name  : {RECEPTOR_NAME}\n")
                f.write(f"Receptor Path  : {RECEPTOR_FILE}\n")
                f.write(f"Sequence       : {clean_seq}\n")
                f.write(f"Fitness        : {best_individual.fitness}\n")
                f.write(f"Binding Energy : {binding_en}\n")
                f.write(f"Hydrophobicity : {weighted_hydro}\n")
                f.write(f"Birthdate      : {birth_date_str}\n")
                f.write("")
                f.write(print_arguments(options, JOB_ID, False))
                f.write("\n")
            
            # Rebuilds the path where Vina saved the temporary file
            unique_id = f"{clean_seq}_{JOB_ID}"       
            base_name     = job_temp_dir / f"p_{unique_id}" / unique_id
            output_dir    = OUTPUT

            try:
                pdb_src   = f"{base_name}.pdb"
                pdbqt_src = f"{base_name}.pdbqt"

                if Path(pdb_src).exists():
                    shutil.copy2(pdb_src, output_dir)
                    print_verbose(f"Copied final PDB: {pdb_src} -> {output_dir}")
                else:
                    print_warning(f"PDB file not found: {pdb_src}")

                if Path(pdbqt_src).exists():
                    shutil.copy2(pdbqt_src, output_dir)
                    
            except Exception as e:
                print_error(f"Error copying final files: {e}")

        finally:
            # To prevent them from being left hanging in the event of a genetic algorithm crash
            statistics_file .close()
            individuals_file.close()

            plot_observer_statistics(observer_file_path = statistics_filepath, plot_folder_directory = OUTPUT)
            plot_energy_vs_hydrophobicity(individuals_file = individuals_filepath, plot_folder_directory = OUTPUT, hydrophobicity_weight = HYDROPHOBICITY_WEIGHT)
            plot_observer_statistics_2(observer_file_path = statistics_filepath, plot_folder_directory = OUTPUT)
            plot_energy_vs_hydrophobicity_2(individuals_file = individuals_filepath, plot_folder_directory = OUTPUT, hydrophobicity_weight = HYDROPHOBICITY_WEIGHT)

            # Removes the entire job's temporary folder at the end
            if job_temp_dir.exists() and not NO_DELETE:
                print_verbose(f"Cleaning up the temporary folder: {job_temp_dir} ...", to_print = VERBOSE, lock = print_lock)
                shutil.rmtree(job_temp_dir)
                print_verbose("Done", to_print = VERBOSE, lock = print_lock)
        # Close the files opened for statistics and individuals

if __name__ == '__main__':
    run_peptide_ga()
