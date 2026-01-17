######################### IMPORTS #########################

import io
import json
import pickle
import random
import shutil
import tarfile
import time

from argparse                 import Namespace
from contextlib               import nullcontext
from inspyred.ec.observers    import file_observer
from multiprocessing.managers import DictProxy
from pathlib                  import Path
from sys                      import stderr, exit
from threading                import Lock
from typing                   import Any


from constants import SAFETY_MARGIN

################# VARIABLES AND CONSTANTS #################

COL_BLACK     : str = "\x1b[{}m".format("0;30")
COL_RED       : str = "\x1b[{}m".format("0;31")
COL_GREEN     : str = "\x1b[{}m".format("0;32")
COL_YELLOW    : str = "\x1b[{}m".format("0;33")
COL_BLUE      : str = "\x1b[{}m".format("0;34")
COL_PURPLE    : str = "\x1b[{}m".format("0;35")
COL_LIGHT_BLUE: str = "\x1b[{}m".format("0;36")
COL_WHITE     : str = "\x1b[{}m".format("0;37")
COL_RESET     : str = "\x1b[0m"

INFO          : str = "<INFO> ----->"
ERROR         : str = "<ERROR> ---->"
VERB          : str = "<VERBOSE> -->"
WARNING       : str = "<WARNING> -->"

######################## FUNCTIONS ########################

def color_text(color_code: str, *text: str, sep: str | None = " ", end: str | None = "",) -> str:
    """
    Format a string with ANSI color codes.

    This utility wraps the provided text between the specified color code 
    and the global reset code (COL_RESET).

    Parameters
    ----------
    color_code : str
        The ANSI escape sequence for the desired color (e.g., COL_RED).
    *text : str
        One or more strings to be colored.
    sep : str, optional
        The separator used to join multiple text arguments. Defaults to " ".
    end : str, optional
        An optional string to append at the very end (after the reset code). 
        Defaults to an empty string.

    Returns
    -------
    str
        The formatted string with color codes applied.
        
    Examples
    --------
    >>> color_text(COL_RED, "Critical", "Error")
    '\\033[31mCritical Error\\033[0m'
    """
    return f"{color_code}{sep.join(text)}{COL_RESET}{end}"

def _base_print(*values: object, color: str | None = None, prefix: str = "", to_print: bool = True, sep: str = " ", end: str = "\n", file: Any = None, flush: bool = False, lock: Lock = None) -> None:
    """
    Internal helper to handle thread-safe printing with optional colors and prefixes.
    """
    if not to_print:
        return

    context = lock if lock else nullcontext()
    
    with context:
        # Let's build the arguments: [Color, Prefix (if any), Values..., Reset]
        output = []
        if color : output.append(color)
        if prefix: output.append(prefix)
        
        # Let's put everything together and print
        print(
            *output, *values, COL_RESET if color else "",
            sep   = sep,
            end   = end,
            file  = file,
            flush = flush
        )

def print_info(*values: object, **kwargs) -> None:
    """
    Print an informational message in light blue.

    Parameters
    ----------
    *values : object
        One or more values or expressions to be printed.
    color : str, optional
        ANSI color escape sequence. Defaults to COL_LIGHT_BLUE.
    prefix : str, optional
        Prefix to display before the message. Defaults to the global INFO constant.
        If set to an empty string (""), no prefix is shown.
    to_print : bool, optional
        If False, the message will not be printed. Defaults to True.
    sep : str, optional
        Separator between objects. Defaults to " ".
    end : str, optional
        String appended after the last value. Defaults to "\\n".
    file : file-like object, optional
        Output stream. Defaults to sys.stdout.
    flush : bool, optional
        Whether to forcibly flush the stream. Defaults to False.
    lock : multiprocessing.Lock, optional
        Lock object to ensure process-safe output.
    """
    kwargs.setdefault('color', COL_LIGHT_BLUE)
    kwargs.setdefault('prefix', INFO)
    _base_print(*values, **kwargs)
    

def print_error(*values: object, code: int = 0, **kwargs) -> None:
    """
    Print an error message in red and optionally exit the program.

    Parameters
    ----------
    *values : object
        One or more values or expressions to be printed.
    code : int, optional
        Exit code for the program. If non-zero, the program terminates 
        with this status. Defaults to 0 (no exit).
    color : str, optional
        ANSI color escape sequence. Defaults to COL_RED.
    prefix : str, optional
        Prefix to display before the message. Defaults to the global ERROR constant.
        If set to an empty string (""), no prefix is shown.
    to_print : bool, optional
        If False, the message will not be printed. Defaults to True.
    sep : str, optional
        Separator between objects. Defaults to " ".
    end : str, optional
        String appended after the last value. Defaults to "\\n".
    file : file-like object, optional
        Output stream. Defaults to sys.stderr.
    flush : bool, optional
        Whether to forcibly flush the stream. Defaults to False.
    lock : multiprocessing.Lock, optional
        Lock object to ensure process-safe output.
    """
    kwargs.setdefault('file', stderr)
    kwargs.setdefault('color', COL_RED)
    kwargs.setdefault('prefix', ERROR)

    _base_print(*values, **kwargs)

    if code != 0:
        exit(code)

def print_verbose(*values: object, to_print: bool = True, **kwargs) -> None:
    """
    Print a warning message, typically in purple.

    Parameters
    ----------
    *values : object
        One or more values or expressions to be printed.
    color : str, optional
        ANSI color escape sequence. Defaults to COL_YELLOW.
    prefix : str, optional
        Prefix to display before the message. Defaults to the global VERBOSE constant.
        If set to an empty string (""), no prefix is shown.
    to_print : bool, optional
        If False, the message will not be printed. Defaults to True.
    sep : str, optional
        Separator between objects. Defaults to " ".
    end : str, optional
        String appended after the last value. Defaults to "\\n".
    file : file-like object, optional
        Output stream. Defaults to sys.stdout.
    flush : bool, optional
        Whether to forcibly flush the stream. Defaults to False.
    lock : multiprocessing.Lock, optional
        Lock object to ensure process-safe output.
    """
    kwargs.setdefault('color', COL_YELLOW)
    kwargs.setdefault('prefix', VERB)

    _base_print(*values, **kwargs)

def print_warning(*values: object, **kwargs) -> None:
    """
    Print a warning message, typically in purple.

    Parameters
    ----------
    *values : object
        One or more values or expressions to be printed.
    color : str, optional
        ANSI color escape sequence. Defaults to COL_PURPLE.
    prefix : str, optional
        Prefix to display before the message. Defaults to the global WARNING constant.
        If set to an empty string (""), no prefix is shown.
    to_print : bool, optional
        If False, the message will not be printed. Defaults to True.
    sep : str, optional
        Separator between objects. Defaults to " ".
    end : str, optional
        String appended after the last value. Defaults to "\\n".
    file : file-like object, optional
        Output stream. Defaults to sys.stdout.
    flush : bool, optional
        Whether to forcibly flush the stream. Defaults to False.
    lock : multiprocessing.Lock, optional
        Lock object to ensure process-safe output.
    """
    kwargs.setdefault('color', COL_PURPLE)
    kwargs.setdefault('prefix', WARNING)
    _base_print(*values, **kwargs)

def time_based_termination(population, num_generations, num_evaluations, args):
    """
    Ferma le generazioni se abbiamo superato la deadline globale (meno il buffer).
    """
    global_deadline = args.get('global_deadline', float('inf'))
    generation_num  = args.get('generation_num', None)
    print_lock      = args.get('print_lock', None)
    
    # If the current time has exceeded the deadline set for Vina,
    # it means the last generation has been 'cut' or we're at the limit.
    # Let's stop evolution to save the data.
    if time.time() >= global_deadline - SAFETY_MARGIN:
        print_warning(f"[TIME LIMIT] Global deadline reached. Evolution halted{f" at generation {generation_num.value} " if generation_num is not None else " "}to save.", lock = print_lock)
        return True
    return False

def generation_tracker_observer(population, num_generations, num_evaluations, args):
    """
    Aggiorna il contatore di generazione condiviso tra i processi.
    """
    generation_num  = args.get('generation_num', None)
    print_lock      = args.get('print_lock', None)
    offset          = args.get('generation_offset', 0)
    
    if generation_num is not None:
        # num_generations is automatically provided by inspyred (0, 1, 2...)
        generation_num.value = num_generations + offset
        if print_lock:
            with print_lock:
                print(f"Rated generation {generation_num.value}.")

def print_arguments(arguments: Namespace, name = "", to_print: bool = True) -> str | None:
    output = list()
    output.append("\n" + "="*40)
    output.append(f"       JOB CONFIGURATION: {name}")
    output.append("="*40)
        
    for key, value in vars(arguments).items():
        output.append(f"{key.ljust(22)}: {value}")
    output.append("="*40 + "\n")

    if to_print:
        print("\n".join(output))
    else:
        return "\n".join(output)

def checkpoint_observer(population, num_generations, num_evaluations, args):
    """
    Salva lo stato corrente dell'evoluzione su file (pickle) usando pathlib.
    """
    global_deadline = args.get('global_deadline', float('inf'))
    if time.time() >= global_deadline - SAFETY_MARGIN:  # 10 second margin to avoid incomplete saves
        print_warning(f"[CHECKPOINT] Global deadline reached, not subject to final checkpoint.", lock = args.get('print_lock'))
        return
    
    # Retrieves the path. If it doesn't exist, use a default.
    # Assume args['checkpoint_file'] is already a string or a Path.
    job_id           = args.get('job_id', None)
    results_dir      = Path(args.get("results_dir", "."))
    filename         = results_dir / f"checkpoint{f'_{job_id}' if job_id else ''}.pkl"
    offset           = args.get('generation_offset', 0)
    cache_proxy      = args.get('multiprocessing_cache', None)
    current_temp_dir = Path(args.get('temp_dir', '.'))
    rand             = args.get('rand', random.Random())

    # Force FLUSH open log files
    # This ensures that the last generated line is actually written to disk before reading it.
    if 'statistics_file' in args and hasattr(args['statistics_file'], 'flush'):
        args['statistics_file'].flush()
    if 'individuals_file' in args and hasattr(args['individuals_file'], 'flush'):
        args['individuals_file'].flush()
    
    stats_path = results_dir / f"ga_observer_{job_id}.csv"
    inds_path  = results_dir / f"ga_individuals_{job_id}.csv"

    csv_stats_content = stats_path.read_text(encoding='utf-8') if stats_path.exists() else ""
    csv_inds_content  = inds_path.read_text(encoding='utf-8')  if inds_path.exists()  else ""
    
    plain_cache = {}
    if cache_proxy is not None:
        plain_cache = dict(cache_proxy)
    
    # Creating the TAR archive in memory (RAM)
    # We will only save the cached sequence files to save space.
    tar_buffer        = io.BytesIO()
    sequences_to_save = list(plain_cache.keys())
    files_added_count = 0
    if current_temp_dir.exists():
        with tarfile.open(fileobj = tar_buffer, mode = 'w:bz2') as tar:
            for seq in sequences_to_save:
                # Expected structure: temp_dir / p_{seq} / ...
                # NOTE: We assume the folder is named p_{seq} or similar.
                # If you use a UUID, the logic below must know how to map seq -> folder.
                # In your case, if you use seq as the folder ID:
                seq_folder_name = f"p_{seq}_{job_id}" if job_id else f"p_{seq}"
                seq_folder_path = current_temp_dir / seq_folder_name
                
                if seq_folder_path.exists():
                    # Adds the entire folder to the tar, keeping the name "p_{seq}/..."
                    # arcname means "save it in the tar with this relative name."
                    tar.add(seq_folder_path, arcname = seq_folder_name)
                    files_added_count += 1
    
    compressed_files_blob = tar_buffer.getvalue()

    safe_args = dict()
    safe_args['job_temp_dir_snapshot'] = str(current_temp_dir)          # Path to the current temporary folder so it can be found on the resume
    safe_args['job_id_snapshot']       = job_id                         # Job ID to rename the files in the temporary folder on the resume
    safe_args['results_dir_snapshot']  = str(results_dir)               # Path to the results folder so it can be found on the resume

    data = {
        "generation"          : num_generations + offset,
        "population"          : population,
        "args"                : safe_args,
        "evaluation_cache"    : plain_cache,
        "temp_files_tar"      : compressed_files_blob,  # Compressed PDB files
        "csv_statistics_dump" : csv_stats_content,      # Statistics CSV text
        "csv_individuals_dump": csv_inds_content,        # Individuals CSV text
        "rand_state"          : rand.getstate(),
    }
    
    temp_filename = filename.with_suffix('.tmp')
    try:
        with temp_filename.open("wb") as f:
            pickle.dump(data, f)
        temp_filename.replace(filename)
        print_info(f"[CHECKPOINT] Saved Gen {num_generations + offset}. Temp Archive: {files_added_count} folders ({len(compressed_files_blob) / (1024 * 1024):.2f} MB).", lock = args.get('print_lock'))
    except Exception as e:
        print_error(f"Checkpoint saving error: {e}", lock = args.get('print_lock'))

def restore_context(checkpoint_path: Path, new_temp_dir: Path, new_job_id: str, new_results_dir: Path) -> tuple[int, list, dict]:
    """
    Carica il pickle, sposta i file temp, unisce i CSV. Restituisce: (generations_done, population_seeds)
    """
    ckpt              = Path(checkpoint_path)

    print_verbose(f"[RESUME] Attempting to restore from checkpoint: {ckpt} ...")
    if not ckpt.exists():
        print_warning(f"[RESUME] Checkpoint file {ckpt} not found. Starting from scratch.")
        return 0, [], {}
    print_verbose(f"[RESUME] Checkpoint {ckpt} found, starting recovery ...")
    
    
    data = None
    try:
        with ckpt.open("rb") as f:
            print_verbose(f"[RESUME] Reading file ...")
            data = pickle.load(f)
            
    except FileNotFoundError:
        print_error(f"[RESUME] CRITICAL ERROR: The file existed a moment ago but now it won't open (Probable .pkl file incompatibility for a different version). Path: {ckpt}")
        print_error(f"[RESUME] Current CWD: {Path('.').resolve()}")
        return 0, [], {}
    except Exception as e:
        print_error(f"[RESUME] Error loading pickle: {e}")
        return 0, [], {}

    if data is None:
        print_warning(f"[RESUME] Empty or invalid checkpoint. Start from scratch.")
        return 0, [], {}
    
    print_verbose(f"[RESUME] Checkpoint loaded successfully.")
    
    gen_done          = data['generation']
    population        = data['population']
    saved_args        = data['args']
    saved_cache       = data.get('evaluation_cache', {})
    temp_archive_blob = data.get('temp_files_tar', None)
    csv_stats_dump    = data.get('csv_statistics_dump', '')
    csv_inds_dump     = data.get('csv_individuals_dump', '')
    rand_state        = data.get('rand_state', None)
    
    old_job_id     = saved_args.get('job_id_snapshot', 'unknown')
    old_temp_dir   = Path(saved_args.get('job_temp_dir_snapshot', ''))
    old_res_dir    = Path(saved_args.get('results_dir_snapshot', '.'))
    new_stats_path = new_results_dir / f"ga_observer_{new_job_id}.csv"
    new_inds_path  = new_results_dir / f"ga_individuals_{new_job_id}.csv"

    print_info(f"[RESUME] Restore state to generation {gen_done}.")
    
    if temp_archive_blob:
        print_info(f"[RESUME] Extract and update temporary file IDs ...")
        try:
            buffer = io.BytesIO(temp_archive_blob)
            with tarfile.open(fileobj = buffer, mode = 'r:bz2') as tar:
                members = tar.getmembers()
                count   = 0
                for member in members:
                    # Example original name: "p_SEQ_OLDID/SEQ_OLDID.pdb"
                    # Replace OLDID with NEWID everywhere in the path.
                    if old_job_id and old_job_id != 'unknown' and old_job_id is not None:
                        member.name = member.name.replace(old_job_id, new_job_id)
                    
                    # Extracting a single file/folder with a new name
                    tar.extract(member, path = new_temp_dir)
                    count += 1
                
            print_info(f"[RESUME] Extraction complete. Files/folders extracted.: {count}")
            
        except Exception as e:
            print_error(f"[RESUME] Archive unzipping error: {e}")
    else:
        # Old-fashioned fallback (if the checkpoint does not have internal storage)
        old_temp_dir_raw = saved_args.get('job_temp_dir_snapshot', '')
        if old_temp_dir_raw:
            old_temp_dir = Path(old_temp_dir_raw).resolve()
            if old_temp_dir.exists():
                print_warning(f"[RESUME] Internal archive is missing. I'm trying to migrate from a physical folder {old_temp_dir} ...")
                for item in old_temp_dir.glob("**/*"):
                    if item.is_file():
                        try:
                            rel_path         = item.relative_to(old_temp_dir)
                            new_rel_path_str = str(rel_path).replace(old_job_id, new_job_id)
                            dest_path        = new_temp_dir / new_rel_path_str
                            dest_path.parent.mkdir(parents = True, exist_ok = True)
                            shutil.copy2(item, dest_path)
                        except: pass
            else:
                print_warning(f"[RESUME] No recoverable temporary files.")

    # CSV Merge (Statistics and Individuals)
    # Look for old files in the saved results folder.
    
    # Internal function to copy
    def append_csv(src, dst):
        if src.exists():
            shutil.copy2(src, dst) # Brutal Copy: Inspyred Will Append After

    def prune_ge_cutoff(filepath: Path, cutoff_gen: int):
        if not filepath.exists(): return
        try:
            with filepath.open('r') as f: 
                lines = f.readlines()
            
            if len(lines) <= 1: return # Header only or empty

            clean_lines = []
            clean_lines.append(lines[0]) 

            removed_count = 0
            # Analyze the data (from row 1 onwards)
            for line in lines[1:]:
                parts = line.split(',')
                # If the line is valid and the generation is < cutoff, we keep it
                if parts[0].strip().isdigit():
                    gen_in_row = int(parts[0])
                    if gen_in_row < cutoff_gen:
                        clean_lines.append(line)
                    else:
                        removed_count += 1
                else:
                    # If we can't read the generation (e.g. corrupted line), we discard it for safety
                    pass

            # Let's rewrite the clean file
            with filepath.open('w') as f:
                f.writelines(clean_lines)
            
            if removed_count > 0:
                print_verbose(f"[RESUME] Pruning {filepath.name}: removed {removed_count} rows (Gen >= {cutoff_gen}).")
                
        except Exception as e:
            print_warning(f"[RESUME] Pruning error {filepath.name}: {e}")
    
    # If the dump exists, we write it
    print_info(f"[RESUME] Merging CSV files of statistics and individuals ...")
    
    if csv_stats_dump:
        try:
            new_stats_path.write_text(csv_stats_dump, encoding = 'utf-8')
            print_info(f"[RESUME] CSV Statistics restored ({len(csv_stats_dump)} bytes)")
        except Exception as e:
            print_warning(f"[RESUME] Error writing CSV stats: {e}")
    else:
        old_stats = old_res_dir     / f"ga_observer_{old_job_id}.csv"
        new_stats = new_results_dir / f"ga_observer_{new_job_id}.csv"
        append_csv(old_stats, new_stats)

    prune_ge_cutoff(new_stats_path, gen_done)

    if csv_inds_dump:
        try:
            new_inds_path.write_text(csv_inds_dump, encoding = 'utf-8')
            print_verbose(f"[RESUME] CSV Individuals restored ({len(csv_inds_dump)} bytes)")
        except Exception as e:
            print_warning(f"[RESUME] CSV writing error inds: {e}")
    else:
        old_inds  = old_res_dir     / f"ga_individuals_{old_job_id}.csv"
        new_inds  = new_results_dir / f"ga_individuals_{new_job_id}.csv"
        append_csv(old_inds, new_inds)

    prune_ge_cutoff(new_inds_path, gen_done)
    
    print_info(f"[RESUME] CSV file merge completed.")
    
    return gen_done, population, saved_cache, rand_state

def custom_file_observer(population, num_generations, num_evaluations, args: dict):
    """
    Crea un file observer personalizzato per la correzione del numero di generazioni in caso di riavvio
    """
    offset      = args.get('generation_offset', 0)
    # Calls the standard observer to save statistics and individuals
    file_observer(population, num_generations + offset, num_evaluations, args)    

# utils.py

def inspect_checkpoint(checkpoint_path: str | Path, show_all: bool = False) -> None:
    """
    Legge e stampa la struttura interna di un file checkpoint (.pkl).
    Allinea dinamicamente le colonne per una lettura pulita.
    
    Parameters
    ----------
    checkpoint_path : str | Path
        Percorso al file .pkl.
    show_all : bool
        Se True, mostra TUTTI gli elementi. Se False, mostra un'anteprima.
    """
    import pickle
    from pathlib import Path

    ckpt = Path(checkpoint_path)

    if not ckpt.exists():
        print_error(f"[INSPECT] File not found: {ckpt}")
        return

    try:
        with ckpt.open("rb") as f:
            print_verbose(f"[INSPECT] Reading file in progress: {ckpt}")
            data = pickle.load(f)

        print(f"\n{'='*60}")
        print(f" CHECKPOINT CONTENT: {ckpt.name}")
        print(f"{'='*60}")
        
        for key, value in data.items():
            print(f"- {key:<25} : {type(value).__name__}")
            
            # List Management
            if isinstance(value, list):
                print(f"\tTotal length: {len(value)}")
                limit = len(value) if show_all else 3
                
                for i, item in enumerate(value[:limit]): 
                    print(f"\t  [{i}] {item}")
                
                if not show_all and len(value) > limit:
                    print(f"\t  ... (other {len(value)-limit} hidden elements)")
            
            # Dictionary Management (WITH DYNAMIC ALIGNMENT)
            elif isinstance(value, dict):
                print(f"\tTotal keys: {len(value)}")
                
                # Let's get the items to display
                items_list = list(value.items())
                limit = len(items_list) if show_all else 5
                display_items = items_list[:limit]
                
                if display_items:
                    # Calculate the length of the longest key in this group to align
                    max_key_len = max(len(str(k)) for k, _ in display_items)
                    
                    for k, v in display_items:
                        # Value formatting for cleanup (rounded floats)
                        if isinstance(v, float):
                            val_str = f"{v:.4f}"
                        elif isinstance(v, (str, int)):
                            val_str = str(v)
                        else:
                            val_str = type(v).__name__

                        # Truncate long strings if we are not in show_all
                        if not show_all and len(val_str) > 50:
                            val_str = val_str[:50] + "..."
                        
                        # PRINT ALIGNED: Use max_key_len for padding
                        print(f"\t  {str(k):<{max_key_len}} -> {val_str}")
                    
                    if not show_all and len(items_list) > limit:
                        print(f"\t  ... (other {len(items_list)-limit} hidden keys)")
            
            # Binary Data Management
            elif isinstance(value, bytes):
                size_mb = len(value) / (1024 * 1024)
                print(f"\t[BINARY DATA] Size: {len(value)} bytes ({size_mb:.2f} MB)")
                print(f"\t(Hidden content)")

            # Other
            else:
                val_str = str(value)
                if not show_all and len(val_str) > 200:
                    print(f"\tValue (truncated): {val_str[:200]} ...")
                else:
                    print(f"\tValue: {val_str}")
            
            print("-" * 60)

    except Exception as e:
        print_error(f"[INSPECT] Error during inspection: {e}")

