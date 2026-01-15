######################### IMPORTS #########################

import time

from argparse   import Namespace
from contextlib import nullcontext
from sys        import stderr, exit
from threading  import Lock
from typing     import Any

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
        # Costruiamo gli argomenti: [Colore, Prefisso (se c'è), Valori..., Reset]
        output = []
        if color : output.append(color)
        if prefix: output.append(prefix)
        
        # Uniamo tutto e stampiamo
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
    
    # Se il tempo attuale ha superato la deadline fissata per Vina,
    # significa che l'ultima generazione è stata 'tagliata' o siamo al limite.
    # Fermiamo l'evoluzione per salvare i dati.
    if time.time() >= global_deadline:
        print_warning(f"[TIME LIMIT] Raggiunta la deadline globale. Arresto evoluzione{f" alla generazione {generation_num.value} " if generation_num is not None else " "}per salvataggio.", lock = print_lock)
        return True
    return False

def generation_tracker_observer(population, num_generations, num_evaluations, args):
    """
    Aggiorna il contatore di generazione condiviso tra i processi.
    """
    generation_num  = args.get('generation_num', None)
    print_lock      = args.get('print_lock', None)
    
    if generation_num is not None:
        # num_generations è fornito automaticamente da inspyred (0, 1, 2...)
        generation_num.value = num_generations
        if print_lock:
            with print_lock:
                print(f"Valutata generazione {"iniziale" if generation_num.value == 0 else generation_num.value}.")

def print_arguments(arguments: Namespace, name = "", to_print: bool = True) -> str | None:
    output = list()
    output.append("\n" + "="*40)
    output.append(f"      CONFIGURAZIONE JOB: {name}")
    output.append("="*40)
        
    for key, value in vars(arguments).items():
        output.append(f"{key.ljust(22)}: {value}")
    output.append("="*40 + "\n")

    if to_print:
        print("\n".join(output))
    else:
        return "\n".join(output)
