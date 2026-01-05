######################### IMPORTS #########################

import os

# from _typeshed import SupportsWrite
from pathlib import Path, PurePath
from sys import stderr, exit
from typing import Any, Union
from typing_extensions import Literal

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

######################## CLASSES ########################

class AutoResolvePath:
    def __init__(self, path: Union["AutoResolvePath", bytes, os.PathLike, Path, str]) -> None:
        if isinstance(path, str):
            if len(path) > 0 and path[0] == "~":
                self._path = Path(f"{Path.home()}/{path[1:]}")
            else:
                self._path = Path(path)
        elif isinstance(path, AutoResolvePath):
            self._path = Path(path._path)
        else:
            self._path = Path(path)
    
    def __getattr__(self, name: str) -> Any:
        return getattr(self._path, name)
    
    def __truediv__(self, other: Union[PurePath, "AutoResolvePath", str]) -> "AutoResolvePath":
        if isinstance(other, PurePath) or isinstance(other, AutoResolvePath):
            return AutoResolvePath(f"{self._path}/{other._path}")
        elif isinstance(other, str):
            return AutoResolvePath(f"{self._path}/{other}")
        else:
            raise TypeError(f"'other' object of type '{type(other)}' can't be concatenated to an 'AutoResolvePath' object")

    def __str__(self) -> str:
        return self._path.__str__()
    
    def __repr__(self) -> str:
        return self._path.__repr__()

######################## FUNCTIONS ########################

def color_text(color_code: str, *text: str, sep: str | None = " ", end: str | None = "",) -> str:
    return f"{color_code}{sep.join(text)}{COL_RESET}{end}"

def print_info(*values: object, to_print: bool = True, sep: str | None = " ", end: str | None = "\n", file = None, flush: Literal[False] = False, prefix: str = INFO) -> None:
    """
    Print info in `*values` to the stream `text` file, separed by `sep` and followed by `end`.

    Parameters
    ----------
    *values : `object`
        One or more values or expressions to print on file or stream.
    
    to_print : `bool`, default `True`, optional
        If `True`, print to file or stream.
    
    sep : `str` | `None`, default ` `, optional
        Separator between objects.

    end : `str` | `None`, default `\\n`, optional
        Character to add at the end.
    
    file : `SupportsWrite[str]` | `None`, default `None`, optional
        File or stream to write output to.
    
    flush : `Literal[False]`, default `False`, optional
        If `True`, flush the buffer immediately.
    """
    if to_print:
        if prefix == "":
            print(
                *(COL_LIGHT_BLUE, *values, COL_RESET),
                sep = sep,
                end = end,
                file = file,
                flush = flush
            ) 
        else:
            print(
                *(COL_LIGHT_BLUE, prefix, *values, COL_RESET),
                sep = sep,
                end = end,
                file = file,
                flush = flush
            )

def print_error(*values: object, code: int = 0, to_print: bool = True, sep: str | None = " ", end: str | None = "\n", file = stderr, flush: Literal[False] = False) -> None:
    """
    Print error in `*values` to the stream `text` file, separed by `sep` and followed by `end`. If `code` different of 0, exit with `code`.

    Parameters
    ----------
    *values : `object`
        One or more values or expressions to print on file or stream.
    
    code : `int`, default `-1`, optional
        Program execution exit code. If equal to `0` then the error won't cause the program to exit.
    
    to_print : `bool`, default `True`, optional
        If `True`, print to file or stream.
    
    sep : `str` | `None`, default ` `, optional
        Separator between objects.

    end : `str` | `None`, default `\\n`, optional
        Character to add at the end.
    
    file : `SupportsWrite[str]` | `None`, default `None`, optional
        File or stream to write output to.
    
    flush : `Literal[False]`, default `False`, optional
        If `True`, flush the buffer immediately.
    """
    if to_print:
        print(
            *(COL_RED, ERROR, *values, COL_RESET),
            sep = sep,
            end = end,
            file = file,
            flush = flush
        )

    if code != 0:
        exit(code)

def print_verbose(*values: object, to_print: bool = True, sep: str | None = " ", end: str | None = "\n", file = None, flush: Literal[False] = False) -> None:
    """
    Print verbose in `*values` to the stream `text` file, separed by `sep` and followed by `end`.

    Parameters
    ----------
    *values : `object`
        One or more values or expressions to print on file or stream.
    
    to_print : `bool`, default `True`, optional
        If `True`, print to file or stream.
    
    sep : `str` | `None`, default ` `, optional
        Separator between objects.

    end : `str` | `None`, default `\\n`, optional
        Character to add at the end.
    
    file : `SupportsWrite[str]` | `None`, default `None`, optional
        File or stream to write output to.
    
    flush : `Literal[False]`, default `False`, optional
        If `True`, flush the buffer immediately.
    """
    if to_print:
        print(
            *(VERB, *values),
            sep = sep,
            end = end,
            file = file,
            flush = flush
        )

def print_warning(*values: object, to_print: bool = True, sep: str | None = " ", end: str | None = "\n", file = None, flush: Literal[False] = False) -> None:
    """
    Print warning in `*values` to the stream `text` file, separed by `sep` and followed by `end`.

    Parameters
    ----------
    *values : `object`
        One or more values or expressions to print on file or stream.
    
    to_print : `bool`, default `True`, optional
        If `True`, print to file or stream.
    
    sep : `str` | `None`, default ` `, optional
        Separator between objects.

    end : `str` | `None`, default `\\n`, optional
        Character to add at the end.
    
    file : `SupportsWrite[str]` | `None`, default `None`, optional
        File or stream to write output to.
    
    flush : `Literal[False]`, default `False`, optional
        If `True`, flush the buffer immediately.
    """
    if to_print:
        print(
            *(COL_PURPLE, WARNING, *values, COL_RESET),
            sep = sep,
            end = end,
            file = file,
            flush = flush
        )
