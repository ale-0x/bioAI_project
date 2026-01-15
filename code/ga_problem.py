# ga_problem.py
import numpy as np
import random
import subprocess
import time
import uuid

from multiprocessing.managers import DictProxy
from openbabel                import openbabel, pybel
from pathlib                  import Path
from rdkit                    import Chem
from rdkit.Chem               import AllChem
from subprocess               import TimeoutExpired
from typing                   import Any, Dict, List, Tuple


import constants as C

from peptide_operators import get_hydrophobicity
from utils      import print_error, print_verbose, print_warning

# --- Funzione di Generazione della Popolazione Iniziale ---

def peptide_generator(random: random.Random, args: Dict[str, Any]) -> str:
    """
    Generatore della sequenza peptidica iniziale per la popolazione.

    Crea una sequenza casuale di aminoacidi di lunghezza fissa 
    (definita da `PEPTIDE_LENGTH` o dagli argomenti).

    Parameters
    ----------
    random : `random.Random`
        L'istanza dell'oggetto casuale fornita dal framework `inspyred`.
    args : `Dict[str, Any]`
        Dizionario di argomenti opzionali. Può contenere la chiave 
        `'peptide_length'` (`int`) per specificare la lunghezza.

    Returns
    -------
    `str`
        Una sequenza peptidica casuale (stringa di codici a una lettera).

    Examples
    --------
    >>> generator = peptide_generator(random_instance, {'peptide_length': 5})
    >>> print(generator)
    'ALTSV' # Esempio di output casuale
    """
    length: int = args.get('peptide_length', C.PEPTIDE_LENGTH)
    return "".join(random.choice(C.AMINO_ACIDS) for _ in range(length))


# --- Funzioni Helper e pipeline di preparazione (RDKit + Meeko) per Struttura e Docking ---

def generate_pdb_rdkit(sequence: str, output_file: str | Path) -> None:
    """
    Genera una conformazione 3D iniziale per un peptide data la sua sequenza.

    Utilizza RDKit per costruire la catena aminoacidica, aggiungere idrogeni
    ed eseguire una minimizzazione energetica rapida (MMFF94).

    Parameters
    ----------
    sequence : `str`
        Sequenza peptidica (es. 'ACDEF').
    output_pdb : `str`
        Percorso del file PDB di output da creare.
    """
    output_file = Path(output_file)

    try:
        mol = Chem.MolFromSequence(sequence)
        if mol is None:
            raise ValueError(f"RDKit non riesce a fare il parse della sequenza: {sequence}")
        mol = Chem.AddHs(mol)
        
        # 1. Embedding (Generazione coordinate iniziali)
        # Usa ETKDGv3 che è più robusto per molecole grandi/cicliche
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        
        if result == -1:
            # Fallback se fallisce
            result = AllChem.EmbedMolecule(mol, useRandomCoords = True)
        
        if result == -1:
            raise ValueError(f"Impossibile generare conformero 3D per {sequence}")
            
        # 2. MINIMIZZAZIONE (Fondamentale!)
        # Rilassa la struttura per evitare atomi sovrapposti
        AllChem.MMFFOptimizeMolecule(mol)
        
        Chem.MolToPDBFile(mol, str(output_file))
    except Exception as e:
        print_error(f"Errore generazione RDKit per {sequence}: {e}")

def prepare_ligand_openbabel(pdb_path: str | Path, pdbqt_output_path: str | Path, center: Tuple[float, float, float]) -> bool:
    """
    Converte e prepara un file ligando (es. PDB) in formato PDBQT utilizzando OpenBabel.

    Questa funzione esegue passaggi critici per il docking:
    1. Legge il file molecolare di input.
    2. Aggiunge gli idrogeni polari (protonazione a pH 7.4).
    3. Calcola le cariche parziali (metodo Gasteiger).
    4. Genera l'albero delle torsioni (necessario per la flessibilità in Vina) e scrive il file .pdbqt.

    Parameters
    ----------
    input_file : `str`
        Il percorso al file del ligando di input (es. 'ligand.pdb').
    
    output_pdbqt_file : `str`
        Il percorso dove salvare il file PDBQT preparato (es. 'ligand.pdbqt').
    
    verbose : `bool`, default `False`, optional
        Se `True`, stampa informazioni dettagliate sul processo di conversione.

    Returns
    -------
    `bool`
        `True` se la preparazione e la scrittura del file sono avvenute con successo,
        `False` in caso di errori (es. file non trovato, errore di parsing).

    Notes
    -----
    Richiede che OpenBabel (bindings Python) sia installato correttamente nell'ambiente.
    L'output PDBQT includerà automaticamente i rami ROOT/BRANCH/TORSDOF gestiti da OpenBabel.
    """
    pdb_path            = Path(pdb_path)
    pdbqt_output_path   = Path(pdbqt_output_path)

    try:
        # 1. Leggi PDB (Generato da RDKit, quindi sicuro)
        mol = next(pybel.readfile("pdb", str(pdb_path)))
        
        # 2. Centratura nella tasca (Box Vina)
        atoms         = [atom.coords for atom in mol] 
        centroid      = np.mean(atoms, axis = 0)
        target_center = np.array([center[0], center[1], center[2]])
        move_v        = target_center - centroid
        mol.OBMol.Translate(openbabel.vector3(move_v[0], move_v[1], move_v[2]))
        
        # 3. Scrittura PDBQT (OpenBabel calcola le cariche Gasteiger automaticamente)
        mol.write("pdbqt", str(pdbqt_output_path), overwrite = True)
        
        return pdbqt_output_path.exists() and pdbqt_output_path.stat().st_size > 0
    except StopIteration:
        print_error(f"Errore OpenBabel: File PDB vuoto o invalido: {pdb_path}")
        return False
    except Exception as e:
        print_error(f"Conversione fallita per {pdb_path}: {e}")
        return False
    

def run_vina_real(vina_exe_path: str | Path, pdbqt_ligand: str | Path, receptor_file: str | Path, center: Tuple[float, float, float], box_size: Tuple[int, int, int], vina_output: str | Path, cpu: int = 1, exhaustiveness: int = C.EXHAUSTIVENESS, time_left: float = float("inf"), verbose: bool = False) -> float:
    """
    Esegue il docking molecolare lanciando il processo AutoDock Vina e restituisce l'energia di legame migliore.

    La funzione costruisce ed esegue il comando da terminale per Vina, specificando
    il recettore, il ligando, la search box e i parametri di precisione.
    Cattura l'output standard per estrarre il punteggio di affinità (energia libera)
    del primo modo di binding (il migliore).

    Parameters
    ----------
    vina_exe_path : `str`
        Percorso dell'eseguibile di AutoDock Vina.
    
    ligand_pdbqt_file : `str`
        Percorso al file del ligando preparato (.pdbqt).
    
    receptor_pdbqt_file : `str`
        Percorso al file del recettore preparato (.pdbqt).
    
    output_file : `str`
        Percorso dove salvare il file di output contenente le pose di docking (.pdbqt).
    
    center : `Tuple[float, float, float]`
        Coordinate (x, y, z) del centro della box di ricerca (in Ångstrom).
    
    box_size : `Tuple[float, float, float]`
        Dimensioni (x, y, z) della box di ricerca (in Ångstrom).
    
    exhaustiveness : `int`, default `8`, optional
        Parametro di esaustività della ricerca globale di Vina (valori più alti = ricerca più accurata ma lenta).
    
    vina_output : `str`
        Output file di AutoDock Vina.
    
    cpu : `int`, default `1`, optional
        Numero di CPU/thread da dedicare a questa singola esecuzione di Vina.
    
    verbose : `bool`, default `False`, optional
        Se `True`, stampa il comando eseguito e l'output grezzo di Vina in caso di errore.

    Returns
    -------
    `float`
        L'energia di legame del miglior modo (in kcal/mol). 
        Restituisce `0.0` (o un valore positivo alto di penalità) se il docking fallisce o non vengono trovati modi.

    Raises
    ------
    RuntimeError
        Se l'eseguibile di Vina non viene trovato o restituisce un codice di errore.
    """
    vina_exe_path = Path(vina_exe_path)
    pdbqt_ligand  = Path(pdbqt_ligand)
    receptor_file = Path(receptor_file)
    vina_output   = Path(vina_output)

    if not receptor_file.exists():
        print_error(f"ERRORE CRITICO: Il file recettore non esiste: {receptor_file}")
        return float('inf')
            
    cmd = [
        vina_exe_path,
        "--receptor"      , str(receptor_file),
        "--ligand"        , str(pdbqt_ligand),
        "--center_x"      , str(center[0]),
        "--center_y"      , str(center[1]),
        "--center_z"      , str(center[2]),
        "--size_x"        , str(box_size[0]),
        "--size_y"        , str(box_size[1]),
        "--size_z"        , str(box_size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--cpu"           , str(cpu),
        "--out"           , str(vina_output),
        # "--score_only" 
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output = True,
            text           = True,
            check          = True,
            timeout        = time_left
        )
        
        for line in result.stdout.splitlines():
            if line not in {
                '#################################################################',
                '# If you used AutoDock Vina in your work, please cite:          #',
                '#                                                               #',
                '# O. Trott, A. J. Olson,                                        #',
                '# AutoDock Vina: improving the speed and accuracy of docking    #',
                '# with a new scoring function, efficient optimization and       #',
                '# multithreading, Journal of Computational Chemistry 31 (2010)  #',
                '# 455-461                                                       #',
                '#                                                               #',
                '# DOI 10.1002/jcc.21334                                         #',
                '#                                                               #',
                '# Please see http://vina.scripps.edu for more information.      #',
                '#################################################################',
                'Reading input ... done.',
                'Setting up the scoring function ... done.',
                'Analyzing the binding site ... done.',
                'Performing search ...               ',
                '',
                '0%   10   20   30   40   50   60   70   80   90   100%',
                '|----|----|----|----|----|----|----|----|----|----|',
                '***************************************************',
                'done.',
                'Refining results ... done.',
            } and "Using random seed:" not in line and "Performing search ..." not in line:
                # print_verbose(line, to_print = verbose)
                # Cerca pattern tipo: "   1        -8.5      0.000      0.000"
                if line.strip().startswith("1"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            return float(parts[1])
                        except ValueError:
                            continue
                
                # Cerca pattern tipo: "Affinity: -8.5"
                if "Affinity:" in line:
                    parts = line.split()
                    if len(parts) >= 2: 
                        # se il parsing è riuscito, ritorna il valore richiesto (energy) 
                        # e inoltre salva l'output per vedere il docking
                        # open(vina_output, "w").write(result.stdout)
                        # close(vina_output)
                        return float(parts[1])
        return float('inf') # Se il parsing fallisce
    except subprocess.CalledProcessError as e:
        # Vina ha restituito un codice di errore
        print_error(f"Vina Error Output: {e.stderr}")
        return float('inf')
    except TimeoutExpired:
        print_warning(f"[TIME LIMIT] Vina interrotto per sequenza '{Path(pdbqt_ligand).stem.split('_')[0]}': Raggiunto il limite globale del Job.")
        return float('inf')
    except Exception as e:
        print_error(f"Errore esecuzione Vina (subprocess): {e}")
        return float('inf')


# --- Funzione di Valutazione Principale ---

def evaluate_peptide_binding(candidates: List[str], args: Dict[str, Any]) -> List[float]:
    """
    Valuta l'affinità di legame (fitness) per una lista di sequenze peptidiche.

    Per ogni candidato, genera la struttura 3D, la prepara per Vina ed esegue
    il calcolo del docking. Restituisce il valore di energia libera stimato.

    Parameters
    ----------
    candidates : `List[str]`
        Lista di sequenze peptidiche (stringhe di aminoacidi).
    args : `Dict[str, Any]`
        Dizionario di parametri che deve contenere 'receptor_file' e i dati della grid box.

    Returns
    -------
    `List[float]`
        Lista dei valori di fitness (kcal/mol). Valori più bassi indicano legami migliori.
    """
    # Lettura parametri con fallback
    job_id                : str              = args.get('job_id', str(uuid.uuid4()))
    verbose               : bool             = args.get('verbose', False)
    vina_exe_path         : Path             = Path(args.get('vina_exe_path', C.VINA_EXE_PATH))
    hydrophobicity_weight : float            = args.get('hydrophobicity_weight', C.HYDROPHOBICITY_WEIGHT)
    receptor_file         : Path             = Path(args.get('receptor_file', C.RECEPTOR_FILE))
    multiprocessing_cache : DictProxy | None = args.get('multiprocessing_cache', None)
    global_deadline       : float            = args.get('global_deadline', float('inf'))
    print_lock            : Any              = args.get('print_lock', None)
    
    cx = args.get('center_x', C.CENTER_X)
    cy = args.get('center_y', C.CENTER_Y)
    cz = args.get('center_z', C.CENTER_Z)
    
    sx = args.get('size_x', C.SIZE_X)
    sy = args.get('size_y', C.SIZE_Y)
    sz = args.get('size_z', C.SIZE_Z)
        
    # Cartella temporanea del JOB (deve esistere, creata da main.py)
    # Se per qualche motivo non c'è, usa /tmp locale con fallback
    base_temp = Path(args.get('temp_dir', f"tmp/bioai_ga_{job_id}"))
    if not base_temp.exists():
        print_verbose(f"[evaluate_peptide_binding] Creazione directory temporanea '{base_temp}' ...", to_print = verbose, lock = print_lock)
        base_temp.mkdir(parents = True, exist_ok = True)
        if base_temp.exists():
            print_verbose(f"[evaluate_peptide_binding] Done.", to_print = verbose, lock = print_lock)
        else:
            print_error(f"[evaluate_peptide_binding] Directory temporanea '{base_temp}' non creata!", code = -1, lock = print_lock)

    fitnesses = []
    for candidate in candidates:
        # Estrai la sequenza se è un oggetto Individual
        seq = candidate.candidate if hasattr(candidate, 'candidate') else candidate

        if multiprocessing_cache is not None and seq in multiprocessing_cache:
            print_verbose(f"[evaluate_peptide_binding] Cache hit per '{seq}'. Recupero fitness da cache.", to_print = verbose, lock = print_lock)
            fitnesses.append(multiprocessing_cache[seq])
            continue

        hydrophobicity = get_hydrophobicity(seq)
        
        # Crea cartella isolata per questo singolo peptide
        peptide_dir   = base_temp / f"p_{seq}_{job_id}"

        print_verbose(f"[evaluate_peptide_binding] Creazione directory peptide '{seq}' in '{peptide_dir}' ...", to_print = verbose, lock = print_lock)
        peptide_dir.mkdir(parents = True, exist_ok = True)
        if peptide_dir.exists():
            print_verbose(f"[evaluate_peptide_binding] Done.", to_print = verbose, lock = print_lock)
        else:
            print_error(f"[evaluate_peptide_binding] Directory peptide '{seq}' in '{peptide_dir}' non creata!", code = -1, lock = print_lock)

        base_name  = peptide_dir / f"{seq}_{job_id}"
        pdb_file   = Path(f"{base_name}.pdb")
        pdbqt_file = Path(f"{base_name}.pdbqt")
        out_file   = Path(f"{base_name}_out.pdbqt")
        
        try:
            time_left    = global_deadline - time.time()
            if time_left <= 60.0:
                # Solleviamo manualmente l'errore per andare nel blocco except
                raise TimeoutExpired(vina_exe_path, 0, f"Tempo globale in esaurimento ({time_left}s <= 60s) prima dell'avvio")
            
            # 1. Generazione Struttura
            print_verbose(f"[evaluate_peptide_binding] Generazione Struttura di '{seq}' in '{pdb_file}' ...", to_print = verbose, lock = print_lock)
            generate_pdb_rdkit(seq, pdb_file)
            print_verbose(f"[evaluate_peptide_binding] Generazione Struttura di '{seq}' in '{pdb_file}' --> Done.", to_print = verbose, lock = print_lock)
            
            # 2. Pybel -> RDKit (Centratura) -> Meeko e poi Docking
            print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' ...", to_print = verbose, lock = print_lock)
            if prepare_ligand_openbabel(pdb_file, pdbqt_file, (cx, cy, cz)):
                print_verbose(f"[evaluate_peptide_binding] Valutazione di Autodock Vina su '{pdbqt_file}' della sequenza '{seq}' ...", to_print = verbose, lock = print_lock)
                energy = run_vina_real(
                    vina_exe_path  = vina_exe_path,
                    pdbqt_ligand   = pdbqt_file,
                    receptor_file  = receptor_file,
                    center         = (cx, cy, cz),
                    box_size       = (sx, sy, sz),
                    vina_output    = out_file,
                    cpu            = 1,             # 1 CPU per processo figlio
                    exhaustiveness = args.get('exhaustiveness', C.EXHAUSTIVENESS),
                    time_left      = time_left,
                    verbose        = verbose
                )
                print_verbose(f"[evaluate_peptide_binding] Valutazione di Autodock Vina su '{pdbqt_file}' della sequenza '{seq}' --> Done", to_print = verbose, lock = print_lock)
                
                # usa la hydrophobicity average come penalità per restringere il campo di soluzioni possibili
                fitnesses.append(energy + (hydrophobicity * hydrophobicity_weight))
                if multiprocessing_cache is not None:
                    multiprocessing_cache[seq] = fitnesses[-1]
                print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' --> Done", to_print = verbose, lock = print_lock)
            else:
                fitnesses.append(float('inf'))               # Penalità per fallimento prep
                print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' --> Fallito", to_print = verbose, lock = print_lock)
            
        except TimeoutExpired:
            print_warning(f"[TIME LIMIT] Valutazione interrotta per sequenza '{seq}': Raggiunto il limite globale del Job.", lock = print_lock)
            fitnesses.append(float('inf'))                   # Penalità massima per timeout
        
        except Exception as e:
            print_error(f"CRITICAL EVAL ERROR of \"{seq}\": {e}", lock = print_lock)
            fitnesses.append(float('inf'))                   # Penalità massima in caso di fallimento
            
    return fitnesses
