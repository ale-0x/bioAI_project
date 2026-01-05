# ga_problem.py
import os
import numpy as np
import random
import re
import shutil
import subprocess
import uuid

from openbabel  import openbabel, pybel
from rdkit      import Chem
from rdkit.Chem import AllChem
from typing     import Any, Dict, List, Tuple
from utils      import print_error, print_verbose


import constants as C

# Hydrophobicity index at pH 7 
# (from Monera, O. D., et al. Relationship of Sidechain Hydrophobicity and Alpha-Helical Propensity on the Stability of the Single-Stranded Amphipathic Alpha-Helix. J Pept Sci. (1995).)

amino_hydrophobicity_ph7 = {
    'L': 97,   # Leucine
    'I': 99,   # Isoleucine
    'F': 100,  # Phenylalanine
    'W': 97,   # Tryptophan
    'V': 76,   # Valine
    'M': 74,   # Methionine
    'C': 49,   # Cysteine
    'Y': 63,   # Tyrosine
    'A': 41,   # Alanine
    'T': 13,   # Threonine
    'E': -31,  # Glutamate
    'H': 8,    # Histidine
    'G': 0,    # Glycine
    'S': -5,   # Serine
    'Q': -10,  # Glutamine
    'D': -55,  # Aspartate
    'R': -14,  # Arginine
    'K': -23,  # Lysine
    'N': -28,  # Asparagine
    'P': 5     # Proline (using the 6.5 value provided)
}

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

def generate_pdb_rdkit(sequence: str, output_file: str):
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
        
        Chem.MolToPDBFile(mol, output_file)
    except Exception as e:
        print_error(f"Errore generazione RDKit per {sequence}: {e}")

def prepare_ligand_openbabel(pdb_path: str, pdbqt_output_path: str, center: Tuple[float, float, float]) -> bool:
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
    try:
        # 1. Leggi PDB (Generato da RDKit, quindi sicuro)
        mol = next(pybel.readfile("pdb", pdb_path))
        
        # 2. Centratura nella tasca (Box Vina)
        atoms         = [atom.coords for atom in mol] 
        centroid      = np.mean(atoms, axis = 0)
        target_center = np.array([center[0], center[1], center[2]])
        move_v        = target_center - centroid
        mol.OBMol.Translate(openbabel.vector3(move_v[0], move_v[1], move_v[2]))
        
        # 3. Scrittura PDBQT (OpenBabel calcola le cariche Gasteiger automaticamente)
        mol.write("pdbqt", pdbqt_output_path, overwrite = True)
        
        return os.path.exists(pdbqt_output_path) and os.path.getsize(pdbqt_output_path) > 0
    except StopIteration:
        print_error(f"Errore OpenBabel: File PDB vuoto o invalido: {pdb_path}")
        return False
    except Exception as e:
        print_error(f"Conversione fallita per {pdb_path}: {e}")
        return False
    

def run_vina_real(pdbqt_ligand: str, receptor_file: str, center: Tuple[float, float, float], box_size: Tuple[int, int, int], cpu: int = 1, exhaustiveness: int = C.EXHAUSTIVENESS, verbose: bool = False) -> float:
    """
    Esegue il docking molecolare lanciando il processo AutoDock Vina e restituisce l'energia di legame migliore.

    La funzione costruisce ed esegue il comando da terminale per Vina, specificando
    il recettore, il ligando, la search box e i parametri di precisione.
    Cattura l'output standard per estrarre il punteggio di affinità (energia libera)
    del primo modo di binding (il migliore).

    Parameters
    ----------
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
    
    cpu : `int`, default `1`, optional
        Numero di CPU/thread da dedicare a questa singola esecuzione di Vina.
    
    vina_exe_path : `str`, default `'vina'`, optional
        Percorso dell'eseguibile di AutoDock Vina.
    
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
    if not os.path.exists(receptor_file):
        print_error(f"ERRORE CRITICO: Il file recettore non esiste: {receptor_file}")
        return 0.0
            
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
        "--cpu"           , "1"
        # "--score_only" 
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output = True,
            text           = True,
            check          = True
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
                '#################################################################'
            }:
                print_verbose(line, to_print = verbose)
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
                        return float(parts[1])
        return 0.0 # Se il parsing fallisce
    except subprocess.CalledProcessError as e:
        # Vina ha restituito un codice di errore
        print_error(f"Vina Error Output: {e.stderr}")
        return 0.0
    except Exception as e:
        print_error(f"Errore esecuzione Vina (subprocess): {e}")
        return 0.0


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
    verbose       = args.get('verbose', False)
    receptor_file = args.get('receptor_file', C.RECEPTOR_FILE)
    
    cx = args.get('center_x', C.CENTER_X)
    cy = args.get('center_y', C.CENTER_Y)
    cz = args.get('center_z', C.CENTER_Z)
    
    sx = args.get('size_x', C.SIZE_X)
    sy = args.get('size_y', C.SIZE_Y)
    sz = args.get('size_z', C.SIZE_Z)

    # Cartella temporanea del JOB (deve esistere, creata da main.py)
    # Se per qualche motivo non c'è, usa /tmp locale con fallback
    base_temp = args.get('temp_dir', f"/tmp/ga_fallback_{uuid.uuid4()}")
    if not os.path.exists(base_temp):
        print_verbose(f"[evaluate_peptide_binding] Creazione directory temporanea '{base_temp}' ...", to_print = verbose)
        os.makedirs(base_temp, exist_ok = True)
        print_verbose(f"[evaluate_peptide_binding] Done.", to_print = verbose)

    fitnesses = []

    for candidate in candidates:
        # Estrai la sequenza se è un oggetto Individual
        seq = candidate.candidate if hasattr(candidate, 'candidate') else candidate

        # Calcola la MEDIA dei valori di idrofobicità per la sequenza
        hydrophobicity = sum(amino_hydrophobicity_ph7.get(aa, 0) for aa in seq) / len(seq)
        
        # Crea cartella isolata per questo singolo peptide
        unique_id     = str(uuid.uuid4())
        unique_folder = f"p_{unique_id}"
        work_dir      = os.path.join(base_temp, unique_folder)
        print_verbose(f"[evaluate_peptide_binding] Creazione directory temporanea '{base_temp}' ...", to_print = verbose)
        os.makedirs(work_dir, exist_ok = True)
        print_verbose(f"[evaluate_peptide_binding] Creazione directory temporanea '{base_temp}' ...", to_print = verbose)

        base_name  = os.path.join(work_dir, f"seq_{unique_id}")
        pdb_file   = f"{base_name}.pdb"
        pdbqt_file = f"{base_name}.pdbqt"
        out_file   = f"{base_name}_out.pdbqt"
        
        try:
            # 1. Generazione Struttura
            print_verbose(f"[evaluate_peptide_binding] Generazione Struttura di '{seq}' in '{pdb_file}' ...", to_print = verbose)
            generate_pdb_rdkit(seq, pdb_file)
            print_verbose(f"[evaluate_peptide_binding] Generazione Struttura di '{seq}' in '{pdb_file}' --> Done.", to_print = verbose)
            
            # 2. Pybel -> RDKit (Centratura) -> Meeko e poi Docking
            print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' ...", to_print = verbose)
            if prepare_ligand_openbabel(pdb_file, pdbqt_file, (cx, cy, cz)):
                print_verbose(f"[evaluate_peptide_binding] Valutazione di Autodock Vina su '{pdbqt_file}' della sequenza '{seq}' ...", to_print = verbose)
                energy = run_vina_real(
                    pdbqt_ligand   = pdbqt_file,
                    receptor_file  = receptor_file,
                    center         = (cx, cy, cz),
                    box_size       = (sx, sy, sz),
                    cpu            = 1,             # 1 CPU per processo figlio
                    exhaustiveness = args.get('exhaustiveness', C.EXHAUSTIVENESS),
                    verbose        = verbose
                )
                print_verbose(f"[evaluate_peptide_binding] Valutazione di Autodock Vina su '{pdbqt_file}' della sequenza '{seq}' --> Done", to_print = verbose)
                # usa la hydrophobicity sum come penalità per restringere il campo di soluzioni possibili
                fitnesses.append(energy + (hydrophobicity*C.HYDROPHOBICITY_WEIGHT))
                print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' --> Done", to_print = verbose)
            else:
                fitnesses.append(0.0)               # Penalità per fallimento prep
                print_verbose(f"[evaluate_peptide_binding] Conversione di '{pdb_file}' in '{pdbqt_file}' e centrato in ({cx}, {cy}, {cz}) della sequenza '{seq}' --> Fallito", to_print = verbose)
            
        except Exception as e:
            print_error(f"CRITICAL EVAL ERROR of \"{seq}\": {e}")
            fitnesses.append(0.0)                   # Penalità massima in caso di fallimento
        finally:
            # Pulizia: rimuovi la cartella del singolo peptide
            if os.path.exists(work_dir) and not args.get("no_delete", False):
                shutil.rmtree(work_dir)
            
    return fitnesses
