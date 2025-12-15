# ga_problem.py
import os
import PeptideBuilder
import random
import subprocess

from typing import List, Dict, Any
from Bio.PDB import PDBIO
from PeptideBuilder import Geometry

from constants import (
    AMINO_ACIDS, PEPTIDE_LENGTH, 
    RECEPTOR_FILE, CENTER_X, CENTER_Y, CENTER_Z, 
    SIZE_X, SIZE_Y, SIZE_Z, EXHAUSTIVENESS,
    TEMP_DOCKING
)

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
    length: int = args.get('peptide_length', PEPTIDE_LENGTH)

    return "".join(random.choice(AMINO_ACIDS) for _ in range(length))



# --- Funzioni Helper per Struttura e Docking ---

def generate_pdb_fast(sequence: str, filename: str) -> None:
    """
    Genera un file PDB (struttura 3D) per la sequenza peptidica.

    Utilizza la libreria `PeptideBuilder` per creare istantaneamente una
    struttura lineare estesa, assumendo una geometria standard per i legami.
    Questo metodo è estremamente veloce ed è cruciale per la velocità del GA.

    Parameters
    ----------
    sequence : `str`
        La sequenza peptidica (stringa di aminoacidi) da modellare.
    filename : `str`
        Il percorso completo dove salvare il file PDB generato.

    Returns
    -------
    Nessuno. Salva la struttura nel file specificato.

    Examples
    --------
    >>> generate_pdb_fast("ALTSV", "temp/test.pdb")
    # Viene creato il file temp/test.pdb
    """
    # Crea una struttura lineare estesa (phi=-180, psi=180)
    structure = PeptideBuilder.make_structure(
        sequence, 
        [180] * len(sequence), 
        [180] * len(sequence)
    )

    # Salva il PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)

def convert_pdb_to_pdbqt(pdb_file: str, pdbqt_file: str) -> None:
    """
    Converte un file PDB in formato PDBQT usando OpenBabel.

    Il formato PDBQT è richiesto da AutoDock Vina e contiene atomi di idrogeno 
    e cariche parziali (qui calcolate con il metodo Gasteiger).

    Parameters
    ----------
    pdb_file : `str`
        Percorso al file PDB di input.
    pdbqt_file : `str`
        Percorso al file PDBQT di output.

    Returns
    -------
    Nessuno. Salva il ligando preparato nel file specificato.

    Notes
    -----
    Richiede che l'eseguibile `obabel` sia installato e nel PATH.
    """
    # Comando: obabel -ipdb file.pdb -opdbqt -O file.pdbqt --partialcharge gasteiger
    cmd = [
        "obabel", "-ipdb", pdb_file, "-opdbqt", "-O", pdbqt_file, "--partialcharge", "gasteiger"
    ]
    subprocess.run(cmd, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)

def run_vina_real(pdbqt_ligand: str) -> float:
    """
    Esegue il Docking molecolare utilizzando AutoDock Vina.

    Docka il ligando (peptide) preparato sulla proteina recettore 
    configurata in `constants.py` e parsa l'energia di legame (fitness).

    Parameters
    ----------
    pdbqt_ligand : `str`
        Percorso al file PDBQT del peptide da dockare.

    Returns
    -------
    `float`
        L'energia di legame Vina stimata in kcal/mol (valore negativo). 
        Ritorna 0.0 in caso di fallimento o se il recettore non è trovato.

    Notes
    -----
    Richiede che l'eseguibile `vina` sia installato e nel PATH. 
    Viene usato un seed fisso (42) per garantire la riproducibilità, 
    nonostante la natura stocastica della ricerca di Vina. 
    """
    if not os.path.exists(RECEPTOR_FILE):
        return 0.0 # Fallback se manca il recettore
        
    out_file = pdbqt_ligand.replace(".pdbqt", "_out.pdbqt")
    
    cmd = [
        "vina",
        "--receptor", RECEPTOR_FILE,
        "--ligand", pdbqt_ligand,
        "--center_x", str(CENTER_X), "--center_y", str(CENTER_Y), "--center_z", str(CENTER_Z),
        "--size_x"  , str(SIZE_X)  , "--size_y"  , str(SIZE_Y)  , "--size_z"  , str(SIZE_Z),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--out", out_file,
        "--cpu", "6"                # 1 CPU per processo (parallelismo gestito da inspyred se necessario)
    ]
    
    try:
        # result = subprocess.run(cmd, capture_output = True, text = True)

        # Rimuovi capture_output=True e il risultato verrà stampato direttamente
        # Usa check=True per sollevare un'eccezione in caso di errore di Vina
        result = subprocess.run(
            cmd, 
            # NON usare capture_output=True, stdout=PIPE o stderr=PIPE
            check=True,  
            text=True  
        )
        
        # Parsing dell'output di Vina per trovare l'affinità migliore
        # L'output contiene linee come: "   1        -8.5      0.000      0.000"
        best_affinity = 0.0
        for line in result.stdout.splitlines():
            if line.strip().startswith("1"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        best_affinity = float(parts[1])
                        break
                    except ValueError:
                        continue
        return best_affinity

    except Exception as e:
        print(f"Vina Error: {e}")
        return 0.0



# --- Funzione di Valutazione Principale ---

def evaluate_peptide_binding(candidates: List[str], args: Dict[str, Any]) -> List[float]:
    """
    Valutazione della Fitness: Calcola l'energia di legame tramite Docking Vina.

    Implementa il workflow completo:
    1. Generazione della struttura 3D (PeptideBuilder).
    2. Preparazione del ligando (OpenBabel PDB -> PDBQT).
    3. Esecuzione del docking (AutoDock Vina) per ottenere l'energia di legame.

    Parameters
    ----------
    candidates : `List[str]`
        Una lista di sequenze peptidiche (stringhe) da valutare.
    args : `Dict[str, Any]`
        Dizionario di argomenti opzionali.

    Returns
    -------
    `List[float]`
        Una lista di energie di legame Vina (in kcal/mol) per ogni candidato.
        Si mira a MINIMIZZARE questo valore (massima affinità).
    """
    fitnesses = []
    
    fitnesses = []
    
    if not os.path.exists(TEMP_DOCKING):
        os.makedirs(TEMP_DOCKING)

    for i, seq in enumerate(candidates):
        # Genera un nome base univoco
        base_name = os.path.join(TEMP_DOCKING, f"cand_{i}_{seq[:5]}_{random.randint(1000,9999)}")
        pdb_file = f"{base_name}.pdb"
        pdbqt_file = f"{base_name}.pdbqt"
        
        try:
            # 1. Generazione Struttura
            generate_pdb_fast(seq, pdb_file)
            
            # 2. Conversione
            convert_pdb_to_pdbqt(pdb_file, pdbqt_file)
            
            # 3. Docking
            energy = run_vina_real(pdbqt_file)
            
            # 4. Pulizia (opzionale, ma consigliata per non riempire il disco)
            if os.path.exists(pdb_file)  : os.remove(pdb_file)
            if os.path.exists(pdbqt_file): os.remove(pdbqt_file)
            
            fitnesses.append(energy)
            
        except Exception as e:
            print(f"Fail during evaluation of {seq}: {e}")
            fitnesses.append(0.0) # Penalità massima in caso di fallimento
            
    return fitnesses
