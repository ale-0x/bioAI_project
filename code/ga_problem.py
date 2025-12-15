# ga_problem.py
import os
import numpy            as np
import PeptideBuilder
import random
import subprocess

from Bio.PDB        import PDBIO
from meeko          import MoleculePreparation, PDBQTWriterLegacy
from PeptideBuilder import Geometry, Structure
from rdkit          import Chem
from rdkit.Chem     import AllChem
from typing         import Any, Dict, List


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


# --- Funzioni Helper e pipeline di preparazione (RDKit + Meeko) per Struttura e Docking ---

def clean_structure_elements(structure: Structure) -> None:
    """
    Utility per assicurare che BioPython scriva gli elementi chimici corretti.
    """
    for atom in structure.get_atoms():
        # Se l'elemento manca, lo intuiamo dal nome (es. "CA" -> "C")
        if not atom.element:
            # Prende il primo carattere del nome (es. "CB" -> "C")
            atom.element = atom.name.strip()[0].upper()
        else:
            # Assicura che sia uppercase (es. "c" -> "C")
            atom.element = atom.element.upper()

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

    # Corregge i simboli degli elementi per OpenBabel
    clean_structure_elements(structure)

    # Salva il PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename)

def prepare_ligand_meeko(pdb_path: str, pdbqt_output_path: str) -> bool:
    """
    Pipeline Aggiornata (Robustezza RDKit + Meeko Fix Tuple):
    1. RDKit: Legge PDB (sanitize=False per evitare errori di valenza).
    2. Meeko: Gestisce il return type che potrebbe essere una tupla.
    """
    try:
        # --- A. RDKit: Lettura "Gentile" ---
        # sanitize=False è CRUCIALE con PeptideBuilder.
        # Evita che RDKit vada in crash se gli atomi sono troppo vicini.
        mol = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
        
        if mol is None:
            return False

        # Tentativo di correzione chimica manuale
        try:
            mol.UpdatePropertyCache(strict=False)
            # FastFindRings serve per definire aromaticità e cicli
            Chem.GetSymmSSSR(mol) 
        except:
            pass # Se fallisce, proviamo a continuare lo stesso

        # Aggiunge idrogeni (se mancano)
        mol = Chem.AddHs(mol, addCoords=True)

        # --- B. RDKit: Centratura ---
        conf = mol.GetConformer()
        coords = conf.GetPositions()
        centroid = np.mean(coords, axis=0)
        
        target_center = np.array([CENTER_X, CENTER_Y, CENTER_Z])
        translation = target_center - centroid
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            new_pos = (pos.x + translation[0], pos.y + translation[1], pos.z + translation[2])
            conf.SetAtomPosition(i, new_pos)

        # --- C. Meeko v0.5: Generazione PDBQT ---
        preparator = MoleculePreparation()
        
        # Prepara la molecola
        mol_setups = preparator.prepare(mol)
        
        if mol_setups:
            # Ottieni l'output grezzo
            pdbqt_data = PDBQTWriterLegacy.write_string(mol_setups[0])
            
            # --- FIX TUPLE ERROR ---
            # Se Meeko restituisce (string, info), prendiamo solo string
            if isinstance(pdbqt_data, tuple):
                pdbqt_string = pdbqt_data[0]
            else:
                pdbqt_string = pdbqt_data
            
            # Scrittura su file
            with open(pdbqt_output_path, "w") as f:
                f.write(pdbqt_string)
            return True
        else:
            return False

    except Exception as e:
        # Stampa l'errore ma non bloccare tutto il programma
        print(f"Error in ligand prep (RDKit/Meeko) for {pdb_path}: {e}")
        return False



# def run_vina_real(pdbqt_ligand: str) -> float:
#     """
#     Esegue il Docking molecolare utilizzando AutoDock Vina.

#     Docka il ligando (peptide) preparato sulla proteina recettore 
#     configurata in `constants.py` e parsa l'energia di legame (fitness).

#     Parameters
#     ----------
#     pdbqt_ligand : `str`
#         Percorso al file PDBQT del peptide da dockare.

#     Returns
#     -------
#     `float`
#         L'energia di legame Vina stimata in kcal/mol (valore negativo). 
#         Ritorna 0.0 in caso di fallimento o se il recettore non è trovato.

#     Notes
#     -----
#     Richiede che l'eseguibile `vina` sia installato e nel PATH. 
#     Viene usato un seed fisso (42) per garantire la riproducibilità, 
#     nonostante la natura stocastica della ricerca di Vina. 
#     """
#     if not os.path.exists(RECEPTOR_FILE):
#         return 0.0 # Fallback se manca il recettore
        
#     out_file = pdbqt_ligand.replace(".pdbqt", "_out.pdbqt")
    
#     cmd = [
#         "vina",
#         "--receptor", RECEPTOR_FILE,
#         "--ligand", pdbqt_ligand,
#         "--center_x", str(CENTER_X), "--center_y", str(CENTER_Y), "--center_z", str(CENTER_Z),
#         "--size_x"  , str(SIZE_X)  , "--size_y"  , str(SIZE_Y)  , "--size_z"  , str(SIZE_Z),
#         "--exhaustiveness", str(EXHAUSTIVENESS),
#         "--out", out_file,
#         "--cpu", "6"                # 1 CPU per processo (parallelismo gestito da inspyred se necessario)
#     ]
    
#     try:
#         # result = subprocess.run(cmd, capture_output = True, text = True)

#         # Rimuovi capture_output=True e il risultato verrà stampato direttamente
#         # Usa check=True per sollevare un'eccezione in caso di errore di Vina
#         result = subprocess.run(
#             cmd, 
#             check  = True,  
#             text   = True,
#             stdout = subprocess.PIPE
#         )
        
#         # Parsing dell'output di Vina per trovare l'affinità migliore
#         # L'output contiene linee come: "   1        -8.5      0.000      0.000"
#         best_affinity = 0.0
#         for line in result.stdout.splitlines():
#             if line.strip().startswith("1"):
#                 parts = line.split()
#                 if len(parts) >= 2:
#                     try:
#                         best_affinity = float(parts[1])
#                         break
#                     except ValueError:
#                         continue
#         return best_affinity
#     except subprocess.CalledProcessError:
#         return 0.0                  # Vina fallisce occasionalmente se la geometria è pessima
#     except Exception as e:
#         print(f"Vina Exception: {e}")
#         return 0.0

def run_vina_real(pdbqt_ligand: str) -> float:
    """Esegue Vina in modalità DEBUG per capire l'errore."""
    if not os.path.exists(RECEPTOR_FILE):
        print(f"ERRORE CRITICO: Il file recettore non esiste: {RECEPTOR_FILE}")
        return 0.0
        
    out_file = pdbqt_ligand.replace(".pdbqt", "_out.pdbqt")
    
    cmd = [
        "vina",
        "--receptor", RECEPTOR_FILE,
        "--ligand", pdbqt_ligand,
        "--center_x", str(CENTER_X), "--center_y", str(CENTER_Y), "--center_z", str(CENTER_Z),
        "--size_x"  , str(SIZE_X)  , "--size_y"  , str(SIZE_Y)  , "--size_z"  , str(SIZE_Z),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--out", out_file,
        "--cpu", "1" 
    ]
    
    try:
        # --- MODIFICA DEBUG: capture_output=True e check=False ---
        # Catturiamo sia stdout che stderr per vederli
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        # Se Vina ha scritto qualcosa su stderr (errori), stampiamolo!
        if result.stderr:
            print(f"\n--- VINA ERROR LOG ({pdbqt_ligand}) ---")
            print(result.stderr)
            print("---------------------------------------")

        # Parsing dell'output
        best_affinity = 0.0
        found = False
        
        for line in result.stdout.splitlines():
            # Cerca la riga che inizia con "   1" (il primo modo)
            if line.strip().startswith("1"):
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        best_affinity = float(parts[1])
                        found = True
                        break
                    except ValueError:
                        continue
        
        if not found:
            print(f"WARNING: Vina non ha prodotto affinity valide per {pdbqt_ligand}")
            # Stampa l'output standard per capire perché
            # print(result.stdout) 

        return best_affinity

    except Exception as e:
        print(f"Vina Exception: {e}")
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
    
    if not os.path.exists(TEMP_DOCKING):
        os.makedirs(TEMP_DOCKING, exist_ok = True)

    for i, seq in enumerate(candidates):
        # ID univoco per evitare collisioni di file
        unique_id  = f"{i}_{random.randint(10000,99999)}"
        base_name  = os.path.join(TEMP_DOCKING, f"seq_{unique_id}")
        pdb_file   = f"{base_name}.pdb"
        pdbqt_file = f"{base_name}.pdbqt"
        out_file   = f"{base_name}_out.pdbqt"
        
        try:
            # 1. Generazione Struttura
            generate_pdb_fast(seq, pdb_file)
            
            # 2. Pybel -> RDKit (Centratura) -> Meeko
            success = prepare_ligand_meeko(pdb_file, pdbqt_file)
            
            # 3. Docking
            if success:
                energy = run_vina_real(pdbqt_file)
                fitnesses.append(energy)
            else:
                fitnesses.append(0.0)               # Penalità per fallimento prep
            
            # # Pulizia file temporanei
            # for f in [pdb_file, pdbqt_file, out_file]:
            #     if os.path.exists(f):
            #         os.remove(f)
            
        except Exception as e:
            print(f"CRITICAL EVAL ERROR of \"{seq}\": {e}")
            fitnesses.append(0.0)                   # Penalità massima in caso di fallimento
            
    return fitnesses
