# ga_problem.py
import os
import numpy            as np
import random
import re
import subprocess

from openbabel      import openbabel, pybel
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

# def generate_pdb_rdkit(seq: str, filename: str):
#     """
#     Genera PDB con RDKit: chimica corretta, idrogeni presenti, geometria rilassata.
#     """
#     # 1. Crea Molecola da Sequenza
#     mol = Chem.MolFromSequence(seq)
    
#     # 2. Aggiunge Idrogeni (Essenziale per il 3D)
#     mol = Chem.AddHs(mol)
    
#     # 3. Genera 3D (Embedding)
#     params = AllChem.ETKDGv3()
#     params.useRandomCoords = True 
#     if AllChem.EmbedMolecule(mol, params) == -1:
#         # Fallback se fallisce il primo tentativo
#         params.useRandomCoords = True
#         AllChem.EmbedMolecule(mol, params)
        
#     # 4. Minimizzazione Energetica (Rilassa la struttura)
#     try:
#         AllChem.MMFFOptimizeMolecule(mol)
#     except:
#         pass # Se MMFF fallisce, usiamo comunque la struttura generata

#     # 5. Salva PDB
#     Chem.MolToPDBFile(mol, filename)

# In ga_problem.py

def generate_pdb_rdkit(sequence: str, output_file: str):
    try:
        mol = Chem.MolFromSequence(sequence)
        mol = Chem.AddHs(mol)
        
        # 1. Embedding (Generazione coordinate iniziali)
        # Usa ETKDGv3 che è più robusto per molecole grandi/cicliche
        params = AllChem.ETKDGv3()
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(mol, params)
        
        if result == -1:
            # Fallback se fallisce
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
        # 2. MINIMIZZAZIONE (Fondamentale!)
        # Rilassa la struttura per evitare atomi sovrapposti
        AllChem.MMFFOptimizeMolecule(mol)
        
        Chem.MolToPDBFile(mol, output_file)
    except Exception as e:
        print(f"Errore generazione RDKit per {sequence}: {e}")
        # Gestisci l'errore o crea un file vuoto

def prepare_ligand_openbabel(pdb_path: str, pdbqt_output_path: str) -> bool:
    try:
        # 1. Leggi PDB (Generato da RDKit, quindi sicuro)
        mol = next(pybel.readfile("pdb", pdb_path))
        
        # 2. Centratura nella tasca (Box Vina)
        atoms = [atom.coords for atom in mol] 
        centroid = np.mean(atoms, axis=0)
        target_center = np.array([CENTER_X, CENTER_Y, CENTER_Z])
        move_v = target_center - centroid
        mol.OBMol.Translate(openbabel.vector3(move_v[0], move_v[1], move_v[2]))
        
        # 3. Scrittura PDBQT (OpenBabel calcola le cariche Gasteiger automaticamente)
        mol.write("pdbqt", pdbqt_output_path, overwrite=True)
        
        return os.path.exists(pdbqt_output_path) and os.path.getsize(pdbqt_output_path) > 0
    except Exception as e:
        print(f"[ERR] Conversione fallita per {pdb_path}: {e}")
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
        "--cpu", "7" 
    ]
    
    try:
        # --- MODIFICA DEBUG: capture_output=True e check=False ---
        # Catturiamo sia stdout che stderr per vederli
        subprocess.run(cmd, capture_output=False, text=True, check=False)
        
        if not os.path.exists(out_file):
            print(f"[ERR] Vina non ha creato il file di output: {out_file}")
            return 0.0

        # Parsing dell'output
        best_affinity = 0.0
        found = False
        
        with open(out_file, 'r') as f:
            for line in f:
                # Cerchiamo la riga: REMARK VINA RESULT: -8.5 0.000 0.000
                if line.startswith("REMARK VINA RESULT:"):
                    parts = line.split()
                    # L'energia è il terzo elemento (indice 3) perché:
                    # parts[0]="REMARK", [1]="VINA", [2]="RESULT:", [3]="-8.5"
                    try:
                        best_affinity = float(parts[3])
                        found = True
                        break # Abbiamo trovato il primo modello (il migliore), usciamo
                    except ValueError:
                        continue
        
        if not found:
            print(f"WARNING: Vina non ha prodotto affinity valide per {pdbqt_ligand}")
            # Stampa l'output standard per capire perché
            # print(result.stdout) 

        return best_affinity

    except subprocess.CalledProcessError:
        print(f"[ERR] Vina è andato in crash sul file {pdbqt_ligand}")
        return 0.0
    except Exception as e:
        print(f"[ERR] Errore generico in run_vina_real: {e}")
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
            generate_pdb_rdkit(seq, pdb_file)
            
            # 2. Pybel -> RDKit (Centratura) -> Meeko
            success = prepare_ligand_openbabel(pdb_file, pdbqt_file)
            
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
