# prepare_receptor.py
import os
import sys
import subprocess

from typing import Optional

from constants import INPUT_PDB_FILE, OUTPUT_PDBQT_FILE


def prepare_receptor_pdbqt(input_pdb_file: str, output_pdbqt_file: str) -> bool:
    """
    Prepara il recettore da PDB a PDBQT usando OpenBabel.

    Questo passaggio è necessario per aggiungere idrogeni e cariche al recettore
    prima di eseguire Vina. Rimuove anche acqua e altri eteroatomi (flag -xr).

    Parameters
    ----------
    input_pdb_file : `str`
        Percorso al file PDB originale (es. '7CAM.pdb').
    output_pdbqt_file : `str`
        Percorso dove salvare il file PDBQT preparato (es. '7CAM.pdbqt').

    Returns
    -------
    `bool`
        True se la preparazione è stata completata con successo, False altrimenti.
    """
    if not os.path.exists(input_pdb_file):
        print(f"ERRORE: File PDB di input non trovato: {input_pdb_file}")
        print("Assicurati di aver scaricato il file .pdb e rinominato correttamente.")
        return False

    # Comando OpenBabel per la preparazione del recettore
    # -xr: Rimuove HETATM (acqua, ligandi cristallizzati, etc.)
    # --partialcharge gasteiger: Calcola le cariche parziali
    cmd = [
        "obabel", 
        "-ipdb", input_pdb_file, 
        "-opdbqt", 
        "-O", output_pdbqt_file, 
        "-xr", 
        "--partialcharge", "gasteiger"
    ]
    
    try:
        print(f"Avvio preparazione: {input_pdb_file} -> {output_pdbqt_file}")
        # Esegue il comando in modo silenzioso
        subprocess.run(cmd, check = True, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
        print(f"Preparazione recettore completata con successo.")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"ERRORE OpenBabel (Codice {e.returncode}). Controlla l'output di errore:")
        print(f"{e.stderr.decode()}")
        return False
        
    except FileNotFoundError:
        print(f"ERRORE: Il comando 'obabel' non è stato trovato.")
        print("Assicurati che OpenBabel sia installato e che l'ambiente Conda sia attivo.")
        return False


if __name__ == '__main__':
    print("--- Utilità di Preparazione del Recettore Vina ---")
    
    # 1. Verifichiamo i prerequisiti (es. file 7cam.pdb)
    if not os.path.exists(INPUT_PDB_FILE):
        print(f"\n[INFO] Necessario: Devi prima scaricare la struttura PDB e salvarla come '{INPUT_PDB_FILE}'.")
        sys.exit(1)
        
    # 2. Eseguiamo la preparazione
    success = prepare_receptor_pdbqt(INPUT_PDB_FILE, OUTPUT_PDBQT_FILE)
    
    if success:
        print("\nPronto per eseguire l'Algoritmo Genetico.")
        print(f"Prossimo passo: Eseguire 'python main.py'")
