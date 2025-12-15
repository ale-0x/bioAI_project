# get_vina_coords.py
import numpy as np
import os
from pymol import cmd

from constants import INPUT_PDB_FILE

# --- 1. Definizione ---

LIGAND_CODE = "3TL"          # Il codice del residuo del tuo ligando
OBJECT_NAME = "2P3D"         # Nome che PyMOL assegna di default al caricamento di 2P3D.pdb

# Assumiamo che il nome dell'oggetto caricato sia "target"
LIGAND_SELECTION = f"resn {LIGAND_CODE} and {OBJECT_NAME}"

# Funzione per calcolare il centro geometrico (centroid) di una selezione
def get_centroid(selection: str) -> np.ndarray:
    """
    Calcola il centro geometrico (centroid) di una selezione PyMOL.
    """
    model = cmd.get_model(selection)
    if not model.atom:
        raise ValueError(f"Selezione vuota: {selection}. Controlla l'ID del ligando.")
        
    coords = np.array([a.coord for a in model.atom])
    return coords.mean(axis=0)

# --- 2. Esecuzione ---
if __name__ == '__main__':
    if not os.path.exists(INPUT_PDB_FILE):
        print(f"ERRORE: File {INPUT_PDB_FILE} non trovato. Scaricalo e posizionalo qui.")
    else:
        try:
            # Carica la struttura PDB
            cmd.load(INPUT_PDB_FILE)
            
            # Calcola il centro del ligando co-cristallizzato
            center_coords = get_centroid(LIGAND_SELECTION)
            
            # --- 3. Risultato ---
            print("\n=======================================================")
            print("Coordinate Vina (Centro Tasca)")
            print("=======================================================")
            print(f"Centro Tasca (X, Y, Z):")
            # Stampa le coordinate con 3 decimali
            print(f"X: {center_coords[0]:.3f}, Y: {center_coords[1]:.3f}, Z: {center_coords[2]:.3f}")
            
            print("\n-------------------------------------------------------")
            print("DA INSERIRE IN constants.py:")
            print(f"CENTER_X, CENTER_Y, CENTER_Z = {center_coords[0]:.3f}, {center_coords[1]:.3f}, {center_coords[2]:.3f}")
            print("-------------------------------------------------------")
            
            # Il Box Size DEVE essere stimato visivamente!
            print("\nBOX SIZE (Dimensione):")
            print("Il box deve essere sufficientemente grande da contenere il peptide (10AA).")
            print("Si consiglia un box di almeno 20x20x20 o 25x25x25 Ångstrom.")
            print("Esempio: SIZE_X, SIZE_Y, SIZE_Z       = 25, 25, 25")
            
            # Pulisci l'ambiente PyMOL
            cmd.delete("all")

        except Exception as e:
            print(f"Errore durante l'analisi PyMOL: {e}")
            print("Assicurati che PyMOL sia installato correttamente e che la selezione LIGAND_SELECTION sia corretta.")
