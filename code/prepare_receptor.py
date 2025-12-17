# prepare_receptor.py
import os
import sys

# Importiamo l'interfaccia Python di OpenBabel
from openbabel import pybel
from openbabel import openbabel

from constants import INPUT_PDB_FILE, OUTPUT_PDBQT_FILE


def prepare_receptor_pdbqt(input_pdb_file: str, output_pdbqt_file: str) -> bool:
    """
    Prepara il recettore da PDB a PDBQT usando OpenBabel.

    Questo passaggio è necessario per aggiungere idrogeni e cariche al recettore
    prima di eseguire Vina. Rimuove anche acqua e altri eteroatomi (flag -xr).

    Steps:
    1. Legge il PDB.
    2. Rimuove l'acqua (HOH) e i sali.
    3. Aggiunge gli idrogeni polari.
    4. Calcola le cariche parziali (Gasteiger) e scrive il PDBQT.

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
    
    try:
        print(f"Lettura PDB: {input_pdb_file}...")
        
        # 1. Lettura del file PDB
        mol = next(pybel.readfile("pdb", input_pdb_file))
        obmol = mol.OBMol # Oggetto C++ sottostante
        
        # 2. Pulizia Manuale (Rimozione HOH e Ioni)
        initial_residues = obmol.NumResidues()
        deleted_count = 0
        solvent_names = {"HOH", "WAT", "CL", "NA", "MG", "K", "SO4", "PO4"}
        
        for i in range(initial_residues - 1, -1, -1):
            res = obmol.GetResidue(i)
            if res.GetName().strip() in solvent_names:
                obmol.DeleteResidue(res)
                deleted_count += 1
                
        if deleted_count > 0:
            print(f"  - Rimossi {deleted_count} residui di solvente/ioni.")

        # 3. FIX KEKULIZZAZIONE (Rigenerazione Topologia)
        print("  - Ricalcolo della connettività chimica (Fix Kekulization)...")
        obmol.DeleteHydrogens() 
        obmol.ConnectTheDots() 
        obmol.PerceiveBondOrders()

        # 4. Aggiunta Idrogeni (pH 7.4) e Cariche
        print("  - Aggiunta idrogeni polari e calcolo cariche...")
        obmol.AddHydrogens(False, True, 7.4) 

        # 5. Scrittura PDBQT (OpenBabel crea un file con ROOT/BRANCH qui)
        print(f"Scrittura PDBQT: {output_pdbqt_file}...")
        mol.write("pdbqt", output_pdbqt_file, overwrite=True)
        
        # --- [NUOVO PASSAGGIO] 6. Post-processing per Vina (Rimuove flessibilità) ---
        # Riapriamo il file appena creato per rimuovere i tag che mandano in crash Vina
        print("  - Post-processing: Rimozione tag flessibilità (ROOT/BRANCH)...")
        
        with open(output_pdbqt_file, 'r') as f:
            lines = f.readlines()
            
        with open(output_pdbqt_file, 'w') as f:
            for line in lines:
                # Scriviamo la riga SOLO se NON è un tag di flessibilità
                if not (
                    line.startswith('ROOT')      or 
                    line.startswith('ENDROOT')   or 
                    line.startswith('BRANCH')    or 
                    line.startswith('ENDBRANCH') or 
                    line.startswith('TORSDOF')
                ):
                    f.write(line)
        # --------------------------------------------------------------------------

        print("Preparazione completata con successo.")
        return True

    except Exception as e:
        print(f"ERRORE CRITICO durante la preparazione: {e}")
        import traceback
        traceback.print_exc()
        return False

def fix_receptor_pdbqt(pdbqt_path):
    """
    Legge un file PDBQT di recettore e rimuove i tag ROOT/BRANCH/TORSDOF
    che causano l'errore 'Unknown or inappropriate tag' in Vina.
    Sovrascrive il file esistente con la versione corretta (rigida).
    """
    with open(pdbqt_path, 'r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        # Ignora i tag che definiscono la flessibilità (tipici dei ligandi)
        if line.startswith('ROOT') or \
           line.startswith('ENDROOT') or \
           line.startswith('BRANCH') or \
           line.startswith('ENDBRANCH') or \
           line.startswith('TORSDOF'):
            continue
        
        # Ignora i REMARK che elencano le torsioni attive (pulizia opzionale ma consigliata)
        if line.startswith('REMARK') and ('active torsions' in line or 'between atoms' in line):
            continue

        new_lines.append(line)

    # Sovrascrive il file con la versione pulita
    with open(pdbqt_path, 'w') as f:
        f.writelines(new_lines)
    
    print(f"--- File recettore corretto per Vina (rimossi tag flessibilità): {pdbqt_path} ---")


if __name__ == '__main__':
    print("--- Utilità di Preparazione del Recettore Vina ---")
    
    # 1. Verifichiamo i prerequisiti (es. file 7cam.pdb)
    if not os.path.exists(INPUT_PDB_FILE):
        print(f"\n[INFO] Necessario: Devi prima scaricare la struttura PDB e salvarla come '{INPUT_PDB_FILE}'.")
        sys.exit(1)
        
    # 2. Eseguiamo la preparazione
    success = prepare_receptor_pdbqt(INPUT_PDB_FILE, OUTPUT_PDBQT_FILE)
    fix_receptor_pdbqt(OUTPUT_PDBQT_FILE)
    
    if success:
        print("\nPronto per eseguire l'Algoritmo Genetico.")
        print(f"Prossimo passo: Eseguire 'python main.py'")
