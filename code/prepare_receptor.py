# prepare_receptor.py
from openbabel import pybel
from pathlib   import Path


from constants import INPUT_PDB_FILE, OUTPUT_PDBQT_FILE
from utils     import print_error


def prepare_receptor_pdbqt(input_pdb_file: str | Path, output_pdbqt_file: str | Path) -> bool:
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
    input_pdb_file    = Path(input_pdb_file)
    output_pdbqt_file = Path(output_pdbqt_file)

    if not input_pdb_file.exists():
        print_error(f"Input PDB file not found: {input_pdb_file.resolve()}. Make sure you have downloaded the .pdb file and renamed it correctly.")
        return False
    
    try:
        print(f"PDB reading: {input_pdb_file} ...")
        
        # Reading the PDB file
        mol = next(pybel.readfile("pdb", str(input_pdb_file)))
        obmol = mol.OBMol
        
        # Manual Cleaning (HOH and Ion Removal)
        initial_residues = obmol.NumResidues()
        deleted_count = 0
        solvent_names = {"HOH", "WAT", "CL", "NA", "MG", "K", "SO4", "PO4"}
        
        for i in range(initial_residues - 1, -1, -1):
            res = obmol.GetResidue(i)
            if res.GetName().strip() in solvent_names:
                obmol.DeleteResidue(res)
                deleted_count += 1
                
        if deleted_count > 0:
            print(f"  - Removed {deleted_count} solvent/ion residues.")

        # FIX KEKULIZATION (Topology Regeneration)
        print("  - Recalculating chemical connectivity (Fix Kekulization)...")
        obmol.DeleteHydrogens() 
        obmol.ConnectTheDots() 
        obmol.PerceiveBondOrders()

        # Added Hydrogens (pH 7.4) and Charges
        print("  - Added polar hydrogens and calculated charges...")
        obmol.AddHydrogens(False, True, 7.4) 

        # PDBQT writing (OpenBabel creates a file with ROOT/BRANCH here)
        print(f"PDBQT Writing: {output_pdbqt_file} ...")
        mol.write("pdbqt", str(output_pdbqt_file), overwrite = True)
        
        # Post-processing for Vina (Remove Flexibility) ---
        # Let's reopen the newly created file to remove the tags that crash Vina.
        print("  - Post-processing: Removing flexibility tags (ROOT/BRANCH)...")
        
        with output_pdbqt_file.open('r') as f:
            lines = f.readlines()
            
        with output_pdbqt_file.open('w') as f:
            for line in lines:
                # We write the line ONLY if it is NOT a flexibility tag
                if not (
                    line.startswith('ROOT')      or 
                    line.startswith('ENDROOT')   or 
                    line.startswith('BRANCH')    or 
                    line.startswith('ENDBRANCH') or 
                    line.startswith('TORSDOF')
                ):
                    f.write(line)
        # --------------------------------------------------------------------------

        print("Preparation completed successfully.")
        return True

    except Exception as e:
        print(f"CRITICAL ERROR during preparation: {e}")
        import traceback
        traceback.print_exc()
        return False

def fix_receptor_pdbqt(pdbqt_path: str | Path) -> None:
    """
    Post-processa il file PDBQT del recettore per rimuovere le definizioni di flessibilità (rigidificazione).

    OpenBabel, durante la conversione da PDB a PDBQT, calcola spesso le torsioni attive e struttura
    il file con un albero di flessibilità (tag `ROOT`, `BRANCH`). Tuttavia, per il docking standard
    (recettore rigido), AutoDock Vina richiede che il file del recettore non contenga questi tag.
    Questa funzione legge il file generato, filtra le righe problematiche e sovrascrive il file originale.

    Parameters
    ----------
    pdbqt_path : `str`
        Il percorso assoluto o relativo al file PDBQT del recettore da correggere.
        Il file viene modificato **in-place** (il contenuto originale viene sovrascritto).

    Returns
    -------
    None
        La funzione non restituisce valori. Stampa un messaggio di conferma al termine dell'operazione.

    Notes
    -----
    Le righe rimosse iniziano con le seguenti keyword:
    - ``ROOT``, ``ENDROOT``
    - ``BRANCH``, ``ENDBRANCH``
    - ``TORSDOF``
    - ``REMARK`` (solo quelli relativi a "active torsions")
    """
    pdbqt_path = Path(pdbqt_path)

    with pdbqt_path.open('r') as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        # Ignore tags that define flexibility (typical of ligands)
        if line.startswith('ROOT') or \
           line.startswith('ENDROOT') or \
           line.startswith('BRANCH') or \
           line.startswith('ENDBRANCH') or \
           line.startswith('TORSDOF'):
            continue
        
        # Ignore the REMARKs that list active twists (cleaning is optional but recommended)
        if line.startswith('REMARK') and ('active torsions' in line or 'between atoms' in line):
            continue

        new_lines.append(line)

    # Overwrites the file with the clean version
    with pdbqt_path.open('w') as f:
        f.writelines(new_lines)
    
    print(f"--- Fixed receptor file for Vina (flexibility tags removed): {pdbqt_path} ---")


if __name__ == '__main__':
    print("--- Vina Receptor Preparation Utility ---")
    
    # 1. Let's check the prerequisites (e.g. 7cam.pdb file)
    if not INPUT_PDB_FILE.exists():
        print_error(f"\n[INFO] Required: You must first download the PDB structure and save it as '{INPUT_PDB_FILE}'.", code = -1)
        
    # 2. We carry out the preparation
    success = prepare_receptor_pdbqt(INPUT_PDB_FILE, OUTPUT_PDBQT_FILE)
    fix_receptor_pdbqt(OUTPUT_PDBQT_FILE)
    
    if success:
        print("\nReady to run the Genetic Algorithm.")
        print(f"Next step: Run 'python main.py'")
