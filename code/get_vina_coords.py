# get_vina_coords.py
import numpy as np

from pymol import cmd


from constants import INPUT_PDB_FILE

LIGAND_CODE = "3TL"          # Your ligand residue code
OBJECT_NAME = "2P3D"         # Name that PyMOL assigns by default when loading 2P3D.pdb

LIGAND_SELECTION = f"resn {LIGAND_CODE} and {OBJECT_NAME}"

# Function to calculate the geometric center (centroid) of a selection
def get_centroid(selection: str) -> np.ndarray:
    """
    Calculates the geometric center (centroid) of an atomic selection in PyMOL.

    Queries the PyMOL API to retrieve atomic coordinates of the specified selection
    and computes the arithmetic mean for X, Y, and Z axes. Used to define the 
    center of the docking search box.

    Parameters
    ----------
    selection : `str`
        The PyMOL selection string (e.g., 'resn LIG').

    Returns
    -------
    `np.ndarray`
        A NumPy array of shape `(3,)` containing the [x, y, z] coordinates 
        of the centroid in Ångstroms.

    Raises
    ------
    ValueError
        If the selection is empty or does not exist in the PyMOL session.

    Examples
    --------
    >>> center = get_centroid("resn 3TL")
    >>> print(center)
    [10.5, -4.2, 33.1]
    """
    model = cmd.get_model(selection)
    if not model.atom:
        raise ValueError(f"Empty selection: {selection}. Check the ligand ID.")
        
    coords = np.array([a.coord for a in model.atom])
    return coords.mean(axis = 0)

if __name__ == '__main__':
    if not INPUT_PDB_FILE.exists():
        print(f"ERROR: File {INPUT_PDB_FILE.resolve()} not found. Download it and place it here.")
    else:
        try:
            # Load the PDB structure
            cmd.load(str(INPUT_PDB_FILE))
            
            # Calculate the center of the co-crystallized ligand
            center_coords = get_centroid(LIGAND_SELECTION)
            
            print("\n=======================================================")
            print("Vina Coordinates (Pocket Center)")
            print("=======================================================")
            print(f"File analyzed: {INPUT_PDB_FILE}")
            print(f"Pocket Center (X, Y, Z):")
            print(f"X: {center_coords[0]:.3f}, Y: {center_coords[1]:.3f}, Z: {center_coords[2]:.3f}")
            
            print("\n-------------------------------------------------------")
            print(f"CENTER_X, CENTER_Y, CENTER_Z = {center_coords[0]:.3f}, {center_coords[1]:.3f}, {center_coords[2]:.3f}")
            print("-------------------------------------------------------")
            
            # The Box Size MUST be estimated visually
            print("\nBOX SIZE (Dimension):")
            print("The box must be large enough to contain the peptide (10AA).")
            print("A box of at least 20x20x20 or 25x25x25 Ångstroms is recommended.")
            print("Example: SIZE_X, SIZE_Y, SIZE_Z       = 25, 25, 25")
            
            # Clean up the PyMOL environment
            cmd.delete("all")

        except Exception as e:
            print(f"Error while parsing PyMOL: {e}")
            print("Make sure PyMOL is installed correctly and that the LIGAND_SELECTION selection is correct.")
