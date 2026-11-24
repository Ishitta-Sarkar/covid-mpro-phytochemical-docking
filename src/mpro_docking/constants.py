"""Protocol constants for the published 6LU7 run.

These values are part of the method, not outputs. The search box is not a
fixed 22 Å cube: ``cubic_search_box`` sets the edge from the co-crystal
ligand plus ``BOX_PADDING`` Å, and it is at least ``MIN_BOX_SIZE`` Å.
"""

PDB_ID = "6LU7"
RECEPTOR_CHAIN = "A"
LIGAND_CHAIN = "C"
COVALENT_PROTEIN_RESIDUE = ("CYS", 145)

EXHAUSTIVENESS = 16
NUM_MODES = 9
ENERGY_RANGE = 3
MIN_MODE_RMSD = 1
SEED = 42
CPU = 4

BOX_PADDING = 4.0
MIN_BOX_SIZE = 22.0
CENTER_DECIMALS = 3

# Baker and Hubbard, Prog. Biophys. Mol. Biol. 1984, 44, 97-179.
HBOND_MAX_HA = 2.5
HBOND_MIN_ANGLE = 120.0
CONTACT_CUTOFF = 4.0

PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
PUBCHEM_3D_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
PUBCHEM_2D_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=2d"
USER_AGENT = "mpro-docking/1.0 (structure fetch for a published docking workflow)"
