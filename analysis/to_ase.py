import MDAnalysis as mda
import numpy as np
from ase import Atoms
from tqdm.auto import tqdm
from MDAnalysis.lib.mdamath import triclinic_vectors   # for the cell matrix
import ase.io

from shiftml.ase import ShiftML

# 2-A.  Load topology + trajectory --------------------
u = mda.Universe("../data/zenodo_1468560/prod1.tpr",              # topology (atom order + bonds)
                 "../data/zenodo_1468560/prod1.trr",               # coordinates / velocities / box
                 convert_units=True)       # MDAnalysis default: Å & ps

# 2-B.  Prepare a list with the chemical symbols once
symbols = [a.element.title() if a.element
           else a.name[0].upper()          # crude fallback: first letter of atom name
           for a in u.atoms]


model = ShiftML("ShiftML3", device="cuda")
print(len(u.trajectory), "frames in MDAnalysis trajectory.")


def ts_to_ase(ts):
    """Convert a single MDAnalysis timestep to an ASE Atoms object."""
    
    pos  = ts.positions.copy()
    cell = triclinic_vectors(ts.dimensions)
    frame = Atoms(symbols, positions=pos, cell=cell, pbc=True)
    
    if ts.has_velocities:
        frame.set_velocities(ts.velocities / 1000.0)  # ps → fs
    
    return frame

cs_isos = []

for n, ts in tqdm(enumerate(u.trajectory[::5])):                    # iterate over every frame
    # positions are already in Å if convert_units=True
    if n == 0:
        ase.io.write("first_frame.xyz", ts_to_ase(ts))  # write first frame to file
    
    frame = ts_to_ase(ts)  # convert MDAnalysis timestep to ASE Atoms object
    
    Yiso = model.get_cs_iso_ensemble(frame) # shape (N_atoms, N_ensemble)
    cs_isos.append(Yiso)

cs_isos = np.array(cs_isos)
np.save("cs_isos.npy", cs_isos)
print("CS isos:", cs_isos.shape)
