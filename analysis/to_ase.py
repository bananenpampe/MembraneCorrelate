import MDAnalysis as mda
import numpy as np
from ase import Atoms
from MDAnalysis.lib.mdamath import triclinic_vectors   # for the cell matrix


from shiftml.ase import ShiftML

# 2-A.  Load topology + trajectory --------------------
u = mda.Universe("../data/zenodo_1468560/prod1.tpr",              # topology (atom order + bonds)
                 "../data/zenodo_1468560/prod1.trr",               # coordinates / velocities / box
                 convert_units=True)       # MDAnalysis default: Å & ps

# 2-B.  Prepare a list with the chemical symbols once
symbols = [a.element.title() if a.element
           else a.name[0].upper()          # crude fallback: first letter of atom name
           for a in u.atoms]

ase_frames = []

print(len(u.trajectory), "frames in MDAnalysis trajectory.")

for ts in u.trajectory[:2]:                    # iterate over every frame
    # positions are already in Å if convert_units=True
    pos  = ts.positions.copy()

    # build the (3×3) cell matrix from GROMACS box lengths & angles
    cell = triclinic_vectors(ts.dimensions)   # also Å

    # assemble an ase.Atoms for this frame
    frame = Atoms(symbols, positions=pos, cell=cell, pbc=True)

    # optional: attach velocities (Å/fs).  MDAnalysis keeps them in Å/ps
    if ts.has_velocities:
        frame.set_velocities(ts.velocities / 1000.0)  # ps → fs

    ase_frames.append(frame)

print(len(ase_frames), "frames loaded from MDAnalysis trajectory.")
print(len(ase_frames[0]))
print(ase_frames[0])
