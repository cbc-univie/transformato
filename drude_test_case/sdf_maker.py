import pandas as pd
import os
from openbabel import openbabel

#where the data directory is located
path="."


df = pd.read_csv('molecules.csv', sep=';', header=None)  

for arg in df[0]:
    print(f"Processing: {arg}")
    mol=arg.lower()
    final_path=os.path.join("data",mol,"waterbox",mol)
    print(final_path)
    os.chdir(final_path)

    with open("solu.pdb", "r") as infile, open("filtered.pdb", "w") as outfile:
        for line in infile:
            split = line.strip().split()
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                #atom_name = split[2]
                if atom_name.startswith("D") or atom_name.startswith("LP"):
                    continue
            outfile.write(line)

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "sdf")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, "filtered.pdb")

    print(mol.NumAtoms())
    print(mol.NumBonds())
    print(mol.NumResidues())

    obConversion.WriteFile(mol, 'solu.sdf')
