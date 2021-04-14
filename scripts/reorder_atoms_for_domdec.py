import sys


class Atom:
    def __init__(
        self, atom_str, idx, atom_name, resname, segname, resid, x, y, z, element
    ) -> None:
        self.atom_str = atom_str
        self.old_idx = idx
        self.new_idx = -1
        self.atom_name = atom_name
        self.resname = resname
        self.segname = segname
        self.resid = 1
        self.x = x
        self.y = y
        self.z = z
        self.element = element
        self.bonded_to = set()
        self.double_bonds = set()

    def set_bond_partners(self, bonded_to: set):
        if not self.bonded_to:
            self.bonded_to = bonded_to
        else:
            self.double_bonds = bonded_to


def read_pdb(file_path):
    """
    Read in pdb file, generate atom objects and save bonds
    """
    atoms = {}
    bonds = []

    tail = []
    head = []
    print("Reading in pdb ...")
    with open(file_path, "r") as f:
        stop = False
        for line in f.readlines():
            # read in atom records
            if line.startswith("HETATM"):
                # HETATM    1  C1  UNK A 900       3.707  -4.447  18.470  1.00  0.00           C
                stop = True
                (
                    _,
                    idx,
                    atom_name,
                    resname,
                    segname,
                    resid,
                    x,
                    y,
                    z,
                    _,
                    _,
                    element,
                ) = line.split()
                atoms[int(idx)] = Atom(
                    line, idx, atom_name, resname, segname, resid, x, y, z, element
                )

            # read in bonds
            elif line.startswith("CONECT"):
                s = line.split()[1:]
                idx = int(s[0])
                neighbors = [int(idx) for idx in s[1:]]
                atoms[idx].set_bond_partners(neighbors)
                bonds.append([int(idx) for idx in s])
            # save header
            elif stop == False:
                head.append(line)
            else:
                tail.append(line)

    return atoms, bonds, head, tail


def write_pdb_file(atoms, bonds, head, tail):
    """ write out pdb file """
    new_atom_idx = 0
    print("Reordering indices and writing out ...")

    with open("new.pdb", "w+") as f:
        for l in head:
            f.write(l)

        for atom_idx in sorted(atoms):
            heavy_atom = atoms[atom_idx]
            if heavy_atom.element != "H":
                new_atom_idx += 1
                heavy_atom.new_idx = new_atom_idx
                f.write(
                    f"HETATM {new_atom_idx:4d}  {heavy_atom.atom_name:3} {heavy_atom.resname:>3} {heavy_atom.segname} {heavy_atom.resid:>3}      {float(heavy_atom.x):6.3f} {float(heavy_atom.y):7.3f} {float(heavy_atom.z):7.3f}  1.00  0.00           {heavy_atom.element}\n"
                )
                # print(heavy_atom.atom_str)
                atom_idx += 1
                for idx in heavy_atom.bonded_to:
                    neighbor = atoms[idx]
                    if neighbor.element == "H":
                        new_atom_idx += 1
                        neighbor.new_idx = new_atom_idx
                        f.write(
                            f"HETATM {new_atom_idx:4d}  {neighbor.atom_name:3} {neighbor.resname:>3} {neighbor.segname} {neighbor.resid:>3}      {float(neighbor.x):6.3f} {float(neighbor.y):7.3f} {float(neighbor.z):7.3f}  1.00  0.00           {neighbor.element}\n"
                        )

                        # print(neighbor.atom_str)

        for b in bonds:
            # CONECT    1    2    4   15
            if len(b) == 2:
                b1, b2 = b
                b1 = atoms[b1].new_idx
                b2 = atoms[b2].new_idx
                f.write(f"CONECT {b1:4d} {b2:4d}\n")
            elif len(b) == 3:
                b1, b2, b3 = b
                b1 = atoms[b1].new_idx
                b2 = atoms[b2].new_idx
                b3 = atoms[b3].new_idx

                f.write(f"CONECT {b1:4d} {b2:4d} {b3:4d}\n")
            elif len(b) == 4:
                b1, b2, b3, b4 = b
                b1 = atoms[b1].new_idx
                b2 = atoms[b2].new_idx
                b3 = atoms[b3].new_idx
                b4 = atoms[b4].new_idx

                f.write(f"CONECT {b1:4d} {b2:4d} {b3:4d} {b4:4d}\n")
            elif len(b) == 5:
                b1, b2, b3, b4, b5 = b
                b1 = atoms[b1].new_idx
                b2 = atoms[b2].new_idx
                b3 = atoms[b3].new_idx
                b4 = atoms[b4].new_idx
                b5 = atoms[b5].new_idx

                f.write(f"CONECT {b1:4d} {b2:4d} {b3:4d} {b4:4d} {b5:4d}\n")
        for l in tail:
            f.write(l)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        raise RuntimeError("pdb file needs to be specified")
    # file_path = "BMI_original.pdb"
    file_path = sys.argv[1]
    print(file_path)
    atoms, bonds, head, tail = read_pdb(file_path)
    write_pdb_file(atoms, bonds, head, tail)