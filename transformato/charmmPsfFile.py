from parmed.utils import tag_molecules
from parmed.constants import DEG_TO_RAD
from parmed.periodic_table import AtomicNum, element_by_mass
from parmed.topologyobjects import (
    Bond,
    Angle,
    Dihedral,
    Improper,
    AcceptorDonor,
    Group,
    Cmap,
    UreyBradley,
    NoUreyBradley,
    Atom,
    DihedralType,
    ImproperType,
    UnassignedAtomType,
    ExtraPoint,
    DrudeAtom,
    LocalCoordinatesFrame,
    DrudeAnisotropy,
)
from parmed.exceptions import CharmmError
from parmed.structure import Structure
from parmed.utils.io import genopen
from parmed.charmm.psf import _ZeroDict, _FileEOF, _resre

import parmed as pm

# from parmed.charmm.psf import _catchindexerror
@pm.charmm.psf._catchindexerror
class CustomCharmmPsfFile(pm.charmm.CharmmPsfFile):
    def __init__(self, psf_name=None):
        """
        Opens and parses a PSF file, then instantiates a CharmmPsfFile
        instance from the data.
        """
        Structure.__init__(self)
        # Bail out if we don't have a filename
        if psf_name is None:
            return
        # Open the PSF and read the first line. It must start with "PSF"
        if isinstance(psf_name, str):
            fileobj = genopen(psf_name, "r")
            own_handle = True
        else:
            fileobj = psf_name
            own_handle = False
        try:
            self.name = psf_name if isinstance(psf_name, str) else ""
            line = fileobj.readline()
            if not line.startswith("PSF"):
                raise CharmmError(
                    f"Unrecognized PSF file. First line is {line.strip()}"
                )
            # Store the flags
            psf_flags = line.split()[1:]
            # Now get all of the sections and store them in a dict
            fileobj.readline()
            # Now get all of the sections
            psfsections = _ZeroDict()
            while True:
                try:
                    sec, ptr, data = self._parse_psf_section(fileobj)
                except _FileEOF:
                    break
                psfsections[sec] = (ptr, data)
            # store the title
            self.title = psfsections["NTITLE"][1]
            # Next is the number of atoms
            natom = self._convert(psfsections["NATOM"][0], int, "natom")
            # Parse all of the atoms
            drude_alpha_thole = []
            is_drude = "DRUDE" in psf_flags
            for i in range(natom):
                words = psfsections["NATOM"][1][i].split()
                atid = int(words[0])
                if atid != i + 1:
                    raise CharmmError("Nonsequential atoms detected!")
                segid = words[1]
                rematch = _resre.match(words[2])
                if not rematch:
                    raise CharmmError(f"Could not interpret residue number {words[2]}")
                resid, inscode = rematch.groups()
                resid = self._convert(resid, int, "residue number")
                resname = words[3]
                name = words[4]
                attype = words[5]
                # Try to convert the atom type to an integer a la CHARMM
                try:
                    attype = int(attype)
                except ValueError:
                    pass
                charge = self._convert(words[6], float, "partial charge")
                mass = self._convert(words[7], float, "atomic mass")
                props = words[8:]
                if is_drude:
                    alpha = self._convert(words[9], float, "alpha")
                    thole = self._convert(words[10], float, "thole")
                    drude_alpha_thole.append((alpha, thole))
                if is_drude and i >= 1 and drude_alpha_thole[-2] != (0, 0):
                    # This assumes that the Drude atom is defined immediately after its parent atom.
                    # This always appears to be the case, but if it proves untrue at some point then
                    # this will need to be refactored to identify Drude atoms after identifying bonds.
                    my_alpha, my_thole = drude_alpha_thole[-2]
                    atom = DrudeAtom(
                        name=name,
                        type=attype,
                        charge=charge,
                        mass=mass,
                        parent=self.atoms[-1],
                        atomic_number=0,
                        alpha=my_alpha,
                        thole=my_thole,
                        drude_type=attype,
                    )
                elif (
                    name.startswith("LP")
                    or (isinstance(attype, str) and attype.startswith("LP"))
                ) and mass == 0:
                    atom = ExtraPoint(
                        name=name,
                        type=attype,
                        charge=charge,
                        mass=mass,
                        atomic_number=0,
                    )
                else:
                    atom = Atom(
                        name=name,
                        type=attype,
                        charge=charge,
                        mass=mass,
                        atomic_number=AtomicNum[element_by_mass(mass)],
                    )
                atom.props = props
                self.add_atom(
                    atom, resname, resid, chain=segid, inscode=inscode, segid=segid
                )
            # Now get the number of bonds
            nbond = self._convert(psfsections["NBOND"][0], int, "number of bonds")
            if len(psfsections["NBOND"][1]) != nbond * 2:
                raise CharmmError(
                    f"Got {len(psfsections['NBOND'][1])} indexes for {nbond} bonds"
                )
            it = iter(psfsections["NBOND"][1])
            for i, j in zip(it, it):
                self.bonds.append(Bond(self.atoms[i - 1], self.atoms[j - 1]))
            # Now get the number of angles and the angle list
            ntheta = self._convert(psfsections["NTHETA"][0], int, "number of angles")
            if len(psfsections["NTHETA"][1]) != ntheta * 3:
                raise CharmmError(
                    f"Got {len(psfsections['NTHETA'][1])} indexes for {ntheta} angles"
                )
            it = iter(psfsections["NTHETA"][1])
            for i, j, k in zip(it, it, it):
                self.angles.append(
                    Angle(self.atoms[i - 1], self.atoms[j - 1], self.atoms[k - 1])
                )
                self.angles[-1].funct = 5  # urey-bradley
            # Now get the number of torsions and the torsion list
            nphi = self._convert(psfsections["NPHI"][0], int, "number of torsions")
            if len(psfsections["NPHI"][1]) != nphi * 4:
                raise CharmmError(
                    f"Got {len(psfsections['NPHI'])} indexes for {nphi} torsions"
                )
            it = iter(psfsections["NPHI"][1])
            for i, j, k, l in zip(it, it, it, it):
                self.dihedrals.append(
                    Dihedral(
                        self.atoms[i - 1],
                        self.atoms[j - 1],
                        self.atoms[k - 1],
                        self.atoms[l - 1],
                    )
                )
            self.dihedrals.split = False
            # Now get the number of improper torsions
            nimphi = self._convert(psfsections["NIMPHI"][0], int, "number of impropers")
            if len(psfsections["NIMPHI"][1]) != nimphi * 4:
                raise CharmmError(
                    f"Got {len(psfsections['NIMPHI'][1])} indexes for {nimphi} impropers"
                )
            it = iter(psfsections["NIMPHI"][1])
            for i, j, k, l in zip(it, it, it, it):
                self.impropers.append(
                    Improper(
                        self.atoms[i - 1],
                        self.atoms[j - 1],
                        self.atoms[k - 1],
                        self.atoms[l - 1],
                    )
                )
            # Now handle the donors (what is this used for??)
            ndon = self._convert(psfsections["NDON"][0], int, "number of donors")
            if len(psfsections["NDON"][1]) != ndon * 2:
                raise CharmmError(
                    f"Got {len(psfsections['NDON'][1])} indexes for {ndon} donors"
                )
            it = iter(psfsections["NDON"][1])
            for i, j in zip(it, it):
                self.donors.append(AcceptorDonor(self.atoms[i - 1], self.atoms[j - 1]))
            # Now handle the acceptors (what is this used for??)
            nacc = self._convert(psfsections["NACC"][0], int, "number of acceptors")
            if len(psfsections["NACC"][1]) != nacc * 2:
                raise CharmmError(
                    f"Got {len(psfsections['NACC'][1])} indexes for {nacc} acceptors"
                )
            it = iter(psfsections["NACC"][1])
            for i, j in zip(it, it):
                self.acceptors.append(
                    AcceptorDonor(self.atoms[i - 1], self.atoms[j - 1])
                )
            # Now get the group sections
            try:
                ngrp, nst2 = psfsections["NGRP NST2"][0]
            except ValueError:  # pragma: no cover
                raise CharmmError("Could not unpack GROUP pointers")  # pragma: no cover
            tmp = psfsections["NGRP NST2"][1]
            self.groups.nst2 = nst2
            # Now handle the groups
            if len(psfsections["NGRP NST2"][1]) != ngrp * 3:
                raise CharmmError(f"Got {len(tmp)} indexes for {ngrp} groups")
            it = iter(psfsections["NGRP NST2"][1])
            for i, j, k in zip(it, it, it):
                self.groups.append(Group(self.atoms[i], j, k))
            # Assign all of the atoms to molecules recursively
            tmp = psfsections["MOLNT"][1]
            tag_molecules(self)
            if len(psfsections["NUMLP NUMLPH"][1]) != 0:
                # We have a CHARMM PSF file; now do NUMLP/NUMLPH sections
                self._process_lonepair_section(psfsections["NUMLP NUMLPH"])
            # Now process the NUMANISO records if this is a drude PSF
            if is_drude:
                self._process_aniso_section(psfsections["NUMANISO"])
            # Now do the CMAPs
            ncrterm = self._convert(
                psfsections["NCRTERM"][0], int, "Number of cross-terms"
            )
            if len(psfsections["NCRTERM"][1]) != ncrterm * 8:
                raise CharmmError(
                    f"Got {len(psfsections['NCRTERM'])} CMAP indexes for {ncrterm} cmap terms"
                )
            it = iter(psfsections["NCRTERM"][1])
            for i, j, k, l, m, n, o, p in zip(it, it, it, it, it, it, it, it):
                self.cmaps.append(
                    Cmap.extended(
                        self.atoms[i - 1],
                        self.atoms[j - 1],
                        self.atoms[k - 1],
                        self.atoms[l - 1],
                        self.atoms[m - 1],
                        self.atoms[n - 1],
                        self.atoms[o - 1],
                        self.atoms[p - 1],
                    )
                )
            self.unchange()
            self.flags = psf_flags
        finally:
            if own_handle:
                fileobj.close()
