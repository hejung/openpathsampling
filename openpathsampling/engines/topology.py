import string
import numpy as np
import pandas as pd
import mdtraj as md
import MDAnalysis as mda
from openpathsampling.netcdfplus import StorableNamedObject

import logging
logger = logging.getLogger(__name__)


def make_martini_bead_elements():
    """
    Register all possible MARTINI beads as Elements with mdtraj.

    NOTE: This needs to be called **before** loading mdtraj topologies
    that contain MARTINI beads.
    Otherwise the openpathsampling.Storage returns None for the topologies.
    """
    for number, symbol, mass, radius in martini_beads:
        # do not care for return value, mdtraj keeps its own element registry
        name = symbol + ' MARTINI bead'
        symbol = symbol + ' MB'
        _ = md.core.element.Element(number, name, symbol, mass, radius)


# NOTE:
# -> ordering as in the martini v2.2 '.itp' file
# -> we start with atomic number 1000 to make sure we have no collisions with
#    real elements
# -> element symbol is the shorthand as used in the MARTINI .itp files with the
#    additional suffix ' MB' (for MARTINI bead) to be able to differentiate
#    between real elements and the MARTINI beads
#    ('Na' and 'Nd' are ambigous otherwise)
# -> element name is the short hand for the bead suffixed by ' MARTINI bead',
#    e.g. `P5 MARTINI bead' or 'Qd MARTINI bead'
# -> element mass is correct, i.e. 72 for normal beads and 45 for rings
#    for now radius is 0.0 because it depends on the interacting bead type
#    TODO: we could use self interaction values for LJ as radius?
martini_beads = [
                 # tuple is (number, symbol, mass, radius)
                 # note that we can construct the name from the symbol
                 # STANDARD types, 4:1 mapping
                 # polar type
                 (1000, 'P5', 72.0, 0.0),
                 (1001, 'P4', 72.0, 0.0),
                 (1002, 'P3', 72.0, 0.0),
                 (1003, 'P2', 72.0, 0.0),
                 (1004, 'P1', 72.0, 0.0),
                 # intermediate polar
                 (1005, 'Nda', 72.0, 0.0),
                 (1006, 'Nd', 72.0, 0.0),  # this could be neodym!
                 (1007, 'Na', 72.0, 0.0),  # this could be sodium!
                 (1008, 'N0', 72.0, 0.0),
                 # apolar
                 (1009, 'C5', 72.0, 0.0),
                 (1010, 'C4', 72.0, 0.0),
                 (1011, 'C3', 72.0, 0.0),
                 (1012, 'C2', 72.0, 0.0),
                 (1013, 'C1', 72.0, 0.0),
                 # charged
                 (1014, 'Qda', 72.0, 0.0),
                 (1015, 'Qd', 72.0, 0.0),
                 (1016, 'Qa', 72.0, 0.0),
                 (1017, 'Q0', 72.0, 0.0),
                 # RING types, 2-3:1 mapping
                 (1018, 'SP5', 45.0, 0.0),
                 (1019, 'SP4', 45.0, 0.0),
                 (1020, 'SP3', 45.0, 0.0),
                 (1021, 'SP2', 45.0, 0.0),
                 (1022, 'SP1', 45.0, 0.0),
                 (1023, 'SNda', 45.0, 0.0),
                 (1024, 'SNd', 45.0, 0.0),
                 (1025, 'SNa', 45.0, 0.0),
                 (1026, 'SN0', 45.0, 0.0),
                 (1027, 'SC5', 45.0, 0.0),
                 (1028, 'SC4', 45.0, 0.0),
                 (1029, 'SC3', 45.0, 0.0),
                 (1030, 'SC2', 45.0, 0.0),
                 (1031, 'SC1', 45.0, 0.0),
                 (1032, 'SQda', 45.0, 0.0),
                 (1033, 'SQd', 45.0, 0.0),
                 (1034, 'SQa', 45.0, 0.0),
                 (1035, 'SQ0', 45.0, 0.0),
                 # AMINO ACIDS (for Q-C interactions inside proteins)
                 (1036, 'AC2', 72.0, 0.0),
                 (1037, 'AC1', 72.0, 0.0),
                 # BIG particle type (prevent freezing of water)
                 (1038, 'BP4', 72.0, 0.0)
                ]


class Topology(StorableNamedObject):
    """
    Topology is the object that contains all information about the structure
    of the system to be simulated.

    Attributes
    ----------
    n_atoms : int
        number of atoms
    n_spatial : int
        number of spatial dimensions, default is 3
    """

    def __init__(self, n_atoms, n_spatial=3):
        super(Topology, self).__init__()
        self.n_atoms = n_atoms
        self.n_spatial = n_spatial


class MDTopology(Topology):
    def __init__(self, mdtraj_topology=None, mdanalysis_topology=None, martini=False):
        if mdtraj_topology is None and mdanalysis_topology is not None:
            mdtraj_topology = self.mdanalysis_to_mdtraj(mdanalysis_topology, martini=martini)
            original = 'MDAnalysis'
        elif mdanalysis_topology is None and mdtraj_topology is not None:
            mdanalysis_topology = self.mdtraj_to_mdanalysis(mdtraj_topology)
            original = 'MDTraj'
        else:
            raise ValueError('Either mdtraj_topology or mdanalysis_topology '
                             + 'need to be given.')

        super(MDTopology, self).__init__(int(mdtraj_topology.n_atoms), 3)
        self.mdtraj = mdtraj_topology
        self.mdanalysis = mdanalysis_topology
        self._original = original
        self.martini = martini

    def to_dict(self):
        out = dict()

        atom_data = []
        for atom in self.mdtraj.atoms:
            if atom.element is None:
                element_symbol = ""
            else:
                element_symbol = atom.element.symbol

            if hasattr(atom, 'segment_id'):
                atom_data.append((
                    atom.serial, atom.name, element_symbol,
                    int(atom.residue.resSeq), atom.residue.name,
                    atom.residue.chain.index, atom.segment_id))
            else:
                atom_data.append((
                    atom.serial, atom.name, element_symbol,
                    int(atom.residue.resSeq), atom.residue.name,
                    atom.residue.chain.index, ''))

        out['atom_columns'] = ["serial", "name", "element", "resSeq",
                               "resName", "chainID", "segmentID"]
        out['atoms'] = atom_data
        out['bonds'] = [(a.index, b.index) for (a, b) in self.mdtraj.bonds]

        return {'mdtraj': out,
                'original': self._original,
                'martini': self.martini
                }

    @classmethod
    def from_dict(cls, dct):
        top_dict = dct['mdtraj']

        atoms = pd.DataFrame(
            top_dict['atoms'], columns=top_dict['atom_columns'])
        bonds = np.array(top_dict['bonds'])

        try:
            md_topology = md.Topology.from_dataframe(atoms, bonds)
            return cls(md_topology)
        except Exception:
            # we try a fix and add multiples of 10000 to the resSeq

            logger.info('Normal reconstruction of topology failed. '
                        'Trying a fix to the 10k residue ID problem.')

            for ci in np.unique(atoms['chainID']):
                chain_atoms = atoms[atoms['chainID'] == ci]
                indices = chain_atoms.index.tolist()

                old_residue_id = 0
                multiplier = 0
                places = []
                for row, res_id in zip(indices, list(chain_atoms['resSeq'])):
                    if res_id < old_residue_id:
                        if multiplier > 0:
                            atoms.loc[places, 'resSeq'] += 10000 * multiplier

                        places = []
                        multiplier += 1

                    if multiplier > 0:
                        places.append(row)

                    old_residue_id = res_id

                if multiplier > 0:
                    atoms.loc[places, 'resSeq'] += 10000 * multiplier

            # this function is really slow! Reads ~ 1000 atoms per second
            md_topology = md.Topology.from_dataframe(atoms, bonds)

            # that we have successfully created the topology using from_df
            # we remove the wrong multipliers
            # this is weird, but reproduces the current behaviour

            for atom in md_topology.atoms:
                atom.residue.resSeq %= 10000

            ret = cls(mdtraj_topology=md_topology, martini=top_dict['martini'])
            ret._original = top_dict['original']
            return ret

    @classmethod
    def mdanalysis_to_mdtraj(cls, mda_topol, standard_replacements=True,
                             martini=False):
        """
        Convert a MDAnalysis topology into a MDTraj topology.

        Parameters:
        -----------
        mda_topol - :class: MDAnalysis.core.topology.Topology,
                    the mdtraj topology to convert
        standard_replacements - bool, whether to try to replace
                                known residue and atom names,
                                default is True
        martini - bool, whether the MDAnalysis.Topology is a
                  MARTINI topology, i.e. contains MARTINI beads,
                  ensures that the correct 'elements' are used,
                  default is False

        """
        PDBTrajectoryFile = md.formats.PDBTrajectoryFile
        Element = md.core.element.Element

        if martini:
            try:
                make_martini_bead_elements()
            except AssertionError:
                # MARTINI beads already registered,
                # mdtraj asserts there is only be one of each 'element'
                pass

        # empty mdtraj topology to populate
        mdt_topol = md.core.topology.Topology()
        if standard_replacements:
            # load the name replacement tables to let mdtraj guess elements
            PDBTrajectoryFile._loadNameReplacementTables()
        # make an 'empty' universe to have easier access to the autogenerated properties
        mda_u = mda.Universe(mda_topol)
        # loop structure orientend @ how mdtraj builds its topology from `.pdb`
        # NOTE: MDAnalysis uses segments as chains...
        atoms_by_idx = {}
        for orig_seg in mda_u.segments:
            c = mdt_topol.add_chain()
            for orig_res in orig_seg.residues:
                resName = orig_res.resname
                resNum = orig_res.resnum
                resSegID = orig_res.segid
                # the residue name replacement can also be useful for MARTINI
                if standard_replacements and resName in PDBTrajectoryFile._residueNameReplacements:
                    resName = PDBTrajectoryFile._residueNameReplacements[resName]
                r = mdt_topol.add_residue(resName, c, resNum, resSegID)
                # the atom name replacements will only work reliably for atomistic systems
                if (standard_replacements and not martini
                    and resName in PDBTrajectoryFile._atomNameReplacements
                    ):
                    atomReplace = PDBTrajectoryFile._atomNameReplacements[resName]
                else:
                    atomReplace = {}
                for orig_at in orig_res.atoms:
                    atName = orig_at.name
                    # this again only works for atomistic systems
                    if standard_replacements:
                        if atName in atomReplace:
                            atName = atomReplace[atName]
                        # let mdtraj guess the element
                        element = PDBTrajectoryFile._guess_element(atName, resName, len(orig_res.atoms))
                    if martini:
                        # get the MARTINI element
                        atType = orig_at.type
                        # make sure we do not add the suffix to atom types that already have it,
                        # i.e. make sure we can convert back and forth from mdtraj to mdanalysis without changing atomTypes
                        if atType.endswith('MB'):
                            element = Element.getBySymbol(atType)
                        else:
                            element = Element.getBySymbol(atType + ' MB')
                    at = mdt_topol.add_atom(atName, element, r, serial=orig_at.index)
                    atoms_by_idx[orig_at.index] = at
        # make bonds
        # NOTE: I think the most sensible thing is to not create any additional bonds
        # as this makes sure the mdtraj topology is equivalent to the mdanalysis one
        # (at least in terms of connectivity as we allow mdtraj to substitute known residueNames and atomNames)
        #mdt_topol.create_standard_bonds()
        #mdt_topol.create_disulfide_bonds()
        for bond in mda_u.bonds:
            # MDAnalysis always gives the lower index atom first as does mdtraj, no need to check for other ordering
            at1 = atoms_by_idx[bond.indices[0]]
            at2 = atoms_by_idx[bond.indices[1]]
            mdt_topol.add_bond(at1, at2)

        return mdt_topol

    @classmethod
    def mdtraj_to_mdanalysis(cls, mdt_topol):
        """
        Convert a MDTraj topology to a MDAnalysis topology.

        Parameters:
        -----------
        mdt_topol - :class:mdtraj.core.topology.Topology,
                    the mdtraj topology to convert

        """
        change_squash = mda.topology.base.change_squash
        mda_attrs = mda.core.topologyattrs
        # oriented @ MDAnalysis PDBParser
        serials = []
        names = []
        chainids = []
        segids = []
        atomtypes = []
        masses = []
        resids = []
        resnames = []
        # NOTE: mdanalysis uses segments instead of chains,
        # but mdtraj has names only for segments
        for i, c in enumerate(mdt_topol.chains):
            for at in c.atoms:
                serials.append(at.serial)
                names.append(at.name)
                # convert integer to letter to have a `name`
                chainids.append(string.ascii_uppercase[i])
                segids.append(at.segment_id)
                atomtypes.append(at.element.symbol)
                masses.append(at.element.mass)
                resnames.append(at.residue.name)
                resids.append(at.residue.resSeq)

        if not any(segids):
            # if all segids are empty strings we take the chainID instead
            segids, chainids = chainids, None

        attrs = []
        # atom attributes
        for vals, Attr, dtype in (
                    (names, mda_attrs.Atomnames, object),
                    (chainids, mda_attrs.ChainIDs, object),
                    (serials, mda_attrs.Atomids, np.int32),
                    (masses, mda_attrs.Masses, np.float64),
                    (atomtypes, mda_attrs.Atomtypes, object),
                                  ):
                if vals is not None:
                    if not any(vals):
                        # we dont have any values so leave out the attribute
                        continue
                    attrs.append(Attr(np.array(vals, dtype=dtype)))

        # residue attributes
        resids = np.array(resids, dtype=np.int32)
        resnames = np.array(resnames, dtype=object)
        resnums = resids.copy()
        segids = np.array(segids, dtype=object)
        # squash arrays to one entry per residue
        residx, (resids, resnames, resnums, segids) = change_squash((resids, resnames, segids), (resids, resnames, resnums, segids))
        attrs.append(mda_attrs.Resnums(resnums))
        attrs.append(mda_attrs.Resids(resids))
        attrs.append(mda_attrs.Resnames(resnames))
        # squash array to one entry per segment/chain
        segidx, (segids,) = change_squash((segids,), (segids,))
        attrs.append(mda_attrs.Segids(segids))

        # bonds
        mapp = {s: i for i, s in enumerate(serials)}  # mapp serials to indices
        bonds = []
        orders = []
        # TODO: we ignore the bond type for now, the reason is that mdtraj uses Singelton classes like 'Amide'
        # and MDAnalysis uses tuples with the bonding atoms elements
        #types = []
        for bond in mdt_topol.bonds:
            # NOTE: from v1.9.1 to v1.9.2 mdtraj changed the bonds from tuples to namedtuples
            # getting the atoms through indices always works, i.e. for tuples and namedtuples
            bonds.append((mapp[bond[0].serial], mapp[bond[1].serial]))
            # but only the namedtuple knows about bond order
            try:
                o = bond.order
            except AttributeError:
                # mdtraj version < 1.9.2
                o = None
            orders.append(int(o) if o is not None else None)
            #types.append(bond.type)
        attrs.append(mda_attrs.Bonds(bonds, order=orders, types=None))

        # finally create the topology
        mda_topol = mda.core.topology.Topology(n_atoms=mdt_topol.n_atoms,
                                               n_res=mdt_topol.n_residues,
                                               n_seg=mdt_topol.n_chains,
                                               attrs=attrs,
                                               atom_resindex=residx,
                                               residue_segindex=segidx
                                               )
        return mda_topol
