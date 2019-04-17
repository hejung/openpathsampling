from openpathsampling.netcdfplus import StorableNamedObject
from mdtraj.core.element import Element


def make_martini_bead_elements():
    """
    Register all possible MARTINI beads as Elements with mdtraj.

    NOTE: This needs to be called **before** loading mdtraj topologies
    that contain MARTINI beads.
    Otherwise openpathsampling returns None for the topologies.
    """
    for number, symbol, mass, radius in martini_beads:
        # do not care for return value, mdtraj keeps its own element registry
        name = symbol + ' MARTINI bead'
        symbol = symbol + ' MB'
        _ = Element(number, name, symbol, mass, radius)


# NOTE:
# -> ordering as in the martini v2.2 '.itp' file
# -> we start with atomic number 1000 to make sure we have no collisions with
#    'real' elements
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
    def __init__(self, mdtraj_topology=None, mdanalysis_topology=None):
        if mdtraj_topology is None and mdanalysis_topology is not None:
            mdtraj_topology = mdanalysis_to_mdtraj(mdanalysis_topology)
        elif mdanalysis_topology is None and mdtraj_topology is not None:
            mdanalysis_topology = mdtraj_to_mdanalysis(mdtraj_topology)
        else:
            raise ValueError('Either mdtraj_topology or mdanalysis_topology '
                             + 'need to be given.')

        super(MDTrajTopology, self).__init__(int(mdtraj_topology.n_atoms), 3)
        self.mdtraj = mdtraj_topology
        self.mdanalysis = mdanalysis_topology

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

        return {'mdtraj': out}

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

            return cls(md_topology)