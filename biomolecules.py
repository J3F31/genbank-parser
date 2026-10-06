import warnings
from typing import List
import global_variables as gv


class Biomolecule:

    def __init__(self, molecule: str, condition: set[str], molecule_type: int):
        # Not doing anything with return bool value from validate molecule.
        # For now just allow incorrect molecules with a warning.
        self.__validate_molecule(molecule, condition)
        self.molecule = molecule
        self.molecule_type = molecule_type

    def __str__(self) -> str:
        return self.molecule

    # Makes value readable from lists of molecules.
    def __repr__(self) -> str:
        return self.molecule

    @property
    def name(self) -> str:
        match self.molecule_type:
            case 0:
                return gv.AMINOACID_DATA[self.molecule]['name']
            case 1:
                return gv.NUCLEOTIDE_DATA[self.molecule]['name']
            case _:
                return ''

    @property
    def weight(self) -> int:
        match self.molecule_type:
            case 0:
                return gv.AMINOACID_DATA[self.molecule]['weight']
            case 1:
                return gv.NUCLEOTIDE_DATA[self.molecule]['weight']
            case _:
                return 0

    def __validate_molecule(self, molecule: str, condition: set[str]) -> bool:
        """Condition for molecule validation."""
        if len(molecule) > 1:
            warnings.warn(f'Molecule {molecule} too long.')
            return False
        if molecule not in condition:
            warnings.warn(f'Molecule {molecule} not valid for {self.__class__.__name__}.')
            return False
        return True


class Nucleotide(Biomolecule):

    __condition: set[str] = set('atcgu')

    def __init__(self, nucleotide: str, molecule_type: int):
        super().__init__(nucleotide, self.__condition, molecule_type)


class Aminoacid(Biomolecule):

    __condition: set[str] = set('ARNDCQEGHILKMFPSTWYV')

    def __init__(self, aminoacid: str, molecule_type: int):
        super().__init__(aminoacid, self.__condition, molecule_type)


class Sequence:

    def __init__(self, sequence: str, sequence_type: int):
        self.__sequence = []
        self.sequence_type = sequence_type
        self.__init_sequence(sequence, self.sequence_type)

    def __str__(self) -> str:
        sequence = ''
        for item in self.__sequence:
            sequence += str(item)
        return sequence

    @property
    def sequence(self):
        sequence = ''
        for item in self.__sequence:
            sequence += str(item)
        return sequence

    @sequence.setter
    def sequence(self, sequence: str):
        self.__sequence = []
        self.__init_sequence(sequence, self.sequence_type)

    @property
    def weight(self) -> int:
        weight = 0
        for item in self.__sequence:
            weight += int(item.weight)
        return weight

    def __init_sequence(self, sequence: str, sequence_type: int):
        if sequence_type not in gv.MOLECULE_TYPES.values():
            warnings.warn('Sequence type not supported')

        match sequence_type:
            # Not sure why this doesn't work, temp literal int values.
            # case gv.SEQUENCE_TYPES['aminoacid']:
            case 0:
                self.__sequence = self.__build_aminoacid_sequence(sequence)
            # case gv.SEQUENCE_TYPES['nucleotide']:
            case 1:
                self.__sequence = self.__build_nucleotide_sequence(sequence)
            case _:
                warnings.warn('Sequence type not supported')

    def __build_nucleotide_sequence(self, sequence: str) -> List[Nucleotide]:
        nucleotide_list: List[Nucleotide] = []
        for item in sequence:
            nucleotide_list.append(Nucleotide(item, self.sequence_type))
        return nucleotide_list

    def __build_aminoacid_sequence(self, sequence: str) -> List[Aminoacid]:
        aminoacid_list: List[Aminoacid] = []
        for item in sequence:
            aminoacid_list.append(Aminoacid(item, self.sequence_type))
        return aminoacid_list
