"""Genbank class"""
import sys
import re
import warnings
import global_variables as gv
from biomolecules import Sequence


class GenBank:
    """Genbank parser"""

    def __init__(self, filename):
        self.__data = {
            'identifiers': {
                'locus_name': '',
                'accession_number': '',
                'version_number': '',
                'gen_info_number': ''
            },
            'molecule_type': '',
            'genbank_division': '',
            'modification_date': '',
            'definition': '',
            'keywords': [],
            'source_organism': {
                'common_name': '',
                'taxonomic_name': '',
                'full_phylogenetic_name': '',
            },
            'origin': {
                'sequence': '',
                'length': 0
            },
            'features': {
                'CDS': []
            }
        }
        self.__filename = filename

        with (open(filename, 'r') as file):
            # Store readlines object (otherwise can't access within iteration)
            self.__file_lines = file.readlines()
            self.__parse_data()

    def __getitem__(self, index):
        return self.__data['features']['CDS'][index]

    def __gb_file_generator(self):
        for index, line in enumerate(self.__file_lines):
            yield index, line

    def __parse_data(self):
        """Parse function for genbank files, called within class constructor."""
        for index, line in self.__gb_file_generator():
            # Parse for variables: locus_name, molecule_type, genbank_division, modification_name, length.
            if re.search("^LOCUS", line):
                self.__data['identifiers']['locus_name'] = line.split()[1]
                self.__data['molecule_type'] = line.split()[4]
                self.__data['genbank_division'] = line.split()[6]
                self.__data['modification_date'] = line.split()[7]
                self.__data['origin']['length'] = line.split()[2]

            # Parse for variables: definition.
            elif re.search("^DEFINITION", line):
                self.__data['definition'] = ' '.join([item for item in line.split() if item != 'DEFINITION'])

            # Parse for variables: accession_number.
            elif re.search("^ACCESSION", line):
                self.__data['identifiers']['accession_number'] = line.split()[1]

            # Parse for variables: version_number, gen_info_number.
            elif re.search("^VERSION", line):
                self.__data['identifiers']['version_number'] = line.split()[1]
                self.__data['identifiers']['gen_info_number'] = line.split()[2].split(':')[1]

            # Parse for variables: keywords.
            elif re.search("^KEYWORDS", line):
                self.__data['keywords'] = [item for item in line.split() if item != 'KEYWORDS']

            # Parse for variables: common_name.
            elif re.search("^SOURCE", line):
                name = re.search("(?<=[(]).*(?=[)])", line)
                self.__data['source_organism']['common_name'] = name.group(0)

            # Parse for variables: full_phylogenetic_name, taxonomic_name.
            elif re.search("ORGANISM", line):
                idx = index
                # Loop for multiple lines
                while True:
                    idx += 1
                    if re.search("REFERENCE", self.__file_lines[idx]):
                        break
                    new_line = self.__file_lines[idx].strip()
                    self.__data['source_organism']['full_phylogenetic_name'] += new_line
                self.__data['source_organism']['taxonomic_name'] = ' '.join(
                    [item for item in line.split() if item != 'ORGANISM'])

            # Parse for variables: sequence.
            elif re.search("^ORIGIN", line):
                idx = index
                # Loop for multiple lines
                while True:
                    idx += 1
                    if idx > len(self.__file_lines) - 1:
                        break
                    new_line = self.__file_lines[idx]
                    new_line = ''.join([item for item in new_line if not item.isdigit()]).strip()
                    new_line = new_line.replace(' ', '')
                    new_line = new_line.replace('//', '')
                    self.__data['origin']['sequence'] += new_line

            # Parse for variables: CDS.
            elif re.search("^FEATURES", line):
                idx = index
                is_cds = False
                cds_list = []
                raw_cds = []
                # Loop for multiple lines
                while True:
                    new_line = self.__file_lines[idx]
                    if re.search("^ORIGIN", new_line):
                        break
                    new_line_cropped = new_line.removeprefix('     ')
                    if is_cds and new_line_cropped[0] == ' ':
                        raw_cds.append(new_line_cropped)
                    if is_cds and new_line_cropped[0] != ' ':
                        is_cds = False
                        cds_list.append(raw_cds)
                    if re.search("^CDS", new_line_cropped):
                        is_cds = True
                        raw_cds = [new_line_cropped]
                    idx += 1

                # Loop over raw_cds to organize variables into a dict.
                for raw_cds in cds_list:
                    cds = {'range': raw_cds[0].replace('CDS', '').strip()}
                    for i, cds_feature in enumerate(raw_cds):
                        feature = cds_feature.strip()
                        if feature[0] == '/':
                            label = re.search("(?<=[/]).*(?=[=])", feature).group(0)
                            if label == 'translation':
                                indx = i
                                translation = ''
                                # Loop for multiple lines
                                while True:
                                    if indx > len(raw_cds) - 1:
                                        break
                                    translation += raw_cds[indx].replace('/translation="', '').replace('"', '').strip()
                                    indx += 1
                                cds['translation'] = translation
                            else:
                                feature = re.search("(?<==).*", feature).group(0).replace('"', '')
                                cds[label] = feature
                    self.__data['features']['CDS'].append(cds)

    def get_data(self):
        """Return genbank file data object."""
        return self.__data

    def get_nucleotides(self, **kwargs) -> Sequence:
        """Return nucleotides Sequence object. Accepts 'codon_range' string input."""
        codon_range = kwargs.get('codon_range')
        if codon_range is None:
            sequence = self.__data['origin']['sequence']
            return Sequence(sequence, gv.MOLECULE_TYPES['nucleotide'])

        # Handle complement
        if 'complement' in codon_range:
            comp_val = codon_range.replace('complement', '')
            comp_val = re.search("(?<=[(]).*(?=[)])", comp_val).group(0)
            sequence = self.get_nucleotides(codon_range=comp_val).sequence
            replacements = str.maketrans({'a': 't', 't': 'a', 'g': 'c', 'c': 'g'})
            comp_sequence = sequence.translate(replacements)
            comp_sequence = comp_sequence[::-1]
            return Sequence(comp_sequence, gv.MOLECULE_TYPES['nucleotide'])

        # Handle join values
        elif 'join' in codon_range:
            join_values = codon_range.replace('join', '')
            join_values = re.search("(?<=[(]).*(?=[)])", join_values).group(0)
            join_values = join_values.split(',')
            join_sequence = ''
            for value in join_values:
                seq = self.get_nucleotides(codon_range=value).sequence
                join_sequence += seq
            return Sequence(join_sequence, gv.MOLECULE_TYPES['nucleotide'])

        # Normal codon range inputs
        else:
            try:
                codon_start_end = codon_range.split('..')
                if len(codon_start_end) != 2:
                    warnings.warn('Incorrect codon range syntax at get_nucleotides function call.\n'
                                  f'{codon_range} does not have 2 values separated by ".."')
                start_index = int(codon_start_end[0].replace('<', '')) - 1
                end_index = int(codon_start_end[1].replace('>', ''))
                if (start_index < 0) or (end_index < 0) or start_index > end_index:
                    warnings.warn('Incorrect codon range syntax at get_nucleotides function call.\n'
                                  f'{codon_range} values should be > 0 and ascending')
                sequence = self.__data['origin']['sequence'][start_index:end_index]
                return Sequence(sequence, gv.MOLECULE_TYPES['nucleotide'])
            except ValueError:
                warnings.warn('Incorrect codon range syntax at get_nucleotides function call.\n'
                              f'Value {codon_range} should be of type "int..int"')

    def get_aminoacids(self, **kwargs) -> Sequence:
        """Return aminoacid Sequence object. Accepts 'codon_range' string input."""
        nucleotides = self.get_nucleotides(**kwargs).sequence
        aminoacids = ''
        for index, nucleotide in enumerate(nucleotides):
            # Evaluate nucleotide string every 3 characters.
            if (index + 1) % 3 == 0:
                codon = nucleotides[0:3]
                nucleotides = nucleotides[3:]
                aminoacids += gv.CODON_TABLE[codon]
                # If input range is not multiple of 3 some nucleotides might be left.
                # These are not evaluated.

        return Sequence(aminoacids, gv.MOLECULE_TYPES['aminoacid'])


gb_files = []
for file in sys.argv:
    if file.endswith('.gb'):
        gb_files.append({'name': file, 'data': None, 'amino': None})

for gb_file in gb_files:
    gb = GenBank(gb_file['name'])
    gb_file['data'] = gb.get_data()
    gb_file['amino'] = gb.get_aminoacids()

for gb_file in gb_files:
    # Header with brief file info
    # All CDS
    # Name of protein
    # Amino sequence
    # Molecular weight
    print('__________')
    print(f'\n{gv.COLORS['RED']}{gb_file['data']['identifiers']['locus_name']}{gv.COLORS['END']}\n')
    print(f'{gv.COLORS['BOLD']}Molecule Type:{gv.COLORS['END']} {gb_file['data']['molecule_type']}')
    print(f'{gv.COLORS['BOLD']}Keywords:{gv.COLORS['END']} {' '.join(gb_file['data']['keywords'])}')
    for index, feature in enumerate(gb_file['data']['features']['CDS']):
        print(f'\n\tCDS #{index}\n')
        if 'protein_id' in feature:
            print(f'{gv.COLORS['BOLD']}Protein Id:{gv.COLORS['END']} {feature['protein_id']}')
        if 'product' in feature:
            print(f'{gv.COLORS['BOLD']}Product:{gv.COLORS['END']} {feature['product']}')
        if 'translation' in feature:
            print(f'{gv.COLORS['BOLD']}Amino Sequence:{gv.COLORS['END']} {feature['translation']}')
            print(f'{gv.COLORS['BOLD']}Sequence Weight:{gv.COLORS['END']} '
                  f'{Sequence(feature['translation'], gv.MOLECULE_TYPES['aminoacid']).weight}')
    print('__________')
