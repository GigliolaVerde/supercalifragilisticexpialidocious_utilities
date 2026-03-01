from abc import ABC, abstractmethod
from typing import Dict, Union, Tuple
import os
from Bio import SeqIO
from Bio.SeqUtils import GC


class BiologicalSequence(ABC):
    """Abstract class for biological sequences and their properties"""

    def __init__(self, monomers: str):
        self._monomers = monomers
        if not self.check_alphabet():
            raise ValueError(f"Invalid nucleotides in sequence: {monomers}")

    def __len__(self) -> int:
        """Returns the length of the sequence"""
        return len(self._monomers)

    def __getitem__(self, index):
        """Get the separate element of the sequence or a slice"""
        return self.__class__(self._monomers[index])

    def __str__(self) -> str:
        return "".join(self._monomers)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self._monomers}')"

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    """The class for nucleic acids and its properties"""

    def __init__(self, nucleotides):
        super().__init__(nucleotides)

    def _get_complement_pairs(self) -> Dict[str, str]:
        """
        Return the complement pairs depending on the type of acid (DNA or RNA)
        """
        raise NotImplementedError("Should be implemented in child classes")

    def complement(self) -> 'NucleicAcidSequence':
        """Return complementary sequence."""
        complement_pairs = self._get_complement_pairs()
        res = map(lambda nucleotide: complement_pairs[nucleotide], self._monomers)
        return self.__class__(''.join(res))

    def reverse(self) -> 'NucleicAcidSequence':
        """Return reversed sequence."""
        return self.__class__(self._monomers[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        """Return reverse complement sequence."""
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences and their properties"""

    def _get_complement_pairs(self) -> Dict[str, str]:
        """Returns complement pairs for DNA"""
        return {
            'A': 'T', 'a': 't',
            'T': 'A', 't': 'a',
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }

    def transcribe(self) -> 'RNASequence':
        """Transcribe DNA to RNA"""
        return RNASequence(self._monomers.replace('T', 'U').replace('t', 'u'))

    def check_alphabet(self) -> bool:
        """Check if the sequence contains only DNA nucleotides"""
        dna_nucl = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
        return set(self._monomers).issubset(dna_nucl)


class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences and their properties"""

    def _get_complement_pairs(self) -> Dict[str, str]:
        """Returns complement pairs for RNA"""
        return {
            'A': 'U', 'a': 'u',
            'U': 'A', 'u': 'a',
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }

    def check_alphabet(self) -> bool:
        """Check if the sequence contains only RNA nucleotides"""
        rna_nucl = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
        return set(self._monomers).issubset(rna_nucl)


class AminoAcidSequence(BiologicalSequence):
    """Class for Amino acid sequences and their properties"""

    def calculate_mol_weight(self) -> float:
        """
        Calculate approximate molecular weight (taking the average amino acid weight as 110 Da)
        """
        return len(self._monomers) * 110.0

    def check_alphabet(self) -> bool:
        """Check if the sequence contains only amino acids"""
        acids = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy')
        return set(self._monomers).issubset(acids)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, Tuple[int, int]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """
    Filter sequences with correct GC-content, length and quality using Biopython.

    Arguments:
    input_fastq: a path to an input file
    output_fastq: a path to an output file
    gc_bounds: a required percentage bounds of GC-pairs
    length_bounds: a required length bounds of sequences
    quality_threshold: a required threshold of quality (Phred score)

    Returns None
    Raises exception if output file name already exist.
    """
    if not os.path.exists("filtered"):
        os.makedirs("filtered", exist_ok=True)

    output_path: str = os.path.join("filtered", output_fastq)
    if os.path.exists(output_path):
        raise FileExistsError("This file name exists...")

    if isinstance(gc_bounds, int):
        gc_min, gc_max = 0, gc_bounds
    else:
        gc_min, gc_max = gc_bounds

    if isinstance(length_bounds, int):
        len_min, len_max = 0, length_bounds
    else:
        len_min, len_max = length_bounds

    filtered_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
        sequence = str(record.seq)
        seq_length = len(sequence)

        if not (len_min <= seq_length <= len_max):
            continue

        gc_percent = GC(record.seq)
        if not (gc_min <= gc_percent <= gc_max):
            continue

        if quality_threshold > 0:
            avg_quality = sum(record.letter_annotations["phred_quality"]) / seq_length
            if avg_quality < quality_threshold:
                continue

        filtered_records.append(record)

    with open(output_path, "w") as output_file:
        SeqIO.write(filtered_records, output_file, "fastq")
