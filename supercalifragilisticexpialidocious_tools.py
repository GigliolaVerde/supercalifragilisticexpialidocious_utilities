from abc import ABC, abstractmethod
from typing import Dict, Union, Tuple
import os
import argparse
import logging
from Bio import SeqIO
from Bio.SeqUtils import GC


def setup_logger(log_file: str = "fastq_filter.log"):
    "Create logger with writing into a file and the console"

    logger = logging.getLogger("FastQFilter")
    logger.setLevel(logging.DEBUG)
    
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger


class BiologicalSequence(ABC):
    """Abstract class for biological sequences and their properties"""

    def __init__(self, monomers: str):
        self._monomers = monomers
        if not self.check_alphabet():
            raise ValueError(f"Invalid nucleotides in sequence: {monomers}")

    def __len__(self) -> int:
        return len(self._monomers)

    def __getitem__(self, index):
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
        raise NotImplementedError("Should be implemented in child classes")

    def complement(self) -> 'NucleicAcidSequence':
        complement_pairs = self._get_complement_pairs()
        res = map(lambda nucleotide: complement_pairs[nucleotide], self._monomers)
        return self.__class__(''.join(res))

    def reverse(self) -> 'NucleicAcidSequence':
        return self.__class__(self._monomers[::-1])

    def reverse_complement(self) -> 'NucleicAcidSequence':
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences and their properties"""

    def _get_complement_pairs(self) -> Dict[str, str]:
        return {
            'A': 'T', 'a': 't',
            'T': 'A', 't': 'a',
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }

    def transcribe(self) -> 'RNASequence':
        return RNASequence(self._monomers.replace('T', 'U').replace('t', 'u'))

    def check_alphabet(self) -> bool:
        dna_nucl = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c'}
        return set(self._monomers).issubset(dna_nucl)


class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences and their properties"""

    def _get_complement_pairs(self) -> Dict[str, str]:
        return {
            'A': 'U', 'a': 'u',
            'U': 'A', 'u': 'a',
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }

    def check_alphabet(self) -> bool:
        rna_nucl = {'A', 'U', 'G', 'C', 'a', 'u', 'g', 'c'}
        return set(self._monomers).issubset(rna_nucl)


class AminoAcidSequence(BiologicalSequence):
    """Class for Amino acid sequences and their properties"""

    def calculate_mol_weight(self) -> float:
        return len(self._monomers) * 110.0

    def check_alphabet(self) -> bool:
        acids = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy')
        return set(self._monomers).issubset(acids)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[int, Tuple[int, int]] = (0, 100),
    length_bounds: Union[int, Tuple[int, int]] = (0, 2**32),
    quality_threshold: int = 0,
    logger: logging.Logger = None,
) -> Tuple[int, int]:
    """
    Filter sequences with correct GC-content, length and quality using Biopython.

    Arguments:
    input_fastq: a path to an input file
    output_fastq: a path to an output file
    gc_bounds: a required percentage bounds of GC-pairs
    length_bounds: a required length bounds of sequences
    quality_threshold: a required threshold of quality (Phred score)
    logger: logger instance for logging

    Returns Tuple[int, int]: (total_reads, filtered_reads)
    Raises exception if output file name already exist.
    """
    if logger is None:
        logger = setup_logger()
    
    logger.info(f"Starting FastQ filtering: input={input_fastq}, output={output_fastq}")
    
    if not os.path.exists("filtered"):
        os.makedirs("filtered", exist_ok=True)
        logger.debug(f"Created directory: filtered")

    output_path: str = os.path.join("filtered", output_fastq)
    if os.path.exists(output_path):
        error_msg = f"Output file already exists: {output_path}"
        logger.error(error_msg)
        raise FileExistsError(error_msg)

    if isinstance(gc_bounds, int):
        gc_min, gc_max = 0, gc_bounds
    else:
        gc_min, gc_max = gc_bounds

    if isinstance(length_bounds, int):
        len_min, len_max = 0, length_bounds
    else:
        len_min, len_max = length_bounds

    filtered_records = []
    total_records = 0

    try:
        for record in SeqIO.parse(input_fastq, "fastq"):
            total_records += 1
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
    except FileNotFoundError:
        logger.error(f"Input file not found: {input_fastq}")
        raise
    except Exception as e:
        logger.error(f"Error processing FastQ file: {str(e)}")
        raise

    with open(output_path, "w") as output_file:
        SeqIO.write(filtered_records, output_file, "fastq")
    
    filtered_count = len(filtered_records)
    logger.info(f"Filtering complete: {filtered_count} sequences kept out of {total_records}")
    logger.info(f"Output saved to: {output_path}")
    
    return total_records, filtered_count


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="FastQ-sequence filter on GC-content, length, and quality",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTQ file path"
    )
    
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTQ file name (saved in 'filtered' directory)"
    )
    
    parser.add_argument(
        "--gc-bounds",
        nargs="+",
        type=int,
        default=[0, 100],
        help="GC content bounds (one value for max or two for min max). Default value: 0-100"
    )
    
    parser.add_argument(
        "--length-bounds",
        nargs="+",
        type=int,
        default=[0, 2**32],
        help="Sequence length bounds (can be one value for max or two for min max). Default value: 0-inf"
    )
    
    parser.add_argument(
        "--quality-threshold",
        type=int,
        default=0,
        help="Minimum average quality threshold (Phred score).Default value: 0 (disabled)"
    )
    
    parser.add_argument(
        "--log-file",
        type=str,
        default="fastq_filter.log",
        help="Log file path"
    )
    
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="verbose output"
    )
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    
    logger = setup_logger(args.log_file)
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)
    
    gc_bounds = args.gc_bounds
    if len(gc_bounds) == 1:
        gc_bounds = gc_bounds[0]
    elif len(gc_bounds) == 2:
        gc_bounds = tuple(gc_bounds)
    else:
        logger.error("Ooops, gc-bounds must have 1 or 2 values")
        raise ValueError("Ooops, gc-bounds must have 1 or 2 values")
    
    length_bounds = args.length_bounds
    if len(length_bounds) == 1:
        length_bounds = length_bounds[0]
    elif len(length_bounds) == 2:
        length_bounds = tuple(length_bounds)
    else:
        logger.error("Ooops, length-bounds must have 1 or 2 values")
        raise ValueError("Ooops, length-bounds must have 1 or 2 values")
    
    total, filtered = filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=args.quality_threshold,
        logger=logger
    )
    
    print(f"\nYour filtering complete!")
    print(f"Total sequences: {total}")
    print(f"Filtered sequences: {filtered}")
    print(f"Output: filtered/{args.output}")
    print(f"Log file: {args.log_file}")

