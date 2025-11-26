from typing import Dict


def is_nucleic_acid(sequence: str) -> bool:
    """Check if sequence is a valid nucleic acid."""
    valid_nucleotides = {'A', 'T', 'G', 'C', 'U'}
    upper_seq = sequence.upper()

    if not set(upper_seq).issubset(valid_nucleotides):
        return False

    if 'T' in upper_seq and 'U' in upper_seq:
        return False

    return True


def transcribe(sequence: str) -> str:
    """Transcribe DNA to RNA."""
    if 'U' not in sequence.upper():
        return sequence.replace('T', 'U').replace('t', 'u')
    else:
        return sequence


def complement(sequence: str) -> str:
    """Return complementary sequence."""
    upper_seq = sequence.upper()

    if 'U' not in upper_seq:
        complement_pairs: Dict[str, str] = {
            'A': 'T', 'a': 't',
            'T': 'A', 't': 'a', 
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }
    else:
        complement_pairs: Dict[str, str] = {
            'A': 'U', 'a': 'u',
            'U': 'A', 'u': 'a',
            'G': 'C', 'g': 'c',
            'C': 'G', 'c': 'g'
        }

    res = map(lambda nucleotide: complement_pairs[nucleotide], sequence)
    return ''.join(res)


def reverse(sequence: str) -> str:
    """Return reversed sequence."""
    return sequence[::-1]


def reverse_complement(sequence: str) -> str:
    """Return reverse complement sequence."""
    return reverse(complement(sequence))
