#!/usr/bin/env python3


from modules.dna_rna_tools import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement
)
from modules.filter_fastq import (
    check_length,
    check_quality,
    check_gc_content
)

def run_dna_rna_tools(*args):
    """Main function for DNA/RNA sequence operations."""
    if len(args) == 0:
        raise ValueError("Where are the arguments???")

    action_name = args[-1]
    sequences = args[:-1]

    if len(sequences) == 0:
        raise ValueError("Where are the sequences???")

    actions = {
        'is_nucleic_acid': is_nucleic_acid,
        'transcribe': transcribe,
        'reverse': reverse,
        'complement': complement,
        'reverse_complement': reverse_complement
    }

    if action_name not in actions:
        raise ValueError("I don't know this action")

    final_results = []
    for sequence in sequences:
        if not isinstance(sequence, str):
            raise TypeError("Sequences must be strings...")

        if action_name == 'is_nucleic_acid':
            result = actions[action_name](sequence)
        else:
            if not is_nucleic_acid(sequence):
                raise ValueError("Some sequences aren't nucleic acids...")
            result = actions[action_name](sequence)

        final_results.append(result)

    return final_results[0] if len(final_results) == 1 else final_results

def filter_fastq(seqs: dict[str, tuple], gc_bounds = (0, 100), length_bounds= (0, 2**32), quality_threshold = 0):
    """ 
    Filter sequences with correct GC-content, length and quality
    
    Arguments:
    seqs: dict[str, tuple]
    gc_bounds: int / tuple
    length_bounds: int / tuple
    quality_threshold: int
    
    Returns dict. 
    Raises exception if seqs is empty or 
    """
    if not seqs:
        raise ValueError("Where are the sequences???")
    filtered_seqs = {}
    for name, (seq, quality) in seqs.items():
        if check_length(seq, length_bounds):
            if check_quality(quality, quality_threshold):
                if check_gc_content(seq, gc_bounds):
                    filtered_seqs[name] = (seq, quality)
    return filtered_seqs
