#!/usr/bin/env python3

import os
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
    check_gc_content,
    iterate_in_fastq
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


def filter_fastq(
    input_fastq: str, 
    output_fastq: str,
    gc_bounds: int | tuple[int, int] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0,
) -> None:
    """ 
    Filter sequences with correct GC-content, length and quality
    
    Arguments:
    input_fastq: a path to an input file 
    output_fastq: a path to an output file
    gc_bounds: a required percentage bounds of GC-pairs 
    length_bounds: a required length bounds of sequences
    quality_threshold: a required threshold of quality according to the Fred scale
    
    Returns None 
    Raises exception if output file name already exist.
    """
    if not os.path.exists("filtered"):
        os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)
    if os.path.exists(output_path):
        raise FileExistsError("This file name exists...")

    with open(output_path, 'w') as file_out:
        for seq_name, seq, plus, qual in iterate_in_fastq(input_fastq):
           if (
            check_length(seq, length_bounds)
            and check_quality(qual, quality_threshold)
            and check_gc_content(seq, gc_bounds)
        ):
            file_out.write(f"@{seq_name}\n")
            file_out.write(f"{seq}\n")
            file_out.write(f"{plus}\n")                                      
            file_out.write(f"{qual}\n")
