import os


def convert_multiline_fasta_to_oneline(
    input_fasta: str, output_fasta: str = None
) -> None:
    """
    Transform multiline fasta into oneline.

    Args:
        input_fasta: a path to an input file
        output_fasta: a path to an output file

    Returns: None
    """
    if output_fasta is None:
        base_name = os.path.splitext(input_fasta)[0]
        output_fasta = f"{base_name}_oneline.fasta"

    if os.path.exists(input_fasta):
        with (
            open(input_fasta, 'r') as file_in,
            open(output_fasta, 'w') as file_out
        ):
            current_header = None
            current_seq = []
            for line in file_in:
                line = line.strip()
                if line.startswith('>'):
                    if current_header is not None:
                        total_seq = ''.join(current_seq)
                        file_out.write(
                            f"{current_header}\n{total_seq}\n"
                        )
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(
                        line.replace(" ", "").replace("\t", "")
                    )
            if current_header is not None and current_seq:
                total_seq = ''.join(current_seq)
                file_out.write(f"{current_header}\n{total_seq}\n")
    else:
        raise FileNotFoundError(
            f"Input FASTA file '{input_fasta}' not found. "
            f"Current directory: {os.getcwd()}"
        )


def parse_blast_output(
    input_file: str,
    output_file: str = None
) -> None:
    """
    Parse blast results and extract
    the names of the best matches for every query.

    Args:
        input_file: a path to an input file
        output_file: a path to an output file

    Returns: None
    """
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_parsed.txt"

    if os.path.exists(input_file):
        with open(input_file, 'r') as file_in:
            content = file_in.read()
        best_matches = []
        sections = content.split('Query #')
        for section in sections[1:]:
            lines = section.split('\n')
            in_results_section = False
            header_found = False

            for line in lines:
                if 'Sequences producing significant alignments:' in line:
                    in_results_section = True
                    continue
                if in_results_section:
                    if not line.strip():
                        continue
                    if 'Description' in line or 'Scientific' in line:
                        header_found = True
                        continue
                    if header_found and line.strip():
                        parts = line.split('  ')
                        description = parts[0].strip()
                        best_matches.append(description)
                        break
        best_matches = sorted(list(set(best_matches)))

        with open(output_file, 'w') as file_out:
            for match in best_matches:
                file_out.write(match + '\n')
    else:
        raise FileNotFoundError(
            f"Input file '{input_file}' not found. "
            f"Current directory: {os.getcwd()}"
        )
