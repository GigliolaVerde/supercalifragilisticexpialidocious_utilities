import os

def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Transform multiline fasta into oneline.
    
    Args:
    input_fasta: a path to an imput file
    output_fasta: a path to an output file
    
    Returns: None
    """
    if output_fasta is None:
        base_name = os.path.splitext(input_fasta)[0]
        output_fasta = f"{base_name}_oneline.fasta"
    try:
        with open(input_fasta, 'r') as file_in, open(output_fasta, 'w') as file_out:
            current_header = None
            current_seq = []  
            for line in file_in:
                line = line.strip()                    
                if line.startswith('>'):  
                    if current_header is not None:
                        total_seq = ''.join(current_seq)
                        file_out.write(f"{current_header}\n{total_seq}\n")
                    current_header = line
                    current_seq = []
                else:
                    current_seq.append(line.replace(" ", "").replace("\t", ""))
            if current_header is not None and current_seq:
                total_seq = ''.join(current_seq)
                file_out.write(f"{current_header}\n{total_seq}\n")
    except Exception as e:
        print(f"The error: {e}")



def parse_blast_output(input_file: str, output_file: str = None):
    """
    Parse_blast results and extract the names of the best matches for every query. 
    Args:
    input_fasta: a path to an imput file
    output_fasta: a path to an output file
    
    Returns: None
    """
    
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_parsed.txt"
    try:
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
                    if 'Description' in line and 'Scientific Name' in line:
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
    except Exception as e:
        print(f"The error: {e}")
