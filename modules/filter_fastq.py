def iterate_in_fastq(input_fastq: str) -> dict[str, tuple[str, str]]:
    """
    Transform fastq-file into a dict.
    Arguments:
    input_fastq: path to the fastq-file.

    Returns dict where key is sequence name 
    and the value is a tuple (sequence, quality).
    """
    with open(input_fastq,'r') as file_in:
        line = file_in.readline().strip()
        while line:
            if line.startswith('@SRX'):
                seq_name = line[1:].split()[0]
                seq = file_in.readline().strip()
                plus = file_in.readline().strip() 
                qual = file_in.readline().strip()
                yield seq_name, seq, plus, qual         
            line = file_in.readline().strip()

def check_gc_content(seq: str, gc_bounds: int | tuple[int, int]) -> bool:
    """
    Check if the sequence length meets the demands.
    Arguments:
    seq: str
    pct_bounds: int / tuple

    Returns bool. 
    Raises exception if GC-content bounds is not a tuple with two numbers. 
    """
    gc_count = (seq.upper().count('G') + seq.upper().count('C'))
    gc_percent = (gc_count / len(seq)) * 100
    if isinstance(gc_bounds, (int, float)):
        return gc_percent < gc_bounds
    elif len(gc_bounds) == 2:
        if all(map(lambda x: isinstance(x, (int, float)), gc_bounds)):
            return gc_bounds[0] <= gc_percent <= gc_bounds[1]            
    else:
        raise ValueError("Check the GC-percent bounds...")

def check_length(seq: str, length_bounds: int | tuple[int, int]) -> bool:
    """
    Check if the sequence length meets the demands.
    Arguments:
    seq: str
    length_bounds: int / tuple

    Returns bool. 
    Raises exception if y is 0
    """
    if isinstance(length_bounds, (int, float)):
        return len(seq) < length_bounds
    if len(length_bounds) == 2:
        if all(map(lambda x: isinstance(x, (int, float)), length_bounds)):
            return length_bounds[0] <= len(seq) <= length_bounds[1]
    else:
        raise ValueError("Check the length bounds...")
        
def check_quality(quality, quality_threshold: int) -> bool:
    """ 
    Check if the sequence quality meets the necessary threshold.
    Arguments:
    quality: int
    quality_threshold: int

    Returns bool. 
    Raises exception if y is 0
    """
    qual_scores = [ord(letter) - 33 for letter in quality]
    avg_qual = sum(qual_scores) / len(qual_scores)
    if isinstance(quality_threshold, (int, float)):
        return avg_qual >= quality_threshold
    else:
        raise ValueError("The quality_threshold isn't a number...")
