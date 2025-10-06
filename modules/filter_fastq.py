def check_gc_content(seq, gc_bounds):
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

def check_length(seq, length_bounds):
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
        
def check_quality(quality, quality_threshold):
    """ 
    Check if the sequence quality meets the necessary threshold.
    Arguments:
    quality: int / float
    quality_threshold: int / float

    Returns bool. 
    Raises exception if y is 0
    """
    qual_scores = [ord(letter) - 33 for letter in quality]
    avg_qual = sum(qual_scores) / len(qual_scores)
    if isinstance(quality_threshold, (int, float)):
        return avg_qual >= quality_threshold
    else:
        raise ValueError("The quality_threshold isn't a number...")