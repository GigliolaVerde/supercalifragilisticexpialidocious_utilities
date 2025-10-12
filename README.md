# Supercalifragilisticexpialidocious_utilities

A set of bioinformatics utilities for working with DNA/RNA sequences and FASTQ files.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/GigliolaVerde/supercalifragilisticexpialidocious_utilities.git
cd supercalifragilisticexpialidocious_utilities
```
## Usage

## Functionality

- Nucleic acid sequence validation
- DNA to RNA transcription 
- Sequence reversal and complement generation
- FASTQ file quality filtering
- Sequence length-, quality- and GC-content-based filtering

#### Available Operations
is_nucleic_acid - check if sequence is DNA/RNA

transcribe - transcribe DNA → RNA

reverse - reverse sequence

complement - complementary sequence

reverse_complement - reverse complementary sequence

##### Usage Examples

```python
from supercalifragilisticexpialidocious_tools import run_dna_rna_tools

# Check if sequence is valid nucleic acid
is_valid = run_dna_rna_tools("ATCGATCG", "is_nucleic_acid")
print(is_valid)

# Transcribe DNA to RNA
rna_sequence = run_dna_rna_tools("ATCGATCG", "transcribe")
print(rna_sequence)

# Generate reverse complement
rev_comp = run_dna_rna_tools("ATCGATCG", "reverse_complement")
print(rev_comp)

# Process multiple sequences
results = run_dna_rna_tools("ATCG", "GCTA", "TTAA", "reverse_complement")
print(results)
```

#### FASTQ Filtering
The filter_fastq function filters reads by:

Sequence length - filter reads by minimum and/or maximum length

Average quality - filter by minimum average Phred quality score

GC-content - filter by GC percentage bounds

##### Usage Examples
```python
from supercalifragilisticexpialidocious_tools import filter_fastq

# Basic filtering
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="filtered.fastq",
    quality_threshold=20
)

# Filter by length
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="length_filtered.fastq",
    length_bounds=(50, 150)
)

# Filter by GC content
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="gc_filtered.fastq",
    gc_bounds=(40, 60)
)

# Combined filtering
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="high_quality.fastq",
    gc_bounds=(45, 55),
    length_bounds=(75, 125),
    quality_threshold=30
)
```
### FASTA file processing
#### Multiline to oneline conversion
Convert multiline FASTA files to single-line format:
```python
from bio_files_processor import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline("multiline.fasta", "oneline.fasta")
```

#### BLAST output parsing
Extract best matches from BLAST output:
```python
from bio_files_processor import parse_blast_output

parse_blast_output("blast_output.txt", "best_matches.txt")
```

## Contact
For questions and suggestions, please create an issue in the repository.

## License
MIT License
