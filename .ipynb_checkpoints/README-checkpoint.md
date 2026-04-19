markdown
# Supercalifragilisticexpialidocious_utilities

A set of bioinformatics utilities for working with biological sequences (DNA, RNA, proteins) and FASTQ files.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/GigliolaVerde/supercalifragilisticexpialidocious_utilities.git
cd supercalifragilisticexpialidocious_utilities
```

## Features
The package provides two main possibilities:

* Sequence manipulation with classes for DNA, RNA, and protein sequences

* FASTQ quality filtering with multiple filtering criteria



The package implements an object-oriented hierarchy for biological sequences:

* BiologicalSequence (abstract base class)

* NucleicAcidSequence (abstract)

* DNASequence

* RNASequence

* AminoAcidSequence

## Basic Usage
```python
from supercalifragilisticexpialidocious_tools import DNASequence, RNASequence, AminoAcidSequence

# Create DNA sequence
dna = DNASequence("ATCGATCG")
print(dna)  # ATCGATCG
print(len(dna))  # 8

# Reverse complement
rev_comp = dna.reverse_complement()
print(rev_comp)  # CGATCGAT

# Transcribe to RNA
rna = dna.transcribe()
print(rna)  # AUCGAUCG

# RNA operations
rna_complement = rna.complement()
print(rna_complement)  # UAGCUAGC

# Protein sequence
protein = AminoAcidSequence("ACDEFGHIK")
print(protein.calculate_mol_weight())  # 990.0 (9 * 110)
```
## Available Methods
### For DNA/RNA sequences:

* complement() - generate complementary sequence

* reverse() - reverse sequence

* reverse_complement() - reverse complementary sequence

* check_alphabet() - validate sequence characters

### For DNA sequences only:

* transcribe() - convert DNA to RNA

### For protein sequences:

* calculate_mol_weight() - approximate molecular weight (110 Da per amino acid)

## FASTQ Filtering
The filter_fastq function filters reads from FASTQ files based on multiple criteria:

* GC-content - filter by GC percentage bounds

* Sequence length - filter by minimum and/or maximum length

* Quality score - filter by minimum average Phred quality score

### Usage Examples
```python
from supercalifragilisticexpialidocious_tools import filter_fastq

# Basic quality filtering
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="high_quality.fastq",
    quality_threshold=20
)

# Filter by length only (sequences <= 150 bp)
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="short_reads.fastq",
    length_bounds=150  # equivalent to (0, 150)
)

# Filter by GC content (40-60% GC)
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="gc_filtered.fastq",
    gc_bounds=(40, 60)
)

# Combined filtering with all parameters
filter_fastq(
    input_fastq="input.fastq",
    output_fastq="stringent_filtered.fastq",
    gc_bounds=(45, 55),
    length_bounds=(75, 125),
    quality_threshold=30
)
```
## Important Notes
* Output files are automatically saved in a filtered/ directory

* The function raises FileExistsError if the output file already exists

* When single integer values are provided for bounds, they're interpreted as upper limits:

* gc_bounds=50 means GC content between 0-50%

* length_bounds=200 means length between 0-200 bp

## Contact
For questions and suggestions, please create an issue in the repository.

### License
MIT License