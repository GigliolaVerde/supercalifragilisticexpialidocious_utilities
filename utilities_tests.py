import pytest
import os
import tempfile
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from supercalifragilisticexpialidocious_tools import (
    DNASequence, RNASequence, AminoAcidSequence,
    filter_fastq, setup_logger
)


class TestSequences:
    """Tests for DNA, RNA and amino acid sequences"""
    
    def test_dna_invalid(self):
    """Test invalid DNA sequence"""
        with pytest.raises(ValueError, match="Invalid nucleotides in sequence"):
            DNASequence("ATGCX")
    
    def test_dna_complement(self):
        """Test DNA complement"""
        dna = DNASequence("ATGC")
        assert str(dna.complement()) == "TACG"
    
    def test_rna_complement(self):
        """Test RNA complement"""
        rna = RNASequence("AUGC")
        assert str(rna.complement()) == "UACG"
    
    def test_amino_weight(self):
        """Test amino acid weight calculation"""
        aa = AminoAcidSequence("ACDEF")
        assert aa.calculate_mol_weight() == 5 * 110.0


class TestFilter:
    """Tests for FastQ filter"""
    
    @pytest.fixture
    def fastq_file(self):
        """Create temporary FASTQ file"""
        records = [
            SeqRecord(Seq("GCGCGCGCGCGCGCGC"), id="seq1",
                     letter_annotations={"phred_quality": [40] * 16}),
            SeqRecord(Seq("ATATATATATATATAT"), id="seq2",
                     letter_annotations={"phred_quality": [20] * 16}),
        ]
        
        temp_dir = tempfile.mkdtemp()
        temp_file = os.path.join(temp_dir, "test.fastq")
        
        with open(temp_file, "w") as f:
            SeqIO.write(records, f, "fastq")
        
        yield temp_file, temp_dir
        shutil.rmtree(temp_dir)
    
    def test_filter_basic(self, fastq_file):
        """Test basic filtering"""
        input_file, temp_dir = fastq_file
        output_file = "out.fastq"
        
        total, filtered = filter_fastq(
            input_fastq=input_file,
            output_fastq=output_file,
            gc_bounds=(50, 100),
            length_bounds=(10, 20),
            quality_threshold=0
        )
        
        assert total == 2
        assert filtered == 1
        
        output_path = os.path.join("filtered", output_file)
        assert os.path.exists(output_path)
        
      
        if os.path.exists(output_path):
            os.remove(output_path)
        if os.path.exists("filtered") and not os.listdir("filtered"):
            os.rmdir("filtered")
    
    def test_filter_quality(self, fastq_file):
        """Test filtering with quality threshold"""
        input_file, temp_dir = fastq_file
        
        total, filtered = filter_fastq(
            input_fastq=input_file,
            output_fastq="out2.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=30
        )
        
        assert total == 2
        assert filtered == 1
        
        output_path = os.path.join("filtered", "out2.fastq")
        if os.path.exists(output_path):
            os.remove(output_path)
    

    
    def test_filter_input_not_found(self):
        """Test error when input not found"""
        with pytest.raises(FileNotFoundError):
            filter_fastq(
                input_fastq="missing.fastq",
                output_fastq="out.fastq",
                gc_bounds=(0, 100),
                length_bounds=(0, 100),
                quality_threshold=0
            )


class TestLogging:
    """Tests for logging"""
    
    def test_logger_creates_file(self):
        """Test logger creates log file"""
        log_file = "test.log"
        logger = setup_logger(log_file)
        
        logger.info("Test message")
        logger.error("Error message")
        
        assert os.path.exists(log_file)
        
        with open(log_file, "r") as f:
            content = f.read()
            assert "Test message" in content
            assert "Error message" in content
        
        if os.path.exists(log_file):
            os.remove(log_file)
    
    def test_filter_logs_messages(self, fastq_file):
        """Test filter produces log messages"""
        input_file, temp_dir = fastq_file
        log_file = "filter.log"
        logger = setup_logger(log_file)
        
        filter_fastq(
            input_fastq=input_file,
            output_fastq="log_out.fastq",
            gc_bounds=(0, 100),
            length_bounds=(0, 100),
            quality_threshold=0,
            logger=logger
        )
        
        assert os.path.exists(log_file)
        
        with open(log_file, "r") as f:
            content = f.read()
            assert "Starting FastQ filtering" in content
            assert "Filtering complete" in content
        
        output_path = os.path.join("filtered", "log_out.fastq")
        if os.path.exists(output_path):
            os.remove(output_path)
        if os.path.exists(log_file):
            os.remove(log_file)
