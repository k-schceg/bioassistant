import pytest
import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
import bioassistant
from contextlib import nullcontext as does_not_raise


@pytest.mark.parametrize(
    ("original_sequence", "expected_sequence"),
    [
        ('ATGC', 'AUGC'),
        ('AtGc', 'AuGc'),
    ],
)
def test_transcribe(original_sequence: str, expected_sequence: str) -> None:
    dna = bioassistant.DNASequence(original_sequence)
    transcribe_result = str(dna.transcribe())
    assert transcribe_result == expected_sequence, 'Wrong sequence transcibation'

@pytest.mark.parametrize(
    ("original_sequence", "expectation"),
    [
        ('ATGC', does_not_raise()),
        ('ASGC', pytest.raises(ValueError)),
        ('AT.GC', pytest.raises(ValueError)),
        ('A3GC', pytest.raises(ValueError)),
    ],
)
def test_transcribe_wrong_input_raise_error(original_sequence, expectation):
    """Test how much I know division."""
    with expectation:
        dna = bioassistant.DNASequence(original_sequence)
        transcribe_result = str(dna.transcribe())

def test_fasta_filter_create_file(tmpdir):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_fastq_filter_input_file.fastq')
    output_file = tmpdir.join('test_fastq_filter_output_file.fastq')
    bioassistant.filter_fastq(
      input_file, 
      str(output_file),
      gc_bounds=(0, 100),
      length_bounds=(0, 25),
      quality_threshold=17
    )
    assert output_file.exists(), "Output .fasta file doesn't exists"

@pytest.mark.parametrize(
    ("bounds", "expectation"),
    [
        ([0, 100], does_not_raise()),
        ([30, 40], does_not_raise()),
        ([-10, 40], pytest.raises(ValueError)),
        ([30, 130], pytest.raises(ValueError)),
        ([40, 10], pytest.raises(ValueError)),
    ],
)
def test_fasta_filter_gc_bounds(tmpdir, bounds, expectation):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_fastq_filter_input_file.fastq')
    output_file = tmpdir.join('test_fastq_filter_output_file.fastq')
    with expectation:
        bioassistant.filter_fastq(
          input_file, 
          str(output_file),
          gc_bounds=bounds,
          length_bounds=(0, 25),
          quality_threshold=17
        )

@pytest.mark.parametrize(
    ("bounds", "expectation"),
    [
        ([0, 2**32], does_not_raise()),
        ([30, 10_000], does_not_raise()),
        ([-10, 40], pytest.raises(ValueError)),
        ([40, 10], pytest.raises(ValueError)),
    ],
)
def test_fasta_filter_length_bounds(tmpdir, bounds, expectation):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_fastq_filter_input_file.fastq')
    output_file = tmpdir.join('test_fastq_filter_output_file.fastq')
    with expectation:
        bioassistant.filter_fastq(
          input_file, 
          str(output_file),
          gc_bounds=(0, 100),
          length_bounds=bounds,
          quality_threshold=17
        )

def test_fasta_filter_output(tmpdir):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_fastq_filter_input_file.fastq')
    output_file = tmpdir.join('out.fastq')
    expected_file = str(pathlib.Path(__file__).parents[0] / 'test_fastq_filter_output_file.fastq')
    bioassistant.filter_fastq(
      input_file, 
      str(output_file),
      gc_bounds=(0, 100),
      length_bounds=(0, 25),
      quality_threshold=17
    )
    assert open(output_file).read() == open(expected_file).read()

def test_blast_parse_create_file(tmpdir):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_blast_parse_input.txt')
    output_file = tmpdir.join('test_output_file.txt')
    bioassistant.parse_blast_output(input_file, str(output_file))
    assert output_file.exists(), "Output blast .txt file doesn't exists"

def test_blast_parse_output(tmpdir):
    input_file = str(pathlib.Path(__file__).parents[0] / 'test_blast_parse_input.txt')
    output_file = tmpdir.join('out.txt')
    expected_file = str(pathlib.Path(__file__).parents[0] / 'test_blast_parse_output.txt')
    bioassistant.parse_blast_output(input_file, str(output_file))
    assert open(output_file).read() == open(expected_file).read()



