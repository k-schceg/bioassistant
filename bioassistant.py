from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Tuple, Union, TypeVar
from abc import ABC, abstractmethod
from bio_files_processor import parse_blast_output
from datetime import datetime
import argparse
import logging
import os


logger = logging.getLogger(__name__)
logging.basicConfig(
    filename=datetime.now().strftime("%Y_%m_%d-%H_%M_%S_%f_bioassistant_log.txt"),
    encoding='utf-8',
    level=logging.NOTSET,
    format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %I:%M:%S %p'
)

T = TypeVar('T', bound='BiologicalSequence')

class BiologicalSequence(ABC):

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index: Union[int, slice]) -> T:
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def _is_valid_alphabet(self) -> bool:
         pass


class NucleicAcidSequence(BiologicalSequence):

    """
    A class representing a nucleic acid sequence.

    Attributes:
        sequence (str): A string representing the nucleic acid sequence.
    
    Methods:
        
        __len__():
            Returns the length of the nucleic acid sequence.
        
        __getitem__(index):
            Returns a new NucleicAcidSequence object containing the elements by index or slices of the sequence.
        
        __str__():
            Returns the string representation of the nucleic acid sequence.
        
        reverse():
            Returns a new NucleicAcidSequence object with the sequence reversed.
        
        complement():
            Returns a new NucleicAcidSequence object with the sequence complemented.
        
        reverse_complement():
            Returns a new NucleicAcidSequence object with the reverse complement of the sequence.

    Exceptions:
        ValueError: Raised if the sequence contains invalid characters.
        NotImplementedError: Raised if the base class `NucleicAcidSequence` is directly instantiated without specifying a subclass.
    """

    def __init__(self, sequence):
        self.sequence = sequence
        if not self._is_valid_alphabet():
            raise ValueError("Object contains invalid characters")
        logger.info('Successfull NucleicAcidSequence initialization')
        
    def _is_valid_alphabet(self):
        logger.info('Check valid alphabet')
        if type(self) == NucleicAcidSequence:
            raise NotImplementedError("Cannot create NucleicAcidSequence object, specify the class of the object")
        flag_check = all(base in self.alphabet for base in self.sequence)
        if flag_check:
            logger.info('Alphabet check passed')
        else:
            logger.warning('Alphabet check not passed!')
        return flag_check

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return type(self)(self.sequence[index.start:index.stop:index.step])
        else:
            return type(self)(self.sequence[index])
    
    def __str__(self):
        return self.sequence

    def reverse(self):
        return type(self)(self.sequence[::-1])

    def complement(self):
        logger.info('Complement sequence constructing')
        complemented_sequence = ''.join(self.complement_dict[base] for base in self.sequence)
        return type(self)(complemented_sequence)

    def reverse_complement(self):
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    
    """
    A class representing a DNA sequence. It inherits from the NucleicAcidSequence class.

    Attributes:
        alphabet (set): A set of valid nucleotide bases for DNA (A, T, G, C).
        complement_dict (dict): A dictionary mapping each DNA base to its complement.
    
    Methods:
        
        transcribe():
            Transcribes the DNA sequence into an RNA sequence.
        
    Exceptions:
        ValueError: Raised if the sequence contains invalid characters not in the DNA alphabet.
    """

    def __init__(self, sequence):
        logger.info('Init DNASequence')
        self.complement_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }
        self.alphabet = set("ATGCatgc")
        super().__init__(sequence)
        logger.info('Successfull DNASequence initialization')
     
    def transcribe(self) -> 'RNASequence':
        logger.info('Transcribing sequence...')
        seq_transcribed = RNASequence(self.sequence.replace('T', 'U').replace("t", "u"))
        logger.info('Success!')
        return seq_transcribed
    
class RNASequence(NucleicAcidSequence):

    """
    A class representing a RNA sequence. It inherits from the NucleicAcidSequence class.

    Attributes:
        alphabet (set): A set of valid nucleotide bases for RNA (A, U, G, C).
        complement_dict (dict): A dictionary mapping each RNA base to its complement.
        
    Exceptions:
        ValueError: Raised if the sequence contains invalid characters not in the RNA alphabet.
    """

    def __init__(self, sequence):
        logger.info('Init RNASequence')
        self.complement_dict = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }
        self.alphabet = set("AUGCaugc")
        super().__init__(sequence)
        logger.info('Successfull RNASequence initialization')

class AminoAcidSequence(BiologicalSequence):

    """
    A class representing an amino acid sequence.

    Attributes:
        sequence (str): A string representing the amino acid sequence.
        alphabet (set): A set of valid amino acid symbols (A, C, D, E, F, G, H, I, K, L, M, 
                         N, P, Q, R, S, T, V, W, Y).

    Methods:
        
        __len__():
            Returns the length of the amino acid sequence.
        
        __getitem__(index):
            Returns the elements by index and make slices of a sequence.
        
        __str__():
            Returns the string representation of the amino acid sequence.

        calculate_mass(): 
            Calculates the molecular mass of the protein represented by the amino acid sequence.
        
        
    Exceptions:
        ValueError: Raised if the sequence contains invalid characters.
    """

    amino_acid_masses = {
        'A': 89.1, 'C': 121.2, 'D': 133.1, 'E': 147.1, 'F': 165.2, 'G': 75.1,
        'H': 155.2, 'I': 131.2, 'K': 146.2, 'L': 131.2, 'M': 149.2, 'N': 132.1,
        'P': 115.1, 'Q': 146.2, 'R': 174.2, 'S': 105.1, 'T': 119.1, 'V': 117.1,
        'W': 204.2, 'Y': 181.2
    }

    def __init__(self, sequence):
        logger.info('Init AminoAcidSequence')
        self.sequence = sequence
        self.alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        if not self._is_valid_alphabet():
            raise ValueError("Object contains invalid characters")
        logger.info('Successfull AminoAcidSequence initialization')

    def _is_valid_alphabet(self):
        logger.info('Check valid alphabet')
        flag_check = all(base in self.alphabet for base in self.sequence)
        if flag_check:
            logger.info('Alphabet check passed')
        else:
            logger.warning('Alphabet check not passed!')
        return flag_check
    
    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return type(self)(self.sequence[index.start:index.stop:index.step])
        else:
            return type(self)(self.sequence[index])

    def __str__(self):
        return self.sequence
    
    def calculate_mass(self) -> float:
        total_mass = 0.0
        for amino_acid in self.sequence:
            total_mass += self.amino_acid_masses.get(amino_acid.upper(), 0.0)
        return total_mass


def filter_fastq(
    input_file: str,
    output_file: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0
) -> None:
    
    """Function filter_fastq

    Args:
        input_file: str
            path to input fastq file,
        output_file: str
            path to output fastq file,
        gc_bounds: Tuple[float, float] = (0, 100),
        length_bounds: Tuple[int, int] = (0, 2**32),
        quality_threshold: float = 0

    Returns: None
        writes filtered sequences to output fastq file
    """
    if not os.path.isfile(input_file):
        logger.error(f"Input file {input_file} doesn't exists!")
        raise ValueError(f"Input file {input_file} doesn't exists!")
    if os.path.isfile(output_file):
        logger.warning(f"Output file {output_file} already exists. Overwriting")
    if isinstance(gc_bounds, float):
        gc_bounds = (gc_bounds, gc_bounds)
        logger.warning(f'GC boundary is a single number => setup gc_bounds = {gc_bounds}')
    if (gc_bounds[0] < 0) or (gc_bounds[1] > 100) or (gc_bounds[0] > gc_bounds[1]):
        raise ValueError('Wrong GC bounds. gc_bounds should be > 0 and < 100 and left bound <= right bound')
    if isinstance(length_bounds, int):
        length_bounds = (length_bounds, length_bounds)
        logger.warning(f'sequence length is a single number => setup length_bounds = {length_bounds}')
    if (length_bounds[0] < 0) or (length_bounds[0] > length_bounds[1]):
        raise ValueError('Wrong length bounds. length_bounds should be > 0 and left bound <= right bound')
        
    min_gc, max_gc = gc_bounds
    min_length, max_length = length_bounds

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        logger.info(f'Input file {input_file} filtering...')
        filtered_records = []
        for record in SeqIO.parse(infile, "fastq"):
        
            if len(record) < min_length or len(record) > max_length:
                continue
            
            avg_quality = sum(record.letter_annotations["phred_quality"]) / len(record)
            if avg_quality < quality_threshold:
                continue

            gc_content = gc_fraction(record.seq) * 100
            if gc_content < min_gc or gc_content > max_gc:
                continue
            
            filtered_records.append(record)
        logger.info(f'Input file filtered')
        logger.info(f'Writing to output file: {output_file}')
        SeqIO.write(filtered_records, outfile, "fastq")
        logger.info(f'Success!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='Bioassistant', 
        epilog="See '<mode> --help' to read about a specific function."
    )
    
    # Argparser with params that are commot to all functions (not used at this version, for future work)
    general_args_parser = argparse.ArgumentParser(add_help=False)
    # general_args_parser.add_argument("--some-general-param", required=True, help="general param help")
    
    subparsers = parser.add_subparsers(dest='mode', help='Bioassistant functions')

    # TRANSCRIBE PARAMS
    transcribe_parser = subparsers.add_parser('transcribe', help='Make DNA transcription', parents=[general_args_parser])
    transcribe_parser.add_argument(
        '--dna',
        required=True,
        help="DNA sequence for transcription"
    )

    # FILTER FASTQ PARAMS
    filter_fastq_parser = subparsers.add_parser(
        'filter-fastq',
        help='Read and filter input .fastq file and write filtered sequences to output .fastq file',
        parents=[general_args_parser]
    )
    filter_fastq_parser.add_argument(
        '--input-file',
        required=True,
        help="path to input .fastq file"
    )
    filter_fastq_parser.add_argument(
        '--output-file',
        required=True,
        help="path to output .fastq file"
    )
    filter_fastq_parser.add_argument(
        '--gc-bounds',
        required=False,
        help="GC proportion in percents. 2 values expected: lower bound and upper bound",
        nargs=2,
        type=int,
        metavar=('lower_bound', 'upper_bound'),
        default=(0, 100)
    )
    filter_fastq_parser.add_argument(
        '--length-bounds',
        required=False,
        help="sequence length bounds",
        nargs=2,
        type=int,
        metavar=('lower_bound', 'upper_bound'),
        default=(0, 2**32)
    )
    filter_fastq_parser.add_argument(
        '--quality-threshold',
        required=False,
        help="sequence length bounds",
        type=int,
        default=0
    )

    # PARSE BLAST PARAMS
    parse_blast_parser = subparsers.add_parser(
        'parse-blast',
        help='Read and parse input blast .txt file and write writes sorted first row from each description to output blast .txt file',
        parents=[general_args_parser]
    )
    parse_blast_parser.add_argument(
        '--input-file',
        dest='input_blast',
        required=True,
        help="path to input blast .txt file"
    )
    parse_blast_parser.add_argument(
        '--output-file',
        dest='output_blast',
        required=True,
        help="path to output blast .txt file"
    )
    parse_blast_parser.add_argument(
        '--mode',
        required=False,
        help="output file writing mode. If argument is not provided, function will raise error if output file already exists",
        choices=['overwrite', 'append'],
        type=str,
        default=None
    )
    
    args = vars(parser.parse_args())
    func_kwargs = dict(args)
    del func_kwargs['mode']
    
    #import pdb; pdb.set_trace()
    if args['mode'] == 'transcribe':
        dna = DNASequence(func_kwargs['dna'])
        print(dna.transcribe())
    elif args['mode'] == 'filter-fastq':
        filter_fastq(**func_kwargs)
    elif args['mode'] == 'parse-blast':
        parse_blast_output(**func_kwargs)    
    
    parser.parse_args()
