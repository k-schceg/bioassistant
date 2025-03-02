from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Tuple, Union, TypeVar
from abc import ABC, abstractmethod

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
        
    def _is_valid_alphabet(self):
        if type(self) == NucleicAcidSequence:
            raise NotImplementedError("Cannot create NucleicAcidSequence object, specify the class of the object")
        return all(base in self.alphabet for base in self.sequence)

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
        self.alphabet = set("ATGC")
        super().__init__(sequence)
     
    def transcribe(self) -> 'RNASequence':
        return RNASequence(self.sequence.replace('T', 'U').replace("t", "u"))
    
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
        self.alphabet = set("AUGC")
        super().__init__(sequence)

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
        self.sequence = sequence
        self.alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        if not self._is_valid_alphabet():
            raise ValueError("Object contains invalid characters")
        
    def _is_valid_alphabet(self):
        return all(base in self.alphabet for base in self.sequence)

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
    
    if isinstance(gc_bounds, float):
        gc_bounds = (gc_bounds, gc_bounds)

    if isinstance(length_bounds, int):
        length_bounds = (length_bounds, length_bounds)

    min_gc, max_gc = gc_bounds
    min_length, max_length = length_bounds

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
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
        
        SeqIO.write(filtered_records, outfile, "fastq")