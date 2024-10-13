import additional_functions.dna_rna_func as drf
import additional_functions.fastq_func as qf
import additional_functions.general_func as gf
from additional_functions.fastq_func import filter_fastq_from_dict
from typing import List, Dict, Tuple, Union


def run_dna_rna_tools(*args: List[str]) -> List[str]:
    """Function run_dna_rna_tools

    Args:
        *seqs: List[str], proc: str
        Seqs: DNA or RNA
        Proc: "transcribe", "reverse", "complement",
              "reverse_complement", "gc_content", "is_palindrome"

    Returns: List[str]
    """
    *seqs, proc = args
    if proc == "transcribe":
        result = gf.apply(seqs, drf.transcribe)
    elif proc == "reverse":
        result = gf.apply(seqs, drf.reverse)
    elif proc == "complement":
        result = gf.apply(seqs, drf.complement)
    elif proc == "reverse_complement":
        result = gf.apply(seqs, drf.reverse_complement)
    elif proc == "gc_content":
        result = gf.apply(seqs, drf.gc_content)
    elif proc == "is_palindrome":
        result = gf.apply(seqs, drf.is_palindrome)
    if len(result) == 1:
        result = result[0]
    return result


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """Function filter_fastq

    Args:
        input_fastq: str
            path to input fastq file,
        output_fastq: str
            path to output fastq file,
        gc_bounds: Tuple[float, float] = (0, 100),
        length_bounds: Tuple[int, int] = (0, 2**32),
        quality_threshold: float = 0

    Returns: None
        writes filtered sequences to output fastq file
    """
    seq_1 = None
    with open(input_fastq, "r") as input_file:
        with open(output_fastq, "a") as output_file:
            for line in input_file:
                # import pdb; pdb.set_trace()
                if line.startswith("@"):
                    name = line
                elif line.startswith("+"):
                    comment = line
                elif seq_1 is None:
                    seq_1 = line
                else:
                    seq_2 = line
                    seqs = {name: (seq_1, seq_2)}
                    seqs_filtred = filter_fastq_from_dict(
                        seqs, gc_bounds, length_bounds, quality_threshold
                    )
                    if len(seqs_filtred) > 0:
                        output_file.write(name + seq_1 + comment + seq_2)
                    seq_1 = None
                    seq_2 = None
