import additional_functions.dna_rna_func as drf
import additional_functions.fastq_func as qf
import additional_functions.general_func as gf
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
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> Dict[str, Tuple[str, str]]:
    """Function filter_fastq

    Args:
        seqs: Dict[str, Tuple[str, str]],
        gc_bounds: Tuple[float, float] = (0, 100),
        length_bounds: Tuple[int, int] = (0, 2**32),
        quality_threshold: float = 0

    Returns: Dict[str, Tuple[str, str]]
    """
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    filtered_fastq = {}
    for name, seq in seqs.items():
        if (drf.gc_content(seq[0]) < gc_bounds[0]) or (
            drf.gc_content(seq[0]) > gc_bounds[1]
        ):
            continue
        if (len(seq[0]) < length_bounds[0]) or (len(seq[0]) > length_bounds[1]):
            continue
        if qf.quality(seq[1]) < quality_threshold:
            continue
        filtered_fastq[name] = seq
    return filtered_fastq
