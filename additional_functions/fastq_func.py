import additional_functions.dna_rna_func as drf
import additional_functions.fastq_func as qf
from typing import List, Dict, Tuple, Union


def quality(seq: str) -> float:
    quality_mean = sum([ord(i) - 33 for i in seq]) / len(seq)
    return quality_mean


def filter_fastq_from_dict(
    seqs: Dict[str, Tuple[str, str]],
    gc_bounds: Union[Tuple[float, float], float] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: float = 0,
) -> Dict[str, Tuple[str, str]]:
    """Function filter_fastq_from_dict

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
