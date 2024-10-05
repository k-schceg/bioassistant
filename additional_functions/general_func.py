from typing import List, Callable


def is_dna(seq: str) -> bool:
    if set(seq).issubset({"A", "T", "G", "C", "a", "t", "g", "c"}) and (len(seq) > 0):
        return True
    else:
        return False


def is_rna(seq: str) -> bool:
    if set(seq).issubset({"A", "U", "G", "C", "a", "u", "g", "c"}) and (len(seq) > 0):
        return True
    else:
        return False


def apply(elements: List[str], func: Callable[[str], str]) -> List[str]:
    result = []
    for arg in elements:
        if (not is_dna(arg)) and (not is_rna(arg)):
            raise ValueError("It's not DNA or RNA!")
        result.append(func(arg))
    return result
