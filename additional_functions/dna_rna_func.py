from additional_functions.general_func import is_dna


def transcribe(dna: str) -> str:
    if not is_dna(dna):
        raise ValueError("Attempt to transcribe RNA!")
    transcript = dna.replace("T", "U").replace("t", "u")
    return transcript


def reverse(seq: str) -> str:
    reverse_seq = seq[::-1]
    return reverse_seq


def complement(seq: str) -> str:
    if is_dna(seq):
        complement_dict = {
            "A": "T",
            "T": "A",
            "G": "C",
            "C": "G",
            "a": "t",
            "t": "a",
            "g": "c",
            "c": "g",
        }
    else:
        complement_dict = {
            "A": "U",
            "U": "A",
            "G": "C",
            "C": "G",
            "a": "u",
            "u": "a",
            "g": "c",
            "c": "g",
        }
    complement_seq = "".join([complement_dict[i] for i in seq])
    return complement_seq


def reverse_complement(seq: str) -> str:
    reverse_complement_value = complement(reverse(seq))
    return reverse_complement_value


def gc_content(seq: str) -> str:
    gc_count = 0
    for letter in seq:
        if letter in ["G", "C", "g", "c"]:
            gc_count += 1
    gc_content_value = round((gc_count / len(seq)) * 100, 3)
    return gc_content_value


def is_palindrome(seq: str) -> bool:
    if reverse_complement(seq) == seq:
        return True
    else:
        return False
