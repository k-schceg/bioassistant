def quality(seq: str) -> float:
    quality_mean = sum([ord(i) - 33 for i in seq]) / len(seq)
    return quality_mean
