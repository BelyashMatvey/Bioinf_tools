from typing import Dict, Union

dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 't': 'a', 'g': 'c', 'c': 'g'}
dna_to_dna = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
rna_to_rna = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
rna_to_dna = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 't': 'a', 'g': 'c', 'c': 'g'}
error_string = 'Invalid input: string must be either DNA or RNA'
error_string_rna = "The string must be DNA. Given an RNA string"
transcr_map_dna = {'A': 'A', 'G': 'G', 'C': 'C', 'T': 'U', 'a': 'a', 'g': 'g', 'c': 'c', 't': 'u'}
transcr_map_rna = {'A': 'A', 'G': 'G', 'C': 'C', 'U': 'T', 'a': 'a', 'g': 'g', 'c': 'c', 'u': 't'}


def reverse(seq: str, type_: int) -> str:
    """
    Function reverse
    Application - returning an expanded sequence

    Аттрибуты
    ----------
    seq : str
        seq - nucleic acid sequence
    type_ : int
        type_ - type of nucleic acid (1 - DNA, 0 - RNA)
    """
    return seq[::-1]


def is_bounded(bounds: tuple[int, int], x: Union[int, float]) -> bool:
    """
    Function is_bounded
    Application - Checks whether the variable x is between bounds[1] and bounds[0] (inclusive).
    Attributes
    ----------
    bounds : tuple[int,int]
        bounds - boundaries for which a check is made to see whether the variable x lies within the given boundaries.
    x: int
        x - variable being tested
    """
    return bounds[0] <= x <= bounds[1]


def complement(seq: str, type_: int) -> str:
    """
    Function complement
    The main application is transcribing the input sequence depending on its type (DNA/RNA).
    To do this, it uses dictionaries of the corresponding remainders - dna_to_dna,rna_to_rna

    Attributes
    ----------
    seq : str
        seq - nucleic acid sequence
    type_ : int
        type_ - type of nucleic acid (1 - DNA, 0 - RNA)"""
    compl_map: Dict[str, str] = dna_to_dna if type_ else rna_to_rna
    return ''.join(compl_map[c] for c in seq)


def reverse_complement(seq: str, type_: int) -> str:
    """Function transcribe
    Application - finding the reverse complementary strand for the input sequence depending on its type (DNA/RNA).
    To do this, it uses functions from the corresponding complement,reverse

    Attributes
    ----------
    seq : str
        seq - nucleic acid sequence
    type_ : int
        type_ - type of nucleic acid (1 - DNA, 0 - RNA)"""
    return reverse(complement(seq, type_), type_)


def transcribe(seq: str, type_: int) -> str:
    """
    Function transcribe
    Application - transcribing the input sequence depending on its type (DNA/RNA).
    To do this, she uses dictionaries of the corresponding remainders - transcr_map_dna,transcr_map_rna

    Attributes
    ----------
    seq : str
        seq - nucleic acid sequence
    type_ : int
        type_ - type of nucleic acid (1 - DNA, 0 - RNA)
    """
    if type_:
        return ''.join(transcr_map_dna[c] for c in seq)
    return ''.join(transcr_map_rna[c] for c in seq)


def transcribe_dna_complement(seq: str, type_: int) -> str:
    """
    Функция transcribe_dna_complement
    Application - transcribing the input DNA sequence, and then constructing and returning its complementary strand.
    The functions complement and transcribe are used for this.
    type_=1 is passed to the transcribe function, since the DNA strand is initially given. type_=0 is passed to the complement function,
    because the RNA strand is transferred there
    Attributes
    ----------
    seq : str
        seq - nucleic acid sequence
    type_ : int
        type_ - type of nucleic acid (1 - DNA, 0 - RNA)
    """
    return complement(transcribe(seq, type_), 1 - type_)
