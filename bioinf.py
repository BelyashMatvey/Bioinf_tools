from typing import Union
from abc import ABC, abstractmethod
import os
import numpy as np

from Bio import SeqIO, SeqUtils


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


def filter_fastq(input_fastq: str, output_fastq: str, gc_bounds: Union[int, tuple[int, int]] = None,
                 length_bounds: Union[int, tuple[int, int]] = None, quality_threshold: int = None) -> None:
    """
    Function filter_fastq
    Application - The function input is the name of the file with data in fastq format. Using the parameters gc_bounds, length_bounds, and
    the quality_threshold function filters the lines of the file and writes to the final file only those that match
    given conditions. To write a file, the function uses the write_fastq module from modules. Writing to a file occurs line by line.
    If the file already exists, it will not be overwritten and an error will be displayed.
    Attributes
    ----------
    input_fastq: str
        Input file name
    output_fastq: str
        Output file name
    gc_bounds: int или tuple[int,int]
        gc_bounds - интервал GC состава (в процентах).
        Если в аргумент передать одно число, то считается, что это верхняя граница.
    length_bounds: int или tuple[int,int]
        Интервал длины для фильтрации. Если в аргумент передать одно число, то считается, что это верхняя граница.
    quality_threshold: int
        Пороговое значение среднего качества рида для фильтрации, по-умолчанию равно 0 (шкала phred33).
    """
    if length_bounds is None:
        length_bounds = (0, 2 ** 32)
    if gc_bounds is None:
        gc_bounds = (0, 100)
    if quality_threshold is None:
        quality_threshold = 0
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    try:
        file = open('filtered/' + output_fastq)
    except IOError:
        pass
    else:
        print('Файл уже существует')
        return
    for seq_record in SeqIO.parse(input_fastq, "fastq"):
        seq_name: str = seq_record.id
        seq: str = seq_record.seq
        quality: str = seq_record.letter_annotations['phred_quality']
        gc_content: float = SeqUtils.GC(seq)
        threshold: int = np.mean(quality)
        if is_bounded(length_bounds, len(seq)) and is_bounded(gc_bounds,
                                                              gc_content) and threshold >= quality_threshold:
            os.makedirs('filtered', exist_ok=True)
            file = open('filtered/' + output_fastq, 'a+')
            file.write(seq_name + '\n' + seq + '\n+\n' + quality + '\n')
            file.close()


# def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
#     """
#     Function run_dna_rna_tools
#     Application - Several string values are supplied to the function input, the last of which is the type of operation applied to the sequences.
#     This function returns the result of the applied operation to the sequence, if it is possible to perform it, otherwise it is printed:
#         1. error_string = 'Invalid input: string must be either DNA or RNA', if the sequence is not DNA or RNA
#         2. error_string_rna = "The string must be DNA. Given an RNA string", if the operation transcribe_dna_complement is attempted
#         to RNA sequence.
#
#     Attributes
#     ----------
#     *args : str или List[str]
#         args - a tuple containing all sequences of input data and the operation applied to them.
#     """
#     operator: str = args[-1]
#     seqs: Union[str, List[str]] = list(args[:-1])
#     output: List[str] = []
#     seqs: List[str] = list(seqs)
#     correct_seqs: List[str] = []
#     type_: List[int] = []
#     for ind in range(len(seqs)):
#         letters: Set[str] = set(seqs[ind])
#         checker = all(ch in "AUTGCautgc" for ch in letters) and not (
#                 any(ch in 'Tt' for ch in letters) and any(ch in 'Uu' for ch in letters))
#         if not checker:
#             correct_seqs.append(utils.error_string)
#             continue
#         correct_seqs.append(seqs[ind])
#         if 'U' in letters or 'u' in letters:
#             type_.append(0)
#             continue
#         type_.append(1)
#     for seq in range(len(correct_seqs)):
#         if correct_seqs[seq] != utils.error_string:
#             if operator == 'transcribe_dna_complement':
#                 if not type_[seq]:
#                     output.append(utils.error_string_rna)
#                     continue
#             if correct_seqs[seq] != utils.error_string:
#                 output.append(operations[operator](correct_seqs[seq], type_[seq]))
#         else:
#             output.append(correct_seqs[seq])
#
#     if len(output) > 1:
#         return output
#     return output[0]


class BiologicalSequence(ABC):
    def __init__(self, sequence: str, alphabet: str):
        if not set(sequence).issubset(set(alphabet)):
            raise TypeError("Invalid character")
        self.sequence = sequence
        self.alphabet = alphabet

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index) -> str:
        return self.sequence[index]

    def __str__(self) -> str:
        return self.sequence

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def validate(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    @abstractmethod
    def __init__(self, sequence: str, alphabet: str, complement_map: dict):
        super().__init__(sequence, alphabet)
        self.type = self.type_sequence()
        self.complement_map = complement_map

    def type_sequence(self) -> int:
        return 0 if 'U' in self.sequence or 'u' in self.sequence else 1

    def complement(self):
        complement_seq = "".join(self.complement_map[base] for base in self.sequence)
        return self.__class__(complement_seq)

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    COMPLEMENT = {'A': "T", "T": "A", 'C': 'G', 'G': 'C'}
    DNA = "ATGCatgc"

    def __init__(self, sequence: str):
        super().__init__(sequence, self.DNA, self.COMPLEMENT)
        self.validate()

    def validate(self):
        if "U" in self.sequence or "u" in self.sequence:
            raise ValueError("Invalid DNA sequence")

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    COMPLEMENT = {"A": "U", "U": "A", "C": "G", "G": "C"}
    RNA = "AUGCaugc"

    def __init__(self, sequence: str):
        super().__init__(sequence, self.RNA, self.COMPLEMENT)
        self.validate()

    def validate(self):
        if "T" in self.sequence or "t" in self.sequence:
            raise ValueError("Invalid RNA sequence")


class AMINOSequence(BiologicalSequence):
    AMINO = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, sequence: str):
        super().__init__(sequence, self.AMINO)
        self.validate()

    def validate(self):
        if not set(self.sequence).issubset(set(self.AMINO)):
            raise ValueError("Invalid Amino Acid sequence")
