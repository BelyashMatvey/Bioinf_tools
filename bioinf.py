from typing import Union
from abc import ABC, abstractmethod
import os
import numpy as np
import click
from Bio import SeqIO, SeqUtils
import logging

logging.basicConfig(
    filename='filter_fastq.log',
    filemode='w',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)


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


@click.command()
@click.argument('input_fastq', type=click.Path(exists=True))
@click.argument('output_fastq', type=str)
@click.option('--gc_bounds', default=None, type=str,
              help="GC content bounds: 'max' или 'min,max'")
@click.option('--length_bounds', default=None, type=str,
              help="Length bounds: 'max' или 'min,max'")
@click.option('--quality_threshold', default=0, type=int,
              help="Минимальное среднее качество рида (по шкале Phred33)")
def filter_fastq(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold) -> None:
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
    logging.info(f"Запуск фильтрации: input={input_fastq},output={output_fastq}, "
                 f"gc_bounds={gc_bounds}, length_bounds={length_bounds}, "
                 f'quality_threshold={quality_threshold}')
    if gc_bounds is None:
        gc_bounds = (0, 100)
    else:
        parts = list(map(int, gc_bounds.split(',')))
        if len(parts) == 1:
            gc_bounds = (0, parts[0])
        elif len(parts) == 2:
            gc_bounds = tuple(parts)
        else:
            logging.error("Неверный формат gc_bounds.")
            raise click.BadParameter("gc_bounds должен быть 'max' или 'min,max'")

    if length_bounds is None:
        length_bounds = (0, 2 ** 32)
    else:
        parts = list(map(int, length_bounds.split(',')))
        if len(parts) == 1:
            length_bounds = (0, parts[0])
        elif len(parts) == 2:
            length_bounds = tuple(parts)
        else:
            raise click.BadParameter("length_bounds должен быть 'max' или 'min,max'")
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
            file.write(f"@{seq_name}\n{str(seq)}\n+\n{''.join(chr(q + 33) for q in quality)}\n")
            file.close()


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
