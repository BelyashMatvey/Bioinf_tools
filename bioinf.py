from typing import Union, List, Set

from modules import write_fastq, utils
from modules.utils import is_bounded

operations = {
    'transcribe': utils.transcribe,
    'reverse': utils.reverse,
    'complement': utils.complement,
    'transcribe_dna_complement': utils.transcribe_dna_complement,
    'reverse_complement': utils.reverse_complement
}


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
    with open(input_fastq) as file:
        seq_name: str = file.readline().split()[0]
        seq: str = file.readline().rstrip('\n')
        quality_name: str = file.readline().rstrip('\n')
        quality: str = file.readline().rstrip('\n')
        while seq_name:
            gc_content: float = (seq.count('G') + seq.count('C')) / len(seq) * 100
            threshold: int = 0
            for char in quality:
                threshold += ord(char) - 33
            threshold /= len(seq)
            if is_bounded(length_bounds, len(seq)) and is_bounded(gc_bounds,
                                                                  gc_content) and threshold >= quality_threshold:
                write_fastq.write_fastq(seq_name + '\n' + seq + '\n+\n' + quality + '\n', output_fastq)
            seq_name_valid: List[str] = file.readline().split()
            if len(seq_name_valid):
                seq_name: str = seq_name_valid[0]
            else:
                break
            seq_name = seq_name_valid[0]
            seq = file.readline().rstrip('\n')
            quality_name: str = file.readline().rstrip('\n')
            quality: str = file.readline().rstrip('\n')


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Function run_dna_rna_tools
    Application - Several string values are supplied to the function input, the last of which is the type of operation applied to the sequences.
    This function returns the result of the applied operation to the sequence, if it is possible to perform it, otherwise it is printed:
        1. error_string = 'Invalid input: string must be either DNA or RNA', if the sequence is not DNA or RNA
        2. error_string_rna = "The string must be DNA. Given an RNA string", if the operation transcribe_dna_complement is attempted
        to RNA sequence.

    Attributes
    ----------
    *args : str или List[str]
        args - a tuple containing all sequences of input data and the operation applied to them.
    """
    operator: str = args[-1]
    seqs: Union[str, List[str]] = list(args[:-1])
    output: List[str] = []
    seqs: List[str] = list(seqs)
    correct_seqs: List[str] = []
    type_: List[int] = []
    for ind in range(len(seqs)):
        letters: Set[str] = set(seqs[ind])
        checker = all(ch in "AUTGCautgc" for ch in letters) and not (
                any(ch in 'Tt' for ch in letters) and any(ch in 'Uu' for ch in letters))
        if not checker:
            correct_seqs.append(utils.error_string)
            continue
        correct_seqs.append(seqs[ind])
        if 'U' in letters or 'u' in letters:
            type_.append(0)
            continue
        type_.append(1)
    for seq in range(len(correct_seqs)):
        if correct_seqs[seq] != utils.error_string:
            if operator == 'transcribe_dna_complement':
                if not type_[seq]:
                    output.append(utils.error_string_rna)
                    continue
            if correct_seqs[seq] != utils.error_string:
                output.append(operations[operator](correct_seqs[seq], type_[seq]))
        else:
            output.append(correct_seqs[seq])

    if len(output) > 1:
        return output
    return output[0]
