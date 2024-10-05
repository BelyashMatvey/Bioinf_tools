from typing import Dict, Union, List, Set

from modules import transcribe, transcribe_dna_complement, reverse, reverse_complement, complement
from modules.const_dicts_and_strs import error_string, error_string_rna


def filter_fastq(seqs: Dict[str, tuple[str, str]], gc_bounds: Union[int, tuple[int, int]] = None,
                 length_bounds: Union[int, tuple[int, int]] = None, quality_threshold: int = None) -> Dict[str, tuple[str, str]]:
    """
    Функция filter_fastq
    Применение - На вход функции подается словарь, где ключ - название последовательности, а значение - список вида
    [последовательность,качество прочтения последовательности]. Используя параметры gc_bounds, length_bounds, и
    quality_threshold функция фильтрует значения словаря и выводит только те элементы словаря, которые подходят под
    заданные условия.
    Аттрибуты
    ----------
    seqs : Dict[str, tuple[str, str]]
        seqs - словарь, состоящий из fastq-сиквенсов. Структура следующая. Ключ - строка, имя последовательности.
        Значение - кортеж из двух строк: последовательность и качество.
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
    result_seqs: Dict[str, tuple[str, str]] = {}
    for name, fastq in seqs.items():
        if isinstance(gc_bounds, int):
            gc_bounds = (0, gc_bounds)
        if isinstance(length_bounds, int):
            length_bounds = (0, length_bounds)
        gc_content: float = (fastq[0].count('G') + fastq[0].count('C')) / len(fastq[0]) * 100
        threshold: int = 0
        for char in fastq[1]:
            threshold += ord(char) - 33
        threshold /= len(fastq[0])
        if length_bounds[0] <= len(fastq[0]) <= length_bounds[1] and gc_bounds[0] <= gc_content <= gc_bounds[1] \
                and threshold >= quality_threshold:
            result_seqs[name] = fastq

    return result_seqs


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Функция run_dna_rna_tools
    Применение - На вход функции подается несколько строковых значений, последнее из которых - тип применяемой операции к последовательностям.
    Данная функция возвращает результат применяемой операции к последовательности, если её возможно совершить, в противном случае выводится:
        1. error_string = 'Некорректный ввод: строка должна быть либо ДНК, либо РНК', если последовательность не является ДНК или РНК
        2. error_string_rna = "Строка должна быть ДНК. Дана строка РНК", если операцию transcribe_dna_complement пытаются применить
        к РНК последовательности.

    Аттрибуты
    ----------
    *args : str или List[str]
        args - кортеж, содержащий в себе все последоватности входных данных и операцию, применяемую к ним.
    """
    operator: str = args[-1]
    seqs: Union[str, List[str]] = list(args[:-1])
    output: List[str] = []
    seqs: List[str] = list(seqs)
    correct_seqs: List[str] = []
    type_: List[int] = []
    for ind in range(len(seqs)):
        letters: Set[str] = set(seqs[ind])
        checker = all(ch in "AUTGCautgc" for ch in letters) and any(ch1 not in letters for ch1 in "TUtu")
        if not checker:
            correct_seqs.append(error_string)
            continue
        correct_seqs.append(seqs[ind])
        if 'U' in letters or 'u' in letters:
            type_.append(0)
            continue
        type_.append(1)
    if operator == "transcribe":
        for seq in range(len(correct_seqs)):
            if correct_seqs[seq] != error_string:
                output.append(transcribe.transcribe(correct_seqs[seq], type_[seq]))
            else:
                output.append(correct_seqs[seq])
    elif operator == 'reverse':
        for seq in range(len(correct_seqs)):
            if correct_seqs[seq] != error_string:
                output.append(reverse.reverse(correct_seqs[seq]))
            else:
                output.append(correct_seqs[seq])
    elif operator == "complement":
        for seq in range(len(correct_seqs)):
            if correct_seqs[seq] != error_string:
                output.append(complement.complement(correct_seqs[seq], type_[seq]))
            else:
                output.append(correct_seqs[seq])
    elif operator == "transcribe_dna_complement":
        for seq in range(len(correct_seqs)):
            if correct_seqs[seq] != error_string:
                if not type_[seq]:
                    output.append(error_string_rna)
                    continue
                output.append(transcribe_dna_complement.transcribe_dna_complement(correct_seqs[seq]))
            else:
                output.append(correct_seqs[seq])
    else:
        for seq in range(len(correct_seqs)):
            if correct_seqs[seq] != error_string:
                output.append(reverse_complement.reverse_complement(correct_seqs[seq], type_[seq]))
            else:
                output.append(correct_seqs[seq])
    if len(output) > 1:
        return output
    return output[0]
