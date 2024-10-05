from typing import Dict

from .const_dicts_and_strs import dna_to_dna, rna_to_rna


def complement(seq: str, type_: int) -> str:
    """Функция complement
    Основное применение - транскрибирование входной последовательности в зависимости от её типа (ДНК/РНК).
    Для этого она использует словари соответствующих остатков из модуля const_dicts_and_strs (dna_to_dna,rna_to_rna)

    Аттрибуты
    ----------
    seq : str
        seq - последовательность нуклеиновой кислоты
    type_ : int
        type_ - тип нуклеиной кислоты (1 - ДНК,0 - РНК)"""
    compl_map: Dict[str, str] = dna_to_dna if type_ else rna_to_rna
    return ''.join(compl_map[c] for c in seq)
