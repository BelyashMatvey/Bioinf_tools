from .const_dicts_and_strs import transcr_map_rna,transcr_map_dna


def transcribe(seq: str, type_: int) -> str:
    """
    Функция transcribe
    Применение - транскрибирование входной последовательности в зависимости от её типа (ДНК/РНК).
    Для этого она использует словари соответствующих остатков из модуля const_dicts_and_strs (transcr_map_dna,transcr_map_rna)

    Аттрибуты
    ----------
    seq : str
        seq - последовательность нуклеиновой кислоты
    type_ : int
        type_ - тип нуклеиной кислоты (1 - ДНК,0 - РНК)
    """
    if type_:
        return ''.join(transcr_map_dna[c] for c in seq)
    return ''.join(transcr_map_rna[c] for c in seq)
