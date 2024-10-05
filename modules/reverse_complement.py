from .complement import complement
from .reverse import reverse


def reverse_complement(seq: str, type_: int) -> str:
    """Функция transcribe
    Применение - нахождение обратной комплементарной цепи для входной последовательности в зависимости от её типа (ДНК/РНК).
    Для этого она использует функции из соответствующих модулей (complement,reverse)

    Аттрибуты
    ----------
    seq : str
        seq - последовательность нуклеиновой кислоты
    type_ : int
        type_ - тип нуклеиной кислоты (1 - ДНК,0 - РНК)"""
    return reverse(complement(seq, type_))
