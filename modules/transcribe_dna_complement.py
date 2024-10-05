from .complement import complement
from .transcribe import transcribe


def transcribe_dna_complement(seq: str) -> str:
    """
    Функция transcribe_dna_complement
    Применение - транскрибирование входной последовательности ДНК, а затем построение и возврат её комплементарной цепи.
    Для этого используются функции из класса complement и transcribe.
    В функцию  transcribe передается type_=1, поскольку у изначально дана цепь ДНК. В фнукцию complement передается type_=0,
    поскольку туда передается цепь РНК.
    Аттрибуты
    ----------
    seq : str
        seq - последовательность нуклеиновой кислоты
    """
    return complement(transcribe(seq, 1), 0)
