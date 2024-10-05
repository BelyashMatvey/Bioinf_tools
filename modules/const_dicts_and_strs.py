dna_to_rna = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 't': 'a', 'g': 'c', 'c': 'g'}
dna_to_dna = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
rna_to_rna = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
rna_to_dna = {'A': 'T', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 't': 'a', 'g': 'c', 'c': 'g'}
error_string = 'Некорректный ввод: строка должна быть либо ДНК, либо РНК'
error_string_rna = "Строка должна быть ДНК. Дана строка РНК"
transcr_map_dna = {'A': 'A', 'G': 'G', 'C': 'C', 'T': 'U', 'a': 'a', 'g': 'g', 'c': 'c', 't': 'u'}
transcr_map_rna = {'A': 'A', 'G': 'G', 'C': 'C', 'U': 'T', 'a': 'a', 'g': 'g', 'c': 'c', 'u': 't'}
