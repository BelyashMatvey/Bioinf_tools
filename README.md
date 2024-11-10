# Module Bioinf_tools
## Annotation
The **Bioinf_tools** module contains all the functions necessary for a novice bioinformatician.
It allows basic conversion of nucleic acid sequences, as well as filtering 
data in _fastq_ format.

## Contents
- [Functions](#main-functions)
- [Data](#data)
- [Authors](#authors)

<a name="main-functions"><h2>Functions</h2></a>
### filter_fastq
**filter_fastq** accepts a file with a sequence in fastq format as input, filters the sequences by parameters 
_gc_bounds_ (minimum and maximum percentage of G and C nucleotides), _length_bounds_ (minimum and maximum
sequence length), _quality_thresholds_ (minimum average read quality).The function creates a filtered folder (if
it doesn’t exist yet), it creates a file called output_file and writes filtered fastq sequences to it. 
If the file already exists, the function does not overwrite it and displays an error. To write to a file, use the module
_write_fastq_ from the _modules_ package.

### run_dna_rna_tools
**run_dna_rna_tools** takes as input several strings of nucleic acids and an operation, 
which needs to be done on them. Returns the result of an operation for each given sequence.
The function uses the _modules_ package containing the _utils_ module for the corresponding operations. If the operation cannot be
a given nucleotide sequence is used or the input string is not a nucleotide sequence -
the corresponding error is displayed.
#### Возможные операции:
 - **complement** - returns the complementary sequence of the given one
 - **reverse** - returns the unwrapped sequence
 - **transcribe** - returns the transcribed sequence
 - **reverse_complement** - returns the reverse complement sequence
 - **transcribe_dna_complement** - returns the complementary sequence for transcription of a given sequence

### parse_blast_output
**parse_blast_output** takes two arguments as input: _input_file_ and _output_file_. The function takes as input a file with 
Blast data with blast format data and for each 'Query' the first row from the 'Description' column is written to the final 
file called output_file.

### convert_multiline_fasta_to_oneline
**convert_multiline_fasta_to_oneline** takes as input the name of the file with fasta format data and each sequence, 
which can be split into multiple lines, is read from the file and written to a new file called output_fastq, 
and all sequences will be written on one line.

<a name="data"><h2>Data</h2></a>
The _modules_ package contains the _utils_ module, which contains dictionaries for operations with nucleic acids and the texts of the corresponding 
errors **const_dicts_and_strs**, which contains the necessary dictionaries for operations on sequences, 
as well as texts of possible errors.

<a name="authors"><h2>Authors</h2></a>
### **Belyakov Matvey** 

