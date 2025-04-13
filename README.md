# Module Bioinf_tools
## Annotation
The **Bioinf_tools** module contains all the functions necessary for a novice bioinformatician.
It allows basic conversion of nucleic acid sequences, as well as filtering 
data in _fastq_ format.

## Contents
- [Functions](#main-functions)
- [Classes](#classes)
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

<a name="classes"><h2>Classes</h2></a>
### BiologicalSequence
**BiologicalSequence** - abstract class, that contains 4 dunder methods and one method:
1. '__len__'(self) - returns the length of sequence
2. '__getitem__'(self, index) - returns character from sequence  or slice
3. '__str__'(self) - returns sequence
4. '__repr__'(self) - returns string representation of class
5. validate - abstract method that returns 

### NucleicAcidSequence
**NucleicAcidSequence** - child class of BiologicalSequence. Contains 4 methods:
1. type_sequence(self) - returns 0 if sequence is RNA and 1 if DNA
2. complement(self) - returns complement sequence of given one.
3. reverse(self) - returns reversed sequence.
4. reverse_complement(self) - returns reversed and complement sequence.
 
### DNASequence
**DNASequence** - child class of NucleicAcidSequence. Contains two attributes and two classes:
1. DNA - contains DNA alphabet (A,T,G,C)
2. COMPLEMENT - Dictionary of complement nucleotides
3. validate - raise ValueError if DNA sequence is incorrect
4. transcribe - returns RNA sequence after transcribing DNA.

### RNASequence
**RNASequence** - child class of NucleicAcidSequence. Contains two attributes and one classes:
1. DNA - contains DNA alphabet (A,U,G,C)
2. COMPLEMENT - Dictionary of complement nucleotides
3. validate - raise ValueError if RNA sequence is incorrect

### AMINOSequence
**AMINOSequence - child class of BiologicalSequence. Contains 1 method and 1 attribute: 
1. AMINO - contains AMINO alphabet 
2. validate - raise ValueError if amino acid  sequence is incorrect

<a name="data"><h2>Data</h2></a>
The _modules_ package contains the _utils_ module, which contains dictionaries for operations with nucleic acids and the texts of the corresponding 
errors **const_dicts_and_strs**, which contains the necessary dictionaries for operations on sequences, 
as well as texts of possible errors.

<a name="authors"><h2>Authors</h2></a>
### **Belyakov Matvey** 

