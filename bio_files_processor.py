from typing import List


def parse_blast_output(input_file: str, output_file: str):
    """
    Function parse_blast_output
    Application - The input of the function is the name of the file with data in blast format and for each Query the first line
    from the Description column is written to the final file with the name output_file.
    Attributes
    ----------
    input_file: str
        Input file name
    output_file: str
        Output file name
    """
    output: List[str] = []
    with open(input_file) as file:
        row: str = file.readline()
        while row:
            if row[0:7] == 'Query #':
                while row[0:11] != 'Description':
                    row = file.readline()
                file.readline()
                row = file.readline().split('\t')[0][:66]
                output.append(row)
            row = file.readline()
    sorted(output, key=lambda x: x.lower())
    with open(output_file, 'w+') as file:
        for row in output:
            file.write(row + '\n')


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> None:
    """
    Function convert_multiline_fasta_to_oneline
    Application - The input of the function is the name of the file with fasta format data and each sequence that
    can be split into several lines, read and written to a new file called output_fastq, while
    all sequences are written on one line.
    Attributes
    ----------
    input_fasta: str
        Input file name
    output_fasta: str
        Output file name
    """
    output: List = []
    with open(input_fasta) as file:
        line: str = file.readline()
        row: List = []
        while line:
            if line[0] == '>':
                row.append(line)
            else:
                row.append(line.rstrip('\n'))
            line = file.readline()
            if len(line) == 0:
                output.append(row)
                break
            if line[0] == '>':
                row.append('\n')
                output.append(row)
                row = []
    with open(output_fasta, 'w+') as f:
        for row in output:
            f.write(''.join(row))


print(parse_blast_output('test.txt', 'fuck.txt'))
