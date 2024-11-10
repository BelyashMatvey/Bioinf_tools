import os


def write_fastq(output_fastq: str, fastq_name: str) -> None:
    """
    Function write_fastq
    Application - The function input is a string that will be written to a fastq format file called fastq_name,
    which will be saved in the filtered folder which will be created if it is not in the current directory.
    Attributes
    ----------
    output_fastq: str
        Line to write to file
    fastq_name: str
        Output file name
    """
    os.makedirs('filtered', exist_ok=True)
    file = open('filtered/' + fastq_name, 'a+')
    file.write(output_fastq)
    file.close()
