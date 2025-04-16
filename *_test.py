import os
import pytest
import tempfile

from anaconda_project.internal.conda_api import result
from click.testing import CliRunner
from Bio import SeqIO
import bioinf
from bioinf import filter_fastq


def create_fastq_file(contents: str) -> str:
    """Вспомогательная функция для создания временного fastq-файла"""
    fd, path = tempfile.mkstemp(suffix=".fastq")
    with os.fdopen(fd, 'w') as tmp:
        tmp.write(contents)
    return path


class TestValidFiltering:
    def test_filter_good_read(self):
        content = "@read1\nATGC\n+\nIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq,
                                   [input_path, "out.fastq", "--gc_bounds", "40,60", "--length_bounds", "2,10",
                                    "--quality_threshold", "30"], catch_exceptions=False)
            assert result.exit_code == 0
            assert os.path.exists("filtered/out.fastq")
            with open("filtered/out.fastq") as f:
                lines = f.readlines()
                assert "@read1" in lines[0]

    def test_filter_low_quality(self):
        content = "@read1\nATGC\n+\n!!!!\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq,
                                   [input_path, "out.fastq", "--quality_threshold", "10"])
            assert result.exit_code == 0
            assert not os.path.exists("filtered/out.fastq")

    def test_gc_bounds_one(self):
        content = "@read1\nGCGA\n+\nIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq,
                                   [input_path, "out.fastq", "--gc_bounds", "90"])
            print("OUTPUT:\n", result.output)
            print("EXCEPTION:\n", result.exception)
            assert result.exit_code == 0
            assert os.path.exists("filtered/out.fastq")

    def test_length_bounds(self):
        content = "@read1\nATGCATGC\n+\nIIIIIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq,
                                   [input_path, "out.fastq", "--length_bounds", "10"])
            assert result.exit_code == 0
            assert os.path.exists("filtered/out.fastq")


class TestErrors:
    def test_existing_out(self):
        content = "@read1\nATGC\n+\nIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            os.makedirs("filtered")
            with open("filtered/out.fastq", "w") as f:
                f.write("place")
            result = runner.invoke(filter_fastq, [input_path, "out.fastq"])
            assert "Файл уже существует" in result.output
            assert result.exit_code == 0

    def test_invalid_gc_format(self):
        content = "@read1\nATGC\n+\nIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq, [input_path, "out.fastq", "--gc_bounds", "10,20,30"])
            assert result.exit_code != 0
            assert "gc_bounds должен быть" in result.output


class TestIO:
    def test_written_correct(self):
        content = "@read1\nATGC\n+\nIIII\n"
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq, [input_path, "result.fastq"])
            assert result.exit_code == 0
            path = "filtered/result.fastq"
            assert os.path.exists(path)
            records = list(SeqIO.parse(path, "fastq"))
            assert len(records) == 1
            assert records[0].id == "read1"

    def test_multiple_reads(self):
        content = (
            "@read1\nATGC\n+\nIIII\n"
            "@read2\nGGGG\n+\nIIII\n"
            "@read3\nAAAA\n+\n!!!!\n"  # низкое качество
        )
        input_path = create_fastq_file(content)
        runner = CliRunner()
        with runner.isolated_filesystem():
            result = runner.invoke(filter_fastq,
                                   [input_path, "filtered.fastq", "--gc_bounds", "0,100", "--quality_threshold", "30"])
            assert result.exit_code == 0
            records = list(SeqIO.parse("filtered/filtered.fastq", "fastq"))
            assert len(records) == 2
