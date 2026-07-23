import pytest
from atlas_helpers.chr_strip import strip_chrom, strip_file, ChrColumnError


def test_strip_chrom():
    assert strip_chrom("chr1") == "1"
    assert strip_chrom("chr22") == "22"
    assert strip_chrom("chrX") == "X"
    assert strip_chrom("chrY") == "Y"
    assert strip_chrom("chrM") == "MT"
    assert strip_chrom("chrMT") == "MT"
    assert strip_chrom("7") == "7"          # already Ensembl -> unchanged (idempotent)


def test_strip_file_chromosome_column(tmp_path):
    p, o = tmp_path / "in.tsv", tmp_path / "out.tsv"
    p.write_text("chromosome\tstart\tend\tcopyNumber\nchr7\t1\t100\t3\nchrMT\t1\t50\t2\n")
    assert strip_file(p, o) == 2
    lines = o.read_text().splitlines()
    assert lines[0] == "chromosome\tstart\tend\tcopyNumber"   # header preserved
    assert lines[1].split("\t")[0] == "7"
    assert lines[2].split("\t")[0] == "MT"                    # chrMT -> MT
    assert lines[1].split("\t")[1:] == ["1", "100", "3"]      # other columns unchanged


def test_only_chromosome_column_changes(tmp_path):
    # a gene column with chr-like text must NOT be altered
    p, o = tmp_path / "g.tsv", tmp_path / "g.out.tsv"
    p.write_text("gene\tchromosome\tlog2\nCHR7orf\tchr7\t0.5\n")
    strip_file(p, o)
    row = o.read_text().splitlines()[1].split("\t")
    assert row == ["CHR7orf", "7", "0.5"]


def test_header_only_passthrough(tmp_path):
    p, o = tmp_path / "h.tsv", tmp_path / "h.out.tsv"
    p.write_text("chromosome\tstart\tend\n")
    assert strip_file(p, o) == 0
    assert o.read_text().splitlines() == ["chromosome\tstart\tend"]


def test_no_chromosome_column_errors(tmp_path):
    p, o = tmp_path / "x.tsv", tmp_path / "y.tsv"
    p.write_text("foo\tbar\n1\t2\n")
    with pytest.raises(ChrColumnError):
        strip_file(p, o)


def test_crlf_terminators_preserved(tmp_path):
    # CRLF line endings must survive verbatim (line-preserving rewrite)
    p, o = tmp_path / "in.tsv", tmp_path / "out.tsv"
    p.write_bytes(b"chromosome\tstart\r\nchr7\t1\r\nchrMT\t2\r\n")
    assert strip_file(p, o) == 2
    data = o.read_bytes()
    assert data == b"chromosome\tstart\r\n7\t1\r\nMT\t2\r\n"
