class ChrColumnError(ValueError):
    ...


_CHROM_HEADERS = ("chromosome", "chrom", "#chrom", "contig")


def strip_chrom(chrom: str) -> str:
    c = chrom.strip()
    if not c.startswith("chr"):
        return c                       # already Ensembl-style / not prefixed
    rest = c[3:]
    return "MT" if rest in ("M", "MT") else rest


def _chrom_index(header: list[str], chrom_column):
    if chrom_column:
        return header.index(chrom_column) if chrom_column in header else None
    for name in _CHROM_HEADERS:
        if name in header:
            return header.index(name)
    return None


def strip_file(in_path, out_path, chrom_column=None) -> int:
    """Rewrite ONLY the chromosome column of a TSV, preserving every other field and the
    original line terminators exactly. Returns the number of data rows rewritten."""
    n = 0
    # newline="" disables universal-newline translation so CRLF/CR terminators survive verbatim.
    with open(in_path, "r", newline="") as fin, open(out_path, "w", newline="") as fout:
        header_line = fin.readline()
        if not header_line:
            return 0
        idx = _chrom_index(header_line.rstrip("\r\n").split("\t"), chrom_column)
        if idx is None:
            raise ChrColumnError(f"no chromosome column in {in_path}: {header_line!r}")
        fout.write(header_line)                  # header passed through verbatim
        for line in fin:
            body = line.rstrip("\r\n")
            term = line[len(body):]              # keep original terminator
            fields = body.split("\t")
            if len(fields) > idx and fields[idx]:
                fields[idx] = strip_chrom(fields[idx])
                n += 1
            fout.write("\t".join(fields) + term)
    return n
