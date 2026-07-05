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
    if chrom_column and chrom_column in header:
        return header.index(chrom_column)
    for name in _CHROM_HEADERS:
        if name in header:
            return header.index(name)
    return None


def strip_file(in_path, out_path, chrom_column=None) -> int:
    """Rewrite ONLY the chromosome column of a TSV, preserving every other field and the
    original line terminators exactly. Returns the number of data rows rewritten."""
    # newline="" disables universal-newline translation so CRLF/CR terminators survive verbatim.
    with open(in_path, "r", newline="") as fh:
        lines = fh.read().splitlines(keepends=True)
    if not lines:
        with open(out_path, "w", newline="") as fh:
            fh.write("")
        return 0
    idx = _chrom_index(lines[0].rstrip("\r\n").split("\t"), chrom_column)
    if idx is None:
        raise ChrColumnError(f"no chromosome column in {in_path}: {lines[0]!r}")
    out, n = [lines[0]], 0                       # header passed through verbatim
    for line in lines[1:]:
        body = line.rstrip("\r\n")
        term = line[len(body):]                  # keep original terminator
        fields = body.split("\t")
        if len(fields) > idx and fields[idx]:
            fields[idx] = strip_chrom(fields[idx])
            n += 1
        out.append("\t".join(fields) + term)
    with open(out_path, "w", newline="") as fh:
        fh.write("".join(out))
    return n
