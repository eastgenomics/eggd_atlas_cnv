import csv
import math
from dataclasses import dataclass
from pathlib import Path


class PurityParseError(ValueError):
    ...


@dataclass(frozen=True)
class PurityFit:
    purity: float
    ploidy: float
    status: str
    sample_sex: str

    def ploidy_int(self) -> int:
        return max(1, round(self.ploidy))


def read_purity_ploidy(path) -> PurityFit:
    with open(Path(path)) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise PurityParseError(f"no data rows in {path}")
    d = rows[0]
    try:
        purity = float(d["purity"])
        ploidy = float(d["ploidy"])
    except (KeyError, ValueError) as e:
        raise PurityParseError(str(e)) from e
    if math.isnan(purity) or math.isnan(ploidy):  # blank/NA/NaN => broken PURPLE run
        raise PurityParseError(f"non-finite purity/ploidy in {path}")
    sex = d.get("gender", "").strip().lower()
    return PurityFit(purity, ploidy, d.get("status", ""), sex)
