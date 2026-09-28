from pathlib import Path
import pytest
from atlas_helpers.purity import read_purity_ploidy, PurityParseError

FIX = Path(__file__).parent / "fixtures"


def test_reads_conditional_fixture():
    fit = read_purity_ploidy(FIX / "purity_conditional.tsv")
    assert fit.purity == pytest.approx(0.31)
    assert fit.ploidy == pytest.approx(3.9)
    assert fit.status == "NORMAL"
    assert fit.sample_sex == "male"       # 'MALE' -> lowercased
    assert fit.ploidy_int() == 4          # round(3.9)


def test_ploidy_int_floor():
    fit = read_purity_ploidy(FIX / "purity_low_ploidy.tsv")  # ploidy 0.4
    assert fit.ploidy_int() == 1          # never below 1


def test_bad_file_raises(tmp_path):
    p = tmp_path / "bad.tsv"
    p.write_text("purity\tploidy\nNA\tNA\n")
    with pytest.raises(PurityParseError):
        read_purity_ploidy(p)


@pytest.mark.parametrize("bad", ["inf", "-inf", "Infinity", "nan"])
def test_non_finite_raises(tmp_path, bad):
    p = tmp_path / "nf.tsv"
    p.write_text(f"purity\tploidy\n{bad}\t2.0\n")
    with pytest.raises(PurityParseError):
        read_purity_ploidy(p)
