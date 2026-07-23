import pytest
from atlas_helpers.ploidy_gate import decide, PloidyMode, PloidyConfigError


def test_none_mode():
    d = decide(max_ploidy=None, threshold=None, cap_value=2)
    assert d.mode is PloidyMode.NONE
    assert d.first_pass_args == []
    assert d.conditional is False


def test_static_mode():
    d = decide(max_ploidy=3, threshold=None, cap_value=2)
    assert d.mode is PloidyMode.STATIC
    assert d.first_pass_args == ["-max_ploidy", "3"]
    assert d.conditional is False


def test_conditional_mode_and_rerun_boundary():
    d = decide(max_ploidy=None, threshold=0.35, cap_value=2)
    assert d.mode is PloidyMode.CONDITIONAL
    assert d.first_pass_args == []
    assert d.conditional is True
    assert d.needs_rerun(0.34) is True     # below threshold -> cap
    assert d.needs_rerun(0.35) is False    # at threshold -> keep
    assert d.needs_rerun(0.36) is False
    assert d.rerun_args() == ["-max_ploidy", "2"]


def test_mutually_exclusive():
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=3, threshold=0.35, cap_value=2)


@pytest.mark.parametrize("bad", [0, -1])
def test_invalid_values(bad):
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=bad, threshold=None, cap_value=2)
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=None, threshold=0.35, cap_value=bad)


@pytest.mark.parametrize("bad", [-0.1, 0.0, 1.5])
def test_invalid_threshold_range(bad):
    with pytest.raises(PloidyConfigError):
        decide(max_ploidy=None, threshold=bad, cap_value=2)
