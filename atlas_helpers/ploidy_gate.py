from dataclasses import dataclass
from enum import Enum


class PloidyConfigError(ValueError):
    ...


class PloidyMode(str, Enum):
    NONE = "none"
    STATIC = "static"
    CONDITIONAL = "conditional"


@dataclass(frozen=True)
class PloidyDecision:
    mode: PloidyMode
    first_pass_args: list[str]
    conditional: bool
    cap_value: int
    threshold: float | None

    def needs_rerun(self, fitted_purity: float) -> bool:
        return (
            self.conditional
            and self.threshold is not None
            and fitted_purity < self.threshold
        )

    def rerun_args(self) -> list[str]:
        return ["-max_ploidy", str(self.cap_value)]


def decide(max_ploidy, threshold, cap_value=2) -> PloidyDecision:
    if max_ploidy is not None and threshold is not None:
        raise PloidyConfigError(
            "max_ploidy and ploidy_cap_purity_threshold are mutually exclusive"
        )
    if max_ploidy is not None:
        if max_ploidy < 1:
            raise PloidyConfigError("max_ploidy must be >= 1")
        return PloidyDecision(
            PloidyMode.STATIC, ["-max_ploidy", str(max_ploidy)], False, cap_value, None
        )
    if threshold is not None:
        if cap_value < 1:
            raise PloidyConfigError("ploidy_cap_value must be >= 1")
        if not (0.0 < float(threshold) <= 1.0):
            raise PloidyConfigError("ploidy_cap_purity_threshold must be in (0, 1]")
        return PloidyDecision(
            PloidyMode.CONDITIONAL, [], True, cap_value, float(threshold)
        )
    return PloidyDecision(PloidyMode.NONE, [], False, cap_value, None)
