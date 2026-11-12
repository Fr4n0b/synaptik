"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import aso, benchmark, biomarker, deposit, design, simulate, translate


def register_all(sub, common) -> None:
    for module in (deposit, design, benchmark, aso, translate, simulate, biomarker):
        module.register(sub, common)
