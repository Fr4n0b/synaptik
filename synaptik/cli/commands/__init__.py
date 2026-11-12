"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import aso, benchmark, biomarker, cite, deposit, design, protocol, simulate, translate


def register_all(sub, common) -> None:
    for module in (deposit, design, benchmark, aso, translate, simulate, biomarker, protocol, cite):
        module.register(sub, common)
