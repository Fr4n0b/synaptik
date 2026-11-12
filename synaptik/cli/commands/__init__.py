"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import aso, benchmark, deposit, design


def register_all(sub, common) -> None:
    for module in (deposit, design, benchmark, aso):
        module.register(sub, common)
