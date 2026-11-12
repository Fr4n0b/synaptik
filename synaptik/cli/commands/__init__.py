"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import benchmark, deposit, design


def register_all(sub, common) -> None:
    for module in (deposit, design, benchmark):
        module.register(sub, common)
