"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import aso, benchmark, biomarker, cite, deposit, design, nof1, protocol, repurpose, simulate, translate, vcf


def register_all(sub, common) -> None:
    for module in (deposit, design, benchmark, aso, translate, simulate, biomarker, protocol, cite, nof1, vcf, repurpose):
        module.register(sub, common)

