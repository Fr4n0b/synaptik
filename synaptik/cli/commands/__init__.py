"""Sub-command modules. New commands live here rather than in cli/main.py."""

from . import aso, audit, benchmark, biomarker, biophysics, cite, clinical, data, deposit, design, nof1, protocol, repurpose, simulate, translate, twin, vcf


def register_all(sub, common) -> None:
    for module in (audit, deposit, design, benchmark, aso, biophysics, data, translate, simulate, biomarker, protocol, cite, nof1, vcf, repurpose, clinical, twin):
        module.register(sub, common)

