"""bioconnect-sprint-py: a transparent, offline ACMG classifier.

The ACMG verdict is what GeneBe gives you and VEP does not. This package supplies it: VEP annotations
-> `vep_to_annotations` -> `classify` (the SQL kernel). Reproducible, no API quota, no liftover.
"""
from .kernel import classify
from .vep_map import vep_to_annotations, REQUIRED_COLS, SO_TO_KERNEL

__all__ = ["classify", "vep_to_annotations", "REQUIRED_COLS", "SO_TO_KERNEL"]
