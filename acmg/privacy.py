"""Privacy defaults for the per-case drivers: which input file is read, and whether clinical text leaves the machine.

Two rules, shared by examples/run_case.py, scripts/run_proband.py and scripts/run_all_probands.py:

1. Input. The drivers default to the DE-IDENTIFIED student file `dataset.parquet`. The curated file
   `dataset_clinical_curated.parquet` carries raw PII (clinician names, DOB, accession IDs) in
   `clinical_indication_text` (see scripts/prepare_bundle.py); reading it needs `--allow-pii-input`.
   The variant calls are identical in both files, so nothing is lost by using the de-identified one.
2. External LLM. HPO extraction defaults to `tool_only` (FastHPOCR + regex negation/scope, fully local).
   Any mode that sends the clinical-indication text to an LLM provider through the `pi` CLI (the LLM HPO modes,
   and the --agent arm) needs `--allow-external-llm`.

Whenever either opt-in is used, a one-line warning is printed to stderr.
"""
from __future__ import annotations
import argparse
import os
import sys

DEID_BUNDLE = "/root/bioconnect/dataset.parquet"
OFFLINE_HPO_MODE = "tool_only"
HPO_MODE_CHOICES = ("tool_only", "augment_select", "candidates_model", "model_only")


def is_pii_bundle(path: str) -> bool:
    """True for the raw curated file (`dataset_clinical_curated.parquet`). A `.deid.` copy is not PII."""
    name = os.path.basename(str(path)).lower()
    return "clinical_curated" in name and ".deid." not in name


def check_bundle(path: str, allow_pii_input: bool) -> str:
    """Refuse the PII file unless explicitly allowed; warn (one line, stderr) when it is used. Returns `path`."""
    if is_pii_bundle(path):
        if not allow_pii_input:
            raise SystemExit(
                f"refusing to read {path}: dataset_clinical_curated.parquet carries raw PII in "
                "clinical_indication_text. Use the de-identified dataset.parquet (the variant calls are identical), "
                "or pass --allow-pii-input if you really mean to read it.")
        print(f"WARNING: reading the PII file {os.path.basename(path)} (--allow-pii-input); "
              "do not share its outputs or send its text to external services.", file=sys.stderr)
    return path


def uses_external_llm(hpo_mode: str) -> bool:
    return hpo_mode != OFFLINE_HPO_MODE


def check_external_llm(hpo_mode: str = OFFLINE_HPO_MODE, allow_external_llm: bool = False, *,
                       agent: bool = False, allow_pii_input: bool = False) -> None:
    """Refuse any path that sends clinical text to an external LLM unless `--allow-external-llm` was passed;
    warn (one line, stderr) when it is used. PII text is never sent out, even with both opt-ins."""
    if hpo_mode not in HPO_MODE_CHOICES:
        raise SystemExit(f"--hpo-mode must be one of {HPO_MODE_CHOICES}")
    what = [f"--hpo-mode {hpo_mode}"] if uses_external_llm(hpo_mode) else []
    if agent:
        what.append("--agent")
    if not what:
        return
    if not allow_external_llm:
        raise SystemExit(
            f"{' and '.join(what)} sends the case's clinical-indication text to an external LLM provider via the "
            f"`pi` CLI. Pass --allow-external-llm to opt in (de-identified text only), or use the offline default "
            f"--hpo-mode {OFFLINE_HPO_MODE}.")
    if allow_pii_input:
        raise SystemExit("--allow-external-llm cannot be combined with --allow-pii-input: "
                         "raw PII text must not be sent to an external LLM.")
    print(f"WARNING: {' and '.join(what)} sends clinical-indication text to an external LLM provider "
          "(--allow-external-llm).", file=sys.stderr)


def add_privacy_args(ap: argparse.ArgumentParser, *, pii_input: bool = True) -> None:
    """Add --hpo-mode / --allow-external-llm (and, if `pii_input`, --allow-pii-input) to a driver's parser."""
    ap.add_argument("--hpo-mode", choices=HPO_MODE_CHOICES, default=OFFLINE_HPO_MODE,
                    help=f"HPO extraction mode (default {OFFLINE_HPO_MODE}: local FastHPOCR, no network). "
                         "The other modes send clinical text to an external LLM and need --allow-external-llm.")
    ap.add_argument("--allow-external-llm", action="store_true",
                    help="opt in to sending de-identified clinical-indication text to an LLM provider via `pi`")
    if pii_input:
        ap.add_argument("--allow-pii-input", action="store_true",
                        help="opt in to reading dataset_clinical_curated.parquet (raw PII); off by default")
