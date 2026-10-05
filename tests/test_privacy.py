"""Privacy defaults (acmg.privacy + the per-case drivers): de-identified input and no external LLM unless the
user opts in. Offline: no network, no real data, no `pi` call (every refusal happens before any work starts)."""
import argparse
import os
import subprocess
import sys

import pytest

from acmg import privacy

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_pii_bundle_detection():
    assert privacy.is_pii_bundle("/root/bioconnect/dataset_clinical_curated.parquet")
    assert privacy.is_pii_bundle("Dataset_Clinical_Curated.parquet")
    assert not privacy.is_pii_bundle("/root/bioconnect/dataset_clinical_curated.deid.parquet")
    assert not privacy.is_pii_bundle(privacy.DEID_BUNDLE)


def test_check_bundle_refuses_pii_by_default(capsys):
    with pytest.raises(SystemExit, match="--allow-pii-input"):
        privacy.check_bundle("x/dataset_clinical_curated.parquet", allow_pii_input=False)
    assert privacy.check_bundle("x/dataset.parquet", allow_pii_input=False) == "x/dataset.parquet"
    assert capsys.readouterr().err == ""                        # de-identified file: no warning


def test_check_bundle_warns_once_when_pii_allowed(capsys):
    privacy.check_bundle("x/dataset_clinical_curated.parquet", allow_pii_input=True)
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("WARNING") and "PII" in err[0]


def test_default_parser_is_offline_and_deidentified():
    ap = argparse.ArgumentParser()
    privacy.add_privacy_args(ap)
    a = ap.parse_args([])
    assert a.hpo_mode == "tool_only" and not a.allow_external_llm and not a.allow_pii_input
    privacy.check_external_llm(a.hpo_mode, a.allow_external_llm)   # default passes silently


@pytest.mark.parametrize("mode", ["augment_select", "candidates_model", "model_only"])
def test_llm_modes_need_opt_in(mode, capsys):
    with pytest.raises(SystemExit, match="--allow-external-llm"):
        privacy.check_external_llm(mode, False)
    privacy.check_external_llm(mode, True)
    err = capsys.readouterr().err.strip().splitlines()
    assert len(err) == 1 and err[0].startswith("WARNING") and "external LLM" in err[0]


def test_agent_needs_opt_in_even_with_offline_hpo():
    with pytest.raises(SystemExit, match="--agent"):
        privacy.check_external_llm("tool_only", False, agent=True)
    privacy.check_external_llm("tool_only", True, agent=True)


def test_pii_text_never_sent_to_external_llm():
    with pytest.raises(SystemExit, match="cannot be combined"):
        privacy.check_external_llm("augment_select", True, allow_pii_input=True)
    privacy.check_external_llm("tool_only", False, allow_pii_input=True)   # PII locally, offline HPO: allowed


def test_unknown_mode_rejected():
    with pytest.raises(SystemExit):
        privacy.check_external_llm("llm_everything", True)


# --- the drivers themselves: refusals happen at argument time, before any data or network is touched ---

def _run(*args):
    env = {k: v for k, v in os.environ.items() if k != "BUNDLE"}
    env["PYTHONPATH"] = ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run([sys.executable, *args], cwd=ROOT, env=env, capture_output=True, text=True, timeout=120)


def test_run_case_refuses_pii_bundle(tmp_path):
    r = _run("examples/run_case.py", "CASE0003", "--bundle", str(tmp_path / "dataset_clinical_curated.parquet"))
    assert r.returncode != 0 and "refusing to read" in r.stderr


def test_run_case_refuses_llm_without_opt_in():
    r = _run("examples/run_case.py", "CASE0003", "--hpo-mode", "augment_select")
    assert r.returncode != 0 and "--allow-external-llm" in r.stderr


def test_run_case_default_reads_deidentified_bundle(tmp_path):
    missing = tmp_path / "nope" / "dataset.parquet"
    r = _run("examples/run_case.py", "CASE0003", "--bundle", str(missing))
    assert r.returncode != 0 and "needs data that is not here yet" in r.stderr   # got past the privacy gates
    assert "WARNING" not in r.stderr
    r = _run("examples/run_case.py", "--help")
    assert r.returncode == 0 and "dataset.parquet" in r.stdout


def test_run_proband_refuses_llm_without_opt_in():
    r = _run("scripts/run_proband.py", "CASE0001", "--hpo-mode", "candidates_model")
    assert r.returncode != 0 and "--allow-external-llm" in r.stderr


def test_run_all_probands_agent_needs_opt_in():
    r = _run("scripts/run_all_probands.py", "--agent")
    assert r.returncode != 0 and "--allow-external-llm" in r.stderr
