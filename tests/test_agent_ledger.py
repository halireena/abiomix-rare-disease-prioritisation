"""Agent run persistence: propose/approve/reject durably record immutable (input, judgment) events to an
append-only ledger, digest-bound and SQL-queryable. Run: `PYTHONPATH=. python3 tests/test_agent_ledger.py`."""
import json
import os
import tempfile
from acmg import agent


def _fake_llm(prompt: str) -> str:
    return "PP4 supported (Supporting). PMID: 12345."


def test_propose_and_approve_are_persisted_and_digest_bound():
    with tempfile.TemporaryDirectory() as d:
        led = agent.AgentLedger(os.path.join(d, "runs.jsonl"))
        p = agent.propose({"gene": "SCN1A", "variant_key": "2-100-A-G"}, _fake_llm, model="test-llm", ledger=led)
        assert p["status"] == "pending" and p["model"] == "test-llm"
        assert led.status_of(p["digest"]) == "proposal"

        a = agent.approve(p, approver="curator", ledger=led)
        assert a["status"] == "approved" and a["approved_digest"] == p["digest"]
        assert led.status_of(p["digest"]) == "approval"

        # both events durably on disk, append-only (2 lines), and reference the same digest
        lines = [json.loads(l) for l in open(led.path) if l.strip()]
        assert len(lines) == 2
        assert [l["event"] for l in lines] == ["proposal", "approval"]
        assert all(l["digest"] == p["digest"] for l in lines)
        # the proposal event carries the FULL run (context + prompt + response + model)
        prop = lines[0]
        assert prop["context"]["gene"] == "SCN1A" and prop["response"].startswith("PP4") and prop["model"] == "test-llm"


def test_records_are_sql_queryable():
    with tempfile.TemporaryDirectory() as d:
        led = agent.AgentLedger(os.path.join(d, "runs.jsonl"))
        p = agent.propose({"gene": "TP53"}, _fake_llm, ledger=led)
        agent.approve(p, approver="c", ledger=led)
        df = led.records()
        assert len(df) == 2 and set(df["event"]) == {"proposal", "approval"}


def test_reject_is_recorded_and_terminal():
    with tempfile.TemporaryDirectory() as d:
        led = agent.AgentLedger(os.path.join(d, "runs.jsonl"))
        p = agent.propose({"gene": "BRCA1"}, _fake_llm, ledger=led)
        r = agent.reject(p, approver="curator", reason="variant-nonspecific abstracts", ledger=led)
        assert r["status"] == "rejected"
        assert led.status_of(p["digest"]) == "rejection"
        # a rejected proposal contributes no points
        assert agent.approved_points(r, 2.0) == 0.0


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print("ok", name)
    print("all passed")
