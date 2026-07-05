"""A GATED agent step — the Python analog of piknit's 'agent chunk', without the pi/R dependency.

Where a sprint pipeline makes ad-hoc LLM calls (HPO validation from clinical text, or the literature-
dependent ACMG criteria PS3/PP4/PM3), this wraps the call so it is REPRODUCIBLE and GATED, which is the
auditability discipline clinical work needs:

  1. the model PROPOSES (it never edits a classification directly);
  2. the proposal is RECORDED with a content digest (an auditable (input, judgment) artifact);
  3. it does NOT affect the score until a human APPROVES that exact digest.

Bring your own LLM: `llm` is any callable `prompt -> text`. Cache it (memoize on the prompt) so a
notebook/Quarto re-render replays byte-for-byte — the render becomes an integration test, like piknit.

The heavyweight, cross-process version of this (a temporal fact store + activation gate) lives in
pi-bio-agent; this is a lightweight, dependency-free version.
"""
from __future__ import annotations
import hashlib
import json
import os
import time
import uuid
from typing import Callable


def _digest(payload: dict) -> str:
    return "sha256:" + hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


class AgentLedger:
    """Durable, APPEND-ONLY record of agent runs — the auditable (input, judgment) ledger clinical work
    requires. One JSON event per line; IMMUTABLE: an approval/rejection is a NEW event referencing the
    proposal's digest, never an overwrite (Datomic-style). Query it as SQL via DuckDB `read_json_auto`. This
    is the file-based analog of pi-bio-agent's temporal fact store; point it at a shared path for a team log."""

    def __init__(self, path: str):
        self.path = path
        d = os.path.dirname(os.path.abspath(path))
        os.makedirs(d, exist_ok=True)

    def append(self, event: dict) -> dict:
        rec = {"event_id": uuid.uuid4().hex, "ts": time.time(), **event}
        with open(self.path, "a") as f:  # append is atomic per line -> safe for concurrent writers
            f.write(json.dumps(rec, default=str, sort_keys=True) + "\n")
        return rec

    def records(self, con=None):
        """All events as a DataFrame (SQL-queryable). Empty frame if nothing logged yet."""
        import duckdb
        if not os.path.exists(self.path) or os.path.getsize(self.path) == 0:
            import pandas as pd
            return pd.DataFrame()
        con = con or duckdb.connect()
        return con.execute(f"SELECT * FROM read_json_auto('{self.path}')").df()

    def status_of(self, digest: str) -> str:
        """Latest event type bound to a proposal digest: 'proposal' | 'approval' | 'rejection' | 'absent'."""
        if not os.path.exists(self.path):
            return "absent"
        evs = [json.loads(l) for l in open(self.path) if l.strip()]
        evs = [e for e in evs if e.get("digest") == digest]
        return evs[-1]["event"] if evs else "absent"


def propose(context: dict, llm: Callable[[str], str], *, kind: str = "literature-criterion",
            render: Callable[[dict], str] | None = None, model: str | None = None,
            ledger: "AgentLedger | None" = None) -> dict:
    """Run one agent proposal over `context` (e.g. a candidate variant + its deterministic criteria).
    Returns a pending, digest-bound artifact; it changes no score by itself. When `ledger` is given, the full
    run (context, prompt, model, response, digest) is durably recorded as an immutable 'proposal' event."""
    prompt = (render or _default_render)(context)
    response = llm(prompt)
    payload = {"kind": kind, "context": context, "response": response}
    artifact = {"digest": _digest(payload), "status": "pending", "prompt": prompt, "model": model, **payload}
    if ledger is not None:
        ledger.append({"event": "proposal", "digest": artifact["digest"], "kind": kind, "model": model,
                       "prompt": prompt, "response": response, "context": context})
    return artifact


def approve(proposal: dict, approver: str, *, ledger: "AgentLedger | None" = None) -> dict:
    """A separate, attested decision, bound to the proposal's exact digest. When `ledger` is given, the
    approval is durably recorded as an immutable event referencing that digest."""
    if not approver:
        raise ValueError("approve: an approver identity is required")
    decided = {**proposal, "status": "approved", "approved_by": approver, "approved_digest": proposal["digest"]}
    if ledger is not None:
        ledger.append({"event": "approval", "digest": proposal["digest"], "approver": approver,
                       "approved_digest": proposal["digest"], "decision": "approved"})
    return decided


def reject(proposal: dict, approver: str, *, reason: str = "", ledger: "AgentLedger | None" = None) -> dict:
    """Terminal REJECT decision, symmetric to approve and equally recorded — a rejection is auditable too."""
    if not approver:
        raise ValueError("reject: an approver identity is required")
    decided = {**proposal, "status": "rejected", "approved_by": approver, "reject_reason": reason}
    if ledger is not None:
        ledger.append({"event": "rejection", "digest": proposal["digest"], "approver": approver,
                       "decision": "rejected", "reason": reason})
    return decided


def approved_points(proposal: dict, points: float, *, ledger: "AgentLedger | None" = None) -> float:
    """`points` count toward the class ONLY when the proposal is approved for its exact digest. When a
    `ledger` is given it is AUTHORITATIVE (the durable, append-only record must carry an 'approval' event for
    this digest) — an in-memory dict with a fabricated status='approved' cannot inject points (pi review #9).
    Without a ledger, falls back to the dict fields (non-authoritative; for offline/testing only)."""
    digest = proposal.get("digest")
    if ledger is not None:
        return points if (digest and ledger.status_of(digest) == "approval") else 0.0
    ok = proposal.get("status") == "approved" and proposal.get("approved_digest") == digest
    return points if ok else 0.0


def _default_render(context: dict) -> str:
    return (
        "You are a clinical-genomics reviewer. Given the structured evidence below, state whether the "
        "named literature-dependent ACMG criterion applies, with a one-line rationale and a citation. "
        "Do NOT assign a final classification.\n\n" + json.dumps(context, indent=2, default=str)
    )
