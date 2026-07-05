"""The literature arm, live: for a 'mixed' candidate VUS, retrieve real literature (Europe PMC) and let
pi PROPOSE an ACMG literature criterion over it — recorded with a digest, GATED behind human approval.

This is what should happen to the 'mixed / adjudicate' candidates from run_case.py: the agent gathers
evidence and proposes; it never moves the class on its own.

Run:  python examples/literature_arm.py
"""
import json
from acmg import evidence, agent, pi

# a 'mixed' candidate from CASE0003 (phenotype-matched VUS). RSID optional (unlocks Ensembl Variation + LitVar2).
CANDIDATE = {"variant_key": "20-49513169-G-T", "gene": "ADNP", "acmg_class": "VUS"}
HPO_LABELS = ["autism", "intellectual disability", "developmental delay", "seizures"]

print(f"[1/3] gather multi-source evidence for {CANDIDATE['gene']} (Ensembl Variation / LitVar2 / MARRVEL / "
      f"Europe PMC / DECIPHER) ...")
context = {**CANDIDATE, "phenotype": HPO_LABELS,
           **evidence.gather(variant_key=CANDIDATE["variant_key"], gene=CANDIDATE["gene"], hpo_labels=HPO_LABELS)}
print(f"      sources: {[k for k in context if k not in ('variant_key','gene','acmg_class','phenotype')]}")


def render(ctx: dict) -> str:
    return (
        "You are a clinical-genomics reviewer applying ACMG/AMP + ClinGen rules. Using ONLY the abstracts "
        "below, state which literature-dependent criterion (PS3 functional / PP4 phenotype-specificity / "
        "PM3 in-trans / PP1 cosegregation / PS4 prevalence) applies to this candidate, at what strength, "
        "with a one-line rationale and the supporting PMID. If none is supported, say so. Do NOT assign a "
        "final classification.\n\n" + json.dumps(ctx, indent=2)
    )


print(f"\n[2/3] pi proposes (recorded + gated) ...")
if not pi.available():
    raise SystemExit("pi CLI not found; wire any prompt->text callable into agent.propose instead")
# full persistence: every run (context, prompt, model, response, digest) + the decision are durably logged
ledger = agent.AgentLedger("agent_runs.jsonl")
proposal = agent.propose(context, llm=pi.as_llm(provider="openai-codex", model="gpt-5.3-codex-spark", timeout=180),
                         kind="literature-criterion", render=render, model="openai-codex/gpt-5.3-codex-spark", ledger=ledger)
print(f"      status: {proposal['status']}   digest: {proposal['digest'][:23]}…")
print("      model proposal:\n" + "\n".join("        " + l for l in proposal["response"].splitlines()[:12]))

print(f"\n[3/3] a curator approves this exact digest (a separate, attested decision) ...")
approved = agent.approve(proposal, approver="curator", ledger=ledger)
print(f"      status: {approved['status']}  approved_by: {approved['approved_by']}  "
      f"digest-bound: {approved['approved_digest'] == proposal['digest']}")
print(f"      persisted -> {ledger.path} ({len(ledger.records())} immutable events; query as SQL via read_json_auto)")
print("\nOnly now would this criterion count toward the class (via the kernel), never automatically.")
