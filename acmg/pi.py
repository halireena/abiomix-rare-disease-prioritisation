"""A small, real wrapper over the `pi` agent CLI — the Abiomix advantage, native in Python.

This is deliberately thin (piknit does the same for R, but in Python there is no reticulate to fight):
`pi` is a subprocess. text / json / rpc output modes. Use it directly for a one-shot agent call, or as
the LLM backend for the GATED step in acmg.agent (`propose(context, llm=pi.as_llm())`) so the model's
judgment is recorded with a digest and a human approves it before it counts.

Requires the `pi` CLI on PATH and a configured provider. Without pi, pass any other `prompt -> text`
callable to acmg.agent instead; the gating discipline is the same.
"""
from __future__ import annotations
import json
import shutil
import subprocess
from typing import Callable


def available() -> bool:
    return shutil.which("pi") is not None


def run(prompt: str, *, provider: str = "openai-codex", model: str = "gpt-5.3-codex-spark",
        system: str | None = None, tools: str | None = None, mode: str = "text",
        extension: str | None = None, timeout: int = 180) -> str | dict:
    """One-shot, non-interactive pi call. `mode='json'` returns the parsed object; else stripped text."""
    if not available():
        raise RuntimeError("pi CLI not found on PATH")
    cmd = ["pi", "--provider", provider, "--model", model, "--mode", mode]
    if system:
        cmd += ["--system-prompt", system]
    if tools:
        cmd += ["-t", tools]
    if extension:
        cmd += ["-e", extension]
    cmd += ["-p", prompt]  # --print (non-interactive) + the prompt
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"pi failed ({r.returncode}): {r.stderr.strip()[:500]}")
    out = r.stdout.strip()
    return json.loads(out) if mode == "json" else out


def as_llm(**kwargs) -> Callable[[str], str]:
    """A `prompt -> text` callable bound to pi, for `acmg.agent.propose(llm=...)`.

    Memoize the returned callable (e.g. functools.lru_cache) if you want a notebook/Quarto re-render to
    replay byte-for-byte instead of re-querying the model."""
    return lambda prompt: run(prompt, **kwargs)
