# Polymer Builder GitHub Preparation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prepare Polymer Builder as a clean, independently publishable private GitHub repository.

**Architecture:** Add standard Python/Git exclusions, preserve the application, chemistry code, UI assets, tests, documentation, and locked dependencies, then remove only reproducible local artifacts.

**Tech Stack:** Python 3.11, uv, Dash, RDKit, OpenBabel, Git

**Spec:** User-approved three-private-repository cleanup arrangement from 2026-10-27.

## Global Constraints

- Preserve source, tests, documentation, assets, `requirements.txt`, and `uv.lock`.
- Do not change application behavior.
- Do not stage virtual environments, caches, secrets, editor settings, or generated outputs.

## Review Focus

- `.gitignore` must cover Python, uv, test, editor, OS, environment, and build artifacts.
- The repository must contain all files needed for `uv sync --locked`.
- Source assets used by the Dash UI must remain tracked.
- The existing test suite must still run after recreating dependencies.
- No credential-shaped file may be staged.

---

### Task 1: Add exclusions and clean the project

**Files:**
- Create: `.gitignore`
- Delete: `.venv/` and any Python/test/tool cache directories

**Interfaces:**
- Consumes: the current source tree
- Produces: a source-only project tree protected from common local artifacts

- [ ] Add the project `.gitignore`.
- [ ] Delete reproducible/generated directories.
- [ ] Verify no required application asset is excluded.

### Task 2: Initialize and verify Git metadata

**Files:**
- Create: `.git/`

**Interfaces:**
- Consumes: the cleaned project tree
- Produces: an initial local Git repository ready for a private GitHub remote

- [ ] Initialize Git with branch `main`.
- [ ] Stage and inspect the complete file list.
- [ ] Run available tests/lint checks after recreating the locked environment if necessary.
- [ ] Create the initial commit.
