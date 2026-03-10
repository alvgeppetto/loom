---
name: lessons
description: 'Record lessons learned during development. Use when: something unexpected happened, a bug was hard to diagnose, a pattern worked well, a tool behaved unexpectedly, knowledge should be preserved for future reference.'
argument-hint: 'Describe the lesson or insight to record'
---

# Lessons Learned Skill

Capture and retrieve knowledge gained during Project LOOM development.

## When to Use

- After debugging a tricky issue
- When discovering a non-obvious behavior
- When a pattern proves effective (or ineffective)
- When tool/library behavior surprises you
- After any "I wish I knew this earlier" moment

## Recording a Lesson

### 1. Create the file

Location: `docs/lessons/YYYY-MM-DD-<short-topic>.md`

### 2. Use this template

```markdown
# <Topic Title>

Date: YYYY-MM-DD
Tags: [relevant, tags, for, search]

## Context

What were you trying to do? What was the situation?

## What Happened

What unexpected behavior or insight occurred?

## The Lesson

What should be remembered for next time?

## Evidence

<code snippets, error messages, or links that support this lesson>
```

### 3. Cross-reference

If related to an ADR, link to it. If related to a TODO item, note which one.

## Retrieving Lessons

Search existing lessons before starting new work:

```bash
grep -r "<keyword>" docs/lessons/
```

## Example Lessons

- `2026-02-27-rust-fm-index-memory.md` — Memory layout considerations for FM-index in Rust
- `2026-02-27-bwt-sentinel-character.md` — Using `$` vs `\0` as BWT sentinel

## Tags Reference

Common tags for categorization:
- `rust`, `python`, `mlx`
- `bwt`, `fm-index`, `indexing`
- `performance`, `memory`, `debugging`
- `testing`, `cucumber`, `ci`
- `llm`, `ollama`, `inference`
