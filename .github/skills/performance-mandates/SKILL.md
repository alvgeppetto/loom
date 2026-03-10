---
name: performance-mandates
description: 'Apply LOOM performance doctrine. Use when implementing, reviewing, benchmarking, or proposing architecture/model/runtime changes. Enforces green compute + commodity hardware + math-first optimization + facade philosophy.'
argument-hint: 'Describe the feature/change and target runtime/hardware'
---

# Performance Mandates Skill (LOOM)

This skill codifies non-negotiable performance policy for Project LOOM.

## Target Hardware (MANDATORY CONTEXT)

- **Build machine:** Mac Studio M3 Ultra — 28 cores (20P + 8E), 96 GB unified RAM
- **ISA:** ARM64 (AArch64) with NEON SIMD. **Never assume x86/AVX/CUDA.**
- **Memory bandwidth:** ~800 GB/s (GPU), ~400 GB/s (CPU) — unified, no PCIe bottleneck
- **BLAS:** Apple Accelerate framework (native). No external BLAS deps.
- **Edge targets:** Raspberry Pi 5 (8 GB), iPads/phones, browsers (WASM)
- **f32::mul_add() is HARMFUL on M3 Ultra** — causes ~6% regression in hot loops
- **Rayon dispatch overhead (~50-100µs)** exceeds benefits for small tensors (T≤1024, D≤128)

See also: `/memories/repo/hardware-profile.md` and `/memories/repo/rust-perf-patterns.md`

## Core Mandates

1. Green compute first
- Minimize operations, memory traffic, and wall-time before scaling model size or hardware.
- Treat energy efficiency as a first-class requirement.

2. Commodity hardware first
- Optimize and benchmark on laptop/CPU/mobile-class hardware first.
- High-end accelerators are optional, not baseline.

3. Reuse efficient architecture first
- Prefer LOOM production paths: FM-index retrieval, Loom-1, WASM sharding.
- Do not introduce heavier alternatives without measured net gain.

4. Math-first optimization doctrine
- Use Hacker's Delight patterns for bit-level/branchless optimization where applicable.
- Use Sakarovitch-style automata/state minimization and transition compression.
- Search for algebraic properties: monoids, groups, symmetry, semiring transforms.
- Use matrix/power methods (binary exponentiation, recurrence transforms) when structure permits.

5. Facade philosophy (Book Ch. 20)
- Replace expensive components that are weakly observable by loss with structurally cheaper surrogates.
- Maximize structural leverage (P/Q) under bounded approximation gap.
- Prefer lossless facades first, then near-lossless facades with ablation evidence.

## Mandatory Evaluation Gates (Before Claiming Improvement)

For any optimization/model/runtime change, report:
- Latency (before/after)
- Throughput (before/after)
- Memory/state footprint (before/after)
- Energy proxy (ops/bytes or equivalent)
- Hardware target used for measurement

A change is not accepted if it improves one metric but violates project priorities:
- worse commodity-hardware latency,
- worse energy/compute efficiency,
- or unacceptable regression against Loom-1 baseline.

## Implementation Checklist

- Identify hottest loop/path first (measure, do not guess).
- Remove branches/allocations from hot loops.
- Extract kernels so compiler can vectorize (noalias-friendly slice signatures).
- Fuse composable element-wise operators (monoid composition).
- Hoist loop-invariant work out of inner loops.
- Validate correctness + benchmark deltas after each optimization step.

## BitLoom Policy

- BitLoom is archived under `research/bitloom/**`.
- Do not treat BitLoom as current production architecture.
- Re-entry requires passing strict latency/efficiency gates on commodity hardware.

## Where These Mandates Came From

- User directives (green compute, commodity hardware, optimize optimize optimize).
- `docs/book/chapters/19-commodity-training.tex`
- `docs/book/chapters/20-facade-architecture.tex`
- `docs/lessons/2026-03-03-algorithmic-training-opts.md`
