# Synthon Enumeration Pipeline
A high-throughput, RDKit-based enumeration pipeline currently for generating single-step analogues by attaching one synthon to a seed at a selected reactive site (“mark” location), using a curated set of reaction templates.

This repository is designed to scale to many seeds and millions of synthons by:
- indexing synthons by mark type for fast retrieval
- pre-indexing reaction templates by compatible mark pairs
- batching outputs for downstream processing (e.g., docking, filtering, scoring)
- supporting both streaming and Parquet materialization

> Note: This repo includes a mention of SyntOn. Detailed SyntOn usage and background are documented in a separate [Readme](SyntOn/README.md) located under the SyntOn directory. It is recommended to read through that first, to get more context on how molecules are processed to annotate chemically significant reactive sites/marks and how they can be reconstructed all using specific reaction templates.

---

## Key Concepts

### Marks / Reactive Sites
Seeds and synthons are annotated with reactive marks using atom-map tags in bracket atoms, e.g.:

- `[C:10]`, `[N:20]`, `[n:20]`, `[S:10]`

A “site” is a specific occurrence of a mark in a seed (not just the mark type). Site identity is deterministic and reproducible. These are chemically determined and can be made using [`mainsynthonsgenerator`](SyntOn/src/SyntOn_BBs.py) from the SyntOn codebase.

### Single-Step Enumeration (Current Scope)
For each seed:
1. Canonicalize seed SMILES (stereochemistry preserved).
2. Identify all reactive sites and assign deterministic a `site_id`.
3. Choose one site to react (random by default, or user-specified).
4. Retrieve all compatible synthons (by mark).
5. Apply only compatible reaction templates (by mark pair) to produce products.

Each product corresponds to a route-level record (even if multiple routes lead to the same final SMILES).

---

## Features

- ✅ Deterministic reactive site IDs on canonicalized seeds
- ✅ Site selection by `site_id` or `MARK@ORDINAL` (e.g., `C:10@0`)
- ✅ Reaction template filtering via a precomputed mark-pair → reaction index
- ✅ Synthon indexing by mark type for fast candidate lookup
- ✅ Output modes:
  - `stream` (generator of batches)
  - `parquet` (sharded Parquet + manifest)
  - `stream+parquet` (yield batches while writing Parquet)
- ✅ Deterministic route_id for stable provenance keys
- ✅ Batched processing for scalability and downstream throughput

---

## Repository Structure

### `standardization.py`
Seed canonicalization / preprocessing

- Parses the input seed SMILES into an RDKit Mol.
- Generates canonical SMILES while preserving stereochemistry.
- Re-parses canonical SMILES to ensure a stable canonical atom order.
- Produces a lightweight `StandardizedSeed` record.

Why it matters: reactive site IDs are defined only on canonicalized seeds, ensuring reproducibility.

---

### `sites.py`
Reactive site inspection + user selection resolution (separate from enumeration)

- Finds atoms with atom-map numbers (`GetAtomMapNum()`).
- Builds normalized `mark_type` strings (`C:10`, `n:20`, etc.).
- Assigns deterministic `site_id` values using canonical atom ranking.
- Computes `mark_occurrence_index` and `site_label` (e.g., `C:10@0`).
- Resolves user constraints:
  - `allowed_sites=[0, 2]`
  - `allowed_site_specs=["C:10@1"]`
- Supports error/warn policies for invalid site requests.

Why it matters: users can target specific reactive sites unambiguously, even if mark types repeat.

---

### `reactions.py`
Reaction template parsing + `ReactionIndex`

- Parses `Setup.xml` reaction definitions:
  - `Labels` → which mark types belong on each reactant side
  - `ReconstructionReaction` → SMARTS used for enumeration
- Pre-indexes templates by mark pair:
  - `(seed_mark, synthon_mark) -> [(reaction_id, order), ...]`
- Caches compiled RDKit reactions for reuse.
- Tracks reactant order:
  - `order=0`: (seed, synthon)
  - `order=1`: (synthon, seed)

Why it matters: avoids trying every reaction template for every candidate.

---

### `synthons.py`
Synthon records + `SynthonIndex`

- Stores minimal synthon metadata:
  - `synthon_id`, `synthon_smiles`, extracted `marks`
- Builds an inverted index:
  - `mark_type -> [synthon_id, ...]`
- Designed for large libraries (millions of synthons).
- RDKit Mol creation is deferred until a synthon is actually tested in a reaction.

Why it matters: fast retrieval of candidates compatible with a seed site.

---

### `output_sinks.py`
Parquet output writer + manifest

- `ParquetSink` writes one Parquet shard per batch.
- Includes `manifest.json` with shard list, counts, and run metadata.
- Uses PyArrow for efficient columnar output and compression.

Why it matters: enables checkpointing, replay, and decoupled downstream pipelines.

---

### `enumeration_single_step.py`
Single-step enumeration engine

- Ties everything together:
  - canonicalize seed
  - inspect sites & resolve constraints
  - choose one site (deterministic RNG)
  - filter compatible marks
  - fetch compatible synthons
  - fetch applicable reaction templates
  - run RDKit reactions
  - emit route-level records in batches
- Supports:
  - `output_mode="stream" | "parquet" | "stream+parquet"`

---

## Installation

### Dependencies
- Python 3.9+
- RDKit
- PyArrow (required for Parquet output modes)

Install dependencies using your preferred environment manager (conda is common for RDKit).

---

## Configuration Inputs

### 1) `Setup.xml` (reaction definitions)
This file defines available reactions and reconstruction SMARTS. The pipeline uses:
- `Labels` to map reaction “handles” to mark types (e.g., `C:10`, `N:20`)
- `ReconstructionReaction` for the actual SMARTS applied during enumeration

### 2) `marks_compatibility` (mark compatibility map)
A dictionary mapping seed mark type → allowed synthon partner mark types, e.g.:

```python
marks_compatibility = {
  "C:10": ["N:20", "O:20"],
  "n:20": ["C:10"],
  # ...
}
```

This is used as a fast first-stage filter before querying reactions or synthons.

---

## Output Records (What You Get)

Each enumerated result is a route-level record (dict-like) containing:

- `route_id` (deterministic hash)
- `seed_id`, `seed_canonical_smiles`
- `seed_site_id`, `seed_site_mark_type`, `seed_site_atom_idx`
- `synthon_id`, `synthon_smiles`, `synthon_mark_type`
- `reaction_id`, `reaction_name`, `reactant_order`
- `product_smiles`, `product_valid`, `failure_reason`
- run metadata: `run_id`, `rng_seed`

> Even if two routes lead to the same final product SMILES, both records are retained (route provenance is preserved).

---

## Quick Start

### Build indices and enumerate
```python
from reactions import ReactionIndex
from synthons import SynthonIndex
from enumeration_single_step import SingleStepEnumerator, SeedSpec

rxn_index = ReactionIndex.from_setup_xml("Setup.xml")
syn_index = SynthonIndex.from_smi_file("synthons.smi")

marks_compatibility = {
    "C:10": ["N:20", "O:20"],
    "N:20": ["C:10"],
    # ...
}

enum = SingleStepEnumerator(
    synthon_index=syn_index,
    reaction_index=rxn_index,
    marks_compatibility=marks_compatibility,
    rng_seed=123,
    invalid_site_policy="error",
)

seeds = [
    SeedSpec(seed_smiles="CC([C:10])N", seed_id="seed1"),
    SeedSpec(seed_smiles="O=C([C:10])c1ccccc1", seed_id="seed2", allowed_site_specs=["C:10@0"]),
]

# Stream
gen, _ = enum.enumerate(seeds, batch_size=1000, output_mode="stream", run_id="run1")
for batch in gen:
    # consume batch (list of records)
    pass
```

### Write Parquet shards
```python
_, summary = enum.enumerate(
    seeds,
    batch_size=500,
    output_mode="parquet",
    out_dir="out_parquet",
    run_id="run2",
)
print(summary)
```

### Stream + Parquet (side-effect writing)
```python
gen, _ = enum.enumerate(
    seeds,
    batch_size=2000,
    output_mode="stream+parquet",
    out_dir="out_parquet2",
    run_id="run3",
)
for batch in gen:
    pass
# manifest.json is written to out_parquet2/ in this example
```

---

## Site Identification & User Selection

To inspect sites (and their IDs) for a seed in case you want to specify during enumeration, you can run the following:

```python
from standardization import canonicalize_seed_smiles
from sites import list_sites_pretty

std, mol = canonicalize_seed_smiles("CC([C:10])N")
print(list_sites_pretty(mol))
```

If a seed has multiple reaction sites, you can restrict enumeration to specific sites via either:
- `allowed_sites=[0]` (site_id)
- `allowed_site_specs=["C:10@0"]` (mark + ordinal occurrence)

If a user specifies an invalid site for a seed, policy determines how to handle it:
- policy `error` raises with a listing of valid sites
- policy `warn_skip_seed` logs and skips that seed
- policy `warn_fallback_all` logs and falls back to all sites and random selection

---

## Reproducibility

Reproducibility and tracability is ensured via:
- canonical seed representation (`seed_canonical_smiles`)
- deterministic site IDs on the canonical seed
- deterministic `route_id` derived from route components
- user-configurable RNG seed for random site selection

The results of every run contains all of the information to understand and trace how it was generated.

---

## Performance Notes and Recommendations (Scaling Guidance)

- Index once, enumerate many: build `SynthonIndex` and `ReactionIndex` once and reuse. Its quick to retrieve by mark type but slower to build and sort many times
- Batch size matters: larger batches (to an extent) reduce overhead (especially for Parquet writing).
- Sharded outputs: for HPC/multi-process workflows, write per-worker output shards to avoid file contention.
- Avoid materializing huge lists: prefer `stream` or `stream+parquet` for large jobs.

---

## SyntOn (Mention Only)

This project interoperates with Synt-On (synthon generation/handling) but does not document it here. This is a proejct by the Cheminfortaics lab in the Univeristy of montreal. Refer to the dedicated SyntOn README for details on synthon generation, preprocessing, and library construction.

---

## Roadmap

Planned extensions (not necessarily implemented yet):

- Multi-step enumeration (iterative growth up to `max_reaction_steps`)
- Constraints on number of reacted sites per seed (`max_sites_reacted`)
- Advanced site targeting (SMARTS-context selection)
- Optional product deduplication strategies (route-preserving vs product-unique)
- Parallel/distributed enumeration strategies (per-seed sharding, per-worker Parquet shards)


---


## Contributing

Contributions are welcome, please fork if you are interested. Suggested areas:
- additional reaction schemas (for synthonization or reconstruction)
- performance profiling and optimization
- improved reaction/template parsing and validation
- support for additional input formats and HPC job orchestration patterns
- expanded test coverage (especially for tricky chemistry edge cases)

---

## License

Add your license information here.
