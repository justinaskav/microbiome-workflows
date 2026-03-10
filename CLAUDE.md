# Microbiome Workflows — Architecture Guide

QIIME2-based amplicon microbiome analysis pipeline in Nextflow DSL2.

## Directory Structure

```
microbiome-workflows/
├── workflows/          # Top-level orchestrators (compose subworkflows)
│   ├── amplicon_wf.nf              # Full per-batch pipeline
│   ├── amplicon_trimming_only_wf.nf # Trim only (for pre-DADA2 merge)
│   └── amplicon_post_trim_wf.nf     # DADA2 + downstream (after merge)
├── subworkflows/       # Composable stages (compose modules, apply config)
│   ├── dada2_wf.nf, filtering_wf.nf, taxonomy_wf.nf, ...
│   └── trimming_wf.nf, diversity_wf.nf, picrust2_wf.nf, ...
├── modules/            # Atomic units (processes + thin workflow wrappers)
│   ├── dada2/, decontam/, diversity/, phylogeny/, picrust2/,
│   │   sequencing/, taxonomy/, quality_control/, ...
│   └── data_processing/, input_processors/, reference_databases/
├── templates/          # Config template for new projects
└── scripts/            # Helper R/Python scripts used by processes
```

**Rule**: Workflows call subworkflows. Subworkflows call modules. Modules contain processes. Don't skip layers (e.g., workflows should not call processes directly).

## Workflow Strategies

Two batch-merging strategies, selected via `params.merge_strategy`:

1. **Post-DADA2 merge** (`amplicon_wf.nf`): Full pipeline per batch, merge feature tables after. Simpler but ASVs may differ between batches.
2. **Pre-DADA2 merge** (`amplicon_trimming_only_wf.nf` + `amplicon_post_trim_wf.nf`): Trim per batch, merge reads, run single DADA2. Better ASV consistency. This is the current active strategy in `main.nf`.

## Parameter Resolution

Parameters follow a 3-tier cascade (highest priority first):

```
batch_option ?: params.section?.key ?: hardcoded_default
```

### How it works:
1. **Config** (`nextflow.config`) defines `params.*` maps: `params.dada2`, `params.taxonomy`, `params.decontam`, `params.phylogeny`, `params.picrust2`, `params.cutadapt`, `params.qc`
2. **Batch config** (`params.batches[]`) can override per-batch: `dada2`, `trimming`, `decontam`, `taxonomy`, `other`
3. **Subworkflows** receive option maps and merge with global params:
   ```groovy
   n_threads: dada2_options.n_threads ?: params.dada2?.n_threads ?: 8
   ```
4. **Processes** get individual values passed from their parent subworkflow/module

### Anti-pattern (DO NOT USE):
```groovy
// BAD: Module-level params declarations create shadow namespaces
params.classifier_confidence = 0.7  // Invisible to config, conflicts with params.taxonomy
```

### Correct pattern for processes accessing global params:
```groovy
// OK: Cross-cutting params accessed directly in processes
cpus params.cutadapt?.cpus ?: 4
conda params.qiime2_conda_env
```

## Batch Configuration

Each batch in `params.batches[]` defines:
- `name`, `metadata_path` — identification
- `input` — data source (data_dir, zip_filename, or accession_ids)
- `primers` — forward/reverse primer sequences
- `trimming` — enabled, use_sampleseq, nextseq_trim
- `dada2` — trunc_len_f/r, trim_left_f/r, and optional overrides
- `decontam` — method override (falls back to `params.decontam`)
- `taxonomy`, `other` — per-batch overrides (falls back to global params)

## Adding New Parameters

1. Add the parameter to `params.*` map in `nextflow.config` with a sensible default
2. Mirror it in `templates/nextflow_template.config` with a comment
3. Wire it through the options map chain at the subworkflow level
4. Use safe navigation (`params.section?.key ?: default`) in fallback expressions
5. Never declare `params.my_param = value` at module level

## Channel Conventions

- Paired reads: `tuple val(sample_id), path([r1, r2])`
- Skipped steps emit `Channel.empty()` — downstream handles gracefully
- Conditional branching happens at subworkflow/workflow level, not in processes
- Use `.mix()` to combine parallel branches, `.collect()` for aggregation

## Known Limitations

- **Primers in non-sampleseq mode**: `trimming_wf.nf` extracts primers from `trimming_options`, but batch config stores primers at the top level (`params.batches[].primers`), not inside `trimming`. In `use_sampleseq` mode this is fine (primers come from barcode_primer_info_ch per sample). For non-sampleseq mode, primers would need to be passed explicitly.
- **Per-batch overrides in pre-DADA2 merge**: `main.nf` uses `params.batches[0].dada2` for the merged DADA2 run. Per-batch `taxonomy` and `other` overrides aren't fully wired — global params are used instead.
- **Two confidence scales**: `params.taxonomy.classification_confidence` (0-1, sklearn) vs DADA2 bootstrap threshold (0-100, hardcoded to 80 in the module). These are intentionally different.
