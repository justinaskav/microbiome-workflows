// Infer Sampling Depth from Feature Table
// Exports a QIIME2 feature table to BIOM, summarizes per-sample counts,
// and computes the sampling depth as a minimum or configurable percentile.

nextflow.enable.dsl = 2

process INFER_SAMPLING_DEPTH {
    tag "infer_sampling_depth"
    conda "${params.qiime2_conda_env}"

    input:
    path feature_table_qza
    val percentile  // 0.0 = minimum, 0.1 = 10th percentile, etc.
    path metadata_file

    output:
    path "sampling_depth.txt", emit: sampling_depth_file

    script:
    """
    # Export QZA to BIOM
    qiime tools export \
        --input-path ${feature_table_qza} \
        --output-path exported

    # Summarize per-sample counts
    biom summarize-table \
        --input-fp exported/feature-table.biom \
        --output-fp biom_summary.txt

    # Parse summary and compute depth at given percentile
    python3 <<'PYEOF'
import re

percentile = ${percentile}
metadata_path = "${metadata_file}"

# Parse BIOM summary for per-sample counts
biom_samples = {}
in_counts_section = False
with open("biom_summary.txt") as f:
    for line in f:
        line = line.strip()
        if line.startswith("Counts/sample detail:"):
            in_counts_section = True
            continue
        if in_counts_section and ": " in line:
            # Split on last ": " to handle sample names that may contain colons
            parts = line.rsplit(": ", 1)
            if len(parts) == 2:
                try:
                    # Remove comma thousand separators before parsing
                    count = float(parts[1].strip().replace(",", ""))
                    biom_samples[parts[0].strip()] = count
                except ValueError:
                    pass

if not biom_samples:
    raise RuntimeError("No per-sample counts found in biom summary")

counts = sorted(biom_samples.values())

# Validate against metadata if provided
if metadata_path:
    metadata_sample_ids = set()
    with open(metadata_path) as f:
        header_skipped = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if not header_skipped:
                header_skipped = True
                continue
            if line.startswith("#"):
                continue
            sample_id = line.split("\\t")[0].strip()
            if sample_id:
                metadata_sample_ids.add(sample_id)

    biom_sample_ids = set(biom_samples.keys())
    missing_from_table = metadata_sample_ids - biom_sample_ids
    if missing_from_table:
        print(f"WARNING: {len(missing_from_table)} sample(s) in metadata but not in feature table: {sorted(missing_from_table)}")
        print(f"  Feature table has {len(biom_sample_ids)} samples, metadata has {len(metadata_sample_ids)} samples")
    extra_in_table = biom_sample_ids - metadata_sample_ids
    if extra_in_table:
        print(f"WARNING: {len(extra_in_table)} sample(s) in feature table but not in metadata: {sorted(extra_in_table)}")

if percentile <= 0.0:
    depth = int(counts[0])
else:
    idx = int(percentile * (len(counts) - 1))
    depth = int(counts[idx])

# Ensure depth is at least 1
depth = max(depth, 1)

with open("sampling_depth.txt", "w") as out:
    out.write(str(depth))

print(f"Inferred sampling depth: {depth} (percentile={percentile}, n_samples={len(counts)}, min={int(counts[0])}, max={int(counts[-1])})")
PYEOF
    """
}

workflow WF_INFER_SAMPLING_DEPTH {
    take:
    feature_table_ch   // path to feature table QZA
    percentile          // val: 0.0 for min, 0.1 for 10th percentile, etc.
    metadata_file       // path to metadata TSV for sample count validation

    main:
    result = INFER_SAMPLING_DEPTH(feature_table_ch, percentile, metadata_file)

    emit:
    sampling_depth_file = result.sampling_depth_file
}
