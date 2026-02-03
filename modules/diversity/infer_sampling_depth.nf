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

counts = []
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
                    counts.append(float(parts[1].strip().replace(",", "")))
                except ValueError:
                    pass

if not counts:
    raise RuntimeError("No per-sample counts found in biom summary")

counts.sort()

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

    main:
    result = INFER_SAMPLING_DEPTH(feature_table_ch, percentile)

    emit:
    sampling_depth_file = result.sampling_depth_file
}
