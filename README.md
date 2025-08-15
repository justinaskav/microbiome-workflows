# Microbiome Workflows

Reusable Nextflow workflows for microbiome analysis projects.

## Overview

This repository contains standardized, modular Nextflow workflows for microbiome research:

- **QIIME2 Amplicon Processing** - DADA2 denoising, taxonomy classification, phylogeny
- **R Microbiome Analysis** - Diversity analysis, differential abundance, visualization  
- **Novogene Data Processing** - Specialized workflows for Novogene amplicon data
- **Quality Control** - Data validation and quality assessment
- **Decontamination** - Remove contaminants using decontam package
- **MaAsLin3 Analysis** - Advanced differential abundance testing
- **Standardized Environments** - Pinned conda environments for reproducibility

## Usage

### In your project's `nextflow.config`:

```groovy
// Include workflows from this repository
includeConfig 'https://github.com/your-username/microbiome-workflows/main/nextflow.config'

// Or for local development:
includeConfig '../microbiome-workflows/nextflow.config'
```

### In your project's main workflow:

```nextflow
include { QIIME2_AMPLICON_PROCESSING } from './microbiome-workflows/workflows/qiime2_processing.nf'
include { NOVOGENE_AMPLICON_WORKFLOW } from './microbiome-workflows/workflows/novogene_amplicon_wf.nf'
include { R_MICROBIOME_ANALYSIS } from './microbiome-workflows/workflows/r_analysis.nf'

workflow {
    // Run QIIME2 processing (all parameters from project config)
    qiime2_results = QIIME2_AMPLICON_PROCESSING(
        params.data_dir,
        params.zip_filename,
        params.metadata,
        params.dada2_trunc_len_f,
        params.dada2_trunc_len_r,
        params.dada2_trim_left_f,
        params.dada2_trim_left_r,
        params.samples_with_additional_barcodes,
        params.use_frequency_decontam,
        params.use_prevalence_decontam,
        params.primer_f,
        params.primer_r,
        params.silva_version,
        params.silva_target,
        params.f_primer,
        params.r_primer
    )
    
    // Run R analysis
    R_MICROBIOME_ANALYSIS(
        qiime2_results.feature_table,
        qiime2_results.taxonomy,
        qiime2_results.phylogeny,
        params.metadata,
        params.formula,
        params.analysis_config
    )
}
```

## Project Structure

```
microbiome-workflows/
├── workflows/              # Main workflow definitions
│   ├── qiime2_processing.nf      # Standard QIIME2 pipeline
│   ├── r_analysis.nf             # R microbiome analysis
│   ├── novogene_amplicon_wf.nf   # Novogene-specific processing
│   ├── graz_amplicon_wf.nf       # Graz lab amplicon workflow
│   └── qc_wf.nf                  # Quality control workflow
├── modules/                # Individual process modules
│   ├── dada2/                    # DADA2 denoising workflows
│   ├── taxonomy/                 # Taxonomy classification
│   ├── phylogeny/                # Phylogenetic analysis
│   ├── diversity/                # Alpha/beta diversity
│   ├── decontam/                 # Contamination removal
│   ├── novogene/                 # Novogene-specific processing
│   ├── silva/                    # SILVA database workflows
│   ├── data_processing/          # Data filtering and processing
│   ├── quality_control/          # Quality assessment
│   ├── r_analysis/               # R statistical analysis
│   │   ├── r_prepare_data.nf     # TSE data preparation
│   │   ├── r_alpha_diversity.nf  # Alpha diversity calculation
│   │   ├── r_beta_diversity.nf   # Beta diversity calculation
│   │   ├── r_ancombc2.nf         # ANCOM-BC2 analysis
│   │   ├── r_maaslin3.nf         # MaAsLin3 analysis
│   │   └── maaslin3_wf.nf        # Legacy MaAsLin3 workflow
│   └── visualization/            # Plotting and reports
├── scripts/                # Supporting R/Python scripts
│   ├── dada2_assign_taxonomy.R
│   ├── format_dada2_taxonomy.R
│   ├── maaslin3.R
│   └── parse_sampleseq_info.py
├── envs/                   # Conda environment definitions
│   ├── qiime2-amplicon.yml
│   └── r-microbiome.yml
├── conf/                   # Configuration profiles
├── tests/                  # Test workflows and data
└── docs/                   # Documentation
```

## Environment Management

This package uses pinned conda environments to ensure reproducibility:

- **qiime2-amplicon.yml** - QIIME2 2025.4 with specific package versions
- **r-microbiome.yml** - R 4.3.x with specific Bioconductor packages

## Contributing

1. Fork this repository
2. Create feature branch: `git checkout -b feature/new-module`
3. Add tests for new functionality
4. Submit pull request

## Version History

- **v1.0.0** - Initial release with QIIME2 and R workflows
- **v1.1.0** - Added differential abundance modules
- **v1.2.0** - Added phylogenetic analysis support
