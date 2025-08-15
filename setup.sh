#!/bin/bash

# Setup script for microbiome-workflows
# This script initializes the environment and validates the setup

set -e

echo "🧬 Setting up microbiome-workflows..."

# Check if conda/mamba is available
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "❌ Error: conda or mamba not found. Please install miniconda/miniforge first."
    exit 1
fi

echo "✅ Using $CONDA_CMD for environment management"

# Create QIIME2 environment
echo "📦 Creating qiime2-amplicon environment..."
$CONDA_CMD env create -f envs/qiime2-amplicon.yml -n qiime2-amplicon || echo "Environment may already exist"

# Create R environment  
echo "📦 Creating r-microbiome environment..."
$CONDA_CMD env create -f envs/r-microbiome.yml -n r-microbiome || echo "Environment may already exist"

# Test QIIME2 environment
echo "🧪 Testing QIIME2 environment..."
$CONDA_CMD run -n qiime2-amplicon qiime --version || echo "QIIME2 test failed - you may need to install QIIME2 manually"

# Test R environment
echo "🧪 Testing R environment..."
$CONDA_CMD run -n r-microbiome R --version || echo "R test failed"

# Check if Nextflow is available
if command -v nextflow &> /dev/null; then
    echo "✅ Nextflow found: $(nextflow -version)"
else
    echo "⚠️  Warning: Nextflow not found. Install with: curl -s https://get.nextflow.io | bash"
fi

echo ""
echo "🎉 Setup complete! You can now:"
echo "   1. Initialize your projects using the microbiome-project-template"
echo "   2. Include workflows from this repository in your projects"
echo "   3. Use standardized environments: qiime2-amplicon and r-microbiome"
echo ""
echo "💡 Next steps:"
echo "   - Copy microbiome-project-template to create new projects"
echo "   - Update project-specific parameters in nextflow.config"
echo "   - Run workflows with: nextflow run workflows/main.nf"
echo ""
echo "📚 Documentation: See README.md for detailed usage instructions"