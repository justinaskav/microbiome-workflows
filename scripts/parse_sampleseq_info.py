#!/usr/bin/env python3
"""
Parse Novogene SampleSeq_info.xls file to extract barcode and primer information.
Based on experiment.py implementation for QIIME2 cutadapt workflow.
"""

import argparse

def parse_sampleseq_file(sampleseq_file, output_file=None):
    """
    Parse SampleSeq_info.xls file and extract barcode/primer information.
    
    Args:
        sampleseq_file (str): Path to SampleSeq_info.xls file
        output_file (str): Output TSV file path
    
    Returns:
        str: Path to output TSV file
    """
    if output_file is None:
        output_file = "sample_barcode_primer_info.tsv"
    
    # Read SampleSeq_info.xls and parse it
    with open(sampleseq_file, 'r') as f:
        lines = f.readlines()
    
    # Check if header exists, if not add it
    if not lines[0].startswith('#SampleID'):
        # Add header
        header = '#SampleID\tBarcodeSequence\tLinkerPrimerSequence\n'
        lines.insert(0, header)
    
    # Write parsed information
    with open(output_file, 'w') as out:
        out.write('SampleID\tBarcode_F\tBarcode_R\tPrimer_F\tPrimer_R\n')
        for line in lines[1:]:  # Skip header
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    sample_id = parts[0]
                    barcodes = parts[1].split(',')
                    primers = parts[2].split(',')
                    if len(barcodes) >= 2 and len(primers) >= 2:
                        barcode_f = barcodes[0]
                        barcode_r = barcodes[1]
                        primer_f = primers[0]
                        primer_r = primers[1]
                        # Use additional barcodes if specified
                        out.write(f'{sample_id}\t{barcode_f}\t{barcode_r}\t{primer_f}\t{primer_r}\n')
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Parse Novogene SampleSeq_info.xls file for cutadapt workflow'
    )
    parser.add_argument(
        'sampleseq_file',
        help='Path to SampleSeq_info.xls file'
    )
    parser.add_argument(
        '-o', '--output',
        default='sample_barcode_primer_info.tsv',
        help='Output TSV file path (default: sample_barcode_primer_info.tsv)'
    )
    
    args = parser.parse_args()
    
    # Parse the SampleSeq file
    output_file = parse_sampleseq_file(
        args.sampleseq_file,
        args.output
    )
    
    print(f"Parsed SampleSeq info written to: {output_file}")
    

if __name__ == "__main__":
    main()
