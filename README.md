# barcode_overlapper

 R script for analyzing and comparing barcodes from GEX (Gene Expression) and ADT (Antibody-Derived Tag) samples in single-cell sequencing data. Here's a breakdown of its main components and functionality:

1. Libraries: The script uses ShortRead for handling FASTQ files, and ggplot2 and reshape2 for data visualization.

2. Functions:
   - `extract_barcodes()`: Extracts barcodes from a FASTQ file using streaming, which is memory-efficient for large files.
   - `calculate_overlap()`: Calculates the overlap between two sets of barcodes.
   - `find_fastq_files()`: Locates R1 FASTQ files in specified directories.

3. File Identification:
   - The script finds GEX and ADT FASTQ files in the specified base directory.

4. Barcode Extraction:
   - Barcodes are extracted from all identified GEX and ADT files.

5. Overlap Calculation:
   - The script calculates the overlap between each pair of GEX and ADT samples, storing the results in a matrix.

6. Visualization:
   - The overlap data is transformed into a long format and visualized as a heatmap using ggplot2.
   - The heatmap shows the barcode overlap between GEX and ADT samples, with color intensity representing the degree of overlap.

7. Output:
   - The heatmap is saved as a PNG file named "barcode_overlap_heatmap.png".

This script is useful for quality control and sample identification in single-cell sequencing experiments, particularly when working with multi-omic data (e.g., combined gene expression and protein abundance measurements). It helps identify potential sample mix-ups or contamination by comparing barcode overlap between different sample types.

The code appears well-structured and uses efficient methods like streaming for handling large files. However, it's worth noting that the subsample size (1e6) and barcode length (16) are hardcoded, which might need adjustment depending on the specific experimental setup.


Expect you folder structure : 
GEX1, GEX2, GEX3, GEX4 ....
ADT1, ADT2, ADT3, ADT4 ....
