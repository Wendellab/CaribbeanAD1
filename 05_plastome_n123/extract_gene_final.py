import sys
from Bio import SeqIO
import os

# Input: list of GenBank files (one per line in a text file)
input_list_file = sys.argv[1]  # File containing paths to GenBank files
output_dir = sys.argv[2]       # Directory to store per-gene FASTA files

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Dictionary to store sequences for each gene, ensuring uniqueness per sample
gene_sequences = {}

# Read the list of GenBank files
with open(input_list_file, "r") as file_list:
    for line in file_list:
        genbank_file = line.strip()
        sample_name = os.path.splitext(os.path.basename(genbank_file))[0]

        # Parse each GenBank file
        for record in SeqIO.parse(genbank_file, "genbank"):
            for feature in record.features:
                # Include gene, CDS, rRNA, and tRNA features
                if feature.type in ["gene", "CDS", "rRNA", "tRNA"]:
                    # Extract gene name, product, or locus_tag as a fallback
                    gene_name = feature.qualifiers.get("gene", [
                        feature.qualifiers.get("product", [
                            feature.qualifiers.get("locus_tag", ["unknown"])[0]
                        ])[0]
                    ])[0]

                    # Skip pseudogenes or partial genes
                    if "pseudo" in feature.qualifiers or "partial" in feature.qualifiers:
                        continue

                    # Extract the gene sequence
                    sequence = str(feature.extract(record.seq))

                    # Create a unique identifier for this gene in this sample
                    unique_id = f"{sample_name}|{gene_name}"

                    # Initialize gene entry if it doesn't exist
                    if gene_name not in gene_sequences:
                        gene_sequences[gene_name] = {}

                    # Store the sequence only if it is unique for the sample
                    if unique_id not in gene_sequences[gene_name]:
                        gene_sequences[gene_name][unique_id] = sequence

# Write each gene to its own FASTA file
for gene_name, sample_sequences in gene_sequences.items():
    output_file = os.path.join(output_dir, f"{gene_name}.fasta")
    with open(output_file, "w") as out_fasta:
        for unique_id, sequence in sample_sequences.items():
            out_fasta.write(f">{unique_id}\n{sequence}\n")
