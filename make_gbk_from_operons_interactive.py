# This script takes operon-mapper output and a nucleotide fasta file, and using user-defined terms, pulls separate operons into .gbk format for usage with programs like Clinker
# This version takes user-inputted filepaths and user-specified genome (or sample or whatever) name.

# Laila Phillips 7/7/25
# Assistance from ChatGPT

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

# ==== USER INPUT PROMPTS ====
genome_name = input("Enter genome name (e.g., lpD04): ").strip()
operon_file = input("Enter path to operon file (list_of_operons_*.txt): ").strip()
protein_file = input("Enter path to protein FASTA file (predicted_protein_sequences_*.txt): ").strip()
genome_fasta = input("Enter path to genome contig FASTA file: ").strip()
keyword = input("Enter keyword to search for (e.g., fimA): ").strip().lower()
output_dir = f"{genome_name}_matched_operons_gbk"
os.makedirs(output_dir, exist_ok=True)

# ==== LOAD INPUT FILES ====

# Load proteins (strip trailing *)
protein_seqs = {
    record.id: str(record.seq).rstrip("*")
    for record in SeqIO.parse(protein_file, "fasta")
}

# Load genome
genome = SeqIO.read(genome_fasta, "fasta")
genome_seq = genome.seq

# Load operon list
df = pd.read_csv(operon_file, sep="\t", comment="#")
df['Operon'] = df['Operon'].ffill()

# Filter operons that contain the keyword
matched_operon_ids = (
    df[df['Function'].fillna("").str.contains(keyword, case=False)]
    ['Operon']
    .unique()
)

# Process each matched operon
for operon_id in matched_operon_ids:
    sub_df = df[df['Operon'] == operon_id]

    # Find region start/end
    region_start = int(sub_df['PosLeft'].min()) - 1
    region_end = int(sub_df['postRight'].max())

    # Create a new SeqRecord
    subseq = genome_seq[region_start:region_end]
    record = SeqRecord(
        subseq,
        id=f"{genome_name}_operon_{operon_id}",
        name=f"{genome_name}_operon_{operon_id}",
        description=f"{genome_name} operon {operon_id} containing '{keyword}'"
    )
    record.annotations["molecule_type"] = "DNA"

    features = []
    for _, row in sub_df.iterrows():
        gene_id = row['IdGene']
        if pd.isna(gene_id): continue
        start = int(row['PosLeft']) - region_start - 1  # adjust to local coords
        end = int(row['postRight']) - region_start
        strand = 1 if row['Strand'] == '+' else -1
        product = str(row['Function']) if pd.notnull(row['Function']) else "hypothetical protein"

        qualifiers = {
            "locus_tag": gene_id,
            "product": product,
        }
        if gene_id in protein_seqs:
            qualifiers["translation"] = protein_seqs[gene_id]

        feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type="CDS",
            qualifiers=qualifiers
        )
        features.append(feature)

    record.features = features

    # Save to .gbk
    out_path = os.path.join(output_dir, f"{genome_name}_operon_{operon_id}_{keyword}.gbk")
    with open(out_path, "w") as f:
        SeqIO.write(record, f, "genbank")

    print(f"Saved: {out_path}")
