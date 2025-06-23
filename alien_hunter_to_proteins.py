# This script utilized ChatGPT for troubleshooting
# Laila E. Phillips 6/23/25
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
import re

def get_input(prompt, default=None):
    user_input = input(f"{prompt} [{default}]: ").strip()
    return user_input if user_input else default

def get_genome_size_from_gbff(gbff_file):
    record = next(SeqIO.parse(gbff_file, "genbank"))
    return len(record.seq)

def parse_alien_embl(alien_file):
    regions = []
    with open(alien_file) as f:
        for line in f:
            line = line.strip()
            match = re.search(r'FT\s+misc_feature\s+(\d+)\.\.(\d+)', line)
            if match:
                start, end = int(match.group(1)), int(match.group(2))
                regions.append((start, end))
    return regions

def get_gene_intervals(feature, genome_length):
    intervals = []
    loc = feature.location
    if isinstance(loc, CompoundLocation):
        for part in loc.parts:
            start = int(part.start) + 1
            end = int(part.end)
            intervals.append((start, end))
    else:
        start = int(loc.start) + 1
        end = int(loc.end)
        intervals.append((start, end))
    adjusted = []
    for s, e in intervals:
        if e < s:
            adjusted.append((s, genome_length))
            adjusted.append((1, e))
        else:
            adjusted.append((s, e))
    return adjusted

def parse_genbank_genes(gbff_file, genome_size):
    genes = []
    for record in SeqIO.parse(gbff_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    gene_id = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                    product = feature.qualifiers.get("product", ["unknown"])[0]
                    protein_id = feature.qualifiers.get("protein_id", ["none"])[0]
                    loc = feature.location
                    g_start = int(loc.start) + 1
                    g_end = int(loc.end)
                    intervals = get_gene_intervals(feature, genome_size)
                    genes.append((gene_id, g_start, g_end, protein_id, product, intervals))
                except Exception as e:
                    print(f"Skipping a CDS due to error: {e}")
    return genes

def check_overlap(start1, end1, start2, end2, genome_size=None):
    if genome_size:
        def interval_set(s, e):
            if e < s:
                return set(range(s, genome_size+1)) | set(range(1, e+1))
            else:
                return set(range(s, e+1))
        set1 = interval_set(start1, end1)
        set2 = interval_set(start2, end2)
        return not set1.isdisjoint(set2)
    else:
        return not (end1 < start2 or end2 < start1)

def match_genes_to_hgt(hgt_regions, genes, genome_size=None):
    matches = []
    for gene in genes:
        gene_id, g_start, g_end, prot_id, product, intervals = gene
        for h_start, h_end in hgt_regions:
            for s, e in intervals:
                if check_overlap(s, e, h_start, h_end, genome_size):
                    matches.append((gene_id, g_start, g_end, h_start, h_end, prot_id, product))
                    break
            else:
                continue
            break
    return matches

def extract_protein_sequences(protein_ids, faa_file, output_fasta):
    records = SeqIO.parse(faa_file, "fasta")
    protein_dict = {}
    for rec in records:
        prot_id = rec.id.split()[0]
        protein_dict[prot_id] = rec

    found = 0
    with open(output_fasta, "w") as out:
        for pid in protein_ids:
            if pid in protein_dict:
                SeqIO.write(protein_dict[pid], out, "fasta")
                found += 1
            else:
                print(f"Protein ID not found in .faa: {pid}")

    print(f"✅ Extracted {found}/{len(protein_ids)} protein sequences to {output_fasta}")

if __name__ == "__main__":
    prefix = get_input("Enter file prefix (e.g. lpD04)", "lpD04")

    alien_file = f"{prefix}.alien"
    gbff_file = f"{prefix}_genomic.gbff"
    faa_file = f"{prefix}_proteins.faa"

    print(f"Using files:\n  Alien_Hunter: {alien_file}\n  GenBank: {gbff_file}\n  Proteins: {faa_file}")

    genome_size = get_genome_size_from_gbff(gbff_file)
    print(f"Detected genome size: {genome_size} bp")

    tsv_output = f"{prefix}_HGT_gene_hits.tsv"
    fasta_output = f"{prefix}_HGT_proteins.faa"

    print("Parsing Alien_Hunter HGT regions...")
    hgt_regions = parse_alien_embl(alien_file)
    print(f"Found {len(hgt_regions)} HGT regions.")

    print("Parsing GenBank genes...")
    genes = parse_genbank_genes(gbff_file, genome_size)
    print(f"Parsed {len(genes)} CDS genes.")

    print("Matching genes to HGT regions...")
    matches = match_genes_to_hgt(hgt_regions, genes, genome_size)
    print(f"Total genes overlapping HGT: {len(matches)}")
    for m in matches[:10]:
        print(m)

    print(f"Writing output table to {tsv_output} ...")
    with open(tsv_output, "w") as out:
        out.write("gene_id\tgene_start\tgene_end\thgt_start\thgt_end\tprotein_id\tproduct\n")
        for row in matches:
            out.write("\t".join(map(str, row)) + "\n")

    protein_ids = [m[5] for m in matches if m[5] != "none"]
    print(f"Extracting protein sequences for {len(protein_ids)} proteins...")
    extract_protein_sequences(protein_ids, faa_file, fasta_output)
