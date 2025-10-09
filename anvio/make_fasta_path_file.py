#!/usr/bin/env python3
import os
import re
from datetime import datetime


# Change this to your directory containing the fasta files
fasta_dir = "/Users/Snacks/Grad_School/Moran_Lab/Projects/Orbaceae/Anvio/anvio_assemblies"

# Get today's date in mmddyy format
date_str = datetime.now().strftime("%m%d%y")
output_file = f"fasta_paths_{date_str}.txt"
# output_file = "fasta_paths.txt" #if you want a simpler version skip the date

# Patterns to remove from filenames
patterns_to_strip = [
    "_annotated_assembly_genomic.fna",
    "_genomic.fna",
    ".fna",
    ".fa"
]

with open(output_file, "w") as out:
    out.write("name\tfasta_path\n")
    for filename in os.listdir(fasta_dir):
        if filename.endswith((".fa", ".fna")):
            # Remove unwanted substrings
            name = filename
            for pattern in patterns_to_strip:
                name = re.sub(re.escape(pattern) + "$", "", name)
            fasta_path = os.path.join(fasta_dir, filename)
            out.write(f"{name}\t{fasta_path}\n")

print(f"✅ File written: {output_file}")
