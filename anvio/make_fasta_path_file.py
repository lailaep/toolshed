import os, re
from datetime import datetime

fasta_dir = "/Users/Snacks/Grad_School/Moran_Lab/Projects/Orbaceae/Anvio/anvio_assemblies/contigs-reformatted"
date_str = datetime.now().strftime("%m%d%y")
output_file = f"fastapaths_clean_{date_str}.txt"

patterns_to_strip = [
    "_annotated_assembly_genomic.fna",
    "_genomic.fna",
    ".fna",
    ".fa"
]

with open(output_file, "w") as out:
    out.write("name\tfasta_path\n")
    for filename in sorted(os.listdir(fasta_dir)):
        if filename.endswith((".fa", ".fna", ".fasta")):
            name = filename
            for pattern in patterns_to_strip:
                name = re.sub(re.escape(pattern) + "$", "", name)
            name = re.sub(r"[^A-Za-z0-9_]", "_", name)
            fasta_path = os.path.join(fasta_dir, filename)
            out.write(f"{name}\t{fasta_path}\n")

print(f"✅ Clean file written: {output_file}")
