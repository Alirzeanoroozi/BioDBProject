import csv

gene_ids = set()
with open("data/part2/biomart_schwannoma_genes.tsv", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        gene_ids.add(row["Gene stable ID"])

with open("data/part4/gene_stable_ids.txt", "w", encoding="utf-8") as out:
    for gid in sorted(gene_ids):
        out.write(gid + "\n")

