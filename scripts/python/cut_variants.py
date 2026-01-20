import pandas as pd

# df = pd.read_csv("data/part2/biomart_variants.tsv", sep="\t")
# # df = df.head(500000)
# df.to_csv("data/part2/biomart_variants_cut.csv", index=False)


# df = pd.read_csv("data/part2/string_interactions.tsv", sep="\t")
# # df = df.head(500000)
# df.to_csv("data/part2/string_interactions.csv", index=False)

# df = pd.read_csv("data/part3/ucsc-table-browser-output.tsv", sep="\t")
# df.to_csv("data/part3/ucsc-table-browser-output_cut.csv", index=False)


# df = pd.read_csv("data/part4/idmapping.tsv", sep="\t")
# df.to_csv("data/part4/idmapping_cut.csv", index=False)


# df = pd.read_csv("data/part6/alphafold_cath.csv")
# df = df.drop(columns=["interactions"])
# df.to_csv("data/part6/alphafold_cath.csv", index=False)

# df = pd.read_csv("data/part2/biomart_schwannoma_genes.tsv", sep="\t")
# df = df.drop(columns=["UniProtKB Gene Name ID"])
# df = df.drop_duplicates()
# df.to_csv("data/part2/biomart_schwannoma_genes_new.csv", index=False)


# df = pd.read_csv("data/part2/biomart_variants.tsv", sep="\t")
# df = df.drop(columns=["Gene Name", "UCSC ID"])
# df = df.drop_duplicates()
# df.to_csv("data/part2/biomart_variants_new.csv", index=False)

df = pd.read_csv("data/part2/biomart_variants_new.csv")
df = df.head(500000)
df.to_csv("data/part2/biomart_variants_new_500000.csv", index=False)