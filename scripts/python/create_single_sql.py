import csv
import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

# Order matters because of foreign key dependencies
TABLE_FILES = [
    "genes.csv",
    "proteins.csv",
    "variants.csv",
    "geo2r.csv",
    "Transcript.csv",
    "alphafoldCath.csv",
    "stringInteractions.csv",
]

# Additional foreign key constraints to strengthen relationships
FK_ADDITIONS = {
    # variants.Gene_stable_ID -> genes.Gene_stable_ID
    "variants": [
        "CONSTRAINT `fk_variants_gene` FOREIGN KEY (`Gene_stable_ID`) "
        "REFERENCES `genes` (`Gene_stable_ID`) ON DELETE NO ACTION ON UPDATE NO ACTION",
    ],
    # geo2r.gene_symbol -> genes.HGNC_symbol
    "geo2r": [
        "CONSTRAINT `fk_geo2r_gene_symbol` FOREIGN KEY (`gene_symbol`) "
        "REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION",
    ],
    # Transcript.name2 (gene symbol) -> genes.HGNC_symbol
    "Transcript": [
        "CONSTRAINT `fk_transcript_gene_symbol` FOREIGN KEY (`name2`) "
        "REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION",
    ],
    # alphafoldcath.uniprot_acc -> proteins.Entry (UniProt accession)
    "alphafoldcath": [
        "CONSTRAINT `fk_alphafoldcath_protein` FOREIGN KEY (`uniprot_acc`) "
        "REFERENCES `proteins` (`Entry`) ON DELETE NO ACTION ON UPDATE NO ACTION",
    ],
    # StringInteractions nodes -> genes.HGNC_symbol
    "StringInteractions": [
        "CONSTRAINT `fk_stringinteractions_node1` FOREIGN KEY (`#node1`) "
        "REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION",
        "CONSTRAINT `fk_stringinteractions_node2` FOREIGN KEY (`node2`) "
        "REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION",
    ],
}


def load_create_statement(csv_path: str) -> tuple[str, str]:
    """
    Read a 2-column CSV (Table, Create Table) and return (table_name, create_sql).
    The CREATE TABLE statement may span multiple lines inside the second column.
    """
    with open(csv_path, newline="") as f:
        reader = csv.reader(f)
        # Skip header
        _ = next(reader, None)
        row = next(reader, None)
        if not row or len(row) < 2:
            raise ValueError(f"Unexpected format in {csv_path}")
        table_name = row[0].strip()
        create_sql = row[1]
        return table_name, create_sql


def add_foreign_keys(table_name: str, create_sql: str) -> str:
    """
    Inject additional foreign key constraints into the CREATE TABLE statement
    just before the final ') ENGINE=InnoDB ...' block.
    """
    additions = FK_ADDITIONS.get(table_name)
    if not additions:
        return create_sql

    marker = ") ENGINE=InnoDB"
    idx = create_sql.rfind(marker)
    if idx == -1:
        # Unexpected format; don't modify
        return create_sql

    before = create_sql[:idx].rstrip()
    after = create_sql[idx:]

    if not before.endswith(")"):
        # Safety check – don't corrupt the statement
        return create_sql

    # Remove the final ')' so we can append more column/constraint lines
    before_without_closing = before[:-1]

    fk_block = ",\n  " + ",\n  ".join(additions) + "\n"

    # Rebuild with the fk_block and put back the closing parenthesis
    rebuilt = before_without_closing + fk_block + ")\n" + after[len(")") :]
    return rebuilt


def main() -> None:
    output_path = os.path.join(BASE_DIR, "get_creation.sql")
    statements: list[str] = []

    for filename in TABLE_FILES:
        csv_path = os.path.join(BASE_DIR, filename)
        table_name, create_sql = load_create_statement(csv_path)
        create_sql = add_foreign_keys(table_name, create_sql)
        create_sql = create_sql.strip()
        if not create_sql.endswith(";"):
            create_sql += ";"
        statements.append(create_sql)

    with open(output_path, "w") as out:
        out.write("\n\n".join(statements) + "\n")


if __name__ == "__main__":
    main()
