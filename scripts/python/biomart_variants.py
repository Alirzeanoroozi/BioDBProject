#!/usr/bin/env python3
from __future__ import annotations

import sys
import time
from typing import List, Tuple
import requests

# ============================================================
# DEBUG DEFAULTS (NO TERMINAL INPUT NEEDED)
# Edit these paths for your machine if needed.
# ============================================================

GENE_FILE = r"string_proteins_gene_stable_id_full.txt"
OUT_TSV   = r"biomart_variants.tsv"
OUT_MISS  = r"missing_gene_ids.txt"

DATASET = "hsapiens_snp"
CHUNK_SIZE = 4

MIRRORS = [
    "https://www.ensembl.org/biomart/martservice",
    "https://uswest.ensembl.org/biomart/martservice",
    "https://useast.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice",
]

# From your BioMart XML view:
FILTER_SOURCE_NAME = "variation_source"
FILTER_GENE_NAME   = "ensembl_gene"
FILTER_SOURCE_VALUE = "dbSNP"

ATTRIBUTES = [
    "ensembl_gene_stable_id",
    "ensembl_gene_name",
    "refsnp_id",
    "chr_name",
    "chrom_start",
    "chrom_end",
    "allele",
    "ucsc_id",
]

LOG_FILE = "run.log"

def returned_gene_ids_from_tsv(tsv_text: str) -> set:
    out = set()
    for ln in tsv_text.splitlines():
        if not ln.strip():
            continue
        if ln.startswith("Gene stable ID"):  # skip header if present
            continue
        out.add(ln.split("\t", 1)[0])  # first column is gene stable id
    return out

def log(msg: str) -> None:
    print(msg)
    with open(LOG_FILE, "a", encoding="utf-8", newline="\n") as f:
        f.write(msg + "\n")

def read_gene_ids(path: str) -> List[str]:
    ids = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.lower().startswith("gene") or s.lower().startswith("ensembl"):
                continue
            if s.startswith("ENSG"):
                ids.append(s)

    seen = set()
    out = []
    for g in ids:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out

def build_query_xml(dataset: str, genes: List[str], header: int = 1) -> str:
    gene_list = ",".join(genes)
    attrs = "\n".join([f'    <Attribute name = "{a}" />' for a in ATTRIBUTES])
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName = "default" formatter = "TSV" header = "{header}" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
  <Dataset name = "{dataset}" interface = "default" >
    <Filter name = "{FILTER_SOURCE_NAME}" value = "{FILTER_SOURCE_VALUE}"/>
    <Filter name = "{FILTER_GENE_NAME}" value = "{gene_list}"/>
{attrs}
  </Dataset>
</Query>
"""

def post_query(mirror: str, xml: str, timeout: int = 180) -> str:
    r = requests.post(mirror, data={"query": xml}, timeout=timeout)
    r.raise_for_status()
    return r.text

def is_biomart_error(txt: str) -> bool:
    low = txt.lower()
    return ("query error" in low) or ("biomart::exception" in low) or ("dataset" in low and "not found" in low)

def query_with_mirrors(xml: str, retries: int = 3) -> Tuple[str, str]:
    last_err = None
    for mirror in MIRRORS:
        for attempt in range(1, retries + 1):
            try:
                txt = post_query(mirror, xml)
                if is_biomart_error(txt):
                    raise RuntimeError(txt.splitlines()[0] if txt else "BioMart error")
                return txt, mirror
            except Exception as e:
                last_err = e
                time.sleep(1.2 * attempt)
    raise RuntimeError(f"All mirrors failed. Last error: {last_err}")

def parse_returned_gene_ids(tsv: str) -> set:
    lines = [ln for ln in tsv.splitlines() if ln.strip()]
    if len(lines) < 2:
        return set()
    header = lines[0].split("\t")
    idx = None
    for key in ("Ensembl Gene Stable ID", "ensembl_gene_stable_id"):
        if key in header:
            idx = header.index(key)
            break
    if idx is None:
        return set()
    out = set()
    for ln in lines[1:]:
        parts = ln.split("\t")
        if len(parts) > idx:
            out.add(parts[idx])
    return out

import os
from pathlib import Path

def write_chunk_file(chunk_dir: str, chunk_index: int, txt: str) -> str:
    """
    Writes one chunk output to a separate TSV file.
    Returns the filepath.
    """
    Path(chunk_dir).mkdir(parents=True, exist_ok=True)
    fp = os.path.join(chunk_dir, f"biomart_variants_chunk_{chunk_index:04d}.tsv")
    with open(fp, "w", encoding="utf-8", newline="\n") as f:
        f.write(txt if txt.endswith("\n") else (txt + "\n"))
    return fp

def merge_chunk_files(chunk_files: List[str], out_tsv: str) -> None:
    """
    Merges chunk TSV files into a single TSV.
    Keeps exactly ONE header (from the first chunk file that has a header).
    If chunks have no header (header=0), it merges as-is.
    """
    # Determine if the first non-empty chunk begins with a header line
    header_line = None
    first_data_seen = False

    with open(out_tsv, "w", encoding="utf-8", newline="\n") as out:
        for fp in chunk_files:
            with open(fp, "r", encoding="utf-8") as f:
                lines = f.read().splitlines()

            # Skip completely empty chunk files
            if not any(ln.strip() for ln in lines):
                continue

            # If we haven't decided on a header yet, detect it
            if header_line is None:
                first_nonempty = next((ln for ln in lines if ln.strip()), "")
                if first_nonempty.startswith("Gene stable ID"):
                    header_line = first_nonempty
                    out.write(header_line + "\n")
                    # write remaining lines excluding header
                    for ln in lines[1:]:
                        if ln.strip():
                            out.write(ln + "\n")
                    first_data_seen = True
                    continue
                else:
                    # No header in this file; just write all lines
                    for ln in lines:
                        if ln.strip():
                            out.write(ln + "\n")
                    first_data_seen = True
                    continue

            # Header already known:
            if lines and lines[0].startswith("Gene stable ID"):
                # Skip the header line in subsequent chunks
                for ln in lines[1:]:
                    if ln.strip():
                        out.write(ln + "\n")
            else:
                for ln in lines:
                    if ln.strip():
                        out.write(ln + "\n")

    if not first_data_seen:
        # Ensure output exists but is empty if nothing was merged
        open(out_tsv, "w", encoding="utf-8").close()


# def run():
#     # reset log
#     open(LOG_FILE, "w", encoding="utf-8").close()

#     log(f"Using dataset={DATASET}")
#     log(f"Gene file: {GENE_FILE}")
#     log(f"Output TSV: {OUT_TSV}")
#     log(f"Chunk size: {CHUNK_SIZE}")

#     gene_ids = read_gene_ids(GENE_FILE)
#     if not gene_ids:
#         log("ERROR: No ENSG IDs found in input file.")
#         sys.exit(1)

#     log(f"Input genes: {len(gene_ids)}")
#     log(f"First 3 gene IDs: {gene_ids[:3]}")

#     # reset outputs
#     open(OUT_TSV, "w", encoding="utf-8").close()

#     header_written = False
#     returned_genes_all = set()

#     for i in range(0, len(gene_ids), CHUNK_SIZE):
#         chunk_genes = gene_ids[i:i + CHUNK_SIZE]

#         # header=1 only for first chunk; then header=0
#         xml = build_query_xml(DATASET, chunk_genes, header=1 if not header_written else 0)

#         # For debugging: write the XML being sent
#         with open("last_query.xml", "w", encoding="utf-8") as f:
#             f.write(xml)

#         txt, used = query_with_mirrors(xml)
#         log(f"Chunk {i//CHUNK_SIZE + 1}: genes {i+1}-{min(i+CHUNK_SIZE, len(gene_ids))} OK (mirror: {used})")
#         log(f"Response length: {len(txt)} chars")

#         with open(OUT_TSV, "a", encoding="utf-8", newline="\n") as f:
#             f.write(txt if txt.endswith("\n") else (txt + "\n"))

#         # track returned genes from first chunk only (has header)
#         if not header_written:
#             returned_genes_all |= parse_returned_gene_ids(txt)
#             header_written = True

#         time.sleep(0.4)

#     # If TSV is empty, log it explicitly
#     try:
#         size = __import__("os").path.getsize(OUT_TSV)
#         log(f"Final TSV size: {size} bytes")
#         if size < 5:
#             log("WARNING: TSV file is essentially empty. Check last_query.xml and run.log.")
#     except Exception as e:
#         log(f"Could not stat output TSV: {e}")

#     missing = [g for g in gene_ids if g not in returned_genes_all]
#     with open(OUT_MISS, "w", encoding="utf-8", newline="\n") as f:
#         f.write("\n".join(missing) + ("\n" if missing else ""))

#     log(f"Genes not observed in first-chunk tracking: {len(missing)} -> {OUT_MISS}")
#     log("DONE.")

def run():
    import os

    # reset log
    open(LOG_FILE, "w", encoding="utf-8").close()

    log(f"Using dataset={DATASET}")
    log(f"Gene file: {GENE_FILE}")
    log(f"Output TSV: {OUT_TSV}")
    log(f"Chunk size: {CHUNK_SIZE}")

    gene_ids = read_gene_ids(GENE_FILE)
    if not gene_ids:
        log("ERROR: No ENSG IDs found in input file.")
        sys.exit(1)

    log(f"Input genes: {len(gene_ids)}")
    log(f"First 3 gene IDs: {gene_ids[:3]}")

    # Directory to store chunk results
    chunk_dir = "chunks"

    # Clean old chunk files (avoid mixing old + new)
    if os.path.isdir(chunk_dir):
        for fn in os.listdir(chunk_dir):
            if fn.startswith("biomart_variants_chunk_") and fn.endswith(".tsv"):
                try:
                    os.remove(os.path.join(chunk_dir, fn))
                except Exception:
                    pass
    else:
        os.makedirs(chunk_dir, exist_ok=True)

    returned_genes_all = set()
    chunk_files: List[str] = []
    failed_chunks: List[int] = []

    # We will request header=1 for EVERY chunk file.
    # Reason: simplifies merge + missing gene detection.
    # We will remove duplicate headers at merge time.
    for idx, i in enumerate(range(0, len(gene_ids), CHUNK_SIZE), start=1):
        chunk_genes = gene_ids[i:i + CHUNK_SIZE]
        xml = build_query_xml(DATASET, chunk_genes, header=1)

        # keep last query for debugging
        with open("last_query.xml", "w", encoding="utf-8") as f:
            f.write(xml)

        try:
            txt, used = query_with_mirrors(xml)
            log(f"Chunk {idx}: genes {i+1}-{min(i+CHUNK_SIZE, len(gene_ids))} OK (mirror: {used})")
            log(f"Response length: {len(txt)} chars")

            # Save this chunk output separately
            fp = write_chunk_file(chunk_dir, idx, txt)
            chunk_files.append(fp)

            # Update returned genes set from the real returned TSV text
            returned_genes_all |= returned_gene_ids_from_tsv(txt)

        except Exception as e:
            failed_chunks.append(idx)
            log(f"ERROR: Chunk {idx} failed. Reason: {e}")

        time.sleep(0.4)

    # Merge chunk files into final output
    if chunk_files:
        merge_chunk_files(chunk_files, OUT_TSV)
        try:
            size = __import__("os").path.getsize(OUT_TSV)
            log(f"Final merged TSV size: {size} bytes")
        except Exception as e:
            log(f"Could not stat merged output TSV: {e}")
    else:
        # no chunks succeeded
        open(OUT_TSV, "w", encoding="utf-8").close()
        log("WARNING: No chunk files produced. Output TSV is empty.")

    # Missing genes = those never observed in returned data
    missing = [g for g in gene_ids if g not in returned_genes_all]
    with open(OUT_MISS, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(missing) + ("\n" if missing else ""))

    log(f"Genes with no returned rows: {len(missing)} -> {OUT_MISS}")
    if failed_chunks:
        log(f"Failed chunks: {failed_chunks} (you can re-run; successful chunks are preserved in {chunk_dir}/)")
    log("DONE.")


if __name__ == "__main__":
    run()
