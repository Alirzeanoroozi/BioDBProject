-- ============================================================================
-- CAMM 535 Final Project - Schwannoma Database
-- Database Creation, Data Import, and SQL Queries for Final Report 2
-- Group 4: Schwannoma
-- ============================================================================

-- ============================================================================
-- DATABASE CREATION
-- ============================================================================

CREATE DATABASE IF NOT EXISTS schwannoma_db CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
USE schwannoma_db;

-- ============================================================================
-- TABLE CREATION
-- ============================================================================

-- Table: Genes
-- Source: BioMart/Ensembl (data/part2/biomart_schwannoma_genes.tsv)
CREATE TABLE IF NOT EXISTS Genes (
    gene_stable_id VARCHAR(20) NOT NULL,
    hgnc_symbol VARCHAR(50),
    gene_name VARCHAR(255),
    chromosome VARCHAR(10),
    gene_start BIGINT,
    gene_end BIGINT,
    uniprot_gene_name_id VARCHAR(50),
    PRIMARY KEY (gene_stable_id),
    UNIQUE KEY idx_hgnc_symbol (hgnc_symbol),
    KEY idx_chromosome (chromosome)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: Variants
-- Source: BioMart dbSNP (data/part2/biomart_variants.tsv)
CREATE TABLE IF NOT EXISTS Variants (
    gene_stable_id VARCHAR(20) NOT NULL,
    gene_name VARCHAR(255),
    refsnp_id VARCHAR(20) NOT NULL,
    chr_name VARCHAR(10),
    chrom_start BIGINT,
    chrom_end BIGINT,
    allele VARCHAR(50),
    ucsc_id VARCHAR(50),
    PRIMARY KEY (gene_stable_id, refsnp_id),
    FOREIGN KEY (gene_stable_id) REFERENCES Genes(gene_stable_id) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    KEY idx_refsnp (refsnp_id),
    KEY idx_chr_pos (chr_name, chrom_start)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: Transcripts
-- Source: UCSC Table Browser (data/part3/ucsc-table-browser-output.tsv)
CREATE TABLE IF NOT EXISTS Transcripts (
    refseq_accession VARCHAR(50) NOT NULL,
    chrom VARCHAR(10),
    strand CHAR(1),
    txStart BIGINT,
    txEnd BIGINT,
    exonCount INT,
    exonStarts TEXT,
    exonEnds TEXT,
    score INT,
    gene_symbol VARCHAR(50),
    PRIMARY KEY (refseq_accession),
    FOREIGN KEY (gene_symbol) REFERENCES Genes(hgnc_symbol) 
        ON DELETE SET NULL ON UPDATE CASCADE,
    KEY idx_gene_symbol (gene_symbol),
    KEY idx_chrom_pos (chrom, txStart)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: Proteins
-- Source: UniProt ID Mapping (data/part4/idmapping.tsv)
CREATE TABLE IF NOT EXISTS Proteins (
    uniprot_id VARCHAR(20) NOT NULL,
    ensembl_gene_id VARCHAR(20),
    entry_name VARCHAR(50),
    length INT,
    protein_names TEXT,
    gene_names VARCHAR(255),
    pdb_ids TEXT,
    alphafolddb_ids TEXT,
    mass INT,
    subcellular_location TEXT,
    function_description TEXT,
    gene_ontology TEXT,
    interpro TEXT,
    pfam TEXT,
    PRIMARY KEY (uniprot_id),
    FOREIGN KEY (ensembl_gene_id) REFERENCES Genes(gene_stable_id) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    KEY idx_ensembl_gene (ensembl_gene_id),
    KEY idx_entry_name (entry_name)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: STRING_Interactions
-- Source: STRING database (data/part2/string_interactions.tsv)
CREATE TABLE IF NOT EXISTS STRING_Interactions (
    node1 VARCHAR(50) NOT NULL,
    node2 VARCHAR(50) NOT NULL,
    node1_string_id VARCHAR(100),
    node2_string_id VARCHAR(100),
    neighborhood_on_chromosome DECIMAL(5,3),
    gene_fusion DECIMAL(5,3),
    phylogenetic_cooccurrence DECIMAL(5,3),
    homology DECIMAL(5,3),
    coexpression DECIMAL(5,3),
    experimentally_determined_interaction DECIMAL(5,3),
    database_annotated DECIMAL(5,3),
    automated_textmining DECIMAL(5,3),
    combined_score DECIMAL(5,3),
    PRIMARY KEY (node1, node2),
    FOREIGN KEY (node1) REFERENCES Genes(hgnc_symbol) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    FOREIGN KEY (node2) REFERENCES Genes(hgnc_symbol) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    KEY idx_node1 (node1),
    KEY idx_node2 (node2),
    KEY idx_combined_score (combined_score)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: AlphaFold_CATH
-- Source: AlphaFold/CATH database (data/part6/alphafold_cath.csv)
CREATE TABLE IF NOT EXISTS AlphaFold_CATH (
    ted_id VARCHAR(100) NOT NULL,
    uniprot_acc VARCHAR(20),
    alphafold VARCHAR(20),
    md5_domain VARCHAR(32),
    consensus_level VARCHAR(20),
    chopping VARCHAR(100),
    nres_domain INT,
    num_segments INT,
    plddt DECIMAL(6,4),
    num_helix_strand_turn INT,
    num_helix INT,
    num_strand INT,
    num_helix_strand INT,
    num_turn INT,
    proteome_id INT,
    cath_label VARCHAR(50),
    cath_assignment_level VARCHAR(10),
    cath_assignment_method VARCHAR(20),
    packing_density DECIMAL(8,3),
    norm_rg DECIMAL(6,3),
    tax_common_name VARCHAR(50),
    tax_scientific_name VARCHAR(100),
    tax_lineage TEXT,
    interactions TEXT,
    PRIMARY KEY (ted_id),
    FOREIGN KEY (uniprot_acc) REFERENCES Proteins(uniprot_id) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    KEY idx_uniprot (uniprot_acc),
    KEY idx_cath_label (cath_label)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- Table: GEO2R_Results
-- Source: GEO2R analysis (data/part5/GSE39645 (GEO2R results).tsv)
CREATE TABLE IF NOT EXISTS GEO2R_Results (
    geo_result_id INT AUTO_INCREMENT,
    probe_id VARCHAR(50),
    adj_p_value DECIMAL(20,10),
    p_value DECIMAL(20,10),
    t_statistic DECIMAL(10,6),
    b_statistic DECIMAL(10,6),
    logFC DECIMAL(10,6),
    gene_symbol VARCHAR(50),
    gene_title VARCHAR(500),
    PRIMARY KEY (geo_result_id),
    FOREIGN KEY (gene_symbol) REFERENCES Genes(hgnc_symbol) 
        ON DELETE CASCADE ON UPDATE CASCADE,
    KEY idx_gene_symbol (gene_symbol),
    KEY idx_adj_p_value (adj_p_value),
    KEY idx_logFC (logFC)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;

-- ============================================================================
-- DATA IMPORT
-- ============================================================================

-- IMPORTANT NOTES:
-- 1. File paths are relative to the MySQL server's data directory or absolute paths
-- 2. For phpMyAdmin: Use the Import tab instead of LOAD DATA INFILE
-- 3. For command line: Ensure files are accessible to MySQL server
-- 4. If LOAD DATA LOCAL INFILE is disabled, enable it with:
--    SET GLOBAL local_infile = 1;
--    Or add to my.cnf: [mysql] local_infile=1
-- 5. Adjust file paths according to your system (absolute paths recommended)

-- Enable local file loading (if needed)
SET GLOBAL local_infile = 1;

-- Import Genes data
-- Source: data/part2/biomart_schwannoma_genes.tsv
LOAD DATA LOCAL INFILE 'data/part2/biomart_schwannoma_genes.tsv'
INTO TABLE Genes
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(gene_stable_id, hgnc_symbol, gene_name, chromosome, gene_start, gene_end, uniprot_gene_name_id);

-- Import Variants data
-- Source: data/part2/biomart_variants.tsv
-- Column order: Gene stable ID, Gene Name, Variant name, Chromosome/scaffold name, 
--               Chromosome/scaffold position start (bp), Chromosome/scaffold position end (bp), 
--               Variant alleles, UCSC ID
LOAD DATA LOCAL INFILE 'data/part2/biomart_variants.tsv'
INTO TABLE Variants
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(gene_stable_id, gene_name, refsnp_id, chr_name, chrom_start, chrom_end, allele, ucsc_id);

-- Import Transcripts data
-- Source: data/part3/ucsc-table-browser-output.tsv
-- Note: First column header starts with #name, but the data rows don't have #
-- The IGNORE 1 ROWS will skip the header row
LOAD DATA LOCAL INFILE 'data/part3/ucsc-table-browser-output.tsv'
INTO TABLE Transcripts
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(refseq_accession, chrom, strand, txStart, txEnd, exonCount, exonStarts, exonEnds, score, gene_symbol);

-- Import Proteins data
-- Source: data/part4/idmapping.tsv
LOAD DATA LOCAL INFILE 'data/part4/idmapping.tsv'
INTO TABLE Proteins
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(ensembl_gene_id, uniprot_id, entry_name, length, protein_names, gene_names, pdb_ids, alphafolddb_ids, mass, subcellular_location, function_description, gene_ontology, interpro, pfam);

-- Import STRING Interactions data
-- Source: data/part2/string_interactions.tsv
-- Note: First row starts with #, so we ignore it
LOAD DATA LOCAL INFILE 'data/part2/string_interactions.tsv'
INTO TABLE STRING_Interactions
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(node1, node2, node1_string_id, node2_string_id, neighborhood_on_chromosome, gene_fusion, phylogenetic_cooccurrence, homology, coexpression, experimentally_determined_interaction, database_annotated, automated_textmining, combined_score);

-- Import AlphaFold_CATH data
-- Source: data/part6/alphafold_cath.csv
LOAD DATA LOCAL INFILE 'data/part6/alphafold_cath.csv'
INTO TABLE AlphaFold_CATH
FIELDS TERMINATED BY ','
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(alphafold, ted_id, uniprot_acc, md5_domain, consensus_level, chopping, nres_domain, num_segments, plddt, num_helix_strand_turn, num_helix, num_strand, num_helix_strand, num_turn, proteome_id, cath_label, cath_assignment_level, cath_assignment_method, packing_density, norm_rg, tax_common_name, tax_scientific_name, tax_lineage, interactions);

-- Import GEO2R Results data
-- Source: data/part5/GSE39645 (GEO2R results).tsv
-- Column order: ID, adj.P.Val, P.Value, t, B, logFC, Gene.symbol, Gene.title
LOAD DATA LOCAL INFILE 'data/part5/GSE39645 (GEO2R results).tsv'
INTO TABLE GEO2R_Results
FIELDS TERMINATED BY '\t'
ENCLOSED BY '"'
ESCAPED BY '\\'
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(probe_id, adj_p_value, p_value, t_statistic, b_statistic, logFC, gene_symbol, gene_title);

-- ============================================================================
-- ALTERNATIVE IMPORT METHOD (for phpMyAdmin or when LOAD DATA INFILE is disabled)
-- ============================================================================

-- If LOAD DATA INFILE is not available, use phpMyAdmin Import feature:
-- 
-- STEP-BY-STEP INSTRUCTIONS FOR phpMyAdmin:
-- 1. Go to phpMyAdmin
-- 2. Select the schwannoma_db database
-- 3. Click on the table name (e.g., "Genes")
-- 4. Click "Import" tab
-- 5. Choose the file and set:
--    - Format: CSV
--    - Columns separated with: Tab (for TSV files) or Comma (for CSV files)
--    - Columns enclosed with: "
--    - Columns escaped with: \
--    - Lines terminated with: Auto
--    - Check "Replace table data with file" (if re-importing)
--    - Check "The first line contains the table column names"
-- 6. Click "Go" to import
--
-- IMPORT ORDER (to respect foreign key constraints):
-- 1. Genes (no dependencies)
-- 2. Variants (depends on Genes)
-- 3. Transcripts (depends on Genes)
-- 4. Proteins (depends on Genes)
-- 5. STRING_Interactions (depends on Genes)
-- 6. AlphaFold_CATH (depends on Proteins)
-- 7. GEO2R_Results (depends on Genes)
--
-- ALTERNATIVE: Use MySQL Workbench or command line with absolute paths:
-- Example: LOAD DATA LOCAL INFILE '/full/path/to/data/part2/biomart_schwannoma_genes.tsv'

-- ============================================================================
-- DATA VALIDATION QUERIES
-- ============================================================================

-- Check record counts after import
SELECT 'Genes' AS table_name, COUNT(*) AS record_count FROM Genes
UNION ALL
SELECT 'Variants', COUNT(*) FROM Variants
UNION ALL
SELECT 'Transcripts', COUNT(*) FROM Transcripts
UNION ALL
SELECT 'Proteins', COUNT(*) FROM Proteins
UNION ALL
SELECT 'STRING_Interactions', COUNT(*) FROM STRING_Interactions
UNION ALL
SELECT 'AlphaFold_CATH', COUNT(*) FROM AlphaFold_CATH
UNION ALL
SELECT 'GEO2R_Results', COUNT(*) FROM GEO2R_Results;

-- Check for missing foreign key references
SELECT 'Genes without matching Variants' AS issue, COUNT(*) AS count
FROM Genes g LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
WHERE v.gene_stable_id IS NULL
UNION ALL
SELECT 'Proteins without matching Genes', COUNT(*)
FROM Proteins p LEFT JOIN Genes g ON p.ensembl_gene_id = g.gene_stable_id
WHERE g.gene_stable_id IS NULL
UNION ALL
SELECT 'Transcripts without matching Genes', COUNT(*)
FROM Transcripts t LEFT JOIN Genes g ON t.gene_symbol = g.hgnc_symbol
WHERE g.hgnc_symbol IS NULL AND t.gene_symbol IS NOT NULL;

-- ============================================================================
-- SQL QUERIES FOR FINAL REPORT 2
-- ============================================================================

-- ============================================================================
-- PART 4: Queries Requiring at Least 2 Tables (9 points)
-- ============================================================================

-- Query 4.1: Retrieve genes with their associated variants and GEO2R expression data
-- This query joins Genes, Variants, and GEO2R_Results tables to show
-- genes that have both genetic variants and differential expression data
SELECT DISTINCT
    g.gene_stable_id,
    g.hgnc_symbol,
    g.gene_name,
    COUNT(DISTINCT v.refsnp_id) AS variant_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    geo.adj_p_value IS NOT NULL
GROUP BY 
    g.gene_stable_id, g.hgnc_symbol, g.gene_name, geo.logFC, geo.adj_p_value, geo.p_value
HAVING 
    COUNT(DISTINCT v.refsnp_id) > 0
ORDER BY 
    variant_count DESC, ABS(geo.logFC) DESC
LIMIT 50;

-- Query 4.2: Find proteins with their structural domains and associated genes
-- This query joins Proteins, AlphaFold_CATH, and Genes tables to show
-- protein structural information linked to gene annotations
SELECT 
    g.hgnc_symbol,
    g.gene_name,
    p.uniprot_id,
    p.entry_name,
    p.length AS protein_length,
    p.pdb_ids,
    COUNT(DISTINCT ac.ted_id) AS domain_count,
    GROUP_CONCAT(DISTINCT ac.cath_label SEPARATOR '; ') AS cath_labels
FROM 
    Genes g
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    LEFT JOIN AlphaFold_CATH ac ON p.uniprot_id = ac.uniprot_acc
GROUP BY 
    g.hgnc_symbol, g.gene_name, p.uniprot_id, p.entry_name, p.length, p.pdb_ids
HAVING 
    COUNT(DISTINCT ac.ted_id) > 0
ORDER BY 
    domain_count DESC, protein_length DESC
LIMIT 50;

-- Query 4.3: Retrieve genes with transcript information and expression data
-- This query joins Genes, Transcripts, and GEO2R_Results to show
-- transcript structure and expression levels for each gene
SELECT 
    g.hgnc_symbol,
    g.gene_name,
    g.chromosome,
    COUNT(DISTINCT t.refseq_accession) AS transcript_count,
    AVG(t.exonCount) AS avg_exon_count,
    MIN(t.exonCount) AS min_exon_count,
    MAX(t.exonCount) AS max_exon_count,
    geo.logFC,
    geo.adj_p_value
FROM 
    Genes g
    LEFT JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    geo.adj_p_value IS NOT NULL
GROUP BY 
    g.hgnc_symbol, g.gene_name, g.chromosome, geo.logFC, geo.adj_p_value
HAVING 
    COUNT(DISTINCT t.refseq_accession) > 0
ORDER BY 
    transcript_count DESC, ABS(geo.logFC) DESC
LIMIT 50;

-- ============================================================================
-- PART 5: Comprehensive Gene Information Query (6 points)
-- Retrieve variation name, PDB IDs, UniProt/SwissProt IDs, RefSeq IDs, 
-- exon counts, and logFC value for a given gene symbol
-- ============================================================================

-- Query 5.1: For gene symbol 'NF2' (Merlin - key schwannoma gene)
SELECT 
    g.hgnc_symbol AS gene_symbol,
    g.gene_name,
    v.refsnp_id AS variation_name,
    p.pdb_ids,
    p.uniprot_id AS uniprot_swissprot_id,
    t.refseq_accession AS refseq_id,
    t.exonCount AS exon_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    LEFT JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    g.hgnc_symbol = 'NF2'
ORDER BY 
    t.exonCount DESC, v.refsnp_id
LIMIT 100;

-- Query 5.2: For gene symbol 'PIK3CB' (PI3K-beta - signaling pathway)
SELECT 
    g.hgnc_symbol AS gene_symbol,
    g.gene_name,
    v.refsnp_id AS variation_name,
    p.pdb_ids,
    p.uniprot_id AS uniprot_swissprot_id,
    t.refseq_accession AS refseq_id,
    t.exonCount AS exon_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    LEFT JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    g.hgnc_symbol = 'PIK3CB'
ORDER BY 
    t.exonCount DESC, v.refsnp_id
LIMIT 100;

-- Query 5.3: For gene symbol 'FGFR2' (Fibroblast growth factor receptor)
SELECT 
    g.hgnc_symbol AS gene_symbol,
    g.gene_name,
    v.refsnp_id AS variation_name,
    p.pdb_ids,
    p.uniprot_id AS uniprot_swissprot_id,
    t.refseq_accession AS refseq_id,
    t.exonCount AS exon_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    LEFT JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    g.hgnc_symbol = 'FGFR2'
ORDER BY 
    t.exonCount DESC, v.refsnp_id
LIMIT 100;

-- Alternative: Generic query template for any gene symbol
-- Replace 'GENE_SYMBOL' with the desired gene symbol
/*
SELECT 
    g.hgnc_symbol AS gene_symbol,
    g.gene_name,
    v.refsnp_id AS variation_name,
    p.pdb_ids,
    p.uniprot_id AS uniprot_swissprot_id,
    t.refseq_accession AS refseq_id,
    t.exonCount AS exon_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    LEFT JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
    LEFT JOIN GEO2R_Results geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    g.hgnc_symbol = 'GENE_SYMBOL'
ORDER BY 
    t.exonCount DESC, v.refsnp_id;
*/

-- ============================================================================
-- PART 6: Genes with PDB IDs and Exon Count >= 4 (3 points)
-- ============================================================================

SELECT DISTINCT
    g.hgnc_symbol,
    g.gene_name,
    g.gene_stable_id,
    p.uniprot_id,
    p.pdb_ids,
    t.refseq_accession,
    t.exonCount
FROM 
    Genes g
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
    INNER JOIN Transcripts t ON g.hgnc_symbol = t.gene_symbol
WHERE 
    p.pdb_ids IS NOT NULL 
    AND p.pdb_ids != ''
    AND t.exonCount >= 4
ORDER BY 
    t.exonCount DESC, g.hgnc_symbol
LIMIT 100;

-- ============================================================================
-- PART 7: GEO2R Significant Genes Analysis (7 points)
-- ============================================================================

-- Part 7(a): Count genes with adjusted p-value < 0.05
SELECT 
    COUNT(DISTINCT gene_symbol) AS significant_gene_count
FROM 
    GEO2R_Results
WHERE 
    adj_p_value < 0.05;

-- Part 7(a) Alternative: If count < 30, use less stringent cutoff (adj_p_value < 0.1)
-- This query shows genes with less stringent significance threshold
SELECT 
    COUNT(DISTINCT gene_symbol) AS significant_gene_count_relaxed
FROM 
    GEO2R_Results
WHERE 
    adj_p_value < 0.1;

-- Part 7(a) Alternative: Top 100 genes by absolute logFC
-- This query retrieves the top 100 most differentially expressed genes
SELECT 
    gene_symbol,
    logFC,
    adj_p_value,
    p_value
FROM 
    GEO2R_Results
WHERE 
    adj_p_value IS NOT NULL
ORDER BY 
    ABS(logFC) DESC
LIMIT 100;

-- Part 7(b): Retrieve UniProt IDs for genes with adj_p_value < 0.05
SELECT DISTINCT
    geo.gene_symbol,
    g.gene_name,
    p.uniprot_id,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    GEO2R_Results geo
    INNER JOIN Genes g ON geo.gene_symbol = g.hgnc_symbol
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
WHERE 
    geo.adj_p_value < 0.05
ORDER BY 
    ABS(geo.logFC) DESC, geo.adj_p_value ASC;

-- Part 7(b) Alternative: If using relaxed cutoff (adj_p_value < 0.1)
SELECT DISTINCT
    geo.gene_symbol,
    g.gene_name,
    p.uniprot_id,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    GEO2R_Results geo
    INNER JOIN Genes g ON geo.gene_symbol = g.hgnc_symbol
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
WHERE 
    geo.adj_p_value < 0.1
ORDER BY 
    ABS(geo.logFC) DESC, geo.adj_p_value ASC
LIMIT 100;

-- Part 7(b) Alternative: Top 100 genes by absolute logFC with UniProt IDs
SELECT DISTINCT
    geo.gene_symbol,
    g.gene_name,
    p.uniprot_id,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    GEO2R_Results geo
    INNER JOIN Genes g ON geo.gene_symbol = g.hgnc_symbol
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensembl_gene_id
WHERE 
    geo.adj_p_value IS NOT NULL
ORDER BY 
    ABS(geo.logFC) DESC
LIMIT 100;

-- ============================================================================
-- ADDITIONAL USEFUL QUERIES
-- ============================================================================

-- Query: Find genes with high interaction scores in STRING network
SELECT 
    si.node1,
    si.node2,
    si.combined_score,
    g1.gene_name AS gene1_name,
    g2.gene_name AS gene2_name
FROM 
    STRING_Interactions si
    INNER JOIN Genes g1 ON si.node1 = g1.hgnc_symbol
    INNER JOIN Genes g2 ON si.node2 = g2.hgnc_symbol
WHERE 
    si.combined_score > 0.7
ORDER BY 
    si.combined_score DESC
LIMIT 50;

-- Query: Genes with variants on specific chromosomes
SELECT 
    g.hgnc_symbol,
    g.gene_name,
    v.chr_name,
    COUNT(DISTINCT v.refsnp_id) AS variant_count
FROM 
    Genes g
    INNER JOIN Variants v ON g.gene_stable_id = v.gene_stable_id
WHERE 
    v.chr_name IN ('22', 'X')  -- Chromosomes commonly associated with schwannoma
GROUP BY 
    g.hgnc_symbol, g.gene_name, v.chr_name
ORDER BY 
    variant_count DESC;

-- Query: Up-regulated vs Down-regulated genes summary
SELECT 
    CASE 
        WHEN logFC > 0 THEN 'Up-regulated'
        WHEN logFC < 0 THEN 'Down-regulated'
        ELSE 'No change'
    END AS regulation_status,
    COUNT(DISTINCT gene_symbol) AS gene_count,
    AVG(logFC) AS avg_logFC,
    AVG(adj_p_value) AS avg_adj_p_value
FROM 
    GEO2R_Results
WHERE 
    adj_p_value < 0.05
GROUP BY 
    CASE 
        WHEN logFC > 0 THEN 'Up-regulated'
        WHEN logFC < 0 THEN 'Down-regulated'
        ELSE 'No change'
    END;

-- Query: Proteins with CATH domain annotations
SELECT 
    p.uniprot_id,
    p.entry_name,
    g.hgnc_symbol,
    COUNT(DISTINCT ac.ted_id) AS domain_count,
    GROUP_CONCAT(DISTINCT ac.cath_label SEPARATOR '; ') AS cath_classifications
FROM 
    Proteins p
    INNER JOIN Genes g ON p.ensembl_gene_id = g.gene_stable_id
    INNER JOIN AlphaFold_CATH ac ON p.uniprot_id = ac.uniprot_acc
WHERE 
    ac.cath_label IS NOT NULL 
    AND ac.cath_label != ''
GROUP BY 
    p.uniprot_id, p.entry_name, g.hgnc_symbol
ORDER BY 
    domain_count DESC
LIMIT 50;

-- ============================================================================
-- END OF SQL QUERIES
-- ============================================================================
