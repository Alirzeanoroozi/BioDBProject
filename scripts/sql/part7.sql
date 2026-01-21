SELECT 
    COUNT(DISTINCT gene_symbol) AS significant_gene_count
FROM
    geo2r
WHERE 
    adj_p_value < 0.05;
    
SELECT 
    DISTINCT gene_symbol AS significant_genes
FROM
    geo2r
WHERE 
    adj_p_value < 0.05;

SELECT DISTINCT
    geo.gene_symbol,
    p.Entry,
    p.Length,
    geo.logFC,
    geo.regulation
FROM 
    geo2r geo
    INNER JOIN Genes g ON geo.gene_symbol = g.HGNC_symbol
    INNER JOIN Proteins p ON g.Gene_stable_ID = p.ensmble_id
WHERE 
    geo.adj_p_value < 0.05;