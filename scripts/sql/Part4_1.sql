-- Query 4.1: Retrieve genes with their associated variants and GEO2R expression data

SELECT DISTINCT
    g.gene_stable_id,
    g.hgnc_symbol,
    COUNT(DISTINCT v.variant_name) AS variant_count,
    geo.logFC,
    geo.adj_p_value,
    geo.p_value
FROM 
    Genes g
    LEFT JOIN variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN geo2r geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    geo.adj_p_value IS NOT NULL
GROUP BY 
    g.gene_stable_id, g.hgnc_symbol, g.gene_name, geo.logFC, geo.adj_p_value, geo.p_value
HAVING 
    COUNT(DISTINCT v.variant_name) > 0
ORDER BY 
    variant_count DESC, ABS(geo.logFC) DESC
LIMIT 50;

