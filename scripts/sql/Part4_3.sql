-- Query 4.3: Retrieve genes with transcript information and expression data

SELECT 
    g.hgnc_symbol,
    g.Chromosome_scaffold,
    COUNT(DISTINCT t.refseq_id) AS transcript_count,
    AVG(t.exonCount) AS avg_exon_count,
    MIN(t.exonCount) AS min_exon_count,
    MAX(t.exonCount) AS max_exon_count,
    geo.logFC,
    geo.adj_p_value
FROM 
    Genes g
    LEFT JOIN Transcript t ON g.hgnc_symbol = t.name2
    LEFT JOIN geo2r geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    geo.adj_p_value IS NOT NULL
GROUP BY 
    g.hgnc_symbol, g.Chromosome_scaffold, geo.logFC, geo.adj_p_value
HAVING 
    COUNT(DISTINCT t.refseq_id) > 0
ORDER BY 
    transcript_count DESC, ABS(geo.logFC) DESC
LIMIT 50;
