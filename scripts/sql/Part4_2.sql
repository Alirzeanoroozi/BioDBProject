-- Query 4.2: Find proteins with their structural domains and associated genes

SELECT 
    g.hgnc_symbol,
    p.Entry,
    p.entry_name,
    p.length AS protein_length,
    p.pdb,
    COUNT(DISTINCT ac.ted_id) AS domain_count,
    GROUP_CONCAT(DISTINCT ac.cath_label SEPARATOR '; ') AS cath_labels
FROM 
    Genes g
    INNER JOIN Proteins p ON g.gene_stable_id = p.ensmble_id
    LEFT JOIN AlphaFoldCATH ac ON p.Entry = ac.Alphafold
GROUP BY 
    g.hgnc_symbol, p.Entry, p.entry_name, p.length, p.pdb
HAVING 
    COUNT(DISTINCT ac.ted_id) > 0
ORDER BY 
    domain_count DESC, protein_length DESC
LIMIT 50;
