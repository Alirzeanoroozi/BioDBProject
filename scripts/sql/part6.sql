SELECT DISTINCT
    g.hgnc_symbol,
    g.gene_stable_id,
    p.Entry,
    p.pdb,
    t.refseq_id,
    t.exonCount
FROM 
    Genes g
    INNER JOIN proteins p ON g.gene_stable_id = p.ensmble_id
    INNER JOIN Transcript t ON g.hgnc_symbol = t.name2
WHERE 
    p.pdb IS NOT NULL 
    AND p.pdb != ''
    AND t.exonCount >= 4
ORDER BY 
    t.exonCount DESC, g.hgnc_symbol
LIMIT 100;
