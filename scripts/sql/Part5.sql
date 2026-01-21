SELECT 
    g.hgnc_symbol AS gene_symbol,
    v.variant_name,
    p.pdb,
    p.Entry AS uniprot_swissprot_id,
    t.refseq_id,
    t.exonCount AS exon_count,
    geo.logFC
FROM 
    genes g
	LEFT JOIN variants v ON g.gene_stable_id = v.gene_stable_id
    LEFT JOIN proteins p ON g.gene_stable_id = p.ensmble_id
    LEFT JOIN Transcript t ON g.hgnc_symbol = t.name2
    LEFT JOIN geo2r geo ON g.hgnc_symbol = geo.gene_symbol
WHERE 
    g.hgnc_symbol = 'DCC' # DCC , NF2, NBEA
ORDER BY 
	t.exonCount DESC
LIMIT 100;