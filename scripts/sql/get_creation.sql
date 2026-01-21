CREATE TABLE `genes` (
  `Gene_stable_ID` varchar(15) NOT NULL,
  `HGNC_symbol` varchar(10) NOT NULL,
  `Gene_name` text DEFAULT NULL,
  `Chromosome/scaffold name` text DEFAULT NULL,
  `Gene start (bp)` int(11) DEFAULT NULL,
  `Gene end (bp)` int(11) DEFAULT NULL,
  PRIMARY KEY (`Gene_stable_ID`,`HGNC_symbol`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `proteins` (
  `ensmble_id` varchar(30) DEFAULT NULL,
  `Entry` text DEFAULT NULL,
  `Entry_Name` text DEFAULT NULL,
  `Length` int(11) DEFAULT NULL,
  `Protein names` text DEFAULT NULL,
  `Gene Names` text DEFAULT NULL,
  `PDB` text DEFAULT NULL,
  `AlphaFoldDB` text DEFAULT NULL,
  `Mass` int(11) DEFAULT NULL,
  `Subcellular location [CC]` text DEFAULT NULL,
  `Function [CC]` text DEFAULT NULL,
  `Gene Ontology (biological process)` text DEFAULT NULL,
  `InterPro` text DEFAULT NULL,
  `Pfam` text DEFAULT NULL,
  KEY `ensembel_idx` (`ensmble_id`),
  CONSTRAINT `ensembel` FOREIGN KEY (`ensmble_id`) REFERENCES `genes` (`Gene_stable_ID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `variants` (
  `Gene_stable_ID` text DEFAULT NULL,
  `Gene_Name` text DEFAULT NULL,
  `Variant_name` text DEFAULT NULL,
  `Chromosome/scaffold name` int(11) DEFAULT NULL,
  `Chromosome/scaffold position start (bp)` int(11) DEFAULT NULL,
  `Chromosome/scaffold position end (bp)` int(11) DEFAULT NULL,
  `Variant alleles` text DEFAULT NULL,
  `UCSC ID` text DEFAULT NULL
  CONSTRAINT `fk_variants_gene` FOREIGN KEY (`Gene_stable_ID`) REFERENCES `genes` (`Gene_stable_ID`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `geo2r` (
  `geo_id` text DEFAULT NULL,
  `gene_symbol` varchar(500) DEFAULT NULL,
  `logFC` double DEFAULT NULL,
  `p_value` double DEFAULT NULL,
  `adj_p_value` double DEFAULT NULL,
  `t_stat` double DEFAULT NULL,
  `regulation` text DEFAULT NULL,
  KEY `gene_symbol__idx` (`gene_symbol`),
  KEY `gene_symbol_new__idx` (`gene_symbol`,
  CONSTRAINT `fk_geo2r_gene_symbol` FOREIGN KEY (`gene_symbol`) REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION
)
 ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `Transcript` (
  `refseq_id` text DEFAULT NULL,
  `chrom` text DEFAULT NULL,
  `strand` text DEFAULT NULL,
  `txStart` int(11) DEFAULT NULL,
  `txEnd` int(11) DEFAULT NULL,
  `exonCount` int(11) DEFAULT NULL,
  `exonStarts` text DEFAULT NULL,
  `exonEnds` text DEFAULT NULL,
  `score` int(11) DEFAULT NULL,
  `name2` text DEFAULT NULL
  CONSTRAINT `fk_transcript_gene_symbol` FOREIGN KEY (`name2`) REFERENCES `genes` (`HGNC_symbol`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `alphafoldcath` (
  `Alphafold` text DEFAULT NULL,
  `ted_id` text DEFAULT NULL,
  `uniprot_acc` text DEFAULT NULL,
  `md5_domain` text DEFAULT NULL,
  `consensus_level` text DEFAULT NULL,
  `chopping` text DEFAULT NULL,
  `nres_domain` int(11) DEFAULT NULL,
  `num_segments` int(11) DEFAULT NULL,
  `plddt` double DEFAULT NULL,
  `num_helix_strand_turn` int(11) DEFAULT NULL,
  `num_helix` int(11) DEFAULT NULL,
  `num_strand` int(11) DEFAULT NULL,
  `num_helix_strand` int(11) DEFAULT NULL,
  `num_turn` int(11) DEFAULT NULL,
  `proteome_id` int(11) DEFAULT NULL,
  `cath_label` text DEFAULT NULL,
  `cath_assignment_level` text DEFAULT NULL,
  `cath_assignment_method` text DEFAULT NULL,
  `packing_density` double DEFAULT NULL,
  `norm_rg` double DEFAULT NULL,
  `tax_common_name` text DEFAULT NULL,
  `tax_scientific_name` text DEFAULT NULL,
  `tax_lineage` text DEFAULT NULL
  CONSTRAINT `fk_alphafoldcath_protein` FOREIGN KEY (`uniprot_acc`) REFERENCES `proteins` (`Entry`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;

CREATE TABLE `StringInteractions` (
  `#node1` text DEFAULT NULL,
  `node2` text DEFAULT NULL,
  `node1_string_id` text DEFAULT NULL,
  `node2_string_id` text DEFAULT NULL,
  `neighborhood_on_chromosome` double DEFAULT NULL,
  `gene_fusion` double DEFAULT NULL,
  `phylogenetic_cooccurrence` double DEFAULT NULL,
  `homology` double DEFAULT NULL,
  `coexpression` double DEFAULT NULL,
  `experimentally_determined_interaction` double DEFAULT NULL,
  `database_annotated` double DEFAULT NULL,
  `automated_textmining` double DEFAULT NULL,
  `combined_score` double DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_general_ci;
