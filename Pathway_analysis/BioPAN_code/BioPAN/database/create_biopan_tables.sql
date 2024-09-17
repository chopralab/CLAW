DROP TABLE IF EXISTS biopan_gene CASCADE;
CREATE TABLE biopan_gene (
    gene_id INTEGER PRIMARY KEY,
    gene_symbol character varying(30)
);
ALTER TABLE biopan_gene OWNER TO biopan;


DROP TABLE IF EXISTS biopan_reaction CASCADE;
CREATE TABLE biopan_reaction (
    reaction_id INTEGER PRIMARY KEY,
    reaction character varying(20),
    class character varying(100),
    sub_class character varying(100),
    species_long character varying(100),
    type character varying(5),
    acyl_add boolean default true,
    compound_require character varying(10)
);
ALTER TABLE biopan_reaction OWNER TO biopan;


DROP TABLE IF EXISTS biopan_reaction_gene;
CREATE TABLE biopan_reaction_gene (
    reaction_id INTEGER NOT NULL,
    gene_id INTEGER REFERENCES biopan_gene(gene_id) ON UPDATE CASCADE ON DELETE CASCADE,
    PRIMARY KEY (reaction_id , gene_id)
);
ALTER TABLE biopan_reaction_gene OWNER TO biopan;


DROP TABLE IF EXISTS biopan_pathway;        
CREATE TABLE biopan_pathway(
    pathway_batch_id INTEGER NOT NULL,
    pathway_id INTEGER NOT NULL,
    sort_id INTEGER NOT NULL,
    reaction_id INTEGER
        REFERENCES biopan_reaction (reaction_id)
        ON UPDATE CASCADE ON DELETE CASCADE,
    
    name character varying(100),
    descr character varying(100),
    PRIMARY KEY (pathway_batch_id, pathway_id, reaction_id)
);
ALTER TABLE biopan_pathway OWNER TO biopan;