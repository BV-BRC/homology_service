CREATE TABLE FeatureType
(
    id char PRIMARY KEY,
    name varchar(32)
);
INSERT INTO FeatureType (id, name) VALUES ('f', 'features'), ('c', 'contigs');

CREATE TABLE DbType
(
    id char PRIMARY KEY,
    name VARCHAR(10),
    suffix VARCHAR(32),
    blast_type VARCHAR(32)
);

INSERT INTO DbType (id, name, suffix, blast_type) VALUES
       ('d', 'DNA features', 'ffn', 'nucl'),
       ('c', 'DNA contigs', 'fna', 'nucl'),
       ('a', 'AA features', 'faa', 'prot');

CREATE TABLE GenomeGroup
(
    id integer primary key,
    name varchar(255) UNIQUE,
    curated BOOLEAN DEFAULT FALSE,
    taxonomy_level varchar(32),
    superkingdom varchar(32)
);

DROP TABLE IF EXISTS BlastDatabase;
CREATE TABLE BlastDatabase
(
    id integer primary key,
    genome_group integer NOT NULL,
    path varchar(255),
    dbtype varchar(32) NOT NULL,
    ftype varchar(32) NOT NULL,
    FOREIGN KEY (genome_group) REFERENCES GenomeGroup(id),
    FOREIGN KEY (dbtype) REFERENCES DbType(id),
    FOREIGN KEY (ftype) REFERENCES FeatureType(id)
);

CREATE TABLE GenomeInDatabase
(
    genome_id VARCHAR(32),
    blast_database integer,
    PRIMARY KEY (genome_id, blast_database),
    FOREIGN KEY (blast_database) REFERENCES BlastDatabase
);

CREATE TABLE TaxonInDatabase
(
    taxon_id INTEGER,
    blast_database integer,
    PRIMARY KEY (taxon_id, blast_database),
    FOREIGN KEY (blast_database) REFERENCES BlastDatabase
);

CREATE TABLE TaxNode
(
    tax_id INTEGER PRIMARY KEY,
    parent INTEGER,
    rank VARCHAR(100)
);

CREATE TABLE TaxonLineage
(
    taxon_id INTEGER,
    lineage_id INTEGER
);
CREATE INDEX tax_idx ON TaxonLineage(taxon_id);
CREATE INDEX lin_idx ON TaxonLineage(lineage_id);
