CREATE TABLE "TCGA_ChiTaRS_combined_fusion_information_on_hg19" (
	"source1" varchar(50) NOT NULL DEFAULT '',
	"source2" varchar(50) NOT NULL DEFAULT '',
	"ctype" varchar(50) NOT NULL DEFAULT '',
	"sample" varchar(50) NOT NULL DEFAULT '',
	"h_gene" varchar(50) NOT NULL DEFAULT '',
	"h_chr" varchar(10) NOT NULL DEFAULT '',
	"h_bp" integer NOT NULL DEFAULT -1,
	"h_strand" char(1) NOT NULL DEFAULT '',
	"t_gene" varchar(50) NOT NULL DEFAULT '',
	"t_chr" varchar(10) NOT NULL DEFAULT '',
	"t_bp" integer NOT NULL DEFAULT -1,
	"t_strand" char(1) NOT NULL DEFAULT ''
);
CREATE INDEX h_gene_index ON TCGA_ChiTaRS_combined_fusion_information_on_hg19(h_gene);
CREATE INDEX t_gene_index ON TCGA_ChiTaRS_combined_fusion_information_on_hg19(t_gene);

CREATE TABLE "TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19" (
	"orf" varchar(50) NOT NULL DEFAULT '',
	"h_enst" varchar(50) NOT NULL DEFAULT '',
	"t_enst" varchar(50) NOT NULL DEFAULT '',
	"source" varchar(50) NOT NULL DEFAULT '',
    "ctype1" varchar(50) NOT NULL DEFAULT '',
    "ctype2" varchar(50) NOT NULL DEFAULT '',
    "ctype3" varchar(50) NOT NULL DEFAULT '',
    "sample" varchar(50) NOT NULL DEFAULT '',
	"h_gene" varchar(50) NOT NULL DEFAULT '',
	"h_chr" varchar(10) NOT NULL DEFAULT '',
	"h_bp" integer NOT NULL DEFAULT -1,
	"h_strand" char(1) NOT NULL DEFAULT '',
    "t_gene" varchar(50) NOT NULL DEFAULT '',
	"t_chr" varchar(10) NOT NULL DEFAULT '',
    "t_bp" integer NOT NULL DEFAULT -1,
	"t_strand" char(1) NOT NULL DEFAULT ''
);
CREATE INDEX h_gene_orf_index ON TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19(h_gene);
CREATE INDEX t_gene_orf_index ON TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19(t_gene);

CREATE TABLE "uniprot_gsymbol" (
	"uniprot_acc" varchar(50) NOT NULL,
	"gene_symbol" varchar(50) NOT NULL DEFAULT ''
);
CREATE INDEX uniprot_acc_symbol_index ON uniprot_gsymbol(uniprot_acc);

CREATE TABLE "fusion_uniprot_related_drugs" (
	"drug_status" varchar(255) NOT NULL DEFAULT '',
	"drug_bank_id" varchar(50) NOT NULL DEFAULT '',
    "drug_name" varchar(50) NOT NULL DEFAULT '',
    "drug_type" varchar(50) NOT NULL DEFAULT '',
    "uniprot_acc" varchar(50) NOT NULL,
    "drug_action" varchar(255) NOT NULL DEFAULT ''
);
CREATE INDEX uniprot_acc_related_index ON fusion_uniprot_related_drugs(uniprot_acc);

CREATE TABLE "fusion_ppi" (
	"h_gene" varchar(50) NOT NULL DEFAULT '',
	"h_gene_interactions" TEXT NOT NULL,
	"t_gene" varchar(50) NOT NULL DEFAULT '',
	"t_gene_interactions" TEXT NOT NULL
);
CREATE INDEX h_gene_ppi_index ON fusion_ppi(h_gene);
CREATE INDEX t_gene_ppi_index ON fusion_ppi(t_gene);

CREATE TABLE "fgene_disease_associations" (
	"code" integer,
	"gene" varchar(50) NOT NULL DEFAULT '',
	"disease_id" varchar(50) NOT NULL,
	"disease_desc" varchar(50) NOT NULL DEFAULT '',
	"disease_prob" float(11),
	"disease_pub" integer,
	"code2" integer,
	"disease_source" varchar(50) NOT NULL DEFAULT ''
);
CREATE INDEX gene_index ON fgene_disease_associations(gene);

.separator "\t"
.import TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19.txt TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19
.import TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt TCGA_ChiTaRS_combined_fusion_information_on_hg19
.import fusion_uniprot_related_drugs.txt fusion_uniprot_related_drugs
.import fusion_ppi.txt fusion_ppi
.import fgene_disease_associations.txt fgene_disease_associations
.import uniprot_gsymbol.txt uniprot_gsymbol
