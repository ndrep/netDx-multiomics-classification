# Load libraries
library(curatedTCGAData);
library(MultiAssayExperiment);
library(TCGAutils);
library(readxl);

# Source code
source("~/tcga_replicateFilter.R");

# Check all available assays
BRCA_assays <- curatedTCGAData(diseaseCode = "BRCA", 
                               assays = "*", version = "2.0.1");

# Informations I can get
# 1 - the release is GRCh37 - hg19
# 2 - for the BRCA, I obtain 17 objects (but they seem old - date in the name
#     is 2016-01-28)
# 3 - See ?curatedTCGA for a minimal explanation of the data
# 4 - From descriptions in point 3, we could use:
#   RNASeqGene (Gene expression values not normalized),
#   RNASeq2Gene (RSEM TPM gene expression values)
#   miRNASeqGene (Gene-level log2 RPM miRNA expression values; raw
#   reads seem not available), 
#   GISTIC_ThresholdedByGene (Gene-level GISTIC2 thresholded discrete copy
#   number values), 
#   GISTIC_Peaks (GISTIC2 thresholded discrete copy number values
#   in recurrent peak regions)
#   RPPAArray (Reverse Phase Protein Array normalized protein
#   expression values)
#   Methylation_methyl450 (Probe-level methylation beta values from Infinium
#   HumanMethylation 450K BeadChip)
#   clinical

# Data download
cnv_thresh <- curatedTCGAData(diseaseCode = "BRCA", 
                              assays = "GISTIC_ThresholdedByGene", 
                              version = "2.0.1", dry.run = FALSE);

#cnv_peaks <- curatedTCGAData(diseaseCode = "BRCA", 
#                             assays = "GISTIC_Peaks", 
#                             version = "2.0.1", dry.run = FALSE);

#methy <- curatedTCGAData(diseaseCode = "BRCA", 
#                         assays = "Methylation_methyl450", 
#                         version = "2.0.1", dry.run = FALSE);

miRNA <- curatedTCGAData(diseaseCode = "BRCA", 
                         assays = "miRNASeqGene", 
                         version = "2.0.1", dry.run = FALSE);

mRNA <- curatedTCGAData(diseaseCode = "BRCA", 
                        assays = "RNASeq2Gene", 
                        version = "2.0.1", dry.run = FALSE);


#clinical <- colData(mRNA)[, getClinicalNames("BRCA")];

protein <- curatedTCGAData(diseaseCode = "BRCA", 
                           assays = "RPPAArray", 
                           version = "2.0.1", dry.run = FALSE);

# Extract data
cnv <- t(assay(cnv_thresh));
colnames(cnv) <- cnv_thresh@ExperimentList@listData$`BRCA_GISTIC_ThresholdedByGene-20160128`@elementMetadata[, 'Gene.Symbol'];
class(cnv) <- "numeric";

mirna <- t(assay(miRNA));
mrna <- t(assay(mRNA));
prot <- t(assay(protein))
#clin <- as.data.frame(clinical);
#met <- t(assay(methy));

# Select only primary solid tumors
# Not necessary for clinical
cnv <- cnv[TCGAsampleSelect(rownames(cnv), c("01")), ];
mirna <- mirna[TCGAsampleSelect(rownames(mirna), c("01")), ];
mrna <- mrna[TCGAsampleSelect(rownames(mrna), c("01")), ];
prot <- prot[TCGAsampleSelect(rownames(prot), c("01")), ];
#met <- met[TCGAsampleSelect(rownames(met), c("01")), ];

# Remove technical replicates (present only in genomic data)
# NOTE: in all cases no technical replicates are found
keep <- tcga_replicateFilter(rownames(cnv), analyte_target=c("DNA"));
keep <- tcga_replicateFilter(rownames(mrna), analyte_target=c("RNA"));
keep <- tcga_replicateFilter(rownames(mirna), analyte_target=c("RNA"));
#keep <- tcga_replicateFilter(rownames(met), analyte_target=c("DNA"));

# Check if barcodes are unique
length(unique(substr(rownames(cnv), 1, 12))) == nrow(cnv)
length(unique(substr(rownames(mirna), 1, 12))) == nrow(mirna)
length(unique(substr(rownames(mrna), 1, 12))) == nrow(mrna)
length(unique(substr(rownames(prot), 1, 12))) == nrow(prot)
#length(unique(substr(rownames(met), 1, 12))) == nrow(met)

# Clean rownames to just sample barcode (strip further details)
# Not necessary for clinical
rownames(cnv) <- substr(rownames(cnv), 1, 12);
rownames(mirna) <- substr(rownames(mirna), 1, 12);
rownames(mrna) <- substr(rownames(mrna), 1, 12);
rownames(prot) <- substr(rownames(prot), 1, 12);
#rownames(met) <- substr(rownames(met), 1, 12);

# Subset to common samples for all data types (except protein)
# in the same order
# NOTE: considering also protein we get 510 samples instead of 616
# Another option is removing methylation and keep protein, but for now
# I'll stick for what we have done in hg38
common.samples <- Reduce(intersect, list(rownames(cnv), 
                                         rownames(mirna), 
                                         rownames(mrna),
                                         rownames(prot)));

cnv <- cnv[common.samples, ];
mirna <- mirna[common.samples, ];
mrna <- mrna[common.samples, ];
prot <- prot[common.samples, ];
#clin <- clin[common.samples, ];
#met <- met[common.samples, ];

# Remove features with |NAs| > 0.20
# NOTE: no need of imputation, after filtering no more NAs are left
prot_filt <- prot[, colSums(is.na(prot))/nrow(prot) < 0.20];

# Add PFI labels from TCGA-CDR
tcga_cdr <- suppressWarnings(read_excel("~/tcga-cdr.xlsx", sheet = "TCGA-CDR"));
labels <- as.data.frame(tcga_cdr[tcga_cdr$type == "BRCA", c("bcr_patient_barcode", "PFI", "DFI")]);
labs <- labels[, -1];
rownames(labs) <- labels[, 1];

labs <- labs[rownames(labs) %in% common.samples, ];
labs_pfi <- as.vector(labs[,'PFI']);
names(labs_pfi) <- rownames(labs);

# Save results as rds and csv
saveRDS(cnv, file = "./dataset/cnv.rds");
saveRDS(mirna, file = "./dataset/mirna.rds");
saveRDS(mrna, file = "./dataset/mrna.rds");
saveRDS(prot_filt, file = "./dataset/proteins.rds");
saveRDS(labs_pfi, file = "./dataset/labels_pfi.rds");


