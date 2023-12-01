#!/path/to/R-4.0.1


## These jobs were run as arrays using SGE_TASK_ID
## The run_coaccessibility.sh script was used to run this script (coaccessibility.R)
## Script calculates pairwise co-accessibility of an ATAC-seq peak at the task index and all ATAC-seq peaks following it until the end of the matrix.

## Load Functions and Packages 
source("functions.R")
source("packages.R")
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(coxme))
suppressPackageStartupMessages(library(kinship2))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(preprocessCore))

args = commandArgs(trailingOnly=TRUE)
task = as.numeric(args[1])

tmm_matrix="tmm_matrix_path"
kinship_matrix="kinship_matrix_path"
out_dir="outdir_path"

### Process Covariates

subject_metadata = read_excel("Table_S1.xlsx")
subject_metadata$sex = ifelse(subject_metadata$Sex == "Male",1,0)
subject_metadata$age = subject_metadata$`Age at enrollment`/mean(subject_metadata$`Age at enrollment`)
subject_covariates = subject_metadata[, c("Subject ID","WGS UUID","sex","age", paste0("PC",1:20))]

atac_metadata = read_excel("Table_S6.xlsx")
## Premerge Sample IDs are redundant and used to map the fastq files on GEO
atac_metadata$`Premerge Sample ID` = NULL
atac_metadata = unique(atac_metadata)
atac_covariates = atac_metadata[, c("Subject ID","Merged Library UUID","iPSC Passage","Reads Passing Filters",
                                    "Mean Fragment Size","Number of Peaks")]

normalize_atac_covariates = function(covariate_label) {
    cov_df = atac_covariates[,c("Subject ID","Merged Library UUID",covariate_label)]
    names(cov_df)[3] = "cov"
    cov_df$normed = cov_df$cov/mean(cov_df$cov)
    names(cov_df)[4] = gsub(" ","_",tolower(covariate_label))
    return(cov_df[,c(1,2,4)])
}

normalized_atac = Reduce(function(...) merge(..., all=T, by=c("Subject ID","Merged Library UUID")), 
    lapply(c("iPSC Passage","Reads Passing Filters","Mean Fragment Size","Number of Peaks"), normalize_atac_covariates))

covariates = merge(subject_covariates,normalized_atac)
colnames(covariates) = gsub(" ","_",tolower(colnames(covariates)))

### Load TMM Matrix 
peak_tmm = add_rownames(fread(tmm_matrix,sep="\t",data.table=F))
inv = transform_standard_normal(peak_tmm)

### Load Kinship MAtrix
kinship_matrix = add_rownames(fread(kinship_matrix_path,
                        data.table=F,sep="\t"))
sparse_kinship = Matrix(as.matrix(kinship_matrix), sparse = TRUE)

### Calculate Co-accessibility
test_mat = inv[ task:nrow(inv),]

calculate_coaccessibility = function(idx) {
    peak1_id = rownames(test_mat[ 1,]) 
    peak2_id = rownames(test_mat[ idx,]) 
    pairtotest = as.data.frame(t(test_mat[c(peak1_id,peak2_id),]))
    pairtotest$merged_library_uuid = rownames(pairtotest)
    pairtotest_covs = merge(pairtotest, covariates)
    
    my_formula = as.formula(paste(peak1_id, 
                                  paste(c(peak2_id,colnames(pairtotest_covs)[6:ncol(pairtotest_covs)], "(1|wgs_uuid)"),
                                        collapse = "+"), sep = "~"))
    model = lmekin(formula = my_formula, data = pairtotest_covs, varlist = sparse_kinship)
    out = extract_coxme_table(model)[2,]
    out$peak1_id = peak1_id
    out$peak2_id = peak2_id
    
    return(out)
}
lmm_coaccessibility = rbindlist(lapply(2:nrow(test_mat),calculate_coaccessibility))

fwrite(lmm_coaccessibility,paste0(out_dir,"/",unique(lmm_coaccessibility$peak1_id),".txt"),
      sep="\t",row.names=F,quote=F)
