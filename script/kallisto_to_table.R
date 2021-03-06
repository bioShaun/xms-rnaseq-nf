suppressMessages(library(argparser))
suppressMessages(library(tximport))
suppressMessages(library(tibble))
suppressMessages(library(omplotr))
suppressMessages(library(edgeR))
options(stringsAsFactors = F)
options(bitmapType = "cairo")

p <- arg_parser("read kallisto quant files generate expression matrix")
p <- add_argument(p, '--kallisto_dir', help = 'kallisto quantification directory')
p <- add_argument(p, '--sample_inf', 
                  help = 'sample information with sample names and group names',
                  default=NULL)
p <- add_argument(p, '--gene2tr', 
                  help = 'gene id and transcript id mapping file',
                  default=NULL)
p <- add_argument(p, '--out_dir',    help = 'diff analyssi output directory')
argv <- parse_args(p)


## read parameters
sample_inf <- argv$sample_inf
kallisto_dir <- argv$kallisto_dir
gene2tr_file <- argv$gene2tr
expression_stat_dir <- argv$out_dir

## directory prepare
dir.create(expression_stat_dir, showWarnings = FALSE)

if (is.na(sample_inf)) {
    sample_ids <- list.files(kallisto_dir)
    samples <- data.frame(condition=sample_ids, sample=sample_ids)
} else {
    samples <- read.delim(sample_inf, stringsAsFactors = F, header = F)
    colnames(samples) <- c('condition', 'sample')
}

files <- file.path(kallisto_dir, samples$sample, "abundance.h5")

names(files) <- samples$sample
print(files)
print(file.exists(files))
if (! is.na(gene2tr_file)) {
    gene2tr <- read.delim(gene2tr_file, header = FALSE)
    colnames(gene2tr) <- c('gene_id', 'transcript_id')
    tx2gene <- gene2tr[,c('transcript_id', 'gene_id')]


    ## normalized expression matrix
    txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
    cts <- txi$counts
    y <- DGEList(cts)
    normfactors <- calcNormFactors(y)
    gene_tpm_matrix <- (txi$abundance)/(normfactors$samples$norm.factors)

    ## output quant table
    out_gene_tpm_matrix <- as.data.frame(gene_tpm_matrix)
    out_gene_tpm_matrix <- round(out_gene_tpm_matrix, 3)
    out_gene_tpm_matrix <- rownames_to_column(out_gene_tpm_matrix, var="Gene_ID")
    out_cts <- as.data.frame(cts)
    out_cts <- round(out_cts, 3)
    out_cts <- rownames_to_column(out_cts, var = 'Gene_ID')
    write.table(out_cts, file = paste(expression_stat_dir, 'Gene.count.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
    write.table(out_gene_tpm_matrix, file = paste(expression_stat_dir, 'Gene.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
}

### transcript level expression matrix
txi.tx <- tximport(files, type = "kallisto", txOut = TRUE)
cts.tx <- txi.tx$counts
y.tx <- DGEList(cts.tx)
normfactors.tx <- calcNormFactors(y.tx)
tx_tpm_matrix <- (txi.tx$abundance)/(normfactors.tx$samples$norm.factors)

### output transcript level quant table
out_tx_cts <- as.data.frame(cts.tx)
out_tx_cts <- round(out_tx_cts, 3)
out_tx_cts <- rownames_to_column(out_tx_cts, var = 'Transcript_ID')
out_tx_tpm_matrix <- as.data.frame(tx_tpm_matrix)
out_tx_tpm_matrix <- round(out_tx_tpm_matrix, 3)
out_tx_tpm_matrix <- rownames_to_column(out_tx_tpm_matrix, var='Transcript_ID')
write.table(out_tx_cts, file = paste(expression_stat_dir, 'Transcript.count.txt', sep = '/'), quote=F, row.names = F, sep = '\t')
write.table(out_tx_tpm_matrix, file = paste(expression_stat_dir, 'Transcript.tpm.txt', sep = '/'), quote = F, row.names = F, sep = '\t')
