#!/usr/bin/env Rscript

# YAMS16 in R and refactored into functions
# By John McCulloch and Jonathan Badger
# Version 1.51 May 28, 2022

# check if packages installed, installs them if needed, and loads them
# by default normal R packages, bioconductor if bioconductor=TRUE
loadInstallPackages <- function(pkgs, bioconductor=FALSE) {
    needinstall <- pkgs[which(pkgs %in% rownames(installed.packages())==FALSE)]
    if (length(needinstall)) {
        print("Installing packages...\n")
        if (bioconductor) {
            if("BiocManager" %in% rownames(installed.packages()) == FALSE) {
                install.packages("BiocManager")
            }
            BiocManager::install(needinstall)
        } else {
            install.packages(needinstall, dependencies=TRUE)
        }
    }
    lapply(pkgs, function(x){
        suppressMessages(library(x, character.only = TRUE))})
}

# convert unit table and a tax table (from e.g. DADA2) to a series of jams files (one per sample)
unit2JAMS <- function(unit_tsv, tax_tsv, output_folder = getwd(), process = "DADA2", amplicon = "16S", verbose = TRUE) {
  unit <- read.csv(unit_tsv, sep="\t", row.names = 1)
  tax <- read.csv(tax_tsv, sep="\t", row.names = 1)
  tax$LKT <- paste0("LKT__",tax$LKT)
  if (!all(row.names(tax)==row.names(unit))) {
    flog.error("row names of unit table and tax table do not match!")
  }
  dir.create(output_folder, showWarnings = verbose)
  proj_info <- data.frame(names=c("Run_info","Sample_name","Run_type","Process"),
                          values=c("Run_value","XXX",amplicon, process))
  unit_ppm <- round(sweep(unit, 2, colSums(unit)/1e6, FUN="/"))
  for(sample in colnames(unit_ppm)) {
    if (verbose) {
      flog.info(paste("Processing sample", sample))
    }
    tax$NumBases <- unit_ppm[,sample]
    tax$PctFromCtg <- 100
    tax$ProbNumGenomes <- 1
    stax <- tax[tax$NumBases > 0,]
    stax <- stax[order(stax$NumBases, decreasing = TRUE),]
    proj_info[2,"values"] <- sample
    jamsdir <- file.path(output_folder, paste0(sample, "_JAMS"))
    jamsfile <- paste0(sample, ".jams")
    dir.create(jamsdir, showWarnings = verbose)
    infofile <- file.path(jamsdir, paste0(sample, "_projinfo.tsv"))
    lktdose <- file.path(jamsdir, paste0(sample, "_LKTdose.tsv"))
    write.table(proj_info, file = infofile, sep = "\t", quote = FALSE, 
                row.names = FALSE, col.names = FALSE)
    write.table(stax, file = lktdose, sep = "\t", quote = FALSE, row.names = FALSE)
    curdir <- getwd()
    setwd(jamsdir)
    system(str_c("tar cfz ", jamsfile, " ", basename(infofile), " ", basename(lktdose)))
    setwd(curdir)
    file.rename(file.path(jamsdir, jamsfile), file.path(output_folder, jamsfile)) 
    unlink(jamsdir, recursive = TRUE)
  }
}

# get path of running script
getScriptPath <- function() {
    cmd.args <- commandArgs()
    m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
    script.dir <- dirname(regmatches(cmd.args, m))
    if(length(script.dir) != 1) 
        stop("can't determine script dir: please call the script with Rscript")
    return(script.dir)
}

# filters fastqs using DADA2
filterFastq <- function(fastqdir, path, threads=2, maxLen = 240, 
                        minLen=180, trimLeft=c(19,20)) {
    fnFs <- sort(list.files(fastqdir, pattern="*_R1_*"))
    sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
    filt_path <- file.path(path, "filtered")
    if (!dir.exists(filt_path)) {
        dir.create(filt_path)
    }

  # filter one by one because sometimes filterAndTrim fails on one sample
  # and is impossible to tell which one :-(
  
    cl <- makeCluster(threads)
    clusterExport(cl, c("str_c","fastqPairedFilter"))
    parLapply(cl, sample.names, function(x) {
        r1 <- Sys.glob(str_c(fastqdir, "/", x,"_*R1_*"))
        r2 <- Sys.glob(str_c(fastqdir, "/", x,"_*R2_*"))
        fn <- c(r1, r2)
        fout <- c(file.path(filt_path, str_c(x, "_F_filt.fastq.gz")),
                  file.path(filt_path, str_c(x, "_R_filt.fastq.gz")))
        if (!file.exists(fout[1])) {
            print(str_c("Filtering ", x, "..."))
            fastqPairedFilter(fn, fout, truncLen=c(maxLen,minLen),
                              trimLeft=trimLeft, maxN=0, OMP=FALSE,
                              maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                              compress=TRUE, verbose=TRUE)
        }
    })
    stopCluster(cl)
}

# learn sequence errors by DADA2
errorLearn <- function(path) {
    filtFs <- Sys.glob(str_c(path,"/","filtered/*_F_*"))
    filtRs <- Sys.glob(str_c(path,"/","filtered/*_R_*"))
    if (!file.exists(file.path(path, "errF.rds"))) {
        flog.info("Learning errors for forward reads...")
        errF <- learnErrors(filtFs, multithread=TRUE)
        saveRDS(errF, file = file.path(path, "errF.rds"))
    }
    if (!file.exists(file.path(path, "errR.rds"))) {
        flog.info("Learning errors for reverse reads...")
        errR <- learnErrors(filtRs, multithread=TRUE)
        saveRDS(errR, file = file.path(path, "errR.rds"))
    }
}

# runs the DADA2 dereplication algorithm
dereplicate <- function(path) {
    if (!file.exists(file.path(path,"seqtab.rds"))) {
        errF <- readRDS(file.path(path, "errF.rds"))
        errR <- readRDS(file.path(path, "errR.rds"))
        filtFs <- sort(Sys.glob(str_c(path,"/","filtered/*_F_*")))
        filtRs <- sort(Sys.glob(str_c(path,"/","filtered/*_R_*")))
        if (!file.exists(file.path(path,"derepFs.rds"))) {
            flog.info("Dereplicating Forward...")
            derepFs <- derepFastq(filtFs, verbose=TRUE)
            saveRDS(derepFs, file.path(path,"derepFs.rds"))
        } else {
            derepFs <- readRDS(file.path(path,"derepFs.rds"))
        }
        if (!file.exists(file.path(path,"derepRs.rds"))) {
            flog.info("Dereplicating Reverse...")
            derepRs <- derepFastq(filtRs, verbose=TRUE)
            saveRDS(derepRs, file.path(path,"derepRs.rds"))
        } else {
            derepRs <- readRDS(file.path(path,"derepRs.rds"))
        }
        if (!file.exists(file.path(path,"dadaFs.rds"))) {
            flog.info("DADA2 Forward...")
            dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
            saveRDS(dadaFs, file.path(path,"dadaFs.rds"))
        } else {
            dadaFs <- readRDS(file.path(path,"dadaFs.rds"))
        }
        if (!file.exists(file.path(path,"dadaRs.rds"))) {
            flog.info("DADA2 Reverse...")
            dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
            saveRDS(dadaRs, file.path(path,"dadaRs.rds"))
        } else {
            dadaRs <- readRDS(file.path(path,"dadaRs.rds"))
        }
        if (!file.exists(file.path(path,"seqtab.rds"))) {
            flog.info("Merging Pairs...")
            sample.names <- unique(rapply(str_split(basename(filtFs),
                                                    pattern="_"), 
                                      function(x){x[1]}))
            mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                                  verbose=TRUE, justConcatenate = FALSE)
            seqtab <- makeSequenceTable(mergers)
            row.names(seqtab) <- sample.names
            saveRDS(seqtab, file.path(path,"seqtab.rds"))
      }
    }
}

# removes chimeras from seqtab by DADA2
removeChimeras <- function(path, method = "consensus") {
    if(!file.exists(file.path(path, "seqtab.nochim.rds"))) {
        flog.info("Removing Chimeras...")
        seqtab <- readRDS(file.path(path, "seqtab.rds"))
        seqtab.nochim <- removeBimeraDenovo(seqtab, method=method, 
                                            multithread=TRUE, 
                                            verbose=TRUE)
        saveRDS(seqtab.nochim, file.path(path, "seqtab.nochim.rds"))
    }
}

# returns the most classified taxon of each taxtable row
infer_LKT <- function(taxtable) {
    apply(taxtable, 1, function(x) {
        undefined = which(x=="Unclassified")
        if (length(undefined)==0) {
            str_c("s__",str_replace(x[6],pattern="g__",""),
                  "_",str_replace(x[7],pattern="s__",""))
        } else {
            str_c(x[min(undefined)-1],"_",x[min(undefined)])
        }
    })
}

# assigns taxonomy using taxonomic database
taxonomy <- function(path, db="silva", yamsdir, species_boot=FALSE) {
    if (!file.exists(file.path(path,"taxa_16S_cons.tsv"))) {
        flog.info("Inferring taxonomy...")
        seqtab.nochim <- readRDS(file.path(path, "seqtab.nochim.rds"))
        taxdb <- Sys.glob(str_c(yamsdir,"/db/", db, "*train_set*"))[1]
        speciesdb <- Sys.glob(str_c(yamsdir,"/db/", db, "*species_assign*"))[1]
        speciesbootdb <- Sys.glob(str_c(yamsdir,"/db/", db,
                                        "*species_comp*"))[1]
        taxa <- as.data.frame(assignTaxonomy(seqtab.nochim, taxdb, 
                                             verbose=TRUE,multithread=TRUE))
        if (species_boot) {
            taxa_Gs <- as.data.frame(
                assignTaxonomy(seqtab.nochim, speciesbootdb,
                               verbose=TRUE,multithread=TRUE,
                               outputBootstraps = TRUE))
            return(taxa_Gs)
        } else {
            taxa_Gs <- as.data.frame(assignSpecies(seqtab.nochim, speciesdb,
                                                   verbose=TRUE))
        }
        taxa$Species <- taxa_Gs$Species
        colnames(taxa)[1] <- "Domain"
        cnames <- colnames(taxa)
        ranks <- c("d__","p__","c__","o__","f__","g__","s__")
        seqs <- row.names(taxa)
        taxa <- as.data.frame(sapply(1:length(ranks),
                                     function(i) taxa[,i] <-
                                                     str_c(ranks[i],taxa[,i])))
        colnames(taxa) <- cnames
        taxa[is.na(taxa)] <- "Unclassified"
        taxa$LKT <- vapply(infer_LKT(taxa), paste, collapse = "|",
                           character(1L))
        lkt_names <- make.names(taxa$LKT,unique = TRUE)
        write.fasta(as.list(seqs), lkt_names, as.string = TRUE, 
                    file.path(path, "rep_set.fa"))
        row.names(taxa) <- lkt_names
        write.table(taxa, file.path(path, "taxa_16S_cons.tsv"), sep="\t", 
                    col.names = NA, quote = FALSE)
        unit_table <- as.data.frame(t(seqtab.nochim))
        row.names(unit_table) <- row.names(taxa)
        write.table(unit_table, file.path(path, "unit_table.tsv"), sep="\t",
                    col.names = NA, quote = FALSE)
        glomByLKT(file.path(path, "unit_table.tsv"),
                  file.path(path, "taxa_16S_cons.tsv"),
                  file.path(path, "rep_set.fa"))
    }
}

glomByLKT <- function(unittable, taxtable, repset=NULL, path=".") {
  unit <- read.csv(unittable, sep="\t", row.names = 1,
                    check.names = FALSE)
                    
  tax <- read.csv(taxtable, sep="\t", row.names = 1, 
                  check.names = FALSE)
  
  unit$LKT <- tax[row.names(unit),]$LKT
  unit_LKT<-aggregate(. ~ LKT, unit, sum)
  row.names(unit_LKT) <- unit_LKT$LKT
  unit_LKT$LKT <- NULL
  unit_LKT <- unit_LKT[order(rowSums(unit_LKT), decreasing = TRUE),]
  unit_LKT <- cbind(Taxon=row.names(unit_LKT), unit_LKT)
  write.table(unit_LKT, file.path(path, "unit_table_LKT.tsv"), row.names=FALSE,
                                        sep="\t", quote=FALSE)
  
  tax <- tax %>% distinct(LKT, .keep_all = TRUE) %>% data.frame(row.names = .$LKT)
  tax <- tax[row.names(unit_LKT),]
  tax <- cbind(Taxon=row.names(unit_LKT), tax)
  write.table(tax, file.path(path, "taxa_16S_cons_LKT.tsv"),
              row.names=FALSE, sep="\t", quote=FALSE)
  if (!is.null(repset)) {
    names <- row.names(tax) %>% gsub("/",".", .) %>%
      gsub("-",".", .) %>% gsub(" ",".", ., fixed=TRUE) %>% gsub("[","", ., fixed=TRUE) %>% gsub("]", "",., fixed=TRUE) %>% gsub("^_Unclassified","X_Unclassified",.)
    seqs <- read.fasta(repset)
    names(seqs) <- names(seqs) %>% gsub("..",".",., fixed=TRUE) %>% gsub("__.","__",., fixed=TRUE)
    seqs <- seqs[names]
    write.fasta(lapply(seqs, str_to_upper), names=names(seqs),
                file.out = file.path(path, "rep_set_LKT.fa"))
  }
}

# main function for parsing arguments and running pipeline
main <- function(args) {
    yamsdir <- getScriptPath()
    adapter <- file.path(yamsdir, "db/adapters/adapter.fasta")
    version <- "1.05 08 Feb 2019"
    if (length(args)==0) {
        args <- c("-h")
    }
    option_list <- list(
        make_option(c("-a", "--adapter"), action="store", default=adapter,
                    type="character",
                    help=str_c("file of Illumina adapters (default: ", 
                               adapter, ")")),
        make_option(c("-i", "--input"), action="store", default=NULL, 
                    type="character",
                    help="path of directory where fastq files are"),
        make_option(c("-p", "--prefix"), action="store", default=NULL, 
                    type="character",
                    help ="filename prefix name for project"),
        make_option(c("-o", "--output"), action="store", default="./", 
                    type="character",
                    help ="path where to output files (default .)"),
        make_option(c("-t", "--threads"), action="store", default=2, 
                    type="numeric",
                    help ="number of threads to use (default 1)"),
        make_option(c("-x", "--maxlength"), action="store",
                    type="numeric", default=240,
                    help="maximum amplicon length (default 240)"),
        make_option(c("-n", "--minlength"), action="store",
                    type="numeric", default=180,
                    help="minimum amplicon length (default 180)")
    )
    opt <- parse_args(OptionParser(option_list=option_list), args)
    
    if (is.null(opt$input)) {
        stop("option -i/--input is required\n")
    }
    if (is.null(opt$prefix)) {
        stop("option -p/--prefix is required\n")
    }
    if (!dir.exists(opt$output)) {
        dir.create(opt$output)
    }
    outprefix <- str_c(opt$output,"/",opt$prefix)
    if (!dir.exists(outprefix)) {
        dir.create(outprefix)
    }
    flog.appender(appender.tee(file.path(outprefix,"YAMS16.log")))
    msg <- loadInstallPackages(c("dada2","phyloseq","metagenomeSeq"), 
                               bioconductor = TRUE)
    filterFastq(opt$input, outprefix, opt$threads, opt$maxlength, opt$minlength)
    errorLearn(outprefix)
    dereplicate(outprefix)
    removeChimeras(outprefix)
    taxonomy(outprefix, yamsdir = yamsdir)
    flog.info("Done")
}

msg <- loadInstallPackages(c("tidyverse", "optparse", "futile.logger",
                             "seqinr", "parallel"))
options(stringsAsFactors = FALSE)
# run main function only if run as script
if (!interactive()) {
    main(commandArgs(trailingOnly = TRUE))
}

