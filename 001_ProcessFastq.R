library(amplicR)
library(readxl)
library(dada2)
library(data.table)
library(ggplot2)
library(ShortRead)

if(dir.exists("C:/Data/Software/Dropbox/Deer/Pop genetics/DELWP led projects/Data"))
dir.in <- "C:/Data/Software/Dropbox/Deer/Pop genetics/DELWP led projects/Data" else
  dir.in <- "E:/Dropbox/Deer/Pop genetics/DELWP led projects/Data"

info_file_R1 <- read_excel(file.path(dir.in, "Deer_AmpliconAndPrimer_Info.xlsx"), 
                           sheet = "Sheet1")

# keep only the ones used
info_file_R1 <- info_file_R1[!is.na(info_file_R1$`In final dataset?`),] 
info_file_R2 <- info_file_R1[,c(1:3, 6)]
names(info_file_R2)[2:3] <- c("Reverse", "Forward") # Reverse the F and R primers for the reverse reads
info_file <- list(info_file_R1, info_file_R2)

orig_fastq_dir <- "Amplicon_Fastqs_Paired_Sep"

# pulling out R1 and R2 names
lfall <- list.files(file.path(dir.in, orig_fastq_dir), 
                    pattern = ".fastq.{,3}$", 
                 full.names = T)


readsDir <- c("R1", "R2")
# Deal with R1 and R2 one at the time
for(Rdir in readsDir) {
  lf <- lfall[grep(Rdir, lfall)]
  dir.create(file.path(dir.in, orig_fastq_dir, Rdir))
  ld <- vector("list", length = length(lf))
  # This loop set the setting file for "deconv" and runs it on each file
  # deconv reads fastq file from 'orig_fastq_dir', separates reads within each 
  # file by contig and write a new fastq files within the relevant folder
  # ("R1" and "R2" respectively and then contig name)
  # reads where the the fwd or rev primers are not found (with 0 missmatches) are dropped
  for(fn in lf) {
    SampleID <- sub(pattern = "_L001_R[1-2]_001.fastq.gz", replacement = "", basename(fn))
    write.csv(cbind(info_file[[which(readsDir == Rdir)]], 
                    Sample_IDs=SampleID,
                    Find="",
                    Rind=""),
                file=file.path(dir.in, orig_fastq_dir, Rdir, 
                               paste0("info_file",SampleID, ".csv")))
    
    ld[[which(fn == lf)]] <- deconv(fn = fn, dir.out = file.path(dir.in, orig_fastq_dir, Rdir),
                                    info.file = file.path(dir.in, orig_fastq_dir, 
                                                Rdir, paste0("info_file",SampleID, ".csv")), 
           Fprimer = "Forward", Rprimer = "Reverse", gene = "Contig")
  }
  # Summarise number of seq processed and dropped
  nProc <- sapply(ld, "[[", 1)
  nFwdFound <- sapply(ld, "[[", 2)
  nRevFound <- sapply(ld, "[[", 3)
  
  s_deconv <- data.frame(Sample_ID=sub(pattern = "_L001_R[1-2]_001.fastq.gz", 
                                       replacement = "", basename(lf)),
                         nProc=nProc, nFwdFound, nRevFound)
  write.csv(s_deconv, file = file.path(dir.in, orig_fastq_dir, Rdir, 
                                       paste0("summary_deconv", Rdir, ".csv")), 
                                       row.names = F)
  }


contigsName <- as.character(info_file_R1$Contig)
path.results <- vector("list", length = length(readsDir))
lf2 <- vector("list", length = length(readsDir))
lrf <- vector("list", length = length(readsDir))
w_lrf <- vector("list", length = length(readsDir))
chk_widths <- vector("list", length = length(readsDir))


for(R in seq_along(readsDir)) {
  # Cross check. read the fastq file to confirm there are reads in there and check their length
path.results[[R]] <- file.path(dir.in, orig_fastq_dir, readsDir[R], info_file[[R]]$Contig)
lf2[[R]] <- sapply(path.results[[R]], function(x) length(list.files(x, pattern = "fastq\\.gz$")))
if(sum(lf2[[R]] == 0)>0) stop()

lrf[[R]] <- lapply(path.results[[R]], readFastq)
w_lrf[[R]] <- lapply(lrf[[R]], width)
chk_widths[[R]] <- lapply(w_lrf[[R]], table)
chk_widths[[R]] <- sapply(chk_widths[[R]], function(x) as.numeric((names(x)[which.max(x)])))
#sum(contigSize == chk_widths) # using the mode because of the large number of discrepancies
}

# check modes are the same 
all.equal(chk_widths[[1]], chk_widths[[2]])
# write chk_widths so don't have to re-run in a new session
save(chk_widths, file=file.path(dir.in, orig_fastq_dir, "chk_width.rds"))
load(file.path(dir.in, orig_fastq_dir, "chk_width.rds"))

# process each contig with data.proc
# Outputs are saved in 'file.path(dir.in, orig_fastq_dir, "Processed",contigsName[g])'
path.results.m <-matrix(unlist(path.results), nrow = 2, byrow = T)
dir.create(file.path(dir.in, orig_fastq_dir, "Processed"))

# Fix name convention 
r1Path <- dirname(path.results[[1]][1])
r2Path <- dirname(path.results[[2]][1])

r1List.files <- list.files(r1Path, pattern = ".fastq.{,3}$", full.names = TRUE, 
                          recursive = TRUE)
r2List.files <- list.files(r2Path, pattern = ".fastq.{,3}$", full.names = TRUE, 
                          recursive = TRUE)
list.files.n <- c(r1List.files, r2List.files)
new.list.files <- gsub(pattern = "_R[1-2].fastq.gz_Ind_", list.files.n, replacement = "_Ind_")
list.files.n <- list.files.n[grep(pattern = "DA2", list.files.n)]
new.list.files <- new.list.files[grep(pattern = "DA2", new.list.files)]

file.rename(list.files.n, new.list.files)

# Changed to retain all length
# allow min overlap 95% (round down) of mode length of fragment
# max 15% mismatch
ldproc <- list()
bp <- 0
for(g in seq_along(contigsName)) {
  minOv <- floor(chk_widths[[1]][g] * 0.95) # USe one since they are the same
  
  fns <- list.files(path=path.results[[1]][g])
  fastqs1 <- fns[grepl(".fastq.{,3}$", fns)]
  fns <- list.files(path=path.results[[2]][g])
  fastqs2 <- fns[grepl(".fastq.{,3}$", fns)]
  if((length(fastqs1) == 0) | (length(fastqs2) == 0)) {
    message(paste("There are no files in", path.results[g],
                  "with either fastq or fastq.gz extension", sep = "\n"))
    next
  } else {
    txt <- capture.output(
      ldproc[[g]] <- data.proc.paired(dir.in=path.results.m[, g], bp=bp, truncQ=2, 
                  dir.out = file.path(dir.in, orig_fastq_dir, "Processed",contigsName[g]),
                                      minOverlap = minOv,
                                      maxMismatch = floor(minOv * 0.15),
                                      trimOverhang = FALSE,
                                verbose=FALSE)
    )
  }
}


save(ldproc, file = file.path(dir.in, orig_fastq_dir, "ldproc.rda"))
# ldproc is a list of lists, where each elements is the output from data.proc, see ?data.proc
# the most important outputs are also saved as .csv within each contig's folder


#### Results ####

# keep only the ones used
info_file_R1 <- info_file_R1[!is.na(info_file_R1$`In final dataset?`),] 

info_file_R2 <- info_file_R1[,c(1:3, 6)]
names(info_file_R2)[2:3] <- c("Reverse", "Forward") # Reverse the F and R primers for the reverse reads
info_file <- list(info_file_R1, info_file_R2)
contigsName <- as.character(info_file_R1$Contig)

load(file = file.path(dir.in, orig_fastq_dir, "ldproc.rda"))
names(ldproc) <- contigsName

##### whole amplicon data #####
dir.out <- file.path(dir.in, orig_fastq_dir, "WholeAmpliconRes")
dir.create(dir.out)
# Extract seq tables
seq_tables <- lapply(ldproc, "[[", 3)
#names(seq_tables) <- contigsName
# get a table wtih n of alleles for each sampels, and each contig
nAlleles <- lapply(seq_tables, function(x) {
  nAlleles <- apply(x, 1, function(x) sum(x>0))
  Sample_ID <- names(nAlleles)
  #Sample_ID <- data.table::tstrsplit(names(nAlleles), "_")[[2]]
  data.frame(Sample_ID, nAlleles)
}
  )

lSampleID <- lapply(nAlleles, function(x) {
  return(x$Sample_ID)
})
snames <- unique(unlist(lSampleID))
nAlleles <- c(list(data.frame(Sample_ID=snames)), nAlleles)
nAlleledf <- plyr::join_all(nAlleles, by="Sample_ID", type="left")

names(nAlleledf)<-c("Sample_ID", names(seq_tables))

# melt df to plot with ggplot
nAlleledt <- melt.data.table(data.table(nAlleledf), id.vars = "Sample_ID", variable.name = "Contig", value.name = "nAlleles")
write.csv(nAlleledt, file = file.path(dir.out, "nAlleledt.csv"), row.names = F)

#heat map of nAlleles: contig X sample
ggplot(nAlleledt, aes(Contig, Sample_ID)) +
  geom_tile(aes(fill=nAlleles)) +
  theme(axis.text.x = element_text(angle = 90))

ggsave(file.path(dir.out, "heatmapnAlleles.png"),
       width = 24, height = 19, units = "cm")

# Distributions of nAlleles for samples with >2 alleles
ggplot(nAlleledt[nAlleles>2,], aes(nAlleles)) +
  geom_histogram(bins=nAlleledt[!is.na(nAlleles),max(nAlleles)]) +
  facet_grid(Contig~.)

ggsave(file.path(dir.out, "barPlotnAlleles.png"),
       width = 19, height = 24, units = "cm")

# Distributions of nAlleles all samples
ggplot(nAlleledt, aes(nAlleles)) +
  geom_histogram(bins=nAlleledt[!is.na(nAlleles),max(nAlleles)]) +
  facet_grid(Contig~.)

ggsave(file.path(dir.out, "barPlotnAllelesAll.png"),
       width = 19, height = 35, units = "cm")

#-----------------------------------------------------------------------------#
#   Manually check some details checked 
#-----------------------------------------------------------------------------#

abund587 <- ldproc[["587"]][[3]]
abundInd <- integer()
for(i in seq_len(nrow(abund587))) {
  abundInd[i] <- which.max(abund587[i,])
}
abundInd # which one is the most abundant seq in each sample

seq587 <- Biostrings::DNAStringSet(ldproc[["587"]][[4]][, "sequence"])
aln587 <- DECIPHER::AlignSeqs(seq587)
aln587

# Compute point distance compared with the first seq
sr <- ShortRead::srdistance(seq587, seq587[1])
sr


seq23911 <- Biostrings::DNAStringSet(ldproc[["23911"]][[4]][, "sequence"])
aln23911 <- DECIPHER::AlignSeqs(seq23911)
aln23911

# Compute point distance compared with the first seq
sr23911 <- ShortRead::srdistance(seq23911, seq23911[1])
sr23911


seq612 <- Biostrings::DNAStringSet(ldproc[["612"]][[4]][, "sequence"])
aln612 <- DECIPHER::AlignSeqs(seq612)
aln612

seq35495 <- Biostrings::DNAStringSet(ldproc[["35495"]][[4]][, "sequence"])
aln35495 <- DECIPHER::AlignSeqs(seq35495)
aln35495
#-----------------------------------------------------------------------------#

# Discarding alleles with less than 10 reads, and retain samples that have no more 
# two alleles 
lgen <- lapply(seq_tables, function(x) {
  #alleles <- colnames(x)
  gen <- apply(x, 1, function(x) {
    x[x<10] <- 0 # if count <10 replace to count zero
    
    newrow <- x[x>0]
    if(length(newrow) == 1) paste(rep(names(newrow), 2), collapse="/") else
      if(length(newrow) == 2) paste(names(newrow), collapse="/") else
        NA
  })
  Sample_ID <- names(gen)
  data.frame(Sample_ID, gen)
})


lSampleID <- lapply(lgen, function(x) {
  return(x$Sample_ID)
})
snames <- unique(unlist(lSampleID))

lgen <- c(list(data.frame(Sample_ID=snames)), lgen)
gen <- plyr::join_all(lgen, by="Sample_ID", type="left")
names(gen)<-c("Sample_ID", names(seq_tables))

library(adegenet)
gi <- df2genind(gen[, -1], sep = "/", ind.names = gen[,1])

write.csv(gen, file = file.path(dir.out, "genotype_table.csv"), 
          row.names = F)
save(gi, file = file.path(dir.out, "genind.rda"))

library(devtools)
install_github("nikostourvas/PopGenUtils")
library("PopGenUtils")
pid <- pid_calc(gi)
save(pid, file = file.path(dir.out, "pid.rda"))
