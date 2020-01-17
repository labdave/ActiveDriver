library(optparse)
option_list <- list ( 
  make_option(c("-o", "--outputdir"), dest = "out_dir", type = "character", help = "output dir name"),
  make_option(c("-r", "--refpath"), dest = "ref_path", type = "character", help = "reference file path"),
  make_option (c("-f","--filelist"),default="blah.txt", help="comma separated list of files (default %default)")
  
)

parser <-OptionParser(option_list=option_list)
arguments <- parse_args (parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args

#myfilelist <- strsplit(opt$filelist, ",")
myfilelist <- unlist(strsplit(opt$filelist, ","))

output.prefix <- opt$out_dir
supp_path = opt$ref_path

print(output.prefix)
print(supp_path)
print(class(myfilelist))
all_vcfs = data.frame()

# Loop through recoded-vcfs converting to ActiveDriver format
for (i in myfilelist) {
  print(paste("Reading", i, "..."))
  x = read.table(i, sep = "\t",header=T, stringsAsFactors = F)
  # Filter based on read-count
  rc = names(x)[(ncol(x)-3)] # (assumes read-count column is 4th from last)
  x = x[which(x[,rc] >= 0.3),] # WARNING: hard-coded!!!
  # Filter based on mutation being non-synonymous SNV
  x = x[x$ExonicFunc.refGene == "nonsynonymous_SNV",]
  # Check if any mutations passed filters
  if (nrow(x) > 0) {
    # Parse out mutation-specific values
    x$gene = NA
    x$hugo = NA
    x$position = NA
    x$wt_residue = NA
    x$mut_residue = NA
    for (i in 1:nrow(x)) {
      my_vect = unlist(strsplit(x[i,"AAChange.refGene"],":"))
      x[i,"hugo"] = my_vect[1]
      x[i,"gene"] = my_vect[2]
      aa_change = gsub("p.","",my_vect[5])
      x[i,"position"] = gsub("[A-Z]","",aa_change)
      x[i,"wt_residue"] = substring(aa_change,0,1)
      x[i,"mut_residue"] = substring(aa_change,nchar(aa_change), nchar(aa_change))
    }
    # Add sample-name column
    sn = substr(rc,2,nchar(rc)) # get sample name (assumes R added an "X" to the column)
    x[,"sample_id"] = sn
    # Add meaningless "cancer_type" column
    x[,"cancer_type"] = "cancer"
    # Create output table
    print(paste("Saving", nrow(x), "variants!"))
    outx = x[,c("gene","cancer_type","sample_id","position","wt_residue","mut_residue","hugo")]
    all_vcfs = rbind(all_vcfs,outx)
  } else {
    print("No variants passed filters.")
  }
}

# Create gene lookup table 
lt = unique(all_vcfs[,c("gene","hugo")])
# Format and save input table
all_vcfs = all_vcfs[,-grep("hugo",names(all_vcfs))]
write.table(all_vcfs,paste(output.prefix, "input.tsv",sep="."),sep="\t",quote = F,row.names = F,col.names = T)

library(ActiveDriver)

# Load support files
sites = read.delim(paste(supp_path,"/PTM_sites.txt",sep=""))
seqs = read_fasta(paste(supp_path,"/refseq_protein_sequences.txt",sep=""))
disorder = read_fasta(paste(supp_path,"/refseq_protein_sequence_disorder.txt",sep=""))

# run ActiveDriver
go = ActiveDriver(seqs, disorder,all_vcfs, sites)

# Get table of pvalues
ptab = merge(lt, go$all_gene_based_fdr, by = "gene", all.x = F, all.y = T)
# Sort by FDR
ptab = ptab[order(ptab$fdr),]
write.table(ptab,paste(output.prefix, "fdr.tsv",sep="."),sep="\t",quote = F,row.names = F,col.names = T)

# Get merged report
mrep = merge(lt, go$merged_report, by = "gene",all.x =F, all.y = F)
write.table(mrep,paste(output.prefix, "merged_report.tsv",sep="."),sep="\t",quote = F,row.names = F,col.names = T)

# Add input and lookup table to output and save as RDS object
go[["input"]] = all_vcfs
go[["gene_dictionary"]] = lt
saveRDS(go, file = paste(output.prefix, "RDS",sep="."))

