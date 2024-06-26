remotes::install_github("josefin-werme/LAVA")
library(LAVA)

input = process.input(input.info.file="input_info.txt",
                      sample.overlap.file=NULL,  #set to NULL, since there is no overlap
                      ref.prefix="g1000_eur",            
                      phenos=c("PD","sodium","potassium")) 

loci = read.loci("test.loci") #preprocessed file obtained from https://github.com/josefin-werme/LAVA, added manually 2 loci from GWAS sumstats files

for (i in 1:nrow(loci)) {
  locus = process.locus(loci[i, ], input)
  if (!is.null(locus)) {
    print(c(locus$chr, locus$start, locus$stop))
    print(run.univ(locus))
  }
}
