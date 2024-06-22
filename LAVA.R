remotes::install_github("josefin-werme/LAVA")
library(LAVA)

input = process.input(input.info.file="C:\\Users\\Derin\\Desktop\\Stat_R\\PD_Sup\\Data\\input_info.txt",
                      sample.overlap.file=NULL,   #set to NULL, since there is no overlap
                      ref.prefix="C:\\Users\\Derin\\Desktop\\Stat_R\\PD_Sup\\Data\\g1000_eur\\g1000_eur",                   
                      phenos=c("PD","sodium","potassium")) 

loci = read.loci("C:\\Users\\Derin\\Desktop\\Stat_R\\PD_Sup\\Data\\test.loci")

for (i in 1:nrow(loci)) {
  locus = process.locus(loci[i, ], input)
  if (!is.null(locus)) {
    print(c(locus$chr, locus$start, locus$stop))
    print(run.univ(locus))
  }
}
