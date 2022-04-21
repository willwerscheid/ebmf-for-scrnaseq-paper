# Pulse-seq dataset from Montoro et al. Run from project home directory.

zz <- download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103354&format=file&file=GSE103354%5FPulseSeq%5FUMI%5Fcounts%2Erds%2Egz",
                    destfile = "./data/pulseseq.rds.gz")
R.utils::gunzip("./data/pulseseq.rds.gz")
