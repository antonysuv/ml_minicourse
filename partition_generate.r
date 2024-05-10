library("optparse")
options(scipen = 100000000)

option_list = list(
  make_option(c("-n", "--nsim"), type="numeric", default=150000, help="Number of TRAIN simulations (default: 150000)", metavar="numeric"),  
  make_option(c("-l", "--len"), type="numeric", default=1000, help="Alignment length for IQTREE partiton file (default: 1000)", metavar="numeric"),
  make_option(c("-f", "--fname"), type="character", default="simulation", help="Name prefix for main file outputs (default: simulation)", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Generate partition file for IQTREE simulatior
get.nexus.part=function(nsim,aln_len,file)
{
    file = paste(file,".part",sep="")  
    write("#nexus\nbegin sets;",file)
    write.table(data.frame('\tcharset',
                           paste("locus",1:nsim,sep="_"),
                           "=",
                           paste("DNA",",",sep=""),
                           paste(seq(1,aln_len*nsim,by=aln_len),
                               "-",
                                 paste(seq(aln_len,aln_len*nsim,by=aln_len),";",sep=""),
                           sep="")),
                           file,
                           append = T,
                           quote = F,
                           row.names = F,
                           col.names = F) 
    write("end;",file,append=T)
}
get.nexus.part(opt$nsim,opt$len,opt$fname)
