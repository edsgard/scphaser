

##Params
sys = 'dna'
sys = 'rs13'


##*###
##Dirs
##*###

if(sys == 'dna'){
    cloud.dir = '/mnt/kauffman/edsgard/cloud/btsync/work/rspd'
}
if(sys == 'rs13'){
    cloud.dir = '/Volumes/Data/cloud/btsync/work/rspd'
}

##IN
in.dirs = c('../nogit/data/marinov', '../nogit/data/borel', '../nogit/data/mousehybrid')
names(in.dirs) = basename(in.dirs)

##OUT
pdf_dir = './ignore/res/combo/pdf'


##*###
##Files
##*###
in.files = file.path(in.dirs, 'ase2nvars.dens.rds')
names(in.files) = names(in.dirs)


main <- function(){

    
    ##*###
    ##Read data
    ##*###
    ase.list = lapply(in.files, readRDS)

    ##combine data
    library('reshape2')
    ase.df = melt(ase.list, id.vars = c('ase', 'cell', 'nvar.dens'))
    colnames(ase.df)[ncol(ase.df)] = 'dataset'

    ##*###
    ##Plot
    ##*###
    library('ggplot2')
    
    alpha = 0.7
    
    gg = ggplot(ase.df, aes_string(x = 'ase', y = 'nvar.dens', colour = 'dataset', fill = 'dataset'))
    gg = gg + geom_line(stat = 'summary', fun.y = 'median')
    ##gg = gg + stat_summary(geom = 'ribbon', fun.data = mean_se)
    ##gg = gg + stat_summary(geom = 'ribbon', fun.ymin = function(x){mean(x) - sd(x)}, fun.ymax = function(x){mean(x) + sd(x)})
    gg = gg + stat_summary(geom = 'ribbon', fun.ymin = function(x){quantile(x, probs = 0.25)}, fun.ymax = function(x){quantile(x, probs = 0.75)}, alpha = alpha)    
    
    ##background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())
    gg = gg + theme(axis.text = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))
    gg = gg + xlab('Allele specific expression')
    gg = gg + ylab('Number of variants (density)')

    ##plot
    j.pdf = file.path(pdf_dir, 'ase2nvars.dens.pdf')
    pdf.w = 4.5
    pdf.h = 3
    pdf(j.pdf, height = pdf.h, width = pdf.w)
    plot(gg)
    dev.off()
    
}
