
###Assess performance under varying parameters
##NB: The paths are relative to the root of the scphaser git repo

##Dirs
out_data_dir = './ignore/res/perf/data'
out_pdf_dir = './ignore/res/perf/pdf'

##Files
perf_rds = file.path(out_data_dir, 'perf.rds')

##Libs
source('./ignore/code/analysis/performance.R')

main <- function(){


    ##*###
    ##Read and prep data
    ##*###
    load('../nogit/data/mousehybrid/mousehybrid.rda')

    ##create acset
    acset = new_acset(mousehybrid$featdata, mousehybrid$refcount, mousehybrid$altcount, mousehybrid$phenodata)
    lapply(acset, dim)
    
    ##subset
    ngenes = 100
    acset = subset_ngenes(acset, ngenes)

    
    ##*###
    ##Specify params
    ##*###
    
    ##phasing method params
    input = c('gt', 'ac')
    method = c('exhaust', 'cluster')
    weigh = c(TRUE, FALSE)    
    phase_params = expand.grid(input, method, weigh, stringsAsFactors = FALSE)
    colnames(phase_params) = c('input', 'method', 'weigh')
        
    ##permutation iterations
    npermiter = 3
    perm_iter = 1:npermiter

    ##min read count for genotyping
    min_acount = c(3, 10)
    
    ##allelic fold-change for mono-allelic call
    fc = c(3, 4)
    
    ##filter vars on min number of cells with mono-allelic call
    nmincells = c(5, 10)


    ##*###
    ##Calculate performance
    ##*###
    ##loop phase_params and get a dataframe of performance with respect to the other params for each phasing parameter set
    n_phasemethods = nrow(phase_params)
    perf_list = list()
    length(perf_list) = n_phasemethods
    phasing_str = sub(' ', '', apply(phase_params, 1, paste, collapse = '.'))
    names(perf_list) = phasing_str
    for(jmethod in 1:n_phasemethods){
        print(jmethod)
        jphase_params = phase_params[jmethod, ]
        perf_list[[jmethod]] = get_perf_paramset(acset, perm_iter, min_acount, fc, nmincells, jphase_params[1, 'input'], jphase_params[1, 'weigh'], jphase_params[1, 'method']);
    }

    ##Dump
    saveRDS(perf_list, file = perf_rds)

    
    ##*###
    ##Plot performance
    ##*###
    ##TODO: matrix of plots for each phasing parameter set

    perf_list = readRDS(perf_rds)
    paramsets = names(perf_list)
    for(jset in paramsets){

        jperf = perf_list[[jset]]
        
        j.pdf = file.path(pdf.dir, paste(jset, '.mcc.pdf', sep = ''))
        pdf(file = j.pdf)
        plot_perf(jperf)
        dev.off()
    }

}

plot_perf <- function(obs.long, your.x, your.colour, y.lab = 'MCC'){

    ##data and mappings
    gg = ggplot(obs.long, aes_string(x = 'min_acount', y = 'mcc', colour = your.colour, group = group.col))

    ##facet
    gg = gg + facet_grid(fc ~ nmincells, scales = 'free_y')
    
    ##mean value as line
    gg = gg + stat_summary(fun.y = median, geom = "line")

    ##95% CI
    gg = gg + stat_summary(data = obs.long, fun.data = mean_cl_boot, geom = "linerange") ## geom = 'errorbar'
    
    ##Background
    gg = gg + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(panel.background = element_blank())

    ##x ticks
    gg = gg + theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), axis.ticks = element_line(colour = 'black'))

    ##gg = gg + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    ##axis lables
    gg = gg + ylab(y.lab) ##+ xlab(x.lab)

    ##colors
    ##gg = gg + scale_color_manual(values = lin2color)

}

get_perf_paramset <- function(acset, perm_iter, min_acount, fc, nmincells, input, weigh, method){

    
    ##specify paramset as all possible combinations of the params
    paramset = expand.grid(perm_iter, min_acount, fc, nmincells, stringsAsFactors = FALSE)
    colnames(paramset) = c('perm_iter', 'min_acount', 'fc', 'nmincells')
    nparamset = nrow(paramset)

    ##loop param set
    perf_list = list()
    length(perf_list) = nparamset
    for(jparam in 1:nparamset){
        jparamset = paramset[jparam, ]
        jperf = filter_get_perf(acset, jparamset[, 'min_acount'], jparamset[, 'fc'], jparamset[, 'nmincells'], input, weigh, method)
        perf_list[[jparam]] = jperf[['var']]
    }
    perf_df = do.call(rbind, perf_list)

    ##bind params and perf
    params2perf_df = cbind(paramset, perf_df)
        
    return(params2perf_df)
}
