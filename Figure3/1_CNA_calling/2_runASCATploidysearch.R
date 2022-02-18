library(copynumber)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
# library(cowplot)
# theme_set(theme_cowplot())

# Plotting settings
run_full_pipeline = TRUE

source("00_general_functions.R")
source("00_runASCATlp.R")

ploid_search = seq(1.5,4,by=0.1)
min_purity   = 0.5
max_lrr      = Inf

blks = unlist(lapply(strsplit(list.files("cna_raw/segments"), "_"), function(i) i[2]))

impure_list = blks[grep("-B[0-9]$", blks)]

patients = unique(unlist(lapply(strsplit(list.files("cna_raw/segments"), split = "_"), 
                                function(i) i[1])))

patients = patients[!patients %in% "LM0041"]

per_patient_solution = list()

for(patient in patients) {

  # Find files
  files = list.files("cna_raw/segments", 
                     pattern = patient, full.names = T)
  
  # Read in data
  segs = lapply(files, function(f) read.table(f)[,1])
  
  # Find files
  files = list.files("cna_raw/bins", 
                     pattern = patient, full.names = T)
  
  # Read in data
  bins = lapply(files, function(f) read.table(f, stringsAsFactors = F))
  
  # Extract names
  sample    = gsub("_dups_valid", "", unlist(lapply(bins, function(b) colnames(b)[5])))
  sample_df = data.frame(do.call(rbind, strsplit(sample, split = "_")), stringsAsFactors = F)
  colnames(sample_df) = c("Project", "Patient", "Block", "Region", "Extraction.Genome")
  sample_df$Paper_name = unlist(lapply(strsplit(basename(files), "_"), function(i) i[2]))
  
  # Gather here the chromosome position
  chr_pos = bins[[1]][,2:4]
  
  ploidy_mean_dist = lapply(ploid_search, function(pld) {
    
    output = list()
  
    ascat = lapply(segs, runASCATlp, fix_ploidy=pld, 
                   min_purity=min_purity, max_lrr = max_lrr, no_fit_psit = 2)
    
    for(s in paste0(sample_df$Block,sample_df$Region)) {
      
      if(s %in% impure_list) {
        
        impure_sample = which(paste0(sample_df$Block,sample_df$Region) == s)
        
        ascat[[impure_sample]] = runASCATlp(segs[[impure_sample]], fix_ploidy=pld, 
                                            min_purity=0.2, max_lrr = max_lrr, no_fit_psit = 2)
      
      }
      
    }
    
    mat = do.call(rbind, lapply(ascat, function(s) s$CN)) / pld
    
    mat = mat[,!apply(do.call(cbind, segs), 1, max) > 1]
    
    dists = as.vector(dist(mat))
    
    output$int_dists = dists
    output$sol_dists = unlist(lapply(ascat, function(i) {
      
      within_dist = i$allsol$dist[1]
      
      # If not fit, you give the measure of fit for a pure sample
      if(is.na(within_dist)) {within_dist = i$fit[rownames(i$fit)==1,]}
      
      return(within_dist)
      
    }))
    
    return(output)
  
  })
  
  inter_sample_dist = unlist(lapply(ploidy_mean_dist, function(i) mean(i$int_dists)))
  CN_solution_dist = unlist(lapply(ploidy_mean_dist, function(i) mean(i$sol_dists, na.rm=T)))
  df = data.frame(inter_sample = inter_sample_dist, in_sample = CN_solution_dist)
  df_scaled = data.frame(apply(df, 2, range01))
  dist_from_zero = apply(df_scaled, 1, function(x) dist(rbind(x, rep(0, times = 2))))
  ploidy_sol_index = which.min(dist_from_zero)
  
  plt.df = df_scaled
  plt.df$Solution = FALSE
  plt.df$Solution[ploidy_sol_index] = TRUE
  
  cols = c("FALSE" = "black", "TRUE" = "red")
  
  plt = ggplot(plt.df, aes(x = in_sample, y = inter_sample, col = Solution)) + geom_point() + 
    annotate("text", x = plt.df$in_sample, y = plt.df$inter_sample+0.03, label = ploid_search) + 
    scale_colour_manual(values = cols) +
    xlab("Within sample mean rel. distance") + ylab("Between samples rel. distance") +
    theme_bw()
  
  if(run_full_pipeline) {
  
    print(plt)

  }
  
  searched_ploidy = ploid_search[ploidy_sol_index]
  
  # Take second option
  if(patient %in% c("LM0032", "LM0082")) {searched_ploidy = ploid_search[order(dist_from_zero)[2]]}
  
  if(run_full_pipeline) {
  
    plot(ploid_search, unlist(lapply(ploidy_mean_dist, function(i) mean(i$int_dists))), type = "l",
         xlab = "Ploidy", ylab = "Mean Intersample Distance")

  }
  
  ascat = lapply(segs, runASCATlp, fix_ploidy=searched_ploidy, 
                 min_purity=min_purity, max_lrr = max_lrr, no_fit_psit = 2)
  
  for(s in paste0(sample_df$Block,sample_df$Region)) {
    
    if(s %in% c(impure_list, "FGR1")) {
      
      impure_sample = which(paste0(sample_df$Block,sample_df$Region) == s)
      
      ascat[[impure_sample]] = runASCATlp(segs[[impure_sample]], fix_ploidy=searched_ploidy, 
                                          min_purity=0.2, max_lrr = max_lrr, no_fit_psit = 2)
      
    }
    
    # Here we also do some manual refits
    if(s %in% c("ABB2", "FGR2")) {
      
      refit_sample = which(paste0(sample_df$Block,sample_df$Region) == s)
      
      if(s %in% impure_list) {refit_purity = 0.2} else {refit_purity = 0.5}
      
      ascat[[refit_sample]] = runASCATlp(segs[[refit_sample]], fix_ploidy=searched_ploidy+0.4, 
                                         min_purity=refit_purity, max_lrr = max_lrr, no_fit_psit = 2)
      
    }
    
    if(s %in% c("FLR1", "FLR2", "FLR3", "ELR1")) {
      
      refit_sample = which(paste0(sample_df$Block,sample_df$Region) == s)
      
      if(s %in% impure_list) {refit_purity = 0.2} else {refit_purity = 0.5}
      
      ascat[[refit_sample]] = runASCATlp(segs[[refit_sample]], fix_ploidy=searched_ploidy-0.4, 
                                         min_purity=refit_purity, max_lrr = max_lrr, no_fit_psit = 2)
      
    }
    
  }
  
  names(ascat) = sample
  
  if(run_full_pipeline) {
  
    for(s in 1:length(ascat)) {
      
      mean_dist = mean(ascat[[s]]$CN - ascat[[s]]$contCN)
      hist(ascat[[s]]$CN - ascat[[s]]$contCN, breaks = 100, xlim = c(-1,1), main = names(ascat)[s], xlab = "Distance")
      abline(v = mean_dist, lty = "dotted")
      
    }
    
  }
  
  for(s in 1:length(ascat)) {if(nrow(ascat[[s]]$allsol)==0) {print(names(ascat)[s])}}
  
  mat = do.call(rbind, lapply(ascat, function(s) s$CN))
  
  # Colour with saturation
  cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
             "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
           rep("#FF4542", times = 100 - 8))
  
  # Name top
  names(cols)[(9:100)+1] = 9:100
  
  plts = lapply(1:length(ascat), function(s) {
    
    sample_ascat = ascat[[s]]
    
    # Make a plotting dataframe
    plt.df = data.frame(genome.bin = 1:nrow(bins[[s]]), 
                        chr = bins[[s]]$chromosome, 
                        Log2ratio = bins[[s]][,5],
                        Call = as.character(sample_ascat$CN),
                        mean_segment = segs[[s]],
                        segment_col = "green")
    
    # Make plot
    p = ggplot(plt.df, aes(x = genome.bin, y = Log2ratio, col = Call)) +
      geom_hline(yintercept = c(-2,-1,0,1,2), lty = c("solid"), lwd = 0.2) +
      geom_point() +
      scale_colour_manual(values = cols) +
      scale_x_continuous(name = "Chromosomes", labels = 1:22, 
                         breaks = as.vector(c(1, cumsum(table(bins[[s]]$chromosome))[-22]) + 
                                              (table(bins[[s]]$chromosome) / 2))) + 
      geom_vline(xintercept = c(1, cumsum(table(bins[[s]]$chromosome))), lty = "dotted") +
      ggtitle(paste0("Low pass calls - ",names(ascat[s])," purity=",
                     sample_ascat$Purity,", psit = ",sample_ascat$PsiT,
                     ", ploidy = ",round(mean(sample_ascat$CN), digits = 2),")")) + 
      scale_y_continuous(limits=c(-2,2), oob=scales::squish) + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 18)) +
      geom_point(aes(y = mean_segment), color="#000000")
    
    print(p)
    
    return(p)
    
  })
  
  names(plts)   = names(ascat)
  rownames(mat) = names(ascat)
  
  if(run_full_pipeline) {
  
    lapply(1:length(ascat), function(i) {
      
      CN_out = data.frame(ascat[[i]]$CN)
      
      colnames(CN_out) = names(ascat)[i]
      
      rownames(CN_out) = bins[[i]]$feature
      
      write.table(CN_out, file = paste0("../2_assessing_diversity/cna_calls/",sample_df[i,"Patient"],"_",
                  sample_df[i,"Paper_name"],"_500kb_GRCh38_cna_calls.txt"),
                  quote = F)
      
    })
    
  }
  
  
  per_patient_solution[[patient]] = searched_ploidy

}

# Make a ploidy list per case
per_case_psit = data.frame(Case = names(per_patient_solution),
                           PsiT = unlist(per_patient_solution))
rownames(per_case_psit) = NULL

if(run_full_pipeline) {

  # Write out the psit 
  write.table(per_case_psit, row.names = F, col.names = T, quote = F, sep = "\t",
              file = "patient_psit_solutions.txt")
  
}
