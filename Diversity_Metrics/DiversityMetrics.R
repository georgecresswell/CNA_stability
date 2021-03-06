###########################################################
# Diversity measures 1, 2, 4=5 (depending on option used) #
###########################################################
diff_state_dist = function(x, y, only_alt_bins = T, absolute_number = F, ploidy = 2) {
  
  frac_genome_alt = length(which(x!=y)) / length(x)
  
  if(only_alt_bins) {
    
    bins_diff = length(which(x!=y))

    c = cbind(x, y)

    bins_ab = length(which(apply(c, 1, function(i) any(i!=ploidy))))

    output = bins_diff / bins_ab
  
  } else {output = frac_genome_alt}

  if(absolute_number) {output = length(which(x!=y))}
  
  return(output)
  
}

#########################################################
# Salpie's divergence measure with ploidy normalisation #
#########################################################
calculateDivergence <- function(to_compare, ploidy = 2) {
    fab <- as.data.frame(table(colSums(to_compare<ploidy)))
    if (length(fab$Freq[fab$Var1==1] > 0)) {
        if (fab$Freq[fab$Var1==1] > 0) {
            score_loss <- fab$Freq[fab$Var1==1]
        } else {
            score_loss <- 0
            }       
        } else {
            score_loss <- 0     
        }
     fab <- as.data.frame(table(colSums(to_compare>ploidy)))
   if (length(fab$Freq[fab$Var1==1] > 0)) {
        if (fab$Freq[fab$Var1==1] > 0) {
            score_gain <- fab$Freq[fab$Var1==1]
        } else {
            score_gain <- 0
            }       
        } else {
            score_gain <- 0     
        }
        anueploidy <- length(to_compare[to_compare<ploidy|to_compare>ploidy])
        loss_anueploidy <- length(to_compare[to_compare<ploidy])
        gain_anueploidy <- length(to_compare[to_compare>ploidy])
    scores <- c(paste(rownames(head(to_compare))[1],rownames(head(to_compare))[2], sep="_"), score_loss, score_gain, anueploidy, loss_anueploidy, gain_anueploidy)

    divergence_score = (score_loss + score_gain) / anueploidy
    
    return(divergence_score)
}

###################################################
# Diversity measures 3 (depending on option used) #
###################################################
genetic_distance = function(x, y, normalise_by_bin_number = T) {
  
  dist = sum(abs(x-y))
  
  if(normalise_by_bin_number) {dist = dist / length(x)}
  
  return(dist)
  
}

##############################################
# L2RSS - a comparison of log2ratio profiles #
##############################################
log2ratio_comparison = function(segs_col_a, segs_col_b, exp_distance = 1848.691, normalise_to_exp = T, min_purity = 0.2) {
  
  # Calculate contineous copy number
  calcCN = function(lrrs, rho, psit, gamma = 1) {
    
    psi = (2*(1 - rho)) + (rho*psit)
    
    n = ((psi*(2^(lrrs/gamma))) - (2 * (1 - rho))) / rho
    
    return(n)
    
  }
  
  # What is our parameter search of purities?
  parameter_comparison = rbind(cbind(seq(min_purity, 0.99, by = 0.01), 1),
                               cbind(1, seq(1, min_purity, by = -0.01)))
  
  # Here we do a search of purity pairs
  search = lapply(1:nrow(parameter_comparison), function(r) {
    
    # Selected parameters for iteration
    rhoA = parameter_comparison[r,1]
    rhoB = parameter_comparison[r,2]
    
    # Continuous copy number calculation 
    CNa  = calcCN(lrrs = segs_col_a, rho = rhoA, psit = 2)
    CNb  = calcCN(lrrs = segs_col_b, rho = rhoB, psit = 2)
    
    # Sum of squared differences (maybe normalise for number of bins?)
    dist = sum((CNa - CNb)^2)
    
    return(dist)
    
  })
  
  # Distance results for parameter comparisons
  res = cbind(parameter_comparison, unlist(search))
  
  # Which has the shortest distance
  R = which.min(res[,3])
  
  # Get the d
  d = res[R,3]
  
  if(normalise_to_exp) {
  
    # Normalise the distance to the cohort (hard coded for now)
    d = d / exp_distance # This number is the median dist in non-same patient comparisons
    if(d>1) {d = 1} # Cap at 1
    
  }
  
  return(d)
  
}

#######################################
# Diversity measure 6 helper function #
#######################################
armCN = function(df, pqs, method = c("median", "mean"), l2r_col = 4, report_NA = F) {
  
  method = match.arg(method)
  
  chrs = unique(df$chromosome)
  
  per_chr = lapply(chrs, function(c) {
    
    chrp = df[df$chromosome==c & pqs=="p",l2r_col]
    chrq = df[df$chromosome==c & pqs=="q",l2r_col]
    
    if(method == "median") {
      
      p = median(chrp, na.rm = T)
      q = median(chrq, na.rm = T)
      
    }
    
    if(method == "mean") {
      
      p = mean(chrp, na.rm = T)
      q = mean(chrq, na.rm = T)
      
    }
    
    out = c(p, q)
    
    names(out) = paste0(c,c("p","q"))
    
    return(out)
    
  })
  
  out = unlist(per_chr)
  
  if(!report_NA) {out = out[!is.na(out)]}
  
  return(out)
  
}

############################################################
# Using biomaRt to calculate number of genes per bin       #
# - may be really slow, might be better to download refseq #
# - only compatible with hg38 right now                    #
############################################################
countGenesPerBin = function(bins, genome = "hg38") {
  
  genome = match.arg(genome)
  
  require("biomaRt")
  
  # Set up BiomaRt
  listMarts(host="www.ensembl.org")
  ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
  filters = listFilters(ensembl)
  
  coords = paste0(bins$chromosome,":",bins$start,":",bins$end)
  
  coords = as.list(coords)
  
  count_entries = lapply(coords, function(b) {
    
    # Get overlapping genes from biomaRt
    results=getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                  filters = c("chromosomal_region", "biotype"),
                  values = list(chromosomal_region=b, biotype="protein_coding"), 
                  mart = ensembl)
    
    nrow(results)
    
  })
  
  count_entries = unlist(count_entries)
  
  out = data.frame(bins[,1:3], gene_number = count_entries)
  
  return(out)
  
}

####################################################
# Using biomaRt to assign each gene to a bin       #
# - only compatible with hg38 right now            #
####################################################
getGeneBinIndex = function(bins_locs, genome = "hg38", chrs = 1:22, saveFile="~/Downloads/human_protein_encoding_genes_ensembl.rds") {
  
  genome = match.arg(genome)
  
  require("biomaRt")
  
  if (!file.exists(saveFile)){
    print("Querying Biomart for protein coding genes")
    ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
    humanProteinCodingGenes = getBM(attributes = c("chromosome_name", "start_position", "end_position", "hgnc_symbol"),
                                    filters = c("biotype"),
                                    values = list(biotype="protein_coding"), 
                                    mart = ensembl)
    saveRDS(humanProteinCodingGenes,file=saveFile)
  } else {
    print("Loading genes from savefile")
    humanProteinCodingGenes = readRDS(saveFile)
  }
  
  # Subset the autosomes
  hPCG_autosomes = humanProteinCodingGenes[humanProteinCodingGenes$chromosome_name %in% chrs,]
  hPCG_autosomes = hPCG_autosomes[order(hPCG_autosomes$chromosome_name, hPCG_autosomes$start_position),]
  rownames(hPCG_autosomes) = NULL
  hPCG_autosomes$mid_point = as.numeric(round(((hPCG_autosomes$end_position - hPCG_autosomes$start_position)/2) + hPCG_autosomes$start_position))
  
  # Run through as a list
  gene_bin_indexes = lapply(1:nrow(hPCG_autosomes), function(i) {
    
    g = hPCG_autosomes[i,]
    
    chr = g[,"chromosome_name"]
    mid = g[,"mid_point"]
    
    which(bins_locs$chr==chr & bins_locs$start<=mid & bins_locs$end>=mid)
    
  })
  
  # Produce the bin indexes
  return(unlist(gene_bin_indexes))
  
}

#######################################
# Breakpoint functions from Salpie :) #
#######################################
convertToBreakpoints = function(cnTable){
  y = cnTable
  y[y > 0] <- 0
  
  for (column in 1:ncol(cnTable)) {
    breakpoints = (which(!!diff(as.numeric(cnTable[,column])))+1) #get indexes
    y[c(breakpoints),column] <- 1
  }
  return(y)
}

calculateRelatednessCn = function(cnTable, pairs, maxgap){
  
  pair_scores <- apply(pairs, 1, function(x){getScoreCN(cnTable, populationBreakpoints, maxgap, as.character(x))})
  
  results <- cbind.data.frame(pairs, pair_scores)
  
  return(results)
}

getScoreCN = function(cnTable, maxgap, pairs){
  
  sample1 <- cnTable[,c(colnames(cnTable) == "Chr" | colnames(cnTable) == "Start" | colnames(cnTable) == "End" | colnames(cnTable) == pairs[1])]
  sample2 <- cnTable[,c(colnames(cnTable) == "Chr" | colnames(cnTable) == "Start" | colnames(cnTable) == "End" | colnames(cnTable) == pairs[2])]
  row_sample1 = apply(sample1, 1, function(row) all(row !=0 ))
  sample1 <- 	sample1[row_sample1,]
  row_sample2 = apply(sample2, 1, function(row) all(row !=0 ))
  sample2 <- 	sample2[row_sample2,]
  
  if (empty(sample1) | empty(sample2)){
    score = 0
  } else {
    #tryCatch creates an empty GRanges object if the list is empty - would error out otherwise
    sample1_granges <- makeGRangesFromDataFrame(sample1[,c("Chr", "Start", "End")], start.field = "Start", end.field = "End")
    sample2_granges <- makeGRangesFromDataFrame(sample2[,c("Chr", "Start", "End")], start.field = "Start", end.field = "End")
    
    hits_start <- suppressWarnings(queryHits(findOverlaps(sample1_granges, sample2_granges, type = "start", maxgap = maxgap)))
    hits_end <- suppressWarnings(queryHits(findOverlaps(sample1_granges, sample2_granges, type = "end", maxgap = maxgap)))
    
    nconcordant_adj <- 2*(length(hits_start)+length(hits_end))
    total_breakpoints <- sum(2*length(sample1_granges)+2*length(sample2_granges))
    
    score = (total_breakpoints-nconcordant_adj)/total_breakpoints
  }
  return(score)
}


###################################
# Average mean length differences #
###################################
get_diffLengths = function(cnTable, pairs, pp = 2, max_size = FALSE) {
  
    # Make a dataframe
    a = cnTable[,c(colnames(cnTable) == "Chr" | colnames(cnTable) == pairs[1])]
    b = cnTable[,c(colnames(cnTable) == "Chr" | colnames(cnTable) == pairs[2])]
    a[,2] = as.numeric(as.character(a[,2]))
    b[,2] = as.numeric(as.character(b[,2]))
    len = NULL
    for (chrs in 1:length(unique(a$Chr))) {
      sub_a = subset(a, a$Chr == unique(a$Chr)[chrs])
      sub_b = subset(b, b$Chr == unique(a$Chr)[chrs])
      pasted = paste0(sub_a[,2], sub_b[,2])
      diff_bins = unlist(lapply(strsplit(rle(pasted)$values, split = ""), function(i) all(!duplicated(i))))
      lengths = rle(pasted)$lengths[diff_bins]
      len = c(len,lengths)
    }
    
    # Catch times when there is no difference and record it as zero
    if(length(len)==0) {len = 0}

    if(max_size) {output = max(len)} else {output = mean(len)}
    
    return(output)
    
}