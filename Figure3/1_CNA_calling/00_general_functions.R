# Functions

# Normalise in range
range01 = function(x){(x-min(x))/(max(x)-min(x))}

# Remove zero function
removeZero = function(x) {x[x!=0]}

# Get tertiles
tertile = function(x) {quantile(x, probs = c(1/3, 2/3))}

# Remove a column's NAs
colNArm = function(x, sub) {x[!is.na(x[,sub]),]}

# Make a summary list
makeSummaryList = function(x) {
  
  res = list()
  length(res) = length(x)
  names(res) = x
  
  return(res)
  
}

# Make location dataframe from chr pos rownames
extractChrRows = function(x) {
  
  m = do.call(rbind, strsplit(x, split = "[:]|-"))
  
  df = data.frame(chr = m[,1],
                  start = as.integer(m[,2]),
                  end = as.integer(m[,3]))
  
  df$chr = factor(df$chr, levels = c(1:22, "X", "Y"))
  
  return(df)
  
}

# Remove generic bam name end
removeBamNameEnding = function(x) {gsub("_dups_valid", "", x)}

# Extract patient name
extractPatientName = function(x) {unlist(lapply(strsplit(x, split = "_"), function(i) i[1]))}

# Get p and q arm assignment
getPQ = function(df, arm_pos) {
  
  psnqs = unlist(lapply(1:nrow(df), function(i) {
    
    r = df[i,]
    
    pos = (r$end - r$start) + r$start
    
    chr = r$chromosome
    
    b = arm_pos[which(arm_pos[,1] == chr),2]
    
    if(pos <  b) {a = "p"}
    if(pos >= b) {a = "q"}
    
    return(a)
    
  }))
  
  return(psnqs)
  
}

# Summarise list to dataframe
list_summary_df = function(l, value_name, summary = c("mean", "median", "sum")) {
  
  summary = match.arg(summary)
  
  if(summary=="mean")   {collapsed_data = unlist(lapply(l, mean))}
  if(summary=="median") {collapsed_data = unlist(lapply(l, median))}
  if(summary=="sum")    {collapsed_data = unlist(lapply(l, sum))}
  
  df = data.frame(Patient = names(collapsed_data), value = collapsed_data)
  
  colnames(df)[2] = value_name
  
  rownames(df) = NULL
  
  return(df)
  
}

# Collapse calls to segments
collapseAsSegments = function(calls, chrs = 1:22) {
  
  do.call(rbind, lapply(chrs, function(chr) {
    
    index = grep(paste0("^", chr,":"), rownames(calls), perl = T)
    
    bin_names = rownames(calls)[index]
    
    chr_cn = calls[index,1]
    
    rle_res = rle(chr_cn)
    
    loc_annot = extractChrRows(bin_names)
    
    data.frame(chr = chr, start = loc_annot$start[c(1, (cumsum(rle_res$lengths)+1)[-length(rle_res$lengths)])], 
               end = loc_annot$end[cumsum(rle_res$lengths)], CN = rle_res$values)
    
  }))
  
}

# Make a function called count
count = function(x) {length(which(x))}

# Make a function that just does a fraciton
fraction = function(x) {length(which(x)) / length(x)}

# Colour with saturation
cols = c(c("0" = "#1981be", "1" = "#56B4E9", "2" = "grey", "3" = "#E69F00", "4" = "#ffc342",
           "5" = "#FFAA42", "6" = "#FF9142", "7" = "#FF7742", "8" = "#FF5E42"), 
         rep("#FF4542", times = 100 - 8))

# Name top
names(cols)[(9:100)+1] = 9:100

dist2d = function(data_point,segment_point1=c(0,1),segment_point2=c(1,0), absolute=F) {
  v1 = segment_point1 - segment_point2
  v2 = data_point - segment_point1
  m = cbind(v1,v2)
  if(absolute) {d = abs(det(m))/sqrt(sum(v1*v1))} else {
    
    d = det(m)/sqrt(sum(v1*v1))
    d = d*-1
    
  }
  return(d)
} 
