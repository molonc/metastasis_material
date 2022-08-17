tree_plot <- function(df = pick, series = "SA919", source = "airtable", resultdir = outputdir){
  resultdir <- paste0(resultdir, source, sep="/")
  if(!dir.exists(resultdir)) dir.create(resultdir)
  
  print(series)
  # series_long <- paste(ifelse(is.na(pat$SA_ID[pat$SA.Identifiers == series]), "", paste(pat$[pat$SA.Identifiers == series]), "_", sep="")), series, sep="")
  xenos <- subset(df, SA.Identifiers == series  )[, c("Transplant_ID", "AtiM.Patient.ID", 
                                                               "Passage", "Transplant.Material")]
  # xenos <- subset(tsp, Patient_ID == series & Grown == 1 & !is.na(Transplant_ID))[, c("Transplant_ID", "Grown")]
  temp = xenos[,  c("AtiM.Patient.ID","Transplant_ID", "Passage", "Transplant.Material")]
  xenos$Xeno <- xenos$Transplant_ID
  if(length(unique(xenos$Xeno)) != length(xenos$Xeno)){
    print("Transplant ID is not unique")
    break
  }
    
  
  l = strsplit(xenos$Transplant_ID, "-")
  
  Transplant_ID_temp = data.frame(plyr::ldply(l, rbind))
  
  Transplant_ID_temp$Transplant_ID= gsub("NA", "", apply(Transplant_ID_temp[,-c(1)], 1, paste0, collapse = ""))
  
  Transplant_ID_temp$Transplant_ID_new = gsub("-NA", "", 
                                              apply(Transplant_ID_temp[, !names(Transplant_ID_temp) %in% c("Transplant_ID")], 
                                                    1, paste0, collapse = "-"))
  
  if(length(unique(Transplant_ID_temp$Transplant_ID_new)) != length(Transplant_ID_temp$Transplant_ID_new)){
    print("Newly generated Transplant ID is not unique")
    break
  }
  # xenos$Transplant_ID <- sapply(strsplit(xenos$Transplant_ID, "-"), tail, 1)
  xenos = merge(xenos, Transplant_ID_temp[, c("Transplant_ID", "Transplant_ID_new")], by.x = "Transplant_ID", by.y = "Transplant_ID_new")
  xenos$pathString <- series
  
  xenos$Transplant_ID = xenos$Transplant_ID.y
  xenos = xenos[, !names(xenos) %in% c('Transplant_ID.y')]
  
  for(i in 1:nrow(xenos)){
    for(j in 1:nchar(xenos$Transplant_ID[i])){
      xenos$pathString[i] <- paste(xenos$pathString[i], substr(xenos$Transplant_ID[i], 1, j), sep="/")
    }
  }
  
  xenos <- merge(xenos, annot, all.x = TRUE)
  xeno_hist <- as.Node(xenos)
  
  
  # Set style of entire tree
  # https://www.rdocumentation.org/packages/data.tree/versions/0.7.5/topics/plot.Node
  # https://graphviz.gitlab.io/_pages/doc/info/attrs.html
  SetNodeStyle(xeno_hist, style = "filled,rounded", fontname = "helvetica", shape = "box", fillcolor = "White", fontcolor = "Black")
  
  for (i in 1:nrow(xenos)) {
    # if (!is.na(xenos$DLP.status[i])) {
    node_str <- paste0("xeno_hist$'", paste(tail(unlist(strsplit(xenos$pathString[i], "/")), -1), collapse = "'$'"), "'")
    node <- eval(parse(text = node_str))
    
    # Magically colour happens here, uses the stats to index into col to select fill color
    status <- xenos$DLP.status[i]
    # pbal_col <- ifelse(is.na(xenos$PBAL.Submission[i]), "Black", "Orange")
    inProj_col <- ifelse(is.na(xenos$Transplant_ID[i]), "Black", "Green")
    # col[status] 
    mixture_stat_col<- ifelse(xenos$Passage[i] =="", "Red", "Grey")
    tenx_shape <- ifelse(is.na(xenos$PBAL.Submission[i]), "box", "ellipse")
    
    SetNodeStyle(node, fillcolor = mixture_stat_col, fontcolor = inProj_col, shape = tenx_shape, inherit = FALSE)  
    # }
  }
  
  export_graph(ToDiagrammeRGraph(xeno_hist), paste(resultdir, series, ".pdf", sep=""))
  write.csv(temp, paste(resultdir, series, ".csv", sep=""))
  
  
  temp_tree_level = xeno_hist$Get('level', traversal = "pre-order")
  temp_tree_level = data.frame(level = temp_tree_level, node = names(temp_tree_level) )
  
  temp_df = xeno_hist$Get(function(x) c(position = x$position, level = x$level), simplify = "regular")
  
  levels = unique(temp_tree_level$level)
  
  # from_to = subset(temp_tree, level == 1 )
  from_to = data.frame(NULL)
  
  
  count = 1
  for(i in head(levels, -1)){
    temp_1 = subset(temp_tree_level, level == i )
    temp_2 = subset(temp_tree_level, level == i +1)
    nrow(temp_2)
    temp = data.frame(NULL)
    for(j in c(1:nrow(temp_1))){
      temptemp = temp_2[grepl(temp_1$node[j], temp_2$node),]
      temp= rbind(temp, cbind(temp_1[rep(j, nrow(temptemp)),], temptemp))
    }
    from_to= rbind(from_to,temp)
    # count = count + nrow(temp_2)
  }
  
  names(from_to) = c("From_Level", "From", "To_Level", "To")
  
  
  from_to_delete = from_to[!xenos$Transplant_ID %in% xenos$Transplant_ID, ]
  from_to_reserve = from_to[(from_to$From %in% xenos$Transplant_ID) |(from_to$To %in% xenos$Transplant_ID), ]
  from_to_reserve = merge(from_to_reserve, xenos[, c("Transplant_ID", "Xeno")], 
                          by.x = "From", by.y = "Transplant_ID", all.x = TRUE)
  from_to_reserve = merge(from_to_reserve, xenos[, c("Transplant_ID", "Xeno")], 
                          by.x = "To", by.y = "Transplant_ID", all.x = TRUE)
  
  from_to_reserve = from_to_reserve[, c("From_Level", "From", "Xeno.x", 
                                        "To_Level", "To", "Xeno.y")]
  
  names(from_to_reserve) = c("From_Level", "From", "From_Transplant_ID", 
                             "To_Level", "To", "To_Transplant_ID")
  write.csv(from_to_reserve, paste(resultdir, "from_to.csv", sep=""))
}