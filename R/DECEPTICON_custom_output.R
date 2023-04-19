#'Integrating analysis strategies and outputting DECEPTICON results
#'
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_custom_output <- function(signature.matrix = signature.matrix, light = FALSE, single.cell = FALSE){
  setwd("./res")
  filename <- list.files(pattern = "_sig_")
  if(single.cell == TRUE) {
    filename = c(filename,"bayes_bayes.txt")
  }
  for (i in 1:length(filename)) {
    var_name <- gsub('.txt', '', filename[i])
    #assign(var_name, read.table(filename[i],header = T, sep = "\t", row.names = 1))
    data_list<- lapply(filename, read.table,header = T, sep = "\t", row.names = 1)
  }
  
  
  #Determining the exported cell types
  for (i in 1:length(signature.matrix)) {
    assign(paste0('sig_','',i,sep=''), read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1))
  }
  if(length(signature.matrix) > 1){
    if(single.cell == FALSE){
      cell_types = BiocGenerics::intersect(colnames(sig_1),colnames(sig_2))
      for (i in 3:length(signature.matrix)) {
        cell_types = BiocGenerics::intersect(cell_types, colnames(get(paste0('sig_','',i,sep=''))))
        if(i>=length(signature.matrix)){
          break;
        }
      }
    } else {
      cell_types = BiocGenerics::intersect(colnames(sig_1),names(data_list[[length(signature.matrix)*10+1]]))
      for (i in 2:length(signature.matrix)) {
        cell_types = BiocGenerics::intersect(cell_types, colnames(get(paste0('sig_','',i,sep=''))))
        if(i>=length(signature.matrix)){
          break;
        }
      }
    }
  } else {
    cell_types = BiocGenerics::intersect(colnames(sig_1),names(data_list[[11]]))
  }

  #choose
  for (j in 1:length(cell_types)) {
    for (i in 1:length(data_list)) {
       assign(sub('....$','',filename[i]),data_list[[i]][,which(colnames(data_list[[i]])==cell_types[j])])
    }
    df = get(sub('....$','',filename[1]))
    df = data.frame(df)
    colnames(df) = sub('....$','',filename[1])
    for (i in 2:length(data_list)) {
          df = cbind(df,get(sub('....$','',filename[i])))
          colnames(df)[i] = sub('....$','',filename[i])
          if(i>=length(signature.matrix)){
            next;
        }
      i = i+length(signature.matrix)-1
    }
    assign(cell_types[j], df)
    assign(paste(cell_types[j], "_cor", sep=""), cor(get(cell_types[j])))
    cell_cor = get(paste(cell_types[j], "_cor", sep=""))
    cell_cor[is.na(cell_cor)]=0
    col = colnames(cell_cor)
    row = rownames(cell_cor)
    if(length(signature.matrix)>1){
      for (i in 1:length(row)) {
        for (k in 1:length(col)) {
          if(substrRight(col[i],4) == substrRight(row[k],4))
            cell_cor[i,k] = 0
        }
      }
    } else {
      diag(cell_cor) = 0
      cell_cor[4,3] = 0
      cell_cor[3,4] = 0
    }
    cell_id= optimal_id(cell_cor,2)
    cell_res= rownames(cell_cor)[cell_id[1]]
    for (i in 2:length(cell_id)) {
      cell_res= c(cell_res,rownames(cell_cor)[cell_id[i]])
    }
    n = 6*length(signature.matrix)+1
    if(length(signature.matrix)>1){
      for (i in n:n+length(signature.matrix)) {
        data_list[[i]] = vegan::decostand(data_list[[i]],MARGIN = 1,"total")
      }
      for (i in 1:length(signature.matrix)) {
        data_list[[i]] = vegan::decostand(data_list[[i]],MARGIN = 1,"total")
      }
      n = 2*length(signature.matrix)+1
      for (i in n:n+length(signature.matrix)) {
        data_list[[i]] = vegan::decostand(data_list[[i]],MARGIN = 1,"total")
      }
    } else {
      data_list[[1]] = vegan::decostand(data_list[[1]],MARGIN = 1,"total")
      n = 6*length(signature.matrix)+1
      data_list[[n]] = vegan::decostand(data_list[[n]],MARGIN = 1,"total")
    }
    df = get(sub('....$','',filename[1]))
    df = data.frame(df)
    colnames(df) = sub('....$','',filename[1])
    for (i in 2:length(data_list)) {
      df = cbind(df,get(sub('....$','',filename[i])))
      colnames(df)[i] = sub('....$','',filename[i])
      if(i>=length(signature.matrix)){
        next;
      }
      i = i+length(signature.matrix)-1
    }
   assign(paste0('res_','',j,sep=''), 
          (data_list[[which(colnames(df) == cell_res[1])]][names(data_list[[which(colnames(df) == cell_res[1])]])==cell_types[j]]*0.25+
           data_list[[which(colnames(df) == cell_res[2])]][names(data_list[[which(colnames(df) == cell_res[2])]])==cell_types[j]]*0.25+
           data_list[[which(colnames(df) == cell_res[3])]][names(data_list[[which(colnames(df) == cell_res[3])]])==cell_types[j]]*0.25+
           data_list[[which(colnames(df) == cell_res[4])]][names(data_list[[which(colnames(df) == cell_res[4])]])==cell_types[j]]*0.25))

  }
  DECEPTICON_res = data.frame(get(paste0('res_','',1,sep='')))
  if(length(cell_types)>1){
    for (i in 2:length(cell_types)) {
      DECEPTICON_res = cbind(DECEPTICON_res ,get(paste0('res_','',i,sep='')))
    }
    colnames(DECEPTICON_res) = cell_types
  }
  
  write.table(DECEPTICON_res, './DECEPTICON.txt', sep='\t', row.names= T, col.names= T, quote=F)
}







