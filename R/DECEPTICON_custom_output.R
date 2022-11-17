#'Integrating analysis strategies and outputting DECEPTICON results
#'
#'@param path The storage path for the results of the analysis strategies
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON__custom_output <- function(path, signature.matrix){
  setwd(path)
  filename <- list.files(pattern = "_sig_")
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
    cell_types = intersect(colnames(sig_1),colnames(sig_2))
    for (i in 3:length(signature.matrix)) {
      cell_types = intersect(cell_types, colnames(get(paste0('sig_','',i,sep=''))))
      if(i>=length(signature.matrix)){
        break;
      }
    }
  } else {
    cell_types = colnames(sig_1)
  }
  
  #choose
  for (j in 1:length(cell_types)) {
    for (i in 1:length(data_list)) {
      if(colnames(data_list[[i]])==cell_types[j]){
       assign(sub('....$','',filename[i]),data_list[[i]][,which(colnames(data_list[[i]])==cell_types[j])])
      }
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
    for (i in 1:length(row)) {
      for (k in 1:length(col)) {
        if(substrRight(col[i],4) == substrRight(row[k],4))
          cell_cor[i,k] = 0
      }
    }
    cell_id= optimal_id(cell_cor,2)
    cell_res= rownames(cell_cor)[cell_id[1]]
    for (i in 2:length(cell_id)) {
      cell_res= c(cell_res,rownames(cell_cor)[cell_id[i]])
    }
    library(vegan)
    library(psych)
    n = 6*length(signature.matrix)+1
    for (i in n:n+length(signature.matrix)) {
      data_list[[i]] = decostand(data_list[[i]],MARGIN = 1,"total")
    }
    for (i in 1:length(signature.matrix)) {
      data_list[[i]] = decostand(data_list[[i]],MARGIN = 1,"total")
    }
   assign(paste0('res_','',j,sep=''), (df[,colnames(df) == cell_res[1]]*0.25+
                                            df[,colnames(df) == cell_res[2]]*0.25+
                                            df[,colnames(df) == cell_res[3]]*0.25+
                                            df[,colnames(df) == cell_res[4]]*0.25)) 
    
  }
  DECEPTICON_res = data.frame(get(paste0('res_','',1,sep=''))) 
  for (i in 2:length(cell_types)) {
    DECEPTICON_res = cbind(DECEPTICON_res ,get(paste0('res_','',i,sep='')))
  }
  colnames(DECEPTICON_res) = cell_types
  res = DECEPTICON_res
  write.table(res, './DECPTICON.txt', sep='\t', row.names= T, col.names= T, quote=F)
}
    
  
  

  
 
  
  