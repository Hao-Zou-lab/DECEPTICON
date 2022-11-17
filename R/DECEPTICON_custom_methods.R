#'Deconvolute using CIBERSORT and custom signature matrix
#'
#'@param mixture_file a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_ciber_custom <- function(mixture_file, signature.matrix){
  source('./CIBERSORT.R')
  ciber_sig_1 = CIBERSORT(signature.matrix[1],mixture_file, perm = 0,QN = TRUE, absolute = F)
  write.table(ciber_sig_1, './res/ciber_sig_1.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(length(signature.matrix) > 1){
    for (i in 2:length(signature.matrix)) {
      assign(paste0('ciber_sig_','',i,sep=''),CIBERSORT(signature.matrix[i],mixture_file, perm = 0,QN = TRUE, absolute = F))
      write.table(get(paste0('ciber_sig_','',i,sep='')), paste("./res/ciber_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
    }
  }
}

#'Deconvolute using CIBERSORT-abs and custom signature matrix
#'
#'@param mixture_file a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_ciber_abs_custom <- function(mixture_file, signature.matrix){
  source('./CIBERSORT.R')
  ciber_abs_sig_1 = CIBERSORT(signature.matrix[1],mixture_file, perm = 0,QN = TRUE, absolute = T)
  write.table(ciber_abs_sig_1, './res/ciber_abs_sig_1.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(length(signature.matrix) > 1){
    for (i in 2:length(signature.matrix)) {
      assign(paste0('ciber_abs_sig_','',i,sep=''),CIBERSORT(signature.matrix[i],mixture_file, perm = 0,QN = TRUE, absolute = T))
      write.table(get(paste0('ciber_abs_sig_','',i,sep='')), paste("./res/ciber_abs_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
    }
  }
}

#'Deconvolute using EPIC and custom signature matrix
#'
#'@param mix.mat a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_epic_custom <- function(mix.mat, signature.matrix){
  library(EPIC)
  mix.mat = read.table(mix.mat, header = T,sep = '\t',row.names = 1)
  sig_1 = read.table(signature.matrix[1],sep = "\t",header = T,row.names = 1)
  epic_sig_1 = EPIC(mix.mat, reference =  list(refProfiles=sig_1,sigGenes = rownames(sig_1)))
  write.table(epic_sig_1[['cellFractions']], './res/EPIC_sig_1.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(length(signature.matrix) > 1){
    for (i in 2:length(signature.matrix)) {
      assign(paste0('sig_','',i,sep=''),read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1))
      assign(paste0('epic_sig_','',i,sep=''),EPIC(mix.mat, reference =  list(refProfiles=get(paste0('sig_','',i,sep='')),sigGenes = rownames(get(paste0('sig_','',i,sep=''))))))
      write.table(get(paste0('epic_sig_','',i,sep=''))[['cellFractions']], paste("./res/EPIC_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
    }
  }
}

#'Deconvolute using DeconRNAseq and custom signature matrix
#'
#'@param dataset a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_decon_custom <- function(dataset, signature.matrix){
  library(DeconRNASeq)
  dataset = read.table(dataset, header = T,sep = '\t',row.names = 1)
  sig_1 = read.table(signature.matrix[1],sep = "\t",header = T,row.names = 1)
  decon_sig_1 = DeconRNASeq(dataset,sig_1,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_sig_1$out.all, './res/Decon_sig_1.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(length(signature.matrix) > 1){
    for (i in 2:length(signature.matrix)) {
      assign(paste0('sig_','',i,sep=''),read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1))
      assign(paste0('decon_sig_','',i,sep=''),DeconRNASeq(dataset,get(paste0('sig_','',i,sep='')),proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE))
      write.table(get(paste0('decon_sig_','',i,sep=''))$out.all, paste("./res/Decon_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
    }
  }
}


#'Deconvolute using MCPcounter and custom signature matrix
#'
#'@param data a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_mcp_custom <- function(data,signature.matrix){
  library(MCPcounter)
  data = read.table(data, header = T,sep = '\t',row.names = 1)
  sig_1 = read.table(signature.matrix[1],sep = "\t",header = T,row.names = 1)
  for (i in 1:length(sig_1)) {
    sig_1 = sig_1[order(-(sig_1[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(sig_1)[1:10])
  }
  HUGO_symbols = cell_1
  for (i in 2:length(sig_1)) {
    HUGO_symbols = c(HUGO_symbols,get(paste0('cell_','',i,sep='')))
  }
  sig_1 = data.frame(HUGO_symbols = HUGO_symbols,
                    Cell_population = rep(colnames(sig_1),each = 10))
  colnames(sig_1) =  c("HUGO symbols", "Cell population")
  mcp_sig_1 = MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=sig_1)
  write.table(t(as.matrix(mcp_sig_1)), './res/MCP_sig_1.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(length(signature.matrix) > 1){
    for (i in 2:length(signature.matrix)) {
      assign(paste0('sig_','',i,sep=''),read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1))
    }
    for (i in 2:length(signature.matrix)) {
      for (j in 1:length(get(paste0('sig_','',i,sep='')))) {
        assign(paste0('sig_','',i,sep=''), get(paste0('sig_','',i,sep=''))[order(-(get(paste0('sig_','',i,sep=''))[,j])),])
        assign(paste0('cell_','',j,sep=''), rownames(get(paste0('sig_','',i,sep='')))[1:10])
      }
      HUGO_symbols = cell_1
      for (n in 2:length(get(paste0('sig_','',i,sep='')))) {
        HUGO_symbols = c(HUGO_symbols,get(paste0('cell_','',n,sep='')))
      }
      assign(paste0('sig_','',i,sep=''), data.frame(HUGO_symbols = HUGO_symbols,
                         Cell_population = rep(colnames(get(paste0('sig_','',i,sep=''))),each = 10)))
      a = get(paste0('sig_','',i,sep=''))
      colnames(a) =  c("HUGO symbols", "Cell population")
      assign(paste0('MCP_sig_','',i,sep=''),MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=a))
      write.table(t(as.matrix(get(paste0('MCP_sig_','',i,sep='')))), paste("./res/MCP_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
    }
  }
}

#'Deconvolute using quanTIseq and custom signature matrix
#'
#'@param mix.mat a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_quan_custom <- function(mix.mat, signature.matrix){
  library(preprocessCore)
  source("./DOCKER_codes.R")
  mix.mat = read.table(mix.mat, header = T,sep = '\t',row.names = 1)
  for (i in 1:length(signature.matrix)) {
    sig.mat = read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1)
    mRNA = rep(1,ncol(sig.mat))
    assign(paste0('quan_sig_','',i,sep=''), quanTlseq(sig.mat,mix.mat,scaling = mRNA,method = "lsei"))
    write.table(get(paste0('quan_sig_','',i,sep='')), paste("./res/quan_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using Bseq-SC and custom signature matrix
#'
#'@param data a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_bseq_custom <- function(data, signature.matrix){
  library(bseqsc)
  for (i in 1:length(signature.matrix)) {
    sig.mat = read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1)
    fdata = rownames(sig.mat)
    pdata = cbind(cellname = colnames(sig.mat),subjects = paste("patient",rep(1,length(sig.mat))))
    eset = getESET(sig.mat,fdata = fdata,pdata = pdata)
    for (j in 1:length(sig.mat)) {
      sig.mat = sig.mat[order(-(sig.mat[,j])),]
      assign(paste0('cell_','',j,sep=''), rownames(sig.mat)[1:10])
    }
    marker = vector("list",length(sig.mat))
    for (k in 1:length(sig.mat)) {
      marker[[k]]=get(paste0('cell_','',k,sep=''))
    }
    names(marker) = colnames(sig.mat)
    B = bseqsc_basis(eset,marker,clusters = "cellname",samples = "subjects",ct.scale = T)
    assign(paste0('bseq_sig_','',i,sep=''), bseqsc_proportions(data,B,verbose = T))
    write.table(t(get(paste0('bseq_sig_','',i,sep=''))$coefficients), paste("./res/bseq_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using SCDC and custom signature matrix
#'
#'@param data a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_scdc_custom <- function(data, signature.matrix){
  for (i in 1:length(signature.matrix)) {
    sig.mat = read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1)
    fdata = rownames(sig.mat)
    pdata = cbind(cellname = colnames(sig.mat),subjects = paste("patient",rep(1,length(sig.mat))))
    eset = getESET(sig.mat,fdata = fdata,pdata = pdata)
    eset.qc = SCDC_qc_ONE(eset,ct.varname = "cellname",sample = "subjects",scsetname = "sig.mat",ct.sub = colnames(sig.mat),qcthreshold = 0.7)
    assign(paste0('scdc_sig_','',i,sep=''), SCDC_prop_ONE(bulk.eset = data ,sc.eset = eset.qc$sc.eset.qc,ct.varname = "cellname",sample = "subjects",
                                    weight.basis = F,ct.sub = colnames(sig.mat)))
    write.table(get(paste0('scdc_sig_','',i,sep=''))$prop.est.mvw, paste("./res/scdc_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using MuSic and custom signature matrix
#'
#'@param data a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_music_custom <- function(data, signature.matrix){
  library(MuSiC)
  library(Biobase)
  for (i in 1:length(signature.matrix)) {
    sig.mat = read.table(signature.matrix[i],sep = "\t",header = T,row.names = 1)
    sig_music = cbind(sig.mat,sig.mat,sig.mat,sig.mat,sig.mat)
    colnames(sig_music) = paste("cell",1:(length(sig.mat)*5))
    fdata = rownames(sig_music)
    pdata = cbind(sampleID = c(rep(1,times=length(sig.mat)),rep(2,times=length(sig.mat)),rep(3,times=length(sig.mat)),rep(4,times=length(sig.mat)),rep(5,times=length(sig.mat))),
                  subjectname = c(rep("num1",times=length(sig.mat)),rep("num2",times=length(sig.mat)),rep("num3",times=length(sig.mat)),rep("num4",times=length(sig.mat)),rep("num5",times=length(sig.mat))),
                  celltypeID = rep(1:length(sig.mat),times=5),
                  celltype = c(colnames(sig.mat),colnames(sig.mat),colnames(sig.mat),colnames(sig.mat),colnames(sig.mat)))
    eset = getESET(sig_music,fdata = fdata,pdata = pdata)
    assign(paste0('music_sig_','',i,sep=''), music_prop(bulk.eset = data,sc.eset = eset,clusters = "celltype",samples = "sampleID"))
    write.table(get(paste0('music_sig_','',i,sep=''))$Est.prop.weighted, paste("./res/music_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using Batman and custom signature matrix
#'
#'@param data a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
DECEPTICON_batman_custom <- function(signature.matrix){
  library(batman)
  for (i in 1:length(signature.matrix)) {
    assign(paste0('batman_sig_','',i,sep=''), batman(runBATMANDir = paste0("./batman/sig_", '', i, sep=""),showPlot = FALSE))
    write.table(t(get(paste0('batman_sig_','',i,sep=''))$beta), paste("./res/batman_sig_", i, ".txt", sep=""), sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'deconvolution for cell type interpretation by consensus
#'
#'@param bulk.samples a m*n matrix with m genes and n samples
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@param RUNpath Working path for storing files
#'@export
DECEPTICON_custom_methods <- function(signature.matrix, bulk.samples, RUNpath){
  bulk.samples = bulk.samples
  signature.matrix = signature.matrix
  RUNpath = RUNpath
  setwd(RUNpath)
  message(paste0("\n",">>> Running ", "CIBERSORT"))
  DECEPTICON_ciber_custom(mixture_file = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "CIBERSORT-abs"))
  DECEPTICON_ciber_abs_custom(mixture_file = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "EPIC"))
  DECEPTICON_epic_custom(mix.mat = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "DeconRNAseq"))
  DECEPTICON_decon_custom(dataset = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "MCPcounter"))
  DECEPTICON_mcp_custom(data = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "quanTIseq"))
  DECEPTICON_quan_custom(mix.mat = bulk.samples, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "Bseq-SC"))
  library(SCDC)
  data = read.table(bulk.samples, header = T,sep = '\t',row.names = 1)
  fdata = rownames(data)
  pdata = cbind(sampleID = c(rep(1:length(data))),subjectname = c(rep("num1",times=length(data))),celltypeID = rep(1:length(data)))
  data = getESET(data,fdata = fdata,pdata = pdata)
  DECEPTICON_bseq_custom(data, signature.matrix = signature.matrix)
  message(paste0("\n",">>> Running ", "SCDC"))
  DECEPTICON_scdc_custom(data, signature.matrix = signature.matrix)
  if(package.version("MuSiC") == c("0.2.0")){
  message(paste0("\n",">>> Running ", "MuSic"))
  DECEPTICON_music_custom(data, signature.matrix = signature.matrix)
  }
  message(paste0("\n",">>> Running ", "Batman, may take a long time to run"))
  DECEPTICON_batman_custom(signature.matrix = signature.matrix)
}
