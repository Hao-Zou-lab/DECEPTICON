#'Deconvolute using CIBERSORT(Robust enumeration of cell subsets from tissue expression profiles)
#'
#'@param mix_file a m*n matrix with m genes and n samples
#'@param light set to TRUE for light version of DECEPTICON
#'@export
DECEPTICON_ciber <- function(mixture_file, light){
  source('./CIBERSORT.R')
  LM22.file = './signature_matrix/LM22.txt'
  TRef.file = './signature_matrix/EPIC_TRef.txt'
  BRef.file = './signature_matrix/EPIC_BRef.txt'
  quan.file = './signature_matrix/TIL10_signature.txt'
  immu.file = './signature_matrix/immu_base.txt'
  ciber_ciber = CIBERSORT(LM22.file,mixture_file,perm = 0,QN = TRUE, absolute = F )
  write.table(ciber_ciber, './res/ciber_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  ciber_quan = CIBERSORT(quan.file,mixture_file,perm = 0,QN = TRUE, absolute = F )
  write.table(ciber_quan, './res/ciber_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  ciber_immu = CIBERSORT(immu.file,mixture_file,perm = 0,QN = TRUE, absolute = F )
  write.table(ciber_immu, './res/ciber_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(light == FALSE){
    ciber_TRef = CIBERSORT(TRef.file,mixture_file,perm = 0,QN = TRUE, absolute = F )
    write.table(ciber_TRef, './res/ciber_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
    ciber_BRef = CIBERSORT(BRef.file,mixture_file,perm = 0,QN = TRUE, absolute = F )
    write.table(ciber_BRef, './res/ciber_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using CIBERSORT-abs
#'
#'@param mix_file a m*n matrix with m genes and n samples
#'@param light set to TRUE for light version of DECEPTICON
#'@export
DECEPTICON_ciber_abs <- function(mixture_file, light){
  source('./CIBERSORT.R')
  LM22.file = './signature_matrix/LM22.txt'
  TRef.file = './signature_matrix/EPIC_TRef.txt'
  BRef.file = './signature_matrix/EPIC_BRef.txt'
  quan.file = './signature_matrix/TIL10_signature.txt'
  immu.file = './signature_matrix/immu_base.txt'
  ciber_ciber = CIBERSORT(LM22.file,mixture_file,perm = 0,QN = TRUE, absolute = T )
  write.table(ciber_ciber, './res/ciber_abs_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  ciber_quan = CIBERSORT(quan.file,mixture_file,perm = 0,QN = TRUE, absolute = T )
  write.table(ciber_quan, './res/ciber_abs_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  ciber_immu = CIBERSORT(immu.file,mixture_file,perm = 0,QN = TRUE, absolute = T )
  write.table(ciber_immu, './res/ciber_abs_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
  if(light == FALSE){
    ciber_TRef = CIBERSORT(TRef.file,mixture_file,perm = 0,QN = TRUE, absolute = T )
    write.table(ciber_TRef, './res/ciber_abs_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
    ciber_BRef = CIBERSORT(BRef.file,mixture_file,perm = 0,QN = TRUE, absolute = T )
    write.table(ciber_BRef, './res/ciber_abs_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  }
}

#'Deconvolute using EPIC(Simultaneous enumeration of cancer an immune cell types from bulk tumor gene expression data)
#'
#'@param mix.mat a m*n matrix with m genes and n samples
#'@export
DECEPTICON_epic <- function(mix.mat){
  library(EPIC)
  mix.mat = read.table(mix.mat, header = T,sep = '\t',row.names = 1)
  epic_BRef <- EPIC::EPIC(mix.mat, reference =list(refProfiles=BRef$refProfiles,sigGenes = BRef$sigGenes,refProfiles.var = BRef$refProfiles.var))
  write.table(epic_BRef[['cellFractions']], './res/EPIC_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  epic_TRef <- EPIC::EPIC(mix.mat, reference =list(refProfiles=TRef$refProfiles,sigGenes = TRef$sigGenes,refProfiles.var = TRef$refProfiles.var))
  write.table(epic_TRef[['cellFractions']], './res/EPIC_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  LM22 = read.table("./signature_matrix/LM22.txt",sep = "\t",header = T,row.names = 1)
  epic_ciber <- EPIC::EPIC(mix.mat, reference =  list(refProfiles=LM22,sigGenes = rownames(LM22)))
  write.table(epic_ciber[['cellFractions']], './res/EPIC_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  quan_base =  read.table("./signature_matrix/TIL10_signature.txt",sep = "\t",header = T,row.names = 1)
  epic_quan <- EPIC::EPIC(mix.mat, reference =  list(refProfiles=quan_base,sigGenes = rownames(quan_base)))
  write.table(epic_quan[['cellFractions']], './res/EPIC_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  immu_base = read.table("./signature_matrix/immu_base.txt",sep = "\t",header = T,row.names = 1)
  epic_immu <- EPIC::EPIC(mix.mat, reference =  list(refProfiles=immu_base,sigGenes = rownames(immu_base)))
  write.table(epic_immu[['cellFractions']], './res/EPIC_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using DeconRNAseq(DeconRNASeq: a statistical framework for deconvolution of heterogeneous tissue samples based on mRNA-Seq data)
#'
#'@param dataset a m*n matrix with m genes and n samples
#'@export
DECEPTICON_decon <- function(dataset){
  library(DeconRNASeq)
  dataset = read.table(dataset, header = T,sep = '\t',row.names = 1)
  LM22 = read.table("./signature_matrix/LM22.txt", header = T,sep = '\t',row.names = 1)
  decon_ciber = DeconRNASeq(dataset,LM22,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_ciber$out.all, './res/Decon_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  TRef = read.table("./signature_matrix/EPIC_TRef.txt", header = T,sep = '\t',row.names = 1)
  decon_TRef = DeconRNASeq(dataset,TRef,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_TRef$out.all, './res/Decon_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  BRef = read.table("./signature_matrix/EPIC_BRef.txt", header = T,sep = '\t',row.names = 1)
  decon_BRef = DeconRNASeq(dataset,BRef,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_BRef$out.all, './res/Decon_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  TIL10 = read.table("./signature_matrix/TIL10_signature.txt", header = T,sep = '\t',row.names = 1)
  decon_quan = DeconRNASeq(dataset,TIL10,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_quan$out.all, './res/Decon_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  immu = read.table("./signature_matrix/immu_base.txt",sep = "\t",header = T,row.names = 1)
  decon_immu = DeconRNASeq(dataset,immu,proportions = NULL,checksig = FALSE,known.prop = FALSE,use.scale  = TRUE,fig  =TRUE)
  write.table(decon_immu$out.all, './res/Decon_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using MCPcounter(Estimating the population abundance of tissue-infiltrating immune and stromal cell populations using gene expression)
#'
#'@param data a m*n matrix with m genes and n samples
#'@export
DECEPTICON_mcp <- function(data){
  data = read.table(data, header = T,sep = '\t',row.names = 1)
  ciber_gene = read.table("./signature_matrix/LM22.txt", header = T,sep = '\t',row.names = 1)
  for (i in 1:length(ciber_gene)) {
    ciber_gene = ciber_gene[order(-(ciber_gene[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(ciber_gene)[1:10])
  }
  ciber_gene = data.frame(HUGO_symbols = c(cell_1,cell_2,cell_3,cell_4,cell_5,cell_6,cell_7,cell_8,cell_9,cell_10,cell_11,
                                          cell_12,cell_13,cell_14,cell_15,cell_16,cell_17,cell_18,cell_19,cell_20,cell_21,cell_22),
                         Cell_population = rep(colnames(ciber_gene),each = 10))
  colnames(ciber_gene) =  c("HUGO symbols", "Cell population")
  mcp_ciber = MCPcounter::MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=ciber_gene)
  write.table(t(as.matrix(mcp_ciber)), './res/MCP_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)

  TRef_gene = read.table('./signature_matrix/EPIC_TRef.txt',header = T,sep = '\t',row.names = 1)
  for (i in 1:length(TRef_gene)) {
    TRef_gene = TRef_gene[order(-(TRef_gene[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(TRef_gene)[1:10])
  }
  TRef_gene = data.frame(HUGO_symbols = c(cell_1,cell_2,cell_3,cell_4,cell_5,cell_6,cell_7),
                          Cell_population = rep(colnames(TRef_gene),each = 10))
  colnames(TRef_gene) =  c("HUGO symbols", "Cell population")
  mcp_TRef = MCPcounter::MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=TRef_gene)
  write.table(t(as.matrix(mcp_TRef)), './res/MCP_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)

  BRef_gene = read.table('./signature_matrix/EPIC_BRef.txt',header = T,sep = '\t',row.names = 1)
  for (i in 1:length(BRef_gene)) {
    BRef_gene = BRef_gene[order(-(BRef_gene[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(BRef_gene)[1:10])
  }
  BRef_gene = data.frame(HUGO_symbols = c(cell_1,cell_2,cell_3,cell_4,cell_5,cell_6),
                         Cell_population = rep(colnames(BRef_gene),each = 10))
  colnames(BRef_gene) =  c("HUGO symbols", "Cell population")
  mcp_BRef = MCPcounter::MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=BRef_gene)
  write.table(t(as.matrix(mcp_BRef)), './res/MCP_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)

  quan_gene = read.table("./signature_matrix/TIL10_signature.txt", header = T,sep = '\t',row.names = 1)
  for (i in 1:length(quan_gene)) {
    quan_gene = quan_gene[order(-(quan_gene[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(quan_gene)[1:10])
  }
  quan_gene = data.frame(HUGO_symbols = c(cell_1,cell_2,cell_3,cell_4,cell_5,cell_6,cell_7,cell_8,cell_9,cell_10),
                         Cell_population = rep(colnames(quan_gene),each = 10))
  colnames(quan_gene) =  c("HUGO symbols", "Cell population")
  mcp_quan = MCPcounter::MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=quan_gene)
  write.table(t(as.matrix(mcp_quan)), './res/MCP_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)

  immu_gene = read.table("./signature_matrix/immu_base.txt",sep = "\t",header = T,row.names = 1)
  for (i in 1:length(immu_gene)) {
    immu_gene = immu_gene[order(-(immu_gene[,i])),]
    assign(paste0('cell_','',i,sep=''), rownames(immu_gene)[1:10])
  }
  immu_gene = data.frame(HUGO_symbols = c(cell_1,cell_2,cell_3,cell_4,cell_5,cell_6,cell_7,cell_8,cell_9,cell_10,cell_11,
                                          cell_12,cell_13,cell_14,cell_15,cell_16,cell_17,cell_18,cell_19,cell_20),
                         Cell_population = rep(colnames(immu_gene),each = 10))
  colnames(immu_gene) =  c("HUGO symbols", "Cell population")
  mcp_immu = MCPcounter::MCPcounter.estimate(data,featuresType="HUGO_symbols",genes=immu_gene)
  write.table(t(as.matrix(mcp_immu)), './res/MCP_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using quanTIseq(Molecular and pharmacological modulators of the tumor immune contexture revealed by deconvolution of RNA-seq data)
#'
#'@param mix.mat a m*n matrix with m genes and n samples
#'@export
DECEPTICON_quan <- function(mix.mat){
  library(preprocessCore)
  source("./DOCKER_codes.R")
  mix.mat = read.table(mix.mat, header = T,sep = '\t',row.names = 1)
  signature_matrix = c(LM22.file = './signature_matrix/LM22.txt',
                       TRef.file = './signature_matrix/EPIC_TRef.txt',
                       BRef.file = './signature_matrix/EPIC_BRef.txt',
                       immu.file = './signature_matrix/immu_base.txt')
  res_id = c("quan_ciber", "quan_TRef", "quan_BRef", "quan_immu")
  for (i in 1:length(signature_matrix)) {
    sig.mat = read.table(signature_matrix[i],sep = "\t",header = T,row.names = 1)
    mRNA = rep(1,ncol(sig.mat))
    assign(res_id[i], quanTlseq(sig.mat,mix.mat,scaling = mRNA,method = "lsei"))
  }
  sig.mat = read.table("./signature_matrix/TIL10_signature.txt", header = T,sep = '\t',row.names = 1)
  mRNA<-read.table("./TIL10_mRNA_scaling.txt",sep="\t",header=FALSE,stringsAsFactors=FALSE)
  colnames(mRNA)<-c("celltype", "scaling")
  mRNA<-as.vector(as.matrix(mRNA$scaling[match(colnames(sig.mat), mRNA$celltype)]))
  mix.mat<-fixMixture(mix.mat, arrays=F)
  lrmgenes<-as.vector(read.table("./TIL10_rmgenes.txt", header=FALSE, sep="\t")[,1])
  n1<-nrow(sig.mat)
  sig.mat<-sig.mat[!rownames(sig.mat) %in% lrmgenes,, drop=FALSE]
  n2<-nrow(sig.mat)
  abgenes<-as.vector(read.table("./TIL10_TCGA_aberrant_immune_genes.txt", header=FALSE, sep="\t")[,1])
  n1<-nrow(sig.mat)
  sig.mat<-sig.mat[!rownames(sig.mat) %in% abgenes,, drop=FALSE]
  n2<-nrow(sig.mat)
  ns<-nrow(sig.mat)
  us<-length(intersect(rownames(sig.mat), rownames(mix.mat)))
  perc<-round(us*100/ns,digits=2)
  quan_quan<-quanTlseq(sig.mat,mix.mat,scaling=mRNA,method="lsei")
  write.table(quan_ciber, './res/quan_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(quan_TRef, './res/quan_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(quan_BRef, './res/quan_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(quan_quan, './res/quan_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(quan_immu, './res/quan_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using Bseq-SC
#'
#'@param data a m*n matrix with m genes and n samples
#'@export
DECEPTICON_bseq <- function(data){
  source("./CIBERSORT.R")
  source("./bseqsc_proportions.R")
  signature_matrix = c(LM22.file = './signature_matrix/LM22.txt',
                       TRef.file = './signature_matrix/EPIC_TRef.txt',
                       BRef.file = './signature_matrix/EPIC_BRef.txt',
                       quan.file = './signature_matrix/TIL10_signature.txt',
                       immu.file = './signature_matrix/immu_base.txt')
  res_id = c("bseq_ciber", "bseq_TRef", "bseq_BRef", "bseq_quan","bseq_immu")
  for (i in 1:length(signature_matrix)) {
    sig.mat = read.table(signature_matrix[i],sep = "\t",header = T,row.names = 1)
    fdata = rownames(sig.mat)
    pdata = cbind(cellname = colnames(sig.mat),subjects = paste("patient",rep(1,length(sig.mat))))
    eset = SCDC::getESET(sig.mat,fdata = fdata,pdata = pdata)
    for (j in 1:length(sig.mat)) {
      sig.mat = sig.mat[order(-(sig.mat[,j])),]
      assign(paste0('cell_','',j,sep=''), rownames(sig.mat)[1:10])
    }
    marker = vector("list",length(sig.mat))
    for (k in 1:length(sig.mat)) {
    marker[[k]]=get(paste0('cell_','',k,sep=''))
    }
    names(marker) = colnames(sig.mat)
    B = bseqsc::bseqsc_basis(eset,marker,clusters = "cellname",samples = "subjects",ct.scale = T)
    assign(res_id[i], bseqsc_proportions(data,B,verbose = T))
  }
  write.table(t(bseq_ciber$coefficients), './res/bseq_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(t(bseq_TRef$coefficients), './res/bseq_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(t(bseq_BRef$coefficients), './res/bseq_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(t(bseq_quan$coefficients), './res/bseq_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(t(bseq_immu$coefficients), './res/bseq_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using SCDC
#'
#'@param data a m*n matrix with m genes and n samples
#'@export
DECEPTICON_scdc <- function(data){
  signature_matrix = c(LM22.file = './signature_matrix/LM22.txt',
                       TRef.file = './signature_matrix/EPIC_TRef.txt',
                       BRef.file = './signature_matrix/EPIC_BRef.txt',
                       quan.file = './signature_matrix/TIL10_signature.txt',
                       immu.file = './signature_matrix/immu_base.txt')
  res_id = c("scdc_ciber", "scdc_TRef", "scdc_BRef", "scdc_quan","scdc_immu")
  for (i in 1:length(signature_matrix)) {
    sig.mat = read.table(signature_matrix[i],sep = "\t",header = T,row.names = 1)
    fdata = rownames(sig.mat)
    pdata = cbind(cellname = colnames(sig.mat),subjects = paste("patient",rep(1,length(sig.mat))))
    eset =  SCDC::getESET(sig.mat,fdata = fdata,pdata = pdata)
    eset.qc = SCDC::SCDC_qc_ONE(eset,ct.varname = "cellname",sample = "subjects",scsetname = "sig.mat",ct.sub = colnames(sig.mat),qcthreshold = 0.7)
    assign(res_id[i], SCDC::SCDC_prop_ONE(bulk.eset = data ,sc.eset = eset.qc$sc.eset.qc,ct.varname = "cellname",sample = "subjects",
                                    weight.basis = F,ct.sub = colnames(sig.mat)))
  }
  write.table(scdc_ciber$prop.est.mvw, './res/scdc_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(scdc_TRef$prop.est.mvw, './res/scdc_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(scdc_BRef$prop.est.mvw, './res/scdc_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(scdc_quan$prop.est.mvw, './res/scdc_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(scdc_immu$prop.est.mvw, './res/scdc_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using MuSic
#'
#'@param data a m*n matrix with m genes and n samples
#'@export
DECEPTICON_music <- function(data){
  source("./music_prop.R")
  source("./music_basis.R")
  signature_matrix = c(LM22.file = './signature_matrix/LM22.txt',
                       TRef.file = './signature_matrix/EPIC_TRef.txt',
                       BRef.file = './signature_matrix/EPIC_BRef.txt',
                       quan.file = './signature_matrix/TIL10_signature.txt',
                       immu.file = './signature_matrix/immu_base.txt')
  res_id = c("music_ciber", "music_TRef", "music_BRef", "music_quan","music_immu")
  for (i in 1:length(signature_matrix)) {
    sig.mat = read.table(signature_matrix[i],sep = "\t",header = T,row.names = 1)
    sig_music = cbind(sig.mat,sig.mat,sig.mat,sig.mat,sig.mat)
    colnames(sig_music) = paste("cell",1:(length(sig.mat)*5))
    fdata = rownames(sig_music)
    pdata = cbind(sampleID = c(rep(1,times=length(sig.mat)),rep(2,times=length(sig.mat)),rep(3,times=length(sig.mat)),rep(4,times=length(sig.mat)),rep(5,times=length(sig.mat))),
                  subjectname = c(rep("num1",times=length(sig.mat)),rep("num2",times=length(sig.mat)),rep("num3",times=length(sig.mat)),rep("num4",times=length(sig.mat)),rep("num5",times=length(sig.mat))),
                  celltypeID = rep(1:length(sig.mat),times=5),
                  celltype = c(colnames(sig.mat),colnames(sig.mat),colnames(sig.mat),colnames(sig.mat),colnames(sig.mat)))
    eset = SCDC::getESET(sig_music,fdata = fdata,pdata = pdata)
    assign(res_id[i], music_prop(bulk.mtx = Biobase::exprs(data),sc.sce = eset,clusters = "celltype",samples = "sampleID"))
  }
  write.table(music_ciber$Est.prop.weighted, './res/music_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(music_TRef$Est.prop.weighted, './res/music_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(music_BRef$Est.prop.weighted, './res/music_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(music_quan$Est.prop.weighted, './res/music_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(music_immu$Est.prop.weighted, './res/music_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'Deconvolute using Batman
#'
#'@param data a m*n matrix with m genes and n samples
#'@export
DECEPTICON_batman <- function(){
  signature_matrix = c(LM22.file = './batman/ciber_base',
                       TRef.file = './batman/TRef_base',
                       BRef.file = './batman/BRef_base',
                       quan.file = './batman/quan_base',
                       immu.file = './batman/immu_base')
  res_id = c("batman_ciber", "batman_TRef", "batman_BRef", "batman_quan","batman_immu")
  for (i in 1:length(signature_matrix)) {
    assign(res_id[i], batman::batman(runBATMANDir = signature_matrix[i],showPlot = FALSE))
  }
  write.table(batman_ciber$beta, './res/batman_ciber.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(batman_TRef$beta, './res/batman_TRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(batman_BRef$beta, './res/batman_BRef.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(batman_quan$beta, './res/batman_quan.txt', sep='\t', row.names= T, col.names= T, quote=F)
  write.table(batman_immu$beta, './res/batman_immu.txt', sep='\t', row.names= T, col.names= T, quote=F)
}


#'deconvolution for cell type interpretation by consensus
#'
#'@param bulk.samples a m*n matrix with m genes and n samples
#'@param RUNpath Working path for storing files
#'@param light set to TRUE for light version of DECEPTICON
#'@export
DECEPTICON_methods <- function(bulk.samples, RUNpath, light){
  bulk.samples = bulk.samples
  RUNpath = RUNpath
  light = light
  setwd(RUNpath)
  message(paste0("\n",">>> Running ", "CIBERSORT"))
  DECEPTICON_ciber(mixture_file = bulk.samples, light = light)
  message(paste0("\n",">>> Running ", "CIBERSORT-abs"))
  DECEPTICON_ciber_abs(mixture_file = bulk.samples, light = light)
  message(paste0("\n",">>> Running ", "EPIC"))
  DECEPTICON_epic(mix.mat = bulk.samples)
  message(paste0("\n",">>> Running ", "DeconRNAseq"))
  DECEPTICON_decon(dataset = bulk.samples)
  message(paste0("\n",">>> Running ", "MCPcounter"))
  DECEPTICON_mcp(data = bulk.samples)
  message(paste0("\n",">>> Running ", "quanTIseq"))
  DECEPTICON_quan(mix.mat = bulk.samples)
  data = read.table(bulk.samples, header = T,sep = '\t',row.names = 1)
  fdata = rownames(data)
  pdata = cbind(sampleID = c(rep(1:length(data))),subjectname = c(rep("num1",times=length(data))),celltypeID = rep(1:length(data)))
  data = SCDC::getESET(data,fdata = fdata,pdata = pdata)
  message(paste0("\n",">>> Running ", "Bseq-SC"))
  DECEPTICON_bseq(data)
  message(paste0("\n",">>> Running ", "SCDC"))
  DECEPTICON_scdc(data)
  message(paste0("\n",">>> Running ", "MuSic"))
  DECEPTICON_music(data)
  if(light == FALSE){
    message(paste0("\n",">>> Running ", "Batman"))
    DECEPTICON_batman()
  }
}


