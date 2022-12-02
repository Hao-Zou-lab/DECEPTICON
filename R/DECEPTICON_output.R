#'Determine whether the expression template of the analysis strategies is the same
#'
#'@param x The name of the analysis strategies
#'@param n The last n digits of the  name of the analysis strategies
#'@export
substrRight = function(x,n){
  substr(x,nchar(x)-n+1, nchar(x))
}

#'Select the analysis strategies with the highest correlation
#'
#'@param matrix a m*m correlation matrix with m analysis strategies
#'@param n The top n groups of analysis strategies with the highest correlation
#'@export
optimal_id<-function(matrix,n){

  res_t=NULL
  for (nc in 1:nrow(matrix)) {
    res_t<-c(res_t,as.numeric(matrix[nc,]))
  }
  idx_optimal=order(res_t,decreasing = T)[1:(2*n)]####top

  info_all<-NULL

  idx_optimal<-idx_optimal[c(seq(1,2*n,by=2))]

  for (i in idx_optimal) {

    ##optimal filter method
    optimal_id1<-ceiling(i/nrow(matrix))
    ##optimal DE method
    if(i%%nrow(matrix)==0){
      optimal_id2<-nrow(matrix)
    }else{
      optimal_id2<-i%%nrow(matrix)
    }
    ##information of the method
    info_temp<-c(optimal_id1,optimal_id2)
    info_all<-rbind(info_all,info_temp)
  }
  return(info_all)
}

#'Integrating analysis strategies and outputting DECEPTICON results
#'
#'@param path The storage path for the results of the analysis strategies
#'@export
DECEPTICON_output <- function(light = FALSE){
  #cibersort
  res1 = read.table('./res/ciber_ciber.txt',header = T, sep = "\t", row.names = 1)
  res4 = read.table('./res/ciber_quan.txt',header = T, sep = "\t", row.names = 1)
  res5 = read.table('./res/ciber_immu.txt',header = T, sep = "\t", row.names = 1)
  if(light == FALSE){
    res2 = read.table('./res/ciber_TRef.txt',header = T, sep = "\t", row.names = 1)
    res3 = read.table('./res/ciber_BRef.txt',header = T, sep = "\t", row.names = 1)
  }



  #epic
  res6 = read.table('./res/EPIC_ciber.txt',header = T, sep = "\t", row.names = 1)
  res7 = read.table('./res/EPIC_TRef.txt',header = T, sep = "\t", row.names = 1)
  res8 = read.table('./res/EPIC_BRef.txt',header = T, sep = "\t", row.names = 1)
  res9 = read.table('./res/EPIC_quan.txt',header = T, sep = "\t", row.names = 1)
  res10 = read.table('./res/EPIC_immu.txt',header = T, sep = "\t", row.names = 1)

  #mcp
  res11 = read.table('./res/MCP_ciber.txt',header = T, sep = "\t", row.names = 1)
  res12 = read.table('./res/MCP_TRef.txt',header = T, sep = "\t", row.names = 1)
  res13 = read.table('./res/MCP_BRef.txt',header = T, sep = "\t", row.names = 1)
  res14 = read.table('./res/MCP_quan.txt',header = T, sep = "\t", row.names = 1)
  res15 = read.table('./res/MCP_immu.txt',header = T, sep = "\t", row.names = 1)

  #quanTlseq
  res16 = read.table('./res/quan_ciber.txt',header = T, sep = "\t", row.names = 1)
  res17 = read.table('./res/quan_TRef.txt',header = T, sep = "\t", row.names = 1)
  res18 = read.table('./res/quan_BRef.txt',header = T, sep = "\t", row.names = 1)
  res19 = read.table('./res/quan_quan.txt',header = T, sep = "\t", row.names = 1)
  res20 = read.table('./res/quan_immu.txt',header = T, sep = "\t", row.names = 1)

  #deconRNAseq
  res21 = read.table('./res/Decon_ciber.txt',header = T, sep = "\t", row.names = 1)
  res22 = read.table('./res/Decon_TRef.txt',header = T, sep = "\t", row.names = 1)
  res23 = read.table('./res/Decon_BRef.txt',header = T, sep = "\t", row.names = 1)
  res24 = read.table('./res/Decon_quan.txt',header = T, sep = "\t", row.names = 1)
  res25 = read.table('./res/Decon_immu.txt',header = T, sep = "\t", row.names = 1)

  #cibersort_abs
  res26 = read.table('./res/ciber_abs_ciber.txt',header = T, sep = "\t", row.names = 1)
  res29 = read.table('./res/ciber_abs_quan.txt',header = T, sep = "\t", row.names = 1)
  res30 = read.table('./res/ciber_abs_immu.txt',header = T, sep = "\t", row.names = 1)
  if(light == FALSE){
    res27 = read.table('./res/ciber_abs_TRef.txt',header = T, sep = "\t", row.names = 1)
    res28 = read.table('./res/ciber_abs_BRef.txt',header = T, sep = "\t", row.names = 1)
  }



  #scdc
  res31 = read.table('./res/scdc_ciber.txt',header = T, sep = "\t", row.names = 1)
  res32 = read.table('./res/scdc_TRef.txt',header = T, sep = "\t", row.names = 1)
  res33 = read.table('./res/scdc_BRef.txt',header = T, sep = "\t", row.names = 1)
  res34 = read.table('./res/scdc_quan.txt',header = T, sep = "\t", row.names = 1)
  res35 = read.table('./res/scdc_immu.txt',header = T, sep = "\t", row.names = 1)

  #music
  res36 = read.table('./res/music_ciber.txt',header = T, sep = "\t", row.names = 1)
  res37 = read.table('./res/music_TRef.txt',header = T, sep = "\t", row.names = 1)
  res38 = read.table('./res/music_BRef.txt',header = T, sep = "\t", row.names = 1)
  res39 = read.table('./res/music_quan.txt',header = T, sep = "\t", row.names = 1)
  res40 = read.table('./res/music_immu.txt',header = T, sep = "\t", row.names = 1)

  #bseq-sc
  res41 = read.table('./res/bseq_ciber.txt',header = T, sep = "\t", row.names = 1)
  res42 = read.table('./res/bseq_TRef.txt',header = T, sep = "\t", row.names = 1)
  res43 = read.table('./res/bseq_BRef.txt',header = T, sep = "\t", row.names = 1)
  res44 = read.table('./res/bseq_quan.txt',header = T, sep = "\t", row.names = 1)
  res45 = read.table('./res/bseq_immu.txt',header = T, sep = "\t", row.names = 1)

  if(light == FALSE){
    #batman
    res46 = data.frame(t(read.table('./res/batman_ciber.txt',header = T, sep = "\t", row.names = 1)))
    res47 = data.frame(t(read.table('./res/batman_TRef.txt',header = T, sep = "\t", row.names = 1)))
    res48 = data.frame(t(read.table('./res/batman_BRef.txt',header = T, sep = "\t", row.names = 1)))
    res49 = data.frame(t(read.table('./res/batman_quan.txt',header = T, sep = "\t", row.names = 1)))
    res50 = data.frame(t(read.table('./res/batman_immu.txt',header = T, sep = "\t", row.names = 1)))
  }


  #choose
    #B cell
    if(light == FALSE){
      b_cell_data = data.frame(ciber_ciber=res1$B.cells.naive+res1$B.cells.memory+res1$Plasma.cells,ciber_TRef = res2$B.cells.naive,ciber_BRef = res3$Bcells,ciber_quan=res4$B.cells,ciber_immu=res5$naive_B_cell+res5$memory_B_cell+res5$plasma_cell,
                               EPIC_ciber=res6$B.cells.naive+res6$B.cells.memory+res6$Plasma.cells,EPIC_TRef = res7$Bcells,EPIC_BRef = res8$Bcells,EPIC_quan=res9$B.cells,EPIC_immu=res10$naive_B_cell+res10$memory_B_cell+res10$plasma_cell,
                               MCP_ciber=res11$B.cells.naive+res11$B.cells.memory+res11$Plasma.cells,MCP_TRef = res12$B.cells.naive,MCP_BRef = res13$Bcells,MCP_quan=res14$B.cells,MCP_immu=res15$naive_B_cell+res15$memory_B_cell+res15$plasma_cell,
                               quan_ciber=res16$B.cells.naive+res16$B.cells.memory+res16$Plasma.cells,quan_TRef = res17$B.cells.naive,quan_BRef = res18$Bcells,quan_quan=res19$B.cells,quan_immu=res20$naive_B_cell+res20$memory_B_cell+res20$plasma_cell,
                               Decon_ciber=res21$B.cells.naive+res21$B.cells.memory+res21$Plasma.cells,Decon_TRef = res22$B.cells.naive,Decon_BRef = res23$Bcells,Decon_quan=res24$B.cells,Decon_immu=res25$naive_B_cell+res25$memory_B_cell+res25$plasma_cell,
                               ciber_abs_ciber=res26$B.cells.naive+res26$B.cells.memory+res26$Plasma.cells,ciber_abs_TRef = res27$B.cells.naive,ciber_abs_BRef = res28$Bcells,ciber_abs_quan=res29$B.cells,ciber_abs_immu=res30$naive_B_cell+res30$memory_B_cell+res30$plasma_cell,
                               scdc_ciber=res31$B.cells.naive+res31$B.cells.memory+res31$Plasma.cells,scdc_TRef = res32$B.cells.naive,scdc_BRef = res33$Bcells,scdc_quan=res34$B.cells,scdc_immu=res35$naive_B_cell+res35$memory_B_cell+res35$plasma_cell,
                               music_ciber=res36$B.cells.naive+res36$B.cells.memory+res36$Plasma.cells,music_TRef = res37$B.cells.naive,music_BRef = res38$Bcells,music_quan=res39$B.cells,music_immu=res40$naive_B_cell+res40$memory_B_cell+res40$plasma_cell,
                               bseq_ciber=res41$B.cells.naive+res41$B.cells.memory+res41$Plasma.cells,bseq_TRef = res42$B.cells.naive,bseq_BRef = res43$Bcells,bseq_quan=res44$B.cells,bseq_immu=res45$naive_B_cell+res45$memory_B_cell+res45$plasma_cell,
                               batman_ciber=res46$B.cells.naive+res46$B.cells.memory,batman_TRef = res47$B.cells.naive,batman_BRef = res48$Bcells,batman_quan=res49$B.cells,batman_immu=res50$naive_B_cell+res50$memory_B_cell+res50$plasma_cell)

      cd4_data = data.frame(ciber_ciber=res1$T.cells.CD4.naive+res1$T.cells.CD4.memory.resting+res1$T.cells.CD4.memory.activated,ciber_TRef = res2$T.cells.CD4.naive,ciber_BRef = res3$CD4_Tcells,ciber_quan=res4$T.cells.CD4,ciber_immu=res5$CD4_positive_alpha_beta_T_cell,
                            EPIC_ciber=res6$T.cells.CD4.naive+res6$T.cells.CD4.memory.resting+res6$T.cells.CD4.memory.activated,EPIC_TRef = res7$CD4_Tcells,EPIC_BRef = res8$CD4_Tcells,EPIC_quan=res9$T.cells.CD4,EPIC_immu=res10$CD4_positive_alpha_beta_T_cell,
                            MCP_ciber=res11$T.cells.CD4.naive+res11$T.cells.CD4.memory.resting+res11$T.cells.CD4.memory.activated,MCP_TRef = res12$T.cells.CD4.naive,MCP_BRef = res13$CD4_Tcells,MCP_quan=res14$T.cells.CD4,MCP_immu=res15$CD4_positive_alpha_beta_T_cell,
                            quan_ciber=res16$T.cells.CD4.naive+res16$T.cells.CD4.memory.resting+res16$T.cells.CD4.memory.activated,quan_TRef = res17$T.cells.CD4.naive,quan_BRef = res18$CD4_Tcells,quan_quan=res19$T.cells.CD4,quan_immu=res20$CD4_positive_alpha_beta_T_cell,
                            Decon_ciber=res21$T.cells.CD4.naive+res21$T.cells.CD4.memory.resting+res21$T.cells.CD4.memory.activated,Decon_TRef = res22$T.cells.CD4.naive,Decon_BRef = res23$CD4_Tcells,Decon_quan=res24$T.cells.CD4,Decon_immu=res25$CD4_positive_alpha_beta_T_cell,
                            ciber_abs_ciber=res26$T.cells.CD4.naive+res26$T.cells.CD4.memory.resting+res26$T.cells.CD4.memory.activated,ciber_abs_TRef = res27$T.cells.CD4.naive,ciber_abs_BRef = res28$CD4_Tcells,ciber_abs_quan=res29$T.cells.CD4,ciber_abs_immu=res30$CD4_positive_alpha_beta_T_cell,
                            scdc_ciber=res31$T.cells.CD4.naive+res31$T.cells.CD4.memory.resting+res31$T.cells.CD4.memory.activated,scdc_TRef = res32$T.cells.CD4.naive,scdc_BRef = res33$CD4_Tcells,scdc_quan=res34$T.cells.CD4,scdc_immu=res35$CD4_positive_alpha_beta_T_cell,
                            music_ciber=res36$T.cells.CD4.naive+res36$T.cells.CD4.memory.resting+res36$T.cells.CD4.memory.activated,music_TRef = res37$T.cells.CD4.naive,music_BRef = res38$CD4_Tcells,music_quan=res39$T.cells.CD4,music_immu=res40$CD4_positive_alpha_beta_T_cell,
                            bseq_ciber=res41$T.cells.CD4.naive+res41$T.cells.CD4.memory.resting+res41$T.cells.CD4.memory.activated,bseq_TRef = res42$T.cells.CD4.naive,bseq_BRef = res43$CD4_Tcells,bseq_quan=res44$T.cells.CD4,bseq_immu=res45$CD4_positive_alpha_beta_T_cell,
                            batman_ciber=res46$T.cells.CD4.naive+res46$T.cells.CD4.memory.resting+res46$T.cells.CD4.memory.activated,batman_TRef = res47$T.cells.CD4.naive,batman_BRef = res48$CD4_Tcells,batman_quan=res49$T.cells.CD4,batman_immu=res50$CD4_positive_alpha_beta_T_cell)

      cd8_data = data.frame(ciber_ciber=res1$T.cells.CD8,ciber_TRef = res2$CD8_T.cells,ciber_BRef = res3$CD8_Tcells,ciber_quan=res4$T.cells.CD8,ciber_immu=res5$CD8_positive_alpha_beta_T_cell,
                            EPIC_ciber=res6$T.cells.CD8,EPIC_TRef = res7$CD8_Tcells,EPIC_BRef = res8$CD8_Tcells,EPIC_quan=res9$T.cells.CD8,EPIC_immu=res10$CD8_positive_alpha_beta_T_cell,
                            MCP_ciber=res11$T.cells.CD8,MCP_TRef = res12$CD8_T.cells,MCP_BRef = res13$CD8_Tcells,MCP_quan=res14$T.cells.CD8,MCP_immu=res15$CD8_positive_alpha_beta_T_cell,
                            quan_ciber=res16$T.cells.CD8,quan_TRef = res17$CD8_T.cells,quan_BRef = res18$CD8_Tcells,quan_quan=res19$T.cells.CD8,quan_immu=res20$CD8_positive_alpha_beta_T_cell,
                            Decon_ciber=res21$T.cells.CD8,Decon_TRef = res22$CD8_T.cells,Decon_BRef = res23$CD8_Tcells,Decon_quan=res24$T.cells.CD8,Decon_immu=res25$CD8_positive_alpha_beta_T_cell,
                            ciber_abs_ciber=res26$T.cells.CD8,ciber_abs_TRef = res27$CD8_T.cells,ciber_abs_BRef = res28$CD8_Tcells,ciber_abs_quan=res29$T.cells.CD8,ciber_abs_immu=res30$CD8_positive_alpha_beta_T_cell,
                            scdc_ciber=res31$T.cells.CD8,scdc_TRef = res32$CD8_T.cells,scdc_BRef = res33$CD8_Tcells,scdc_quan=res34$T.cells.CD8,scdc_immu=res35$CD8_positive_alpha_beta_T_cell,
                            music_ciber=res36$T.cells.CD8,music_TRef = res37$CD8_T.cells,music_BRef = res38$CD8_Tcells,music_quan=res39$T.cells.CD8,music_immu=res40$CD8_positive_alpha_beta_T_cell,
                            bseq_ciber=res41$T.cells.CD8,bseq_TRef = res42$CD8_T.cells,bseq_BRef = res43$CD8_Tcells,bseq_quan=res44$T.cells.CD8,bseq_immu=res45$CD8_positive_alpha_beta_T_cell,
                            batman_ciber=res46$T.cells.CD8,batman_TRef = res47$CD8_T.cells,batman_BRef = res48$CD8_Tcells,batman_quan=res49$T.cells.CD8,batman_immu=res50$CD8_positive_alpha_beta_T_cell)

      treg_data = data.frame(ciber_ciber=res1$T.cells.regulatory..Tregs.,ciber_TRef = 0,ciber_BRef = 0,ciber_quan=res4$Tregs,ciber_immu=0,
                             EPIC_ciber=res6$T.cells.regulatory..Tregs.,EPIC_TRef = 0,EPIC_BRef = 0,EPIC_quan=res9$Tregs,EPIC_immu=0,
                             MCP_ciber=res11$T.cells.regulatory..Tregs.,MCP_TRef = 0,MCP_BRef = 0,MCP_quan=res14$Tregs,MCP_immu=0,
                             quan_ciber=res16$T.cells.regulatory..Tregs.,quan_TRef = 0,quan_BRef = 0,quan_quan=res19$Tregs,quan_immu=0,
                             Decon_ciber=res21$T.cells.regulatory..Tregs.,Decon_TRef = 0,Decon_BRef = 0,Decon_quan=res24$Tregs,Decon_immu=0,
                             ciber_abs_ciber=res26$T.cells.regulatory..Tregs.,ciber_abs_TRef = 0,ciber_abs_BRef = 0,ciber_abs_quan=res29$Tregs,ciber_abs_immu=0,
                             scdc_ciber=res31$T.cells.regulatory..Tregs.,scdc_TRef = 0,scdc_BRef = 0,scdc_quan=res34$Tregs,scdc_immu=0,
                             music_ciber=res36$T.cells.regulatory..Tregs.,music_TRef = 0,music_BRef = 0,music_quan=res39$Tregs,music_immu=0,
                             bseq_ciber=res41$T.cells.regulatory..Tregs.,bseq_TRef = 0,bseq_BRef = 0,bseq_quan=res44$Tregs,bseq_immu=0,
                             batman_ciber=res46$T.cells.regulatory..Tregs.,batman_TRef = 0,batman_BRef = 0,batman_quan=res49$Tregs,batman_immu=0)

      nk_data = data.frame(ciber_ciber=res1$NK.cells.resting+res1$NK.cells.activated,ciber_TRef = res2$NK.cells.activated,ciber_BRef = res3$NKcells,ciber_quan=res4$NK.cells,ciber_immu=res5$CD56bright_natural_killer_cell+res5$CD56dim_natural_killer_cell,
                           EPIC_ciber=res6$NK.cells.resting+res6$NK.cells.activated,EPIC_TRef = res7$NKcells,EPIC_BRef = res8$NKcells,EPIC_quan=res9$NK.cells,EPIC_immu=res10$CD56bright_natural_killer_cell+res10$CD56dim_natural_killer_cell,
                           MCP_ciber=res11$NK.cells.resting+res11$NK.cells.activated,MCP_TRef = res12$NK.cells.activated,MCP_BRef = res13$NKcells,MCP_quan=res14$NK.cells,MCP_immu=res15$CD56bright_natural_killer_cell+res15$CD56dim_natural_killer_cell,
                           quan_ciber=res16$NK.cells.resting+res16$NK.cells.activated,quan_TRef = res17$NK.cells.activated,quan_BRef = res18$NKcells,quan_quan=res19$NK.cells,quan_immu=res20$CD56bright_natural_killer_cell+res20$CD56dim_natural_killer_cell,
                           Decon_ciber=res21$NK.cells.resting+res21$NK.cells.activated,Decon_TRef = res22$NK.cells.activated,Decon_BRef = res23$NKcells,Decon_quan=res24$NK.cells,Decon_immu=res25$CD56bright_natural_killer_cell+res25$CD56dim_natural_killer_cell,
                           ciber_abs_ciber=res26$NK.cells.resting+res26$NK.cells.activated,ciber_abs_TRef = res27$NK.cells.activated,ciber_abs_BRef = res28$NKcells,ciber_abs_quan=res29$NK.cells,ciber_abs_immu=res30$CD56bright_natural_killer_cell+res30$CD56dim_natural_killer_cell,
                           scdc_ciber=res31$NK.cells.resting+res31$NK.cells.activated,scdc_TRef = res32$NK.cells.activated,scdc_BRef = res33$NKcells,scdc_quan=res34$NK.cells,scdc_immu=res35$CD56bright_natural_killer_cell+res35$CD56dim_natural_killer_cell,
                           music_ciber=res36$NK.cells.resting+res36$NK.cells.activated,music_TRef = res37$NK.cells.activated,music_BRef = res38$NKcells,music_quan=res39$NK.cells,music_immu=res40$CD56bright_natural_killer_cell+res40$CD56dim_natural_killer_cell,
                           bseq_ciber=res41$NK.cells.resting+res41$NK.cells.activated,bseq_TRef = res42$NK.cells.activated,bseq_BRef = res43$NKcells,bseq_quan=res44$NK.cells,bseq_immu=res45$CD56bright_natural_killer_cell+res45$CD56dim_natural_killer_cell,
                           batman_ciber=res46$NK.cells.resting+res46$NK.cells.activated,batman_TRef = res47$NK.cells.activated,batman_BRef = res48$NKcells,batman_quan=res49$NK.cells,batman_immu=res50$CD56bright_natural_killer_cell+res50$CD56dim_natural_killer_cell)

      mm_data = data.frame(ciber_ciber=res1$Monocytes+res1$Macrophages.M0+res1$Macrophages.M1+res1$Macrophages.M2,ciber_TRef = res2$Macrophages.M0,ciber_BRef = res3$Monocytes,ciber_quan=res4$Macrophages.M1+res4$Macrophages.M2+res4$Monocytes,ciber_immu=res5$CD14_positive_monocyte+res5$CD16_positive_monocyte+res5$macrophage_m0+res5$macrophage_m1+res5$macrophage_m2,
                           EPIC_ciber=res6$Monocytes+res6$Macrophages.M0+res6$Macrophages.M1+res6$Macrophages.M2,EPIC_TRef = res7$Macrophages,EPIC_BRef = res8$Monocytes,EPIC_quan=res9$Macrophages.M1+res9$Macrophages.M2+res9$Monocytes,EPIC_immu=res10$CD14_positive_monocyte+res10$CD16_positive_monocyte+res10$macrophage_m0+res10$macrophage_m1+res10$macrophage_m2,
                           MCP_ciber=res11$Monocytes+res11$Macrophages.M0+res11$Macrophages.M1+res11$Macrophages.M2,MCP_TRef = res12$Macrophages.M0,MCP_BRef = res13$Monocytes,MCP_quan=res14$Macrophages.M1+res14$Macrophages.M2+res14$Monocytes,MCP_immu=res15$CD14_positive_monocyte+res15$CD16_positive_monocyte+res15$macrophage_m0+res15$macrophage_m1+res15$macrophage_m2,
                           quan_ciber=res16$Monocytes+res16$Macrophages.M0+res16$Macrophages.M1+res16$Macrophages.M2,quan_TRef = res17$Macrophages.M0,quan_BRef = res18$Monocytes,quan_quan=res19$Macrophages.M1+res19$Macrophages.M2+res19$Monocytes,quan_immu=res20$CD14_positive_monocyte+res20$CD16_positive_monocyte+res20$macrophage_m0+res20$macrophage_m1+res20$macrophage_m2,
                           Decon_ciber=res21$Monocytes+res21$Macrophages.M0+res21$Macrophages.M1+res21$Macrophages.M2,Decon_TRef = res22$Macrophages.M0,Decon_BRef = res23$Monocytes,Decon_quan=res24$Macrophages.M1+res24$Macrophages.M2+res24$Monocytes,Decon_immu=res25$CD14_positive_monocyte+res25$CD16_positive_monocyte+res25$macrophage_m0+res25$macrophage_m1+res25$macrophage_m2,
                           ciber_abs_ciber=res26$Monocytes+res26$Macrophages.M0+res26$Macrophages.M1+res26$Macrophages.M2,ciber_abs_TRef = res27$Macrophages.M0,ciber_abs_BRef = res28$Monocytes,ciber_abs_quan=res29$Macrophages.M1+res29$Macrophages.M2+res29$Monocytes,ciber_abs_immu=res30$CD14_positive_monocyte+res30$CD16_positive_monocyte+res30$macrophage_m0+res30$macrophage_m1+res30$macrophage_m2,
                           scdc_ciber=res31$Monocytes+res31$Macrophages.M0+res31$Macrophages.M1+res31$Macrophages.M2,scdc_TRef = res32$Macrophages.M0,scdc_BRef = res33$Monocytes,scdc_quan=res34$Macrophages.M1+res34$Macrophages.M2+res34$Monocytes,scdc_immu=res35$CD14_positive_monocyte+res35$CD16_positive_monocyte+res35$macrophage_m0+res35$macrophage_m1+res35$macrophage_m2,
                           music_ciber=res36$Monocytes+res36$Macrophages.M0+res36$Macrophages.M1+res36$Macrophages.M2,music_TRef = res37$Macrophages.M0,music_BRef = res38$Monocytes,music_quan=res39$Macrophages.M1+res39$Macrophages.M2+res39$Monocytes,music_immu=res40$CD14_positive_monocyte+res40$CD16_positive_monocyte+res40$macrophage_m0+res40$macrophage_m1+res40$macrophage_m2,
                           bseq_ciber=res41$Monocytes+res41$Macrophages.M0+res41$Macrophages.M1+res41$Macrophages.M2,bseq_TRef = res42$Macrophages.M0,bseq_BRef = res43$Monocytes,bseq_quan=res44$Macrophages.M1+res44$Macrophages.M2+res44$Monocytes,bseq_immu=res45$CD14_positive_monocyte+res45$CD16_positive_monocyte+res45$macrophage_m0+res45$macrophage_m1+res45$macrophage_m2,
                           batman_ciber=res46$Monocytes+res46$Macrophages.M0+res46$Macrophages.M1+res46$Macrophages.M2,batman_TRef = res47$Macrophages.M0,batman_BRef = res48$Monocytes,batman_quan=res49$Macrophages.M1+res49$Macrophages.M2+res49$Monocytes,batman_immu=res50$CD14_positive_monocyte+res50$CD16_positive_monocyte+res50$macrophage_m0+res50$macrophage_m1+res50$macrophage_m2)
      } else {
        b_cell_data = data.frame(ciber_ciber=res1$B.cells.naive+res1$B.cells.memory+res1$Plasma.cells,ciber_quan=res4$B.cells,ciber_immu=res5$naive_B_cell+res5$memory_B_cell+res5$plasma_cell,
                                 EPIC_ciber=res6$B.cells.naive+res6$B.cells.memory+res6$Plasma.cells,EPIC_TRef = res7$Bcells,EPIC_BRef = res8$Bcells,EPIC_quan=res9$B.cells,EPIC_immu=res10$naive_B_cell+res10$memory_B_cell+res10$plasma_cell,
                                 MCP_ciber=res11$B.cells.naive+res11$B.cells.memory+res11$Plasma.cells,MCP_TRef = res12$B.cells.naive,MCP_BRef = res13$Bcells,MCP_quan=res14$B.cells,MCP_immu=res15$naive_B_cell+res15$memory_B_cell+res15$plasma_cell,
                                 quan_ciber=res16$B.cells.naive+res16$B.cells.memory+res16$Plasma.cells,quan_TRef = res17$B.cells.naive,quan_BRef = res18$Bcells,quan_quan=res19$B.cells,quan_immu=res20$naive_B_cell+res20$memory_B_cell+res20$plasma_cell,
                                 Decon_ciber=res21$B.cells.naive+res21$B.cells.memory+res21$Plasma.cells,Decon_TRef = res22$B.cells.naive,Decon_BRef = res23$Bcells,Decon_quan=res24$B.cells,Decon_immu=res25$naive_B_cell+res25$memory_B_cell+res25$plasma_cell,
                                 ciber_abs_ciber=res26$B.cells.naive+res26$B.cells.memory+res26$Plasma.cells,ciber_abs_quan=res29$B.cells,ciber_abs_immu=res30$naive_B_cell+res30$memory_B_cell+res30$plasma_cell,
                                 scdc_ciber=res31$B.cells.naive+res31$B.cells.memory+res31$Plasma.cells,scdc_TRef = res32$B.cells.naive,scdc_BRef = res33$Bcells,scdc_quan=res34$B.cells,scdc_immu=res35$naive_B_cell+res35$memory_B_cell+res35$plasma_cell,
                                 music_ciber=res36$B.cells.naive+res36$B.cells.memory+res36$Plasma.cells,music_TRef = res37$B.cells.naive,music_BRef = res38$Bcells,music_quan=res39$B.cells,music_immu=res40$naive_B_cell+res40$memory_B_cell+res40$plasma_cell,
                                 bseq_ciber=res41$B.cells.naive+res41$B.cells.memory+res41$Plasma.cells,bseq_TRef = res42$B.cells.naive,bseq_BRef = res43$Bcells,bseq_quan=res44$B.cells,bseq_immu=res45$naive_B_cell+res45$memory_B_cell+res45$plasma_cell)


        cd4_data = data.frame(ciber_ciber=res1$T.cells.CD4.naive+res1$T.cells.CD4.memory.resting+res1$T.cells.CD4.memory.activated,ciber_quan=res4$T.cells.CD4,ciber_immu=res5$CD4_positive_alpha_beta_T_cell,
                              EPIC_ciber=res6$T.cells.CD4.naive+res6$T.cells.CD4.memory.resting+res6$T.cells.CD4.memory.activated,EPIC_TRef = res7$CD4_Tcells,EPIC_BRef = res8$CD4_Tcells,EPIC_quan=res9$T.cells.CD4,EPIC_immu=res10$CD4_positive_alpha_beta_T_cell,
                              MCP_ciber=res11$T.cells.CD4.naive+res11$T.cells.CD4.memory.resting+res11$T.cells.CD4.memory.activated,MCP_TRef = res12$T.cells.CD4.naive,MCP_BRef = res13$CD4_Tcells,MCP_quan=res14$T.cells.CD4,MCP_immu=res15$CD4_positive_alpha_beta_T_cell,
                              quan_ciber=res16$T.cells.CD4.naive+res16$T.cells.CD4.memory.resting+res16$T.cells.CD4.memory.activated,quan_TRef = res17$T.cells.CD4.naive,quan_BRef = res18$CD4_Tcells,quan_quan=res19$T.cells.CD4,quan_immu=res20$CD4_positive_alpha_beta_T_cell,
                              Decon_ciber=res21$T.cells.CD4.naive+res21$T.cells.CD4.memory.resting+res21$T.cells.CD4.memory.activated,Decon_TRef = res22$T.cells.CD4.naive,Decon_BRef = res23$CD4_Tcells,Decon_quan=res24$T.cells.CD4,Decon_immu=res25$CD4_positive_alpha_beta_T_cell,
                              ciber_abs_ciber=res26$T.cells.CD4.naive+res26$T.cells.CD4.memory.resting+res26$T.cells.CD4.memory.activated,ciber_abs_quan=res29$T.cells.CD4,ciber_abs_immu=res30$CD4_positive_alpha_beta_T_cell,
                              scdc_ciber=res31$T.cells.CD4.naive+res31$T.cells.CD4.memory.resting+res31$T.cells.CD4.memory.activated,scdc_TRef = res32$T.cells.CD4.naive,scdc_BRef = res33$CD4_Tcells,scdc_quan=res34$T.cells.CD4,scdc_immu=res35$CD4_positive_alpha_beta_T_cell,
                              music_ciber=res36$T.cells.CD4.naive+res36$T.cells.CD4.memory.resting+res36$T.cells.CD4.memory.activated,music_TRef = res37$T.cells.CD4.naive,music_BRef = res38$CD4_Tcells,music_quan=res39$T.cells.CD4,music_immu=res40$CD4_positive_alpha_beta_T_cell,
                              bseq_ciber=res41$T.cells.CD4.naive+res41$T.cells.CD4.memory.resting+res41$T.cells.CD4.memory.activated,bseq_TRef = res42$T.cells.CD4.naive,bseq_BRef = res43$CD4_Tcells,bseq_quan=res44$T.cells.CD4,bseq_immu=res45$CD4_positive_alpha_beta_T_cell)

        cd8_data = data.frame(ciber_ciber=res1$T.cells.CD8,ciber_quan=res4$T.cells.CD8,ciber_immu=res5$CD8_positive_alpha_beta_T_cell,
                              EPIC_ciber=res6$T.cells.CD8,EPIC_TRef = res7$CD8_Tcells,EPIC_BRef = res8$CD8_Tcells,EPIC_quan=res9$T.cells.CD8,EPIC_immu=res10$CD8_positive_alpha_beta_T_cell,
                              MCP_ciber=res11$T.cells.CD8,MCP_TRef = res12$CD8_T.cells,MCP_BRef = res13$CD8_Tcells,MCP_quan=res14$T.cells.CD8,MCP_immu=res15$CD8_positive_alpha_beta_T_cell,
                              quan_ciber=res16$T.cells.CD8,quan_TRef = res17$CD8_T.cells,quan_BRef = res18$CD8_Tcells,quan_quan=res19$T.cells.CD8,quan_immu=res20$CD8_positive_alpha_beta_T_cell,
                              Decon_ciber=res21$T.cells.CD8,Decon_TRef = res22$CD8_T.cells,Decon_BRef = res23$CD8_Tcells,Decon_quan=res24$T.cells.CD8,Decon_immu=res25$CD8_positive_alpha_beta_T_cell,
                              ciber_abs_ciber=res26$T.cells.CD8,ciber_abs_quan=res29$T.cells.CD8,ciber_abs_immu=res30$CD8_positive_alpha_beta_T_cell,
                              scdc_ciber=res31$T.cells.CD8,scdc_TRef = res32$CD8_T.cells,scdc_BRef = res33$CD8_Tcells,scdc_quan=res34$T.cells.CD8,scdc_immu=res35$CD8_positive_alpha_beta_T_cell,
                              music_ciber=res36$T.cells.CD8,music_TRef = res37$CD8_T.cells,music_BRef = res38$CD8_Tcells,music_quan=res39$T.cells.CD8,music_immu=res40$CD8_positive_alpha_beta_T_cell,
                              bseq_ciber=res41$T.cells.CD8,bseq_TRef = res42$CD8_T.cells,bseq_BRef = res43$CD8_Tcells,bseq_quan=res44$T.cells.CD8,bseq_immu=res45$CD8_positive_alpha_beta_T_cell)

        treg_data = data.frame(ciber_ciber=res1$T.cells.regulatory..Tregs.,ciber_quan=res4$Tregs,ciber_immu=0,
                               EPIC_ciber=res6$T.cells.regulatory..Tregs.,EPIC_TRef = 0,EPIC_BRef = 0,EPIC_quan=res9$Tregs,EPIC_immu=0,
                               MCP_ciber=res11$T.cells.regulatory..Tregs.,MCP_TRef = 0,MCP_BRef = 0,MCP_quan=res14$Tregs,MCP_immu=0,
                               quan_ciber=res16$T.cells.regulatory..Tregs.,quan_TRef = 0,quan_BRef = 0,quan_quan=res19$Tregs,quan_immu=0,
                               Decon_ciber=res21$T.cells.regulatory..Tregs.,Decon_TRef = 0,Decon_BRef = 0,Decon_quan=res24$Tregs,Decon_immu=0,
                               ciber_abs_ciber=res26$T.cells.regulatory..Tregs.,ciber_abs_quan=res29$Tregs,ciber_abs_immu=0,
                               scdc_ciber=res31$T.cells.regulatory..Tregs.,scdc_TRef = 0,scdc_BRef = 0,scdc_quan=res34$Tregs,scdc_immu=0,
                               music_ciber=res36$T.cells.regulatory..Tregs.,music_TRef = 0,music_BRef = 0,music_quan=res39$Tregs,music_immu=0,
                               bseq_ciber=res41$T.cells.regulatory..Tregs.,bseq_TRef = 0,bseq_BRef = 0,bseq_quan=res44$Tregs,bseq_immu=0)

        nk_data = data.frame(ciber_ciber=res1$NK.cells.resting+res1$NK.cells.activated,ciber_quan=res4$NK.cells,ciber_immu=res5$CD56bright_natural_killer_cell+res5$CD56dim_natural_killer_cell,
                             EPIC_ciber=res6$NK.cells.resting+res6$NK.cells.activated,EPIC_TRef = res7$NKcells,EPIC_BRef = res8$NKcells,EPIC_quan=res9$NK.cells,EPIC_immu=res10$CD56bright_natural_killer_cell+res10$CD56dim_natural_killer_cell,
                             MCP_ciber=res11$NK.cells.resting+res11$NK.cells.activated,MCP_TRef = res12$NK.cells.activated,MCP_BRef = res13$NKcells,MCP_quan=res14$NK.cells,MCP_immu=res15$CD56bright_natural_killer_cell+res15$CD56dim_natural_killer_cell,
                             quan_ciber=res16$NK.cells.resting+res16$NK.cells.activated,quan_TRef = res17$NK.cells.activated,quan_BRef = res18$NKcells,quan_quan=res19$NK.cells,quan_immu=res20$CD56bright_natural_killer_cell+res20$CD56dim_natural_killer_cell,
                             Decon_ciber=res21$NK.cells.resting+res21$NK.cells.activated,Decon_TRef = res22$NK.cells.activated,Decon_BRef = res23$NKcells,Decon_quan=res24$NK.cells,Decon_immu=res25$CD56bright_natural_killer_cell+res25$CD56dim_natural_killer_cell,
                             ciber_abs_ciber=res26$NK.cells.resting+res26$NK.cells.activated,ciber_abs_quan=res29$NK.cells,ciber_abs_immu=res30$CD56bright_natural_killer_cell+res30$CD56dim_natural_killer_cell,
                             scdc_ciber=res31$NK.cells.resting+res31$NK.cells.activated,scdc_TRef = res32$NK.cells.activated,scdc_BRef = res33$NKcells,scdc_quan=res34$NK.cells,scdc_immu=res35$CD56bright_natural_killer_cell+res35$CD56dim_natural_killer_cell,
                             music_ciber=res36$NK.cells.resting+res36$NK.cells.activated,music_TRef = res37$NK.cells.activated,music_BRef = res38$NKcells,music_quan=res39$NK.cells,music_immu=res40$CD56bright_natural_killer_cell+res40$CD56dim_natural_killer_cell,
                             bseq_ciber=res41$NK.cells.resting+res41$NK.cells.activated,bseq_TRef = res42$NK.cells.activated,bseq_BRef = res43$NKcells,bseq_quan=res44$NK.cells,bseq_immu=res45$CD56bright_natural_killer_cell+res45$CD56dim_natural_killer_cell)

        mm_data = data.frame(ciber_ciber=res1$Monocytes+res1$Macrophages.M0+res1$Macrophages.M1+res1$Macrophages.M2,ciber_quan=res4$Macrophages.M1+res4$Macrophages.M2+res4$Monocytes,ciber_immu=res5$CD14_positive_monocyte+res5$CD16_positive_monocyte+res5$macrophage_m0+res5$macrophage_m1+res5$macrophage_m2,
                             EPIC_ciber=res6$Monocytes+res6$Macrophages.M0+res6$Macrophages.M1+res6$Macrophages.M2,EPIC_TRef = res7$Macrophages,EPIC_BRef = res8$Monocytes,EPIC_quan=res9$Macrophages.M1+res9$Macrophages.M2+res9$Monocytes,EPIC_immu=res10$CD14_positive_monocyte+res10$CD16_positive_monocyte+res10$macrophage_m0+res10$macrophage_m1+res10$macrophage_m2,
                             MCP_ciber=res11$Monocytes+res11$Macrophages.M0+res11$Macrophages.M1+res11$Macrophages.M2,MCP_TRef = res12$Macrophages.M0,MCP_BRef = res13$Monocytes,MCP_quan=res14$Macrophages.M1+res14$Macrophages.M2+res14$Monocytes,MCP_immu=res15$CD14_positive_monocyte+res15$CD16_positive_monocyte+res15$macrophage_m0+res15$macrophage_m1+res15$macrophage_m2,
                             quan_ciber=res16$Monocytes+res16$Macrophages.M0+res16$Macrophages.M1+res16$Macrophages.M2,quan_TRef = res17$Macrophages.M0,quan_BRef = res18$Monocytes,quan_quan=res19$Macrophages.M1+res19$Macrophages.M2+res19$Monocytes,quan_immu=res20$CD14_positive_monocyte+res20$CD16_positive_monocyte+res20$macrophage_m0+res20$macrophage_m1+res20$macrophage_m2,
                             Decon_ciber=res21$Monocytes+res21$Macrophages.M0+res21$Macrophages.M1+res21$Macrophages.M2,Decon_TRef = res22$Macrophages.M0,Decon_BRef = res23$Monocytes,Decon_quan=res24$Macrophages.M1+res24$Macrophages.M2+res24$Monocytes,Decon_immu=res25$CD14_positive_monocyte+res25$CD16_positive_monocyte+res25$macrophage_m0+res25$macrophage_m1+res25$macrophage_m2,
                             ciber_abs_ciber=res26$Monocytes+res26$Macrophages.M0+res26$Macrophages.M1+res26$Macrophages.M2,ciber_abs_quan=res29$Macrophages.M1+res29$Macrophages.M2+res29$Monocytes,ciber_abs_immu=res30$CD14_positive_monocyte+res30$CD16_positive_monocyte+res30$macrophage_m0+res30$macrophage_m1+res30$macrophage_m2,
                             scdc_ciber=res31$Monocytes+res31$Macrophages.M0+res31$Macrophages.M1+res31$Macrophages.M2,scdc_TRef = res32$Macrophages.M0,scdc_BRef = res33$Monocytes,scdc_quan=res34$Macrophages.M1+res34$Macrophages.M2+res34$Monocytes,scdc_immu=res35$CD14_positive_monocyte+res35$CD16_positive_monocyte+res35$macrophage_m0+res35$macrophage_m1+res35$macrophage_m2,
                             music_ciber=res36$Monocytes+res36$Macrophages.M0+res36$Macrophages.M1+res36$Macrophages.M2,music_TRef = res37$Macrophages.M0,music_BRef = res38$Monocytes,music_quan=res39$Macrophages.M1+res39$Macrophages.M2+res39$Monocytes,music_immu=res40$CD14_positive_monocyte+res40$CD16_positive_monocyte+res40$macrophage_m0+res40$macrophage_m1+res40$macrophage_m2,
                             bseq_ciber=res41$Monocytes+res41$Macrophages.M0+res41$Macrophages.M1+res41$Macrophages.M2,bseq_TRef = res42$Macrophages.M0,bseq_BRef = res43$Monocytes,bseq_quan=res44$Macrophages.M1+res44$Macrophages.M2+res44$Monocytes,bseq_immu=res45$CD14_positive_monocyte+res45$CD16_positive_monocyte+res45$macrophage_m0+res45$macrophage_m1+res45$macrophage_m2)
        }

  #b_cell
  b_cell_cor = cor(b_cell_data)
  b_cell_cor[is.na(b_cell_cor)]=0
  col = colnames(b_cell_cor)
  row = rownames(b_cell_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        b_cell_cor[i,j] = 0
    }
  }
  b_id= optimal_id(b_cell_cor,2)
  b_res= rownames(b_cell_cor)[b_id[1]]
  for (i in 2:length(b_id)) {
    b_res= c(b_res,rownames(b_cell_cor)[b_id[i]])
  }


  #cd4
  cd4_cor = cor(cd4_data)
  cd4_cor[is.na(cd4_cor)]=0
  col = colnames(cd4_cor)
  row = rownames(cd4_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        cd4_cor[i,j] = 0
    }
  }
  cd4_id = optimal_id(cd4_cor,2)
  cd4_res = rownames(cd4_cor)[cd4_id[1]]
  for (i in 2:length(cd4_id)) {
    cd4_res= c(cd4_res,rownames(cd4_cor)[cd4_id[i]])
  }

  #cd8
  cd8_cor = cor(cd8_data)
  cd8_cor[is.na(cd8_cor)]=0
  col = colnames(cd8_cor)
  row = rownames(cd8_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        cd8_cor[i,j] = 0
    }
  }
  cd8_id = optimal_id(cd8_cor,2)
  cd8_res = rownames(cd8_cor)[cd8_id[1]]
  for (i in 2:length(cd8_id)) {
    cd8_res = c(cd8_res,rownames(cd8_cor)[cd8_id[i]])
  }

  #treg
  treg_cor = cor(treg_data)
  treg_cor[is.na(treg_cor)]=0
  col = colnames(treg_cor)
  row = rownames(treg_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        treg_cor[i,j] = 0
    }
  }
  treg_id = optimal_id(treg_cor,2)
  treg_res = rownames(treg_cor)[treg_id[1]]
  for (i in 2:length(treg_id)) {
    treg_res = c(treg_res,rownames(treg_cor)[treg_id[i]])
  }

  #nk
  nk_cor = cor(nk_data)
  nk_cor[is.na(nk_cor)]=0
  col = colnames(nk_cor)
  row = rownames(nk_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        nk_cor[i,j] = 0
    }
  }
  nk_id = optimal_id(nk_cor,2)
  nk_res = rownames(nk_cor)[nk_id[1]]
  for (i in 2:length(nk_id)) {
    nk_res = c(nk_res,rownames(nk_cor)[nk_id[i]])
  }

  #mm
  mm_cor = cor(mm_data)
  mm_cor[is.na(mm_cor)]=0
  col = colnames(mm_cor)
  row = rownames(mm_cor)
  for (i in 1:length(row)) {
    for (j in 1:length(col)) {
      if(substrRight(col[i],4) == substrRight(row[j],4))
        mm_cor[i,j] = 0
    }
  }
  mm_id = optimal_id(mm_cor,2)
  mm_res = rownames(mm_cor)[mm_id[1]]
  for (i in 2:length(mm_id)) {
    mm_res = c(mm_res,rownames(mm_cor)[mm_id[i]])
  }

  #Assign different weights
  res11 = vegan::decostand(res11,MARGIN = 1,"total")
  res12 = vegan::decostand(res12,MARGIN = 1,"total")
  res13 = vegan::decostand(res13,MARGIN = 1,"total")
  res14 = vegan::decostand(res14,MARGIN = 1,"total")
  res15 = vegan::decostand(res15,MARGIN = 1,"total")
  if(light == FALSE){
    res46 = vegan::decostand(res46,MARGIN = 1,"total")
    res47 = vegan::decostand(res47,MARGIN = 1,"total")
    res48 = vegan::decostand(res48,MARGIN = 1,"total")
    res49 = vegan::decostand(res49,MARGIN = 1,"total")
    res50 = vegan::decostand(res50,MARGIN = 1,"total")
  }
  #B cell
  if(light == FALSE){
    b_cell_data = data.frame(ciber_ciber=res1$B.cells.naive+res1$B.cells.memory+res1$Plasma.cells,ciber_TRef = res2$B.cells.naive,ciber_BRef = res3$Bcells,ciber_quan=res4$B.cells,ciber_immu=res5$naive_B_cell+res5$memory_B_cell+res5$plasma_cell,
                             EPIC_ciber=res6$B.cells.naive+res6$B.cells.memory+res6$Plasma.cells,EPIC_TRef = res7$Bcells,EPIC_BRef = res8$Bcells,EPIC_quan=res9$B.cells,EPIC_immu=res10$naive_B_cell+res10$memory_B_cell+res10$plasma_cell,
                             MCP_ciber=res11$B.cells.naive+res11$B.cells.memory+res11$Plasma.cells,MCP_TRef = res12$B.cells.naive,MCP_BRef = res13$Bcells,MCP_quan=res14$B.cells,MCP_immu=res15$naive_B_cell+res15$memory_B_cell+res15$plasma_cell,
                             quan_ciber=res16$B.cells.naive+res16$B.cells.memory+res16$Plasma.cells,quan_TRef = res17$B.cells.naive,quan_BRef = res18$Bcells,quan_quan=res19$B.cells,quan_immu=res20$naive_B_cell+res20$memory_B_cell+res20$plasma_cell,
                             Decon_ciber=res21$B.cells.naive+res21$B.cells.memory+res21$Plasma.cells,Decon_TRef = res22$B.cells.naive,Decon_BRef = res23$Bcells,Decon_quan=res24$B.cells,Decon_immu=res25$naive_B_cell+res25$memory_B_cell+res25$plasma_cell,
                             ciber_abs_ciber=res26$B.cells.naive+res26$B.cells.memory+res26$Plasma.cells,ciber_abs_TRef = res27$B.cells.naive,ciber_abs_BRef = res28$Bcells,ciber_abs_quan=res29$B.cells,ciber_abs_immu=res30$naive_B_cell+res30$memory_B_cell+res30$plasma_cell,
                             scdc_ciber=res31$B.cells.naive+res31$B.cells.memory+res31$Plasma.cells,scdc_TRef = res32$B.cells.naive,scdc_BRef = res33$Bcells,scdc_quan=res34$B.cells,scdc_immu=res35$naive_B_cell+res35$memory_B_cell+res35$plasma_cell,
                             music_ciber=res36$B.cells.naive+res36$B.cells.memory+res36$Plasma.cells,music_TRef = res37$B.cells.naive,music_BRef = res38$Bcells,music_quan=res39$B.cells,music_immu=res40$naive_B_cell+res40$memory_B_cell+res40$plasma_cell,
                             bseq_ciber=res41$B.cells.naive+res41$B.cells.memory+res41$Plasma.cells,bseq_TRef = res42$B.cells.naive,bseq_BRef = res43$Bcells,bseq_quan=res44$B.cells,bseq_immu=res45$naive_B_cell+res45$memory_B_cell+res45$plasma_cell,
                             batman_ciber=res46$B.cells.naive+res46$B.cells.memory,batman_TRef = res47$B.cells.naive,batman_BRef = res48$Bcells,batman_quan=res49$B.cells,batman_immu=res50$naive_B_cell+res50$memory_B_cell+res50$plasma_cell)

    cd4_data = data.frame(ciber_ciber=res1$T.cells.CD4.naive+res1$T.cells.CD4.memory.resting+res1$T.cells.CD4.memory.activated,ciber_TRef = res2$T.cells.CD4.naive,ciber_BRef = res3$CD4_Tcells,ciber_quan=res4$T.cells.CD4,ciber_immu=res5$CD4_positive_alpha_beta_T_cell,
                          EPIC_ciber=res6$T.cells.CD4.naive+res6$T.cells.CD4.memory.resting+res6$T.cells.CD4.memory.activated,EPIC_TRef = res7$CD4_Tcells,EPIC_BRef = res8$CD4_Tcells,EPIC_quan=res9$T.cells.CD4,EPIC_immu=res10$CD4_positive_alpha_beta_T_cell,
                          MCP_ciber=res11$T.cells.CD4.naive+res11$T.cells.CD4.memory.resting+res11$T.cells.CD4.memory.activated,MCP_TRef = res12$T.cells.CD4.naive,MCP_BRef = res13$CD4_Tcells,MCP_quan=res14$T.cells.CD4,MCP_immu=res15$CD4_positive_alpha_beta_T_cell,
                          quan_ciber=res16$T.cells.CD4.naive+res16$T.cells.CD4.memory.resting+res16$T.cells.CD4.memory.activated,quan_TRef = res17$T.cells.CD4.naive,quan_BRef = res18$CD4_Tcells,quan_quan=res19$T.cells.CD4,quan_immu=res20$CD4_positive_alpha_beta_T_cell,
                          Decon_ciber=res21$T.cells.CD4.naive+res21$T.cells.CD4.memory.resting+res21$T.cells.CD4.memory.activated,Decon_TRef = res22$T.cells.CD4.naive,Decon_BRef = res23$CD4_Tcells,Decon_quan=res24$T.cells.CD4,Decon_immu=res25$CD4_positive_alpha_beta_T_cell,
                          ciber_abs_ciber=res26$T.cells.CD4.naive+res26$T.cells.CD4.memory.resting+res26$T.cells.CD4.memory.activated,ciber_abs_TRef = res27$T.cells.CD4.naive,ciber_abs_BRef = res28$CD4_Tcells,ciber_abs_quan=res29$T.cells.CD4,ciber_abs_immu=res30$CD4_positive_alpha_beta_T_cell,
                          scdc_ciber=res31$T.cells.CD4.naive+res31$T.cells.CD4.memory.resting+res31$T.cells.CD4.memory.activated,scdc_TRef = res32$T.cells.CD4.naive,scdc_BRef = res33$CD4_Tcells,scdc_quan=res34$T.cells.CD4,scdc_immu=res35$CD4_positive_alpha_beta_T_cell,
                          music_ciber=res36$T.cells.CD4.naive+res36$T.cells.CD4.memory.resting+res36$T.cells.CD4.memory.activated,music_TRef = res37$T.cells.CD4.naive,music_BRef = res38$CD4_Tcells,music_quan=res39$T.cells.CD4,music_immu=res40$CD4_positive_alpha_beta_T_cell,
                          bseq_ciber=res41$T.cells.CD4.naive+res41$T.cells.CD4.memory.resting+res41$T.cells.CD4.memory.activated,bseq_TRef = res42$T.cells.CD4.naive,bseq_BRef = res43$CD4_Tcells,bseq_quan=res44$T.cells.CD4,bseq_immu=res45$CD4_positive_alpha_beta_T_cell,
                          batman_ciber=res46$T.cells.CD4.naive+res46$T.cells.CD4.memory.resting+res46$T.cells.CD4.memory.activated,batman_TRef = res47$T.cells.CD4.naive,batman_BRef = res48$CD4_Tcells,batman_quan=res49$T.cells.CD4,batman_immu=res50$CD4_positive_alpha_beta_T_cell)

    cd8_data = data.frame(ciber_ciber=res1$T.cells.CD8,ciber_TRef = res2$CD8_T.cells,ciber_BRef = res3$CD8_Tcells,ciber_quan=res4$T.cells.CD8,ciber_immu=res5$CD8_positive_alpha_beta_T_cell,
                          EPIC_ciber=res6$T.cells.CD8,EPIC_TRef = res7$CD8_Tcells,EPIC_BRef = res8$CD8_Tcells,EPIC_quan=res9$T.cells.CD8,EPIC_immu=res10$CD8_positive_alpha_beta_T_cell,
                          MCP_ciber=res11$T.cells.CD8,MCP_TRef = res12$CD8_T.cells,MCP_BRef = res13$CD8_Tcells,MCP_quan=res14$T.cells.CD8,MCP_immu=res15$CD8_positive_alpha_beta_T_cell,
                          quan_ciber=res16$T.cells.CD8,quan_TRef = res17$CD8_T.cells,quan_BRef = res18$CD8_Tcells,quan_quan=res19$T.cells.CD8,quan_immu=res20$CD8_positive_alpha_beta_T_cell,
                          Decon_ciber=res21$T.cells.CD8,Decon_TRef = res22$CD8_T.cells,Decon_BRef = res23$CD8_Tcells,Decon_quan=res24$T.cells.CD8,Decon_immu=res25$CD8_positive_alpha_beta_T_cell,
                          ciber_abs_ciber=res26$T.cells.CD8,ciber_abs_TRef = res27$CD8_T.cells,ciber_abs_BRef = res28$CD8_Tcells,ciber_abs_quan=res29$T.cells.CD8,ciber_abs_immu=res30$CD8_positive_alpha_beta_T_cell,
                          scdc_ciber=res31$T.cells.CD8,scdc_TRef = res32$CD8_T.cells,scdc_BRef = res33$CD8_Tcells,scdc_quan=res34$T.cells.CD8,scdc_immu=res35$CD8_positive_alpha_beta_T_cell,
                          music_ciber=res36$T.cells.CD8,music_TRef = res37$CD8_T.cells,music_BRef = res38$CD8_Tcells,music_quan=res39$T.cells.CD8,music_immu=res40$CD8_positive_alpha_beta_T_cell,
                          bseq_ciber=res41$T.cells.CD8,bseq_TRef = res42$CD8_T.cells,bseq_BRef = res43$CD8_Tcells,bseq_quan=res44$T.cells.CD8,bseq_immu=res45$CD8_positive_alpha_beta_T_cell,
                          batman_ciber=res46$T.cells.CD8,batman_TRef = res47$CD8_T.cells,batman_BRef = res48$CD8_Tcells,batman_quan=res49$T.cells.CD8,batman_immu=res50$CD8_positive_alpha_beta_T_cell)

    treg_data = data.frame(ciber_ciber=res1$T.cells.regulatory..Tregs.,ciber_TRef = 0,ciber_BRef = 0,ciber_quan=res4$Tregs,ciber_immu=0,
                           EPIC_ciber=res6$T.cells.regulatory..Tregs.,EPIC_TRef = 0,EPIC_BRef = 0,EPIC_quan=res9$Tregs,EPIC_immu=0,
                           MCP_ciber=res11$T.cells.regulatory..Tregs.,MCP_TRef = 0,MCP_BRef = 0,MCP_quan=res14$Tregs,MCP_immu=0,
                           quan_ciber=res16$T.cells.regulatory..Tregs.,quan_TRef = 0,quan_BRef = 0,quan_quan=res19$Tregs,quan_immu=0,
                           Decon_ciber=res21$T.cells.regulatory..Tregs.,Decon_TRef = 0,Decon_BRef = 0,Decon_quan=res24$Tregs,Decon_immu=0,
                           ciber_abs_ciber=res26$T.cells.regulatory..Tregs.,ciber_abs_TRef = 0,ciber_abs_BRef = 0,ciber_abs_quan=res29$Tregs,ciber_abs_immu=0,
                           scdc_ciber=res31$T.cells.regulatory..Tregs.,scdc_TRef = 0,scdc_BRef = 0,scdc_quan=res34$Tregs,scdc_immu=0,
                           music_ciber=res36$T.cells.regulatory..Tregs.,music_TRef = 0,music_BRef = 0,music_quan=res39$Tregs,music_immu=0,
                           bseq_ciber=res41$T.cells.regulatory..Tregs.,bseq_TRef = 0,bseq_BRef = 0,bseq_quan=res44$Tregs,bseq_immu=0,
                           batman_ciber=res46$T.cells.regulatory..Tregs.,batman_TRef = 0,batman_BRef = 0,batman_quan=res49$Tregs,batman_immu=0)

    nk_data = data.frame(ciber_ciber=res1$NK.cells.resting+res1$NK.cells.activated,ciber_TRef = res2$NK.cells.activated,ciber_BRef = res3$NKcells,ciber_quan=res4$NK.cells,ciber_immu=res5$CD56bright_natural_killer_cell+res5$CD56dim_natural_killer_cell,
                         EPIC_ciber=res6$NK.cells.resting+res6$NK.cells.activated,EPIC_TRef = res7$NKcells,EPIC_BRef = res8$NKcells,EPIC_quan=res9$NK.cells,EPIC_immu=res10$CD56bright_natural_killer_cell+res10$CD56dim_natural_killer_cell,
                         MCP_ciber=res11$NK.cells.resting+res11$NK.cells.activated,MCP_TRef = res12$NK.cells.activated,MCP_BRef = res13$NKcells,MCP_quan=res14$NK.cells,MCP_immu=res15$CD56bright_natural_killer_cell+res15$CD56dim_natural_killer_cell,
                         quan_ciber=res16$NK.cells.resting+res16$NK.cells.activated,quan_TRef = res17$NK.cells.activated,quan_BRef = res18$NKcells,quan_quan=res19$NK.cells,quan_immu=res20$CD56bright_natural_killer_cell+res20$CD56dim_natural_killer_cell,
                         Decon_ciber=res21$NK.cells.resting+res21$NK.cells.activated,Decon_TRef = res22$NK.cells.activated,Decon_BRef = res23$NKcells,Decon_quan=res24$NK.cells,Decon_immu=res25$CD56bright_natural_killer_cell+res25$CD56dim_natural_killer_cell,
                         ciber_abs_ciber=res26$NK.cells.resting+res26$NK.cells.activated,ciber_abs_TRef = res27$NK.cells.activated,ciber_abs_BRef = res28$NKcells,ciber_abs_quan=res29$NK.cells,ciber_abs_immu=res30$CD56bright_natural_killer_cell+res30$CD56dim_natural_killer_cell,
                         scdc_ciber=res31$NK.cells.resting+res31$NK.cells.activated,scdc_TRef = res32$NK.cells.activated,scdc_BRef = res33$NKcells,scdc_quan=res34$NK.cells,scdc_immu=res35$CD56bright_natural_killer_cell+res35$CD56dim_natural_killer_cell,
                         music_ciber=res36$NK.cells.resting+res36$NK.cells.activated,music_TRef = res37$NK.cells.activated,music_BRef = res38$NKcells,music_quan=res39$NK.cells,music_immu=res40$CD56bright_natural_killer_cell+res40$CD56dim_natural_killer_cell,
                         bseq_ciber=res41$NK.cells.resting+res41$NK.cells.activated,bseq_TRef = res42$NK.cells.activated,bseq_BRef = res43$NKcells,bseq_quan=res44$NK.cells,bseq_immu=res45$CD56bright_natural_killer_cell+res45$CD56dim_natural_killer_cell,
                         batman_ciber=res46$NK.cells.resting+res46$NK.cells.activated,batman_TRef = res47$NK.cells.activated,batman_BRef = res48$NKcells,batman_quan=res49$NK.cells,batman_immu=res50$CD56bright_natural_killer_cell+res50$CD56dim_natural_killer_cell)

    mm_data = data.frame(ciber_ciber=res1$Monocytes+res1$Macrophages.M0+res1$Macrophages.M1+res1$Macrophages.M2,ciber_TRef = res2$Macrophages.M0,ciber_BRef = res3$Monocytes,ciber_quan=res4$Macrophages.M1+res4$Macrophages.M2+res4$Monocytes,ciber_immu=res5$CD14_positive_monocyte+res5$CD16_positive_monocyte+res5$macrophage_m0+res5$macrophage_m1+res5$macrophage_m2,
                         EPIC_ciber=res6$Monocytes+res6$Macrophages.M0+res6$Macrophages.M1+res6$Macrophages.M2,EPIC_TRef = res7$Macrophages,EPIC_BRef = res8$Monocytes,EPIC_quan=res9$Macrophages.M1+res9$Macrophages.M2+res9$Monocytes,EPIC_immu=res10$CD14_positive_monocyte+res10$CD16_positive_monocyte+res10$macrophage_m0+res10$macrophage_m1+res10$macrophage_m2,
                         MCP_ciber=res11$Monocytes+res11$Macrophages.M0+res11$Macrophages.M1+res11$Macrophages.M2,MCP_TRef = res12$Macrophages.M0,MCP_BRef = res13$Monocytes,MCP_quan=res14$Macrophages.M1+res14$Macrophages.M2+res14$Monocytes,MCP_immu=res15$CD14_positive_monocyte+res15$CD16_positive_monocyte+res15$macrophage_m0+res15$macrophage_m1+res15$macrophage_m2,
                         quan_ciber=res16$Monocytes+res16$Macrophages.M0+res16$Macrophages.M1+res16$Macrophages.M2,quan_TRef = res17$Macrophages.M0,quan_BRef = res18$Monocytes,quan_quan=res19$Macrophages.M1+res19$Macrophages.M2+res19$Monocytes,quan_immu=res20$CD14_positive_monocyte+res20$CD16_positive_monocyte+res20$macrophage_m0+res20$macrophage_m1+res20$macrophage_m2,
                         Decon_ciber=res21$Monocytes+res21$Macrophages.M0+res21$Macrophages.M1+res21$Macrophages.M2,Decon_TRef = res22$Macrophages.M0,Decon_BRef = res23$Monocytes,Decon_quan=res24$Macrophages.M1+res24$Macrophages.M2+res24$Monocytes,Decon_immu=res25$CD14_positive_monocyte+res25$CD16_positive_monocyte+res25$macrophage_m0+res25$macrophage_m1+res25$macrophage_m2,
                         ciber_abs_ciber=res26$Monocytes+res26$Macrophages.M0+res26$Macrophages.M1+res26$Macrophages.M2,ciber_abs_TRef = res27$Macrophages.M0,ciber_abs_BRef = res28$Monocytes,ciber_abs_quan=res29$Macrophages.M1+res29$Macrophages.M2+res29$Monocytes,ciber_abs_immu=res30$CD14_positive_monocyte+res30$CD16_positive_monocyte+res30$macrophage_m0+res30$macrophage_m1+res30$macrophage_m2,
                         scdc_ciber=res31$Monocytes+res31$Macrophages.M0+res31$Macrophages.M1+res31$Macrophages.M2,scdc_TRef = res32$Macrophages.M0,scdc_BRef = res33$Monocytes,scdc_quan=res34$Macrophages.M1+res34$Macrophages.M2+res34$Monocytes,scdc_immu=res35$CD14_positive_monocyte+res35$CD16_positive_monocyte+res35$macrophage_m0+res35$macrophage_m1+res35$macrophage_m2,
                         music_ciber=res36$Monocytes+res36$Macrophages.M0+res36$Macrophages.M1+res36$Macrophages.M2,music_TRef = res37$Macrophages.M0,music_BRef = res38$Monocytes,music_quan=res39$Macrophages.M1+res39$Macrophages.M2+res39$Monocytes,music_immu=res40$CD14_positive_monocyte+res40$CD16_positive_monocyte+res40$macrophage_m0+res40$macrophage_m1+res40$macrophage_m2,
                         bseq_ciber=res41$Monocytes+res41$Macrophages.M0+res41$Macrophages.M1+res41$Macrophages.M2,bseq_TRef = res42$Macrophages.M0,bseq_BRef = res43$Monocytes,bseq_quan=res44$Macrophages.M1+res44$Macrophages.M2+res44$Monocytes,bseq_immu=res45$CD14_positive_monocyte+res45$CD16_positive_monocyte+res45$macrophage_m0+res45$macrophage_m1+res45$macrophage_m2,
                         batman_ciber=res46$Monocytes+res46$Macrophages.M0+res46$Macrophages.M1+res46$Macrophages.M2,batman_TRef = res47$Macrophages.M0,batman_BRef = res48$Monocytes,batman_quan=res49$Macrophages.M1+res49$Macrophages.M2+res49$Monocytes,batman_immu=res50$CD14_positive_monocyte+res50$CD16_positive_monocyte+res50$macrophage_m0+res50$macrophage_m1+res50$macrophage_m2)
  } else {
    b_cell_data = data.frame(ciber_ciber=res1$B.cells.naive+res1$B.cells.memory+res1$Plasma.cells,ciber_quan=res4$B.cells,ciber_immu=res5$naive_B_cell+res5$memory_B_cell+res5$plasma_cell,
                             EPIC_ciber=res6$B.cells.naive+res6$B.cells.memory+res6$Plasma.cells,EPIC_TRef = res7$Bcells,EPIC_BRef = res8$Bcells,EPIC_quan=res9$B.cells,EPIC_immu=res10$naive_B_cell+res10$memory_B_cell+res10$plasma_cell,
                             MCP_ciber=res11$B.cells.naive+res11$B.cells.memory+res11$Plasma.cells,MCP_TRef = res12$B.cells.naive,MCP_BRef = res13$Bcells,MCP_quan=res14$B.cells,MCP_immu=res15$naive_B_cell+res15$memory_B_cell+res15$plasma_cell,
                             quan_ciber=res16$B.cells.naive+res16$B.cells.memory+res16$Plasma.cells,quan_TRef = res17$B.cells.naive,quan_BRef = res18$Bcells,quan_quan=res19$B.cells,quan_immu=res20$naive_B_cell+res20$memory_B_cell+res20$plasma_cell,
                             Decon_ciber=res21$B.cells.naive+res21$B.cells.memory+res21$Plasma.cells,Decon_TRef = res22$B.cells.naive,Decon_BRef = res23$Bcells,Decon_quan=res24$B.cells,Decon_immu=res25$naive_B_cell+res25$memory_B_cell+res25$plasma_cell,
                             ciber_abs_ciber=res26$B.cells.naive+res26$B.cells.memory+res26$Plasma.cells,ciber_abs_quan=res29$B.cells,ciber_abs_immu=res30$naive_B_cell+res30$memory_B_cell+res30$plasma_cell,
                             scdc_ciber=res31$B.cells.naive+res31$B.cells.memory+res31$Plasma.cells,scdc_TRef = res32$B.cells.naive,scdc_BRef = res33$Bcells,scdc_quan=res34$B.cells,scdc_immu=res35$naive_B_cell+res35$memory_B_cell+res35$plasma_cell,
                             music_ciber=res36$B.cells.naive+res36$B.cells.memory+res36$Plasma.cells,music_TRef = res37$B.cells.naive,music_BRef = res38$Bcells,music_quan=res39$B.cells,music_immu=res40$naive_B_cell+res40$memory_B_cell+res40$plasma_cell,
                             bseq_ciber=res41$B.cells.naive+res41$B.cells.memory+res41$Plasma.cells,bseq_TRef = res42$B.cells.naive,bseq_BRef = res43$Bcells,bseq_quan=res44$B.cells,bseq_immu=res45$naive_B_cell+res45$memory_B_cell+res45$plasma_cell)


    cd4_data = data.frame(ciber_ciber=res1$T.cells.CD4.naive+res1$T.cells.CD4.memory.resting+res1$T.cells.CD4.memory.activated,ciber_quan=res4$T.cells.CD4,ciber_immu=res5$CD4_positive_alpha_beta_T_cell,
                          EPIC_ciber=res6$T.cells.CD4.naive+res6$T.cells.CD4.memory.resting+res6$T.cells.CD4.memory.activated,EPIC_TRef = res7$CD4_Tcells,EPIC_BRef = res8$CD4_Tcells,EPIC_quan=res9$T.cells.CD4,EPIC_immu=res10$CD4_positive_alpha_beta_T_cell,
                          MCP_ciber=res11$T.cells.CD4.naive+res11$T.cells.CD4.memory.resting+res11$T.cells.CD4.memory.activated,MCP_TRef = res12$T.cells.CD4.naive,MCP_BRef = res13$CD4_Tcells,MCP_quan=res14$T.cells.CD4,MCP_immu=res15$CD4_positive_alpha_beta_T_cell,
                          quan_ciber=res16$T.cells.CD4.naive+res16$T.cells.CD4.memory.resting+res16$T.cells.CD4.memory.activated,quan_TRef = res17$T.cells.CD4.naive,quan_BRef = res18$CD4_Tcells,quan_quan=res19$T.cells.CD4,quan_immu=res20$CD4_positive_alpha_beta_T_cell,
                          Decon_ciber=res21$T.cells.CD4.naive+res21$T.cells.CD4.memory.resting+res21$T.cells.CD4.memory.activated,Decon_TRef = res22$T.cells.CD4.naive,Decon_BRef = res23$CD4_Tcells,Decon_quan=res24$T.cells.CD4,Decon_immu=res25$CD4_positive_alpha_beta_T_cell,
                          ciber_abs_ciber=res26$T.cells.CD4.naive+res26$T.cells.CD4.memory.resting+res26$T.cells.CD4.memory.activated,ciber_abs_quan=res29$T.cells.CD4,ciber_abs_immu=res30$CD4_positive_alpha_beta_T_cell,
                          scdc_ciber=res31$T.cells.CD4.naive+res31$T.cells.CD4.memory.resting+res31$T.cells.CD4.memory.activated,scdc_TRef = res32$T.cells.CD4.naive,scdc_BRef = res33$CD4_Tcells,scdc_quan=res34$T.cells.CD4,scdc_immu=res35$CD4_positive_alpha_beta_T_cell,
                          music_ciber=res36$T.cells.CD4.naive+res36$T.cells.CD4.memory.resting+res36$T.cells.CD4.memory.activated,music_TRef = res37$T.cells.CD4.naive,music_BRef = res38$CD4_Tcells,music_quan=res39$T.cells.CD4,music_immu=res40$CD4_positive_alpha_beta_T_cell,
                          bseq_ciber=res41$T.cells.CD4.naive+res41$T.cells.CD4.memory.resting+res41$T.cells.CD4.memory.activated,bseq_TRef = res42$T.cells.CD4.naive,bseq_BRef = res43$CD4_Tcells,bseq_quan=res44$T.cells.CD4,bseq_immu=res45$CD4_positive_alpha_beta_T_cell)

    cd8_data = data.frame(ciber_ciber=res1$T.cells.CD8,ciber_quan=res4$T.cells.CD8,ciber_immu=res5$CD8_positive_alpha_beta_T_cell,
                          EPIC_ciber=res6$T.cells.CD8,EPIC_TRef = res7$CD8_Tcells,EPIC_BRef = res8$CD8_Tcells,EPIC_quan=res9$T.cells.CD8,EPIC_immu=res10$CD8_positive_alpha_beta_T_cell,
                          MCP_ciber=res11$T.cells.CD8,MCP_TRef = res12$CD8_T.cells,MCP_BRef = res13$CD8_Tcells,MCP_quan=res14$T.cells.CD8,MCP_immu=res15$CD8_positive_alpha_beta_T_cell,
                          quan_ciber=res16$T.cells.CD8,quan_TRef = res17$CD8_T.cells,quan_BRef = res18$CD8_Tcells,quan_quan=res19$T.cells.CD8,quan_immu=res20$CD8_positive_alpha_beta_T_cell,
                          Decon_ciber=res21$T.cells.CD8,Decon_TRef = res22$CD8_T.cells,Decon_BRef = res23$CD8_Tcells,Decon_quan=res24$T.cells.CD8,Decon_immu=res25$CD8_positive_alpha_beta_T_cell,
                          ciber_abs_ciber=res26$T.cells.CD8,ciber_abs_quan=res29$T.cells.CD8,ciber_abs_immu=res30$CD8_positive_alpha_beta_T_cell,
                          scdc_ciber=res31$T.cells.CD8,scdc_TRef = res32$CD8_T.cells,scdc_BRef = res33$CD8_Tcells,scdc_quan=res34$T.cells.CD8,scdc_immu=res35$CD8_positive_alpha_beta_T_cell,
                          music_ciber=res36$T.cells.CD8,music_TRef = res37$CD8_T.cells,music_BRef = res38$CD8_Tcells,music_quan=res39$T.cells.CD8,music_immu=res40$CD8_positive_alpha_beta_T_cell,
                          bseq_ciber=res41$T.cells.CD8,bseq_TRef = res42$CD8_T.cells,bseq_BRef = res43$CD8_Tcells,bseq_quan=res44$T.cells.CD8,bseq_immu=res45$CD8_positive_alpha_beta_T_cell)

    treg_data = data.frame(ciber_ciber=res1$T.cells.regulatory..Tregs.,ciber_quan=res4$Tregs,ciber_immu=0,
                           EPIC_ciber=res6$T.cells.regulatory..Tregs.,EPIC_TRef = 0,EPIC_BRef = 0,EPIC_quan=res9$Tregs,EPIC_immu=0,
                           MCP_ciber=res11$T.cells.regulatory..Tregs.,MCP_TRef = 0,MCP_BRef = 0,MCP_quan=res14$Tregs,MCP_immu=0,
                           quan_ciber=res16$T.cells.regulatory..Tregs.,quan_TRef = 0,quan_BRef = 0,quan_quan=res19$Tregs,quan_immu=0,
                           Decon_ciber=res21$T.cells.regulatory..Tregs.,Decon_TRef = 0,Decon_BRef = 0,Decon_quan=res24$Tregs,Decon_immu=0,
                           ciber_abs_ciber=res26$T.cells.regulatory..Tregs.,ciber_abs_quan=res29$Tregs,ciber_abs_immu=0,
                           scdc_ciber=res31$T.cells.regulatory..Tregs.,scdc_TRef = 0,scdc_BRef = 0,scdc_quan=res34$Tregs,scdc_immu=0,
                           music_ciber=res36$T.cells.regulatory..Tregs.,music_TRef = 0,music_BRef = 0,music_quan=res39$Tregs,music_immu=0,
                           bseq_ciber=res41$T.cells.regulatory..Tregs.,bseq_TRef = 0,bseq_BRef = 0,bseq_quan=res44$Tregs,bseq_immu=0)

    nk_data = data.frame(ciber_ciber=res1$NK.cells.resting+res1$NK.cells.activated,ciber_quan=res4$NK.cells,ciber_immu=res5$CD56bright_natural_killer_cell+res5$CD56dim_natural_killer_cell,
                         EPIC_ciber=res6$NK.cells.resting+res6$NK.cells.activated,EPIC_TRef = res7$NKcells,EPIC_BRef = res8$NKcells,EPIC_quan=res9$NK.cells,EPIC_immu=res10$CD56bright_natural_killer_cell+res10$CD56dim_natural_killer_cell,
                         MCP_ciber=res11$NK.cells.resting+res11$NK.cells.activated,MCP_TRef = res12$NK.cells.activated,MCP_BRef = res13$NKcells,MCP_quan=res14$NK.cells,MCP_immu=res15$CD56bright_natural_killer_cell+res15$CD56dim_natural_killer_cell,
                         quan_ciber=res16$NK.cells.resting+res16$NK.cells.activated,quan_TRef = res17$NK.cells.activated,quan_BRef = res18$NKcells,quan_quan=res19$NK.cells,quan_immu=res20$CD56bright_natural_killer_cell+res20$CD56dim_natural_killer_cell,
                         Decon_ciber=res21$NK.cells.resting+res21$NK.cells.activated,Decon_TRef = res22$NK.cells.activated,Decon_BRef = res23$NKcells,Decon_quan=res24$NK.cells,Decon_immu=res25$CD56bright_natural_killer_cell+res25$CD56dim_natural_killer_cell,
                         ciber_abs_ciber=res26$NK.cells.resting+res26$NK.cells.activated,ciber_abs_quan=res29$NK.cells,ciber_abs_immu=res30$CD56bright_natural_killer_cell+res30$CD56dim_natural_killer_cell,
                         scdc_ciber=res31$NK.cells.resting+res31$NK.cells.activated,scdc_TRef = res32$NK.cells.activated,scdc_BRef = res33$NKcells,scdc_quan=res34$NK.cells,scdc_immu=res35$CD56bright_natural_killer_cell+res35$CD56dim_natural_killer_cell,
                         music_ciber=res36$NK.cells.resting+res36$NK.cells.activated,music_TRef = res37$NK.cells.activated,music_BRef = res38$NKcells,music_quan=res39$NK.cells,music_immu=res40$CD56bright_natural_killer_cell+res40$CD56dim_natural_killer_cell,
                         bseq_ciber=res41$NK.cells.resting+res41$NK.cells.activated,bseq_TRef = res42$NK.cells.activated,bseq_BRef = res43$NKcells,bseq_quan=res44$NK.cells,bseq_immu=res45$CD56bright_natural_killer_cell+res45$CD56dim_natural_killer_cell)

    mm_data = data.frame(ciber_ciber=res1$Monocytes+res1$Macrophages.M0+res1$Macrophages.M1+res1$Macrophages.M2,ciber_quan=res4$Macrophages.M1+res4$Macrophages.M2+res4$Monocytes,ciber_immu=res5$CD14_positive_monocyte+res5$CD16_positive_monocyte+res5$macrophage_m0+res5$macrophage_m1+res5$macrophage_m2,
                         EPIC_ciber=res6$Monocytes+res6$Macrophages.M0+res6$Macrophages.M1+res6$Macrophages.M2,EPIC_TRef = res7$Macrophages,EPIC_BRef = res8$Monocytes,EPIC_quan=res9$Macrophages.M1+res9$Macrophages.M2+res9$Monocytes,EPIC_immu=res10$CD14_positive_monocyte+res10$CD16_positive_monocyte+res10$macrophage_m0+res10$macrophage_m1+res10$macrophage_m2,
                         MCP_ciber=res11$Monocytes+res11$Macrophages.M0+res11$Macrophages.M1+res11$Macrophages.M2,MCP_TRef = res12$Macrophages.M0,MCP_BRef = res13$Monocytes,MCP_quan=res14$Macrophages.M1+res14$Macrophages.M2+res14$Monocytes,MCP_immu=res15$CD14_positive_monocyte+res15$CD16_positive_monocyte+res15$macrophage_m0+res15$macrophage_m1+res15$macrophage_m2,
                         quan_ciber=res16$Monocytes+res16$Macrophages.M0+res16$Macrophages.M1+res16$Macrophages.M2,quan_TRef = res17$Macrophages.M0,quan_BRef = res18$Monocytes,quan_quan=res19$Macrophages.M1+res19$Macrophages.M2+res19$Monocytes,quan_immu=res20$CD14_positive_monocyte+res20$CD16_positive_monocyte+res20$macrophage_m0+res20$macrophage_m1+res20$macrophage_m2,
                         Decon_ciber=res21$Monocytes+res21$Macrophages.M0+res21$Macrophages.M1+res21$Macrophages.M2,Decon_TRef = res22$Macrophages.M0,Decon_BRef = res23$Monocytes,Decon_quan=res24$Macrophages.M1+res24$Macrophages.M2+res24$Monocytes,Decon_immu=res25$CD14_positive_monocyte+res25$CD16_positive_monocyte+res25$macrophage_m0+res25$macrophage_m1+res25$macrophage_m2,
                         ciber_abs_ciber=res26$Monocytes+res26$Macrophages.M0+res26$Macrophages.M1+res26$Macrophages.M2,ciber_abs_quan=res29$Macrophages.M1+res29$Macrophages.M2+res29$Monocytes,ciber_abs_immu=res30$CD14_positive_monocyte+res30$CD16_positive_monocyte+res30$macrophage_m0+res30$macrophage_m1+res30$macrophage_m2,
                         scdc_ciber=res31$Monocytes+res31$Macrophages.M0+res31$Macrophages.M1+res31$Macrophages.M2,scdc_TRef = res32$Macrophages.M0,scdc_BRef = res33$Monocytes,scdc_quan=res34$Macrophages.M1+res34$Macrophages.M2+res34$Monocytes,scdc_immu=res35$CD14_positive_monocyte+res35$CD16_positive_monocyte+res35$macrophage_m0+res35$macrophage_m1+res35$macrophage_m2,
                         music_ciber=res36$Monocytes+res36$Macrophages.M0+res36$Macrophages.M1+res36$Macrophages.M2,music_TRef = res37$Macrophages.M0,music_BRef = res38$Monocytes,music_quan=res39$Macrophages.M1+res39$Macrophages.M2+res39$Monocytes,music_immu=res40$CD14_positive_monocyte+res40$CD16_positive_monocyte+res40$macrophage_m0+res40$macrophage_m1+res40$macrophage_m2,
                         bseq_ciber=res41$Monocytes+res41$Macrophages.M0+res41$Macrophages.M1+res41$Macrophages.M2,bseq_TRef = res42$Macrophages.M0,bseq_BRef = res43$Monocytes,bseq_quan=res44$Macrophages.M1+res44$Macrophages.M2+res44$Monocytes,bseq_immu=res45$CD14_positive_monocyte+res45$CD16_positive_monocyte+res45$macrophage_m0+res45$macrophage_m1+res45$macrophage_m2)
  }
  
  DECEPTICON_res = data.frame(B_cell =(b_cell_data[,colnames(b_cell_data) == b_res[1]]*0.25+ b_cell_data[,colnames(b_cell_data) == b_res[2]]*0.25+ b_cell_data[,colnames(b_cell_data) == b_res[3]]*0.25+ b_cell_data[,colnames(b_cell_data) == b_res[4]]*0.25),
                            cd4 = (cd4_data[,colnames(cd4_data) == cd4_res[1]]*0.25+ cd4_data[,colnames(cd4_data) == cd4_res[2]]*0.25+ cd4_data[,colnames(cd4_data) == cd4_res[3]]*0.25+ cd4_data[,colnames(cd4_data) == cd4_res[4]]*0.25),
                            cd8 = (cd8_data[,colnames(cd8_data) == cd8_res[1]]*0.25+ cd8_data[,colnames(cd8_data) == cd8_res[2]]*0.25+ cd8_data[,colnames(cd8_data) == cd8_res[3]]*0.25+ cd8_data[,colnames(cd8_data) == cd8_res[4]]*0.25),
                            treg = (treg_data[,colnames(treg_data) == treg_res[1]]*0.25+ treg_data[,colnames(treg_data) == treg_res[2]]*0.25+ treg_data[,colnames(treg_data) == treg_res[3]]*0.25+ treg_data[,colnames(treg_data) == treg_res[4]]*0.25),
                            nk = (nk_data[,colnames(nk_data) == nk_res[1]]*0.25+ nk_data[,colnames(nk_data) == nk_res[2]]*0.25+ nk_data[,colnames(nk_data) == nk_res[3]]*0.25+ nk_data[,colnames(nk_data) == nk_res[4]]*0.25),
                            mm = (mm_data[,colnames(mm_data) == mm_res[1]]*0.25+ mm_data[,colnames(mm_data) == mm_res[2]]*0.25+ mm_data[,colnames(mm_data) == mm_res[3]]*0.25+ mm_data[,colnames(mm_data) == mm_res[4]]*0.25))
  res = DECEPTICON_res
  write.table(res, './res/DECPTICON.txt', sep='\t', row.names= T, col.names= T, quote=F)
}

#'DEconvolution for CEll Proportion esTImation by CONsensus
#'
#'@param bulk.samples a m*n matrix with m genes and n samples
#'@param RUNpath Working path for storing files
#'@param light set to TRUE for light version of DECEPTICON
#'@param custom.signature set to TRUE for custom signature matrix
#'@param signature.matrix a i*l matrix with i genes and l cell types.
#'@export
run_DECEPTICON <- function(bulk.samples, RUNpath, light = FALSE, custom.signature = FALSE, signature.matrix = NULL){
  bulk.samples = bulk.samples
  RUNpath = RUNpath
  light = light
  dir.create("res")
  if(custom.signature == FALSE){
    DECEPTICON_methods(bulk.samples, RUNpath, light)
    DECEPTICON_output(light)
  } else {
    DECEPTICON_custom_methods(bulk.samples = bulk.samples, RUNpath = RUNpath, signature.matrix = signature.matrix)
    DECEPTICON__custom_output(signature.matrix = signature.matrix)
  }
  message("Deconvolution sucessful!")
}
