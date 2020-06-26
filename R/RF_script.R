```
rm(list=ls())
options(stringsAsFactors = FALSE)
library(DESeq2)
library(RWeka)
library(pROC)
library(ROCR)
library(randomForest)
out.path <-"/Users/wenbinye/Documents/study/methyaltion/result/AGE_matched_16VPC_CPC_with_age_control/RF"
if(!dir.exists(out.path)){
  dir.create(out.path)
}
out.figure <- "Metastatic_Localized.pptx"
p_theme<- theme_bw()+theme(
  text=element_text(family="Helvetica"),
  axis.title.x =element_text(color="black",size=16,family="Helvetica") ,
  axis.title.y =element_text(color="black",size=16,family="Helvetica") ,
  axis.text.x =element_text(color="black",size=14,family="Helvetica") ,
  axis.text.y =element_text(color="black",size=14,family="Helvetica") ,
  legend.text =element_text(color="black",size=12,family="Helvetica"),
  legend.title=element_text(color="black",size=14,family="Helvetica"),
  legend.background = element_blank(),
  panel.border = element_blank(),panel.grid.major = element_blank(),
  panel.grid.minor =element_blank(),axis.line=element_line(colour = "black",size=0.4))


plot.vocano <- function(de.data=NULL,p_theme=NULL){
  p<-ggplot(de.data,aes(x=log2FoldChange,y=-log10(padj),col="#999999"))+
    geom_point(color="#999999",alpha=0.5)+
    geom_point(data=subset(de.data,padj<0.05 & abs(log2FoldChange)>0),
               color="#4DAF4A",alpha=0.5)+
    geom_point(data=subset(de.data,padj<0.05 & abs(log2FoldChange)>1),
               color="#377EB8",alpha=0.5)+
    geom_point(data=subset(de.data,padj<0.05 & abs(log2FoldChange)>2),
               color="#E41A1C",alpha=0.5)+p_theme
  
  p1<- p+geom_hline(yintercept = -log10(0.05), linetype="dotted")+
    geom_vline(xintercept = c(-2,-1, 1,2), linetype="dotted") +
    labs(x=bquote(~Log[2]~fold~change),y=bquote(~-Log[10]~q~value))+
    p_theme 
  return(p1)
}

setwd(out.path)


#====================================================================
#' @title M Mechine learing
#' @param pe.vpc DESeq2 result of age-matched 16 VPC cases versus 16 CPCGENE
#' @param pe_hyper  RPKM table based on 26630 DMRs
pe.vpc <-  readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/PE_control_Age/PE_VPC16_vs_CPC16_removeAge_5_result_deseq.Rdata")
pe_hyper<- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/AGE_matched_16VPC_CPC_with_age_control/RF/Prostate_PE_25630_DMRs_VPC_CPC_153.Rdata")
hyper.index <- rownames(subset(pe.vpc,status=="hyper methylated"))
hypo.index <- rownames(subset(pe.vpc,status=="hypo methylated"))

setdiff(as.character(hyper.index),as.character(rownames(pe_hyper)))
setdiff(as.character(hypo.index),as.character(rownames(pe_hyper)))
index <- which(as.character(rownames(pe_hyper)) %in% as.character(hyper.index))
rownames(pe_hyper)[index] <- paste0("hyper_",rownames(pe_hyper)[index])
rownames(pe_hyper)[-index] <- paste0("hypo_",rownames(pe_hyper)[-index])
table(as.character(sapply(strsplit(rownames(pe_hyper),"\\_"),"[[",1)))


sample.names <- colnames(pe_hyper)
group <- gsub("\\d+","",sample.names)
group[grep("Norm|BL",sample.names)] <- "Control"
group[grep("sorted_V",sample.names)] <- "VPC"
group[grep("sorted_J|G",sample.names)] <- "CPCGENE"
group[grep("sorted_C",sample.names)] <- "Benign"
group[grep("sorted_R",sample.names)] <- "Barrier"
group[grep("bam_Counts.wig.txt",group)] <- "Tony"
table(group)
colnames(pe_hyper) <- paste0(group,"_", as.character(sapply(strsplit(colnames(pe_hyper),"\\_"),"[[",2)))

control.id <- which(group=="Control")
group <- group[-control.id]
pe_hyper <- pe_hyper[,-control.id]


pe_hyper <- log2(pe_hyper+1)

data.rf <- pe_hyper[,grep("VPC|CPCGENE",colnames(pe_hyper))]
data.vd <- pe_hyper[,-grep("VPC|CPCGENE",colnames(pe_hyper))]


sample.id <- as.character(sapply(strsplit(colnames(data.rf),"\\_"),"[[",2))
keep.id <-c("G19","J7","J3","J1","G15","J15","G12","J6","J4","G10","J16",
            "J10","J12","G5","G11","J5",
            "V25","V9","V18","V28",
            "V10-33","V37","V12","V13","V32",
            "V17-33","V36","V21","V38","V10","V27-33","V33-33")
data.rf.train <- data.rf[,which(sample.id %in% keep.id)]
data.rf.remain <- data.rf[,-which(sample.id %in% keep.id)]


data.rf.remain.cpc <- data.rf.remain[,grep("CPCGENE",colnames(data.rf.remain))]
data.rf.remain.vpc <- data.rf.remain[,grep("VPC",colnames(data.rf.remain))]
#===========================================================================
#' @title Split train and test
result.accuracy <- list(NULL)
result.temp<-list(NULL)
result.AUC <- list(NULL)


for (itertion in 1:50){
  print("#####itertion time#########")
  print(itertion)
  ind <- sample(2,ncol(data.rf.remain.cpc),replace=TRUE,prob=c(0.5,0.5))
  trainLC <- t(data.rf.remain.cpc)[ind==1,]
  testLC <- t(data.rf.remain.cpc)[ind==2,]
  
  ind <- sample(2,ncol(data.rf.remain.vpc),replace=TRUE,prob=c(0.7,0.3))
  trainMT <- t(data.rf.remain.vpc)[ind==1,]
  testMT<- t(data.rf.remain.vpc)[ind==2,]
  
  
  trainFeature <- rbind(t(data.rf.train),trainLC,trainMT)
  testFeature <- rbind(testLC,testMT)
  
  
  print(unique(as.character(sapply(strsplit(rownames(trainFeature),"\\_"),"[[",2))))
  print("######################")
  print(unique(as.character(sapply(strsplit(rownames(testFeature),"\\_"),"[[",2))))
  print("######################")
  print(table(as.character(sapply(strsplit(rownames(trainFeature),"\\_"),"[[",1))))
  print(table(as.character(sapply(strsplit(rownames(testFeature),"\\_"),"[[",1))))
  
  train <- data.frame(trainFeature,class= as.character(sapply(strsplit(rownames(trainFeature),"\\_"),"[[",1)) )
  test <- data.frame(testFeature,class= as.character(sapply(strsplit(rownames(testFeature),"\\_"),"[[",1)) )
  
  #table(train$class)
  #table(test$class)
  
  train$class <- as.factor(train$class)
  test$class <- as.factor(test$class)
  
  train.hyper <- train[,grep("hyper",colnames(train))]
  train.hypo <- train[,grep("hypo",colnames(train))]
  train.hyper$class <- train$class
  train.hypo$class <- train$class
  
  #options(expression = 5e5)
  infor.hyper <- InfoGainAttributeEval(class ~.,data=train.hyper[,c(1:10000,ncol(train.hyper))])
  infor.hyper2 <- InfoGainAttributeEval(class ~.,data=train.hyper[,c(10001:ncol(train.hyper))])
  infor.hypo <- InfoGainAttributeEval(class ~.,data=train.hypo)
  infor.hyper <- c(infor.hyper,infor.hyper2)
  
  df.hyer <- data.frame(feature=names(infor.hyper),information=infor.hyper)
  df.hyer <- df.hyer[order(df.hyer$information,decreasing = TRUE),]
  df.hyper.top <- df.hyer[1:150,]
  
  df.hypo <- data.frame(feature=names(infor.hypo),information=infor.hypo)
  df.hypo <- df.hypo[order(df.hypo$information,decreasing = TRUE),]
  df.hypo.top <- df.hypo[1:150,]
  
  df.top <- rbind(df.hyper.top,df.hypo.top)
  
  index <- match(df.top$feature,colnames(train))
  
  trainSet <- train[,index]
  testSet <- test[,index]
  trainSet$class <- train$class
  testSet$class <- test$class
  
  #dim(importance(model))
  #OBB <- model$err.rate[500,1]
  model <- randomForest(class ~ ., data = trainSet, proximity = TRUE, type = "classification", importance = TRUE)
  saveRDS(model,file = paste0("RF_model_",itertion,".Rdata"))
  pred <- predict(model, testSet[, -ncol(testSet)])
  proba <- as.data.frame(predict(model, testSet[, -ncol(testSet)], type = "prob"))
  ROC <- plotROC(proba, testSet$class)
  #setwd(out.path)
  #graph2ppt(file=out.figure, append=TRUE)
  result <- as.matrix(table(pred = pred, true = testSet$class))
  accuracy <- sum(diag(result))/sum(result)
  tp <- diag(result)
  fp <- rowSums(result) - diag(result)
  fn <- colSums(result) - diag(result)
  tn <- sum(result) - tp - fn - fp
  sp <- tn/(tn + fp)
  precision <- diag(result)/rowSums(result)
  recall <- (diag(result)/colSums(result))
  f1 <- 2 * precision * recall/(precision + recall)
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroSP <- mean(sp)
  temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
  result.temp[[itertion]] <- temp
  result.accuracy[[itertion]] <- accuracy
  print("################AUC")
  print(ROC$AUC)
  result.AUC[[itertion]] <- ROC$AUC
}



auc.table <- data.frame(NULL)
for(i in 1:length(result.AUC)){
  auc <- result.AUC[[i]][1]
  accuracy <-result.accuracy[[i]] 
  precision <- result.temp[[i]]$precision[3]
  sp <-result.temp[[i]]$sp[3]
  recall <- result.temp[[i]]$recall[3]
  f1 <- result.temp[[i]]$f1[3]
  temp <- data.frame(auc,accuracy,precision,sp ,recall,f1)
  auc.table <- rbind(auc.table ,temp )
}
auc.table.p<-reshape2::melt(auc.table)
auc.table.p$variable <- factor(auc.table.p$variable,levels = c("auc","accuracy","precision","sp","recall","f1"),
                               labels = c("AUC","Accurary","Precision","Specificity","Recall","F-score"))

library(wesanderson)
wes_palette("Darjeeling1")
colors <- brewer.pal(8,"Dark2")
ggplot(auc.table.p,aes(x=variable,y=value,color=variable))+geom_boxplot()+
   scale_color_manual(values=c(wes_palette("Darjeeling1"),colors[3]))+
  geom_jitter(position = position_jitter(0.2))+p_theme

ggplot(auc.table.p,aes(x=variable,y=value,color=variable))+geom_boxplot()+
  scale_color_manual(values=c("#0072B2",colors[c(1:3,5:6)]))+
  geom_jitter(position = position_jitter(0.2))+p_theme+guides(color=FALSE)+
  labs(x=NULL,y="Values")
  
ggplot(subset(auc.table.p,variable!="Specificity"),aes(x=variable,y=value,color=variable))+geom_boxplot()+
  scale_color_manual(values=c(colors[c(1:3,5:6)]))+
  geom_jitter(position = position_jitter(0.2))+p_theme+guides(color=FALSE)+
  labs(x=NULL,y="Value")
graph2ppt(file=out.figure,append=TRUE,width = 6, height = 4)






data.tony <- data.vd[,grep("Tony",colnames(data.vd))]
tony.accuracy <- list(NULL)
tony.temp<-list(NULL)
tony.AUC <- list(NULL)
proba.table <- data.frame(index=1:11)
for (itertion in 1:50){
  setwd(out.path)
  model <- readRDS(paste0("RF_model_",itertion,".Rdata"))
  index <- match(rownames(importance(model)), rownames(data.tony))
  data.test <- data.tony[index,]
  data.test <- t(data.test)
  
  data.test.cpc <- data.rf.remain.cpc[index,1:2]
  data.test.cpc <- t( data.test.cpc)
  
  pred <- predict(model, rbind(data.test,data.test.cpc))
  proba <- as.data.frame(predict(model, rbind(data.test,data.test.cpc), type = "prob"))
  proba.table <- cbind(proba.table,proba)
  ROC <- plotROC(proba, rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc))))
  
  result <- as.matrix(table(pred = pred, rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc)))))
  accuracy <- sum(diag(result))/sum(result)
  tp <- diag(result)
  fp <- rowSums(result) - diag(result)
  fn <- colSums(result) - diag(result)
  tn <- sum(result) - tp - fn - fp
  sp <- tn/(tn + fp)
  precision <- diag(result)/rowSums(result)
  recall <- (diag(result)/colSums(result))
  f1 <- 2 * precision * recall/(precision + recall)
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroSP <- mean(sp)
  temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
  tony.temp[[itertion]] <- temp
  tony.accuracy[[itertion]] <- accuracy
  print("################AUC")
  print(ROC$AUC)
  tony.AUC[[itertion]] <- ROC$AUC
}


auc.table <- data.frame(NULL)
for(i in 1:length(result.AUC)){
  auc <- tony.AUC[[i]][2]
  accuracy <-tony.accuracy[[i]] 
  precision <- tony.temp[[i]]$precision[2]
  sp <-tony.temp[[i]]$sp[2]
  recall <- tony.temp[[i]]$recall[2]
  f1 <- tony.temp[[i]]$f1[2]
  temp <- data.frame(auc,accuracy,precision,sp ,recall,f1)
  auc.table <- rbind(auc.table ,temp )
}
auc.table.p<-reshape2::melt(auc.table)
auc.table.p$variable <- factor(auc.table.p$variable,levels = c("auc","accuracy","precision","sp","recall","f1"),
                               labels = c("AUC","Accurary","Precision","Specificity","Recall","F-score"))
ggplot(subset(auc.table.p,variable!="Specificity"),aes(x=variable,y=value,color=variable))+geom_boxplot()+
  scale_color_manual(values=c(colors[c(1:3,5:6)]))+geom_point()+
  p_theme+guides(color=FALSE)+
  labs(x=NULL,y="Value")+ylim(0.9,1)
graph2ppt(file=out.figure,append=TRUE,width = 6, height = 4)

#geom_jitter(position = position_jitter(0.8))

proba.table$index <- NULL 
proba.table.cpc <- proba.table[,seq(from=1,to=ncol(proba.table),by=2)]
proba.table.vpc <- proba.table[,seq(from=2,to=ncol(proba.table),by=2)]
proba.table.cpc.pro <- rowSums(proba.table.cpc)/50
proba.table.vpc.pro <- rowSums(proba.table.vpc)/50
proba.table.average <-data.frame(CPCGENE=proba.table.cpc.pro,VPC=proba.table.vpc.pro) 
label <- rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc)))

averageROC <-plotROC(proba.table.average,label)
g <- ggplot(subset(averageROC$ROC,type=="VPC"), aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
  scale_color_manual(breaks = c("CPCGENE", "VPC"),
                     values = c("#1B9E77" ,"#00AFBB")) +
  labs(x = "False positive rate", y = "Ture positive rate",title="Tony (n=9)") +
  theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()


g <- g + annotate("text", x = 0.35, y = 0.75, label = paste("AUC", averageROC$AUC[which(names(averageROC$AUC) == "VPC")]), color = "#00AFBB", size = 6) +p_theme+
     guides(color = FALSE)+geom_segment(aes(x=0,y=0,xend=1,yend=1),color="gray")+theme(plot.title = element_text(hjust=0.5))
print(g)
graph2ppt(file=out.figure,append=TRUE,width = 4, height = 4)






data.barrier <- data.vd[,grep("Barrier",colnames(data.vd))]
barrier.accuracy <- list(NULL)
barrier.temp<-list(NULL)
barrier.AUC <- list(NULL)
proba.table <- data.frame(index=1:16)
for (itertion in 1:50){
  setwd(out.path)
  model <- readRDS(paste0("RF_model_",itertion,".Rdata"))
  index <- match(rownames(importance(model)), rownames(data.barrier))
  data.test <- data.barrier[index,]
  data.test <- t(data.test)
  
  data.test.cpc <- data.rf.remain.cpc[index,1:2]
  data.test.cpc <- t( data.test.cpc)
  
  pred <- predict(model, rbind(data.test,data.test.cpc))
  proba <- as.data.frame(predict(model, rbind(data.test,data.test.cpc), type = "prob"))
  proba.table <- cbind(proba.table,proba)
  ROC <- plotROC(proba, rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc))))
  
  result <- as.matrix(table(pred = pred, rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc)))))
  accuracy <- sum(diag(result))/sum(result)
  tp <- diag(result)
  fp <- rowSums(result) - diag(result)
  fn <- colSums(result) - diag(result)
  tn <- sum(result) - tp - fn - fp
  sp <- tn/(tn + fp)
  precision <- diag(result)/rowSums(result)
  recall <- (diag(result)/colSums(result))
  f1 <- 2 * precision * recall/(precision + recall)
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroSP <- mean(sp)
  temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
  barrier.temp[[itertion]] <- temp
  barrier.accuracy[[itertion]] <- accuracy
  print("################AUC")
  print(ROC$AUC)
  barrier.AUC[[itertion]] <- ROC$AUC
}


auc.table <- data.frame(NULL)
for(i in 1:length(result.AUC)){
  auc <- barrier.AUC[[i]][2]
  accuracy <-barrier.accuracy[[i]] 
  precision <- barrier.temp[[i]]$precision[2]
  sp <-barrier.temp[[i]]$sp[2]
  recall <- barrier.temp[[i]]$recall[2]
  f1 <- barrier.temp[[i]]$f1[2]
  temp <- data.frame(auc,accuracy,precision,sp ,recall,f1)
  auc.table <- rbind(auc.table ,temp )
}
auc.table.p<-reshape2::melt(auc.table)
auc.table.p$variable <- factor(auc.table.p$variable,levels = c("auc","accuracy","precision","sp","recall","f1"),
                               labels = c("AUC","Accurary","Precision","Specificity","Recall","F-score"))
ggplot(subset(auc.table.p,variable!="Specificity"),aes(x=variable,y=value,color=variable))+geom_boxplot()+
  scale_color_manual(values=c(colors[c(1:3,5:6)]))+geom_point()+
  p_theme+guides(color=FALSE)+
  labs(x=NULL,y="Value")
graph2ppt(file=out.figure,append=TRUE,width = 6, height = 4)

#geom_jitter(position = position_jitter(0.8))

proba.table$index <- NULL 
proba.table.cpc <- proba.table[,seq(from=1,to=ncol(proba.table),by=2)]
proba.table.vpc <- proba.table[,seq(from=2,to=ncol(proba.table),by=2)]
proba.table.cpc.pro <- rowSums(proba.table.cpc)/50
proba.table.vpc.pro <- rowSums(proba.table.vpc)/50
proba.table.average <-data.frame(CPCGENE=proba.table.cpc.pro,VPC=proba.table.vpc.pro) 
label <- rep(c("VPC","CPCGENE"),c(nrow(data.test),nrow(data.test.cpc)))

averageROC <-plotROC(proba.table.average,label)
g <- ggplot(subset(averageROC$ROC,type=="VPC"), aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
  scale_color_manual(breaks = c("CPCGENE", "VPC"),
                     values = c("#1B9E77" ,"#00AFBB")) +
  labs(x = "False positive rate", y = "Ture positive rate",title="Barrier (n=14)") +
  theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()


g <- g + annotate("text", x = 0.35, y = 0.75, label = paste("AUC", averageROC$AUC[which(names(averageROC$AUC) == "VPC")]), color = "#00AFBB", size = 6) +p_theme+
  guides(color = FALSE)+geom_segment(aes(x=0,y=0,xend=1,yend=1),color="gray")+theme(plot.title = element_text(hjust=0.5))
print(g)
graph2ppt(file=out.figure,append=TRUE,width = 4, height = 4)



#==============================
data.barrier <- data.vd
barrier.accuracy <- list(NULL)
barrier.temp<-list(NULL)
barrier.AUC <- list(NULL)
proba.table <- data.frame(index=colnames(data.vd))
for (itertion in 1:50){
  setwd(out.path)
  model <- readRDS(paste0("RF_model_",itertion,".Rdata"))
  index <- match(rownames(importance(model)), rownames(data.barrier))
  data.test <- data.barrier[index,]
  data.test <- t(data.test)

  
  pred <- predict(model, data.test)
  proba <- as.data.frame(predict(model,data.test, type = "prob"))
  proba.table <- cbind(proba.table,proba)
  class <- as.character(sapply(strsplit(colnames(data.vd),"\\_"),"[[",1))
  class <- factor(class,levels=c("Benign","Barrier","Tony"),labels=c("CPCGENE","VPC","VPC"))
  ROC <- plotROC(proba,  class)

  result <- as.matrix(table(pred = pred,  class))
  accuracy <- sum(diag(result))/sum(result)
  tp <- diag(result)
  fp <- rowSums(result) - diag(result)
  fn <- colSums(result) - diag(result)
  tn <- sum(result) - tp - fn - fp
  sp <- tn/(tn + fp)
  precision <- diag(result)/rowSums(result)
  recall <- (diag(result)/colSums(result))
  f1 <- 2 * precision * recall/(precision + recall)
  macroPrecision <- mean(precision)
  macroRecall <- mean(recall)
  macroF1 <- mean(f1)
  macroSP <- mean(sp)
  temp <- rbind(data.frame(precision, sp, recall, f1), mean = c(macroPrecision, macroSP, macroRecall, macroF1))
  barrier.temp[[itertion]] <- temp
  barrier.accuracy[[itertion]] <- accuracy
  print("################AUC")
  print(ROC$AUC)
  barrier.AUC[[itertion]] <- ROC$AUC
}

auc.table <- data.frame(NULL)
for(i in 1:length(result.AUC)){
  auc <- barrier.AUC[[i]][2]
  accuracy <-barrier.accuracy[[i]] 
  precision <- barrier.temp[[i]]$precision[3]
  sp <-barrier.temp[[i]]$sp[3]
  recall <- barrier.temp[[i]]$recall[3]
  f1 <- barrier.temp[[i]]$f1[3]
  temp <- data.frame(auc,accuracy,precision,sp ,recall,f1)
  auc.table <- rbind(auc.table ,temp )
}
auc.table.p<-reshape2::melt(auc.table)
auc.table.p$variable <- factor(auc.table.p$variable,levels = c("auc","accuracy","precision","sp","recall","f1"),
                               labels = c("AUC","Accurary","Precision","Specificity","Recall","F-score"))
ggplot(subset(auc.table.p,variable!="Specificity"),aes(x=variable,y=value,color=variable))+geom_boxplot()+
  scale_color_manual(values=c(colors[c(1:3,5:6)]))+geom_point()+
  p_theme+guides(color=FALSE)+
  labs(x=NULL,y="Value")
graph2ppt(file=out.figure,append=TRUE,width = 6, height = 4)

#geom_jitter(position = position_jitter(0.8))

proba.table$index <- NULL 
proba.table.cpc <- proba.table[,seq(from=1,to=ncol(proba.table),by=2)]
proba.table.vpc <- proba.table[,seq(from=2,to=ncol(proba.table),by=2)]
proba.table.cpc.pro <- rowSums(proba.table.cpc)/50
proba.table.vpc.pro <- rowSums(proba.table.vpc)/50
proba.table.average <-data.frame(CPCGENE=proba.table.cpc.pro,VPC=proba.table.vpc.pro) 
label <-class

averageROC <-plotROC(proba.table.average,label)

g <- ggplot(averageROC$ROC, aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
  scale_color_manual(breaks = c("CPCGENE", "VPC"),
                     values = c("#D95F02" ,"#00AFBB")) +
  labs(x = "False positive rate", y = "Ture positive rate",title="Barrier (n=14)/ Tony (n=9)/ Benign(n=18)") +
  theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()


g <- g + annotate("text", x = 0.35, y = 0.75, label = paste("AUC", averageROC$AUC[which(names(averageROC$AUC) == "VPC")]), color = "#00AFBB", size = 6) +p_theme+
  guides(color = FALSE)+geom_segment(aes(x=0,y=0,xend=1,yend=1),color="gray")+theme(plot.title = element_text(hjust=0.5))
print(g)

graph2ppt(file=out.figure,append=TRUE,width = 4, height = 4)



topfeature.total <- c()
for (itertion in 1:50){
  setwd(out.path)
  model <- readRDS(paste0("RF_model_",itertion,".Rdata"))
  temp <- importance(model)
  topfeature <- rownames(temp)
  topfeature.total <- c(topfeature.total,topfeature )
  topfeature.total <- unique(topfeature.total)
}

MeanDecreaseAccuracy <- data.frame(feature=topfeature.total)
MeanDecreaseGini <- data.frame(feature=topfeature.total)
for (itertion in 1:50){
  setwd(out.path)
  model <- readRDS(paste0("RF_model_",itertion,".Rdata"))
  temp <- importance(model)
  temp <- as.data.frame(temp)
  index <- match(MeanDecreaseAccuracy$feature,rownames(temp))
  MeanDecreaseAccuracy[[paste0("N0_",itertion)]] <- temp$MeanDecreaseAccuracy[index]
  MeanDecreaseGini[[paste0("N0_",itertion)]] <- temp$MeanDecreaseGini[index]
}
rownames( MeanDecreaseAccuracy) <- MeanDecreaseAccuracy$feature
rownames( MeanDecreaseGini) <- MeanDecreaseGini$feature
MeanDecreaseAccuracy$feature<- NULL
MeanDecreaseGini$feature<- NULL
MeanDecreaseAccuracy.mean <- apply(MeanDecreaseAccuracy,1,function(x){ return(mean(as.numeric(x),na.rm=TRUE) ) })

MeanDecreaseGini.mean <- apply(MeanDecreaseGini,1,function(x){ return(mean(as.numeric(x),na.rm=TRUE) ) })

MeanDecreaseAccuracy.mean <- MeanDecreaseAccuracy.mean[order(MeanDecreaseAccuracy.mean,decreasing = TRUE)]
MeanDecreaseGini.mean <- MeanDecreaseGini.mean[order(MeanDecreaseGini.mean,decreasing = TRUE)]


txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gene_region <-  genes(txdb,columns="gene_id")
gene_region$type <- "gene"
library(clusterProfiler)
gene.df <- bitr(gene_region$gene_id, fromType = "ENTREZID",
                toType = c("ENSEMBL", "ENTREZID","SYMBOL"),
                OrgDb = org.Hs.eg.db)


index <- match(as.character(gene_region$gene_id),as.character(gene.df$ENTREZID))
gene_region$symbol <-  gene.df$SYMBOL[index]
gene_region$ensembl <-  gene.df$ENSEMBL[index]

load("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/annotation_Granges_genome.Rdata")
MeanDecreaseGini.names <-gsub("hyper\\_|hypo\\_","",names(MeanDecreaseGini.mean))
ov.DecreaseGini <- distanceToNearest(GRanges.genome[as.numeric(MeanDecreaseGini.names)],gene_region)

index <- which(ov.DecreaseGini@elementMetadata$distance<=5000)
ov.DecreaseGini <- ov.DecreaseGini[index]
ov.DecreaseGini@to

MeanDecreaseGini.feature <- data.frame(value=MeanDecreaseGini.mean,index=MeanDecreaseGini.names)
MeanDecreaseGini.feature$Nearestgene <-NA
MeanDecreaseGini.feature$Nearestgene[ov.DecreaseGini@from] <- gene_region$symbol[ov.DecreaseGini@to]


MeanDecreaseAccuracy.feature <- data.frame(value=MeanDecreaseAccuracy.mean,
                                           index=gsub("hyper\\_|hypo\\_","",names(MeanDecreaseAccuracy.mean)))
MeanDecreaseAccuracy.feature$Nearestgene <-NA
index <- match(MeanDecreaseAccuracy.feature$index,MeanDecreaseGini.feature$index)
MeanDecreaseAccuracy.feature$Nearestgene <- MeanDecreaseGini.feature$Nearestgene[index]


gini.hyper <-MeanDecreaseGini.feature[grep("hyper",rownames(MeanDecreaseGini.feature)),]
gini.hypo <-MeanDecreaseGini.feature[grep("hypo",rownames(MeanDecreaseGini.feature)),]

gini.hypo<-gini.hypo[-which(is.na(gini.hypo$Nearestgene)),]
gini.hyper<-gini.hyper[-which(is.na(gini.hyper$Nearestgene)),]

gini.hyper$Nearestgene <- factor(gini.hyper$Nearestgene,levels =  rev(unique(gini.hyper$Nearestgene)))

ggplot(gini.hyper[1:10,],aes(x=Nearestgene,y=value))+
  geom_bar(stat="identity",fill="#00AFBB")+coord_flip()+
  labs(x = NULL, y = "MeanDecreaseGini",title="Top10 Hyper-Features")+p_theme+
  theme(plot.title = element_text(hjust=0.5))
graph2ppt(file=out.figure,append=TRUE,width = 4, height = 4)

gini.hypo$Nearestgene <- factor(gini.hypo$Nearestgene,levels =  rev(unique(gini.hypo$Nearestgene)))
ggplot(gini.hypo[1:10,],aes(x=Nearestgene,y=value))+
  geom_bar(stat="identity",fill="#00AFBB")+coord_flip()+
  labs(x = NULL, y = "MeanDecreaseGini",title="Top10 Hypo-Features")+p_theme+
  theme(plot.title = element_text(hjust=0.5))
graph2ppt(file=out.figure,append=TRUE,width = 4, height = 4)




plotROC <- function(probalility, labels) {
  # ROC <- plotROC(proba, testSet$class)
  name <- colnames(probalility)
  Plotdata <- data.frame(x = 1, y = 2, type = "aa")
  Auc <- c()
  for (i in 1:length(name)) {
    label <- as.character(labels)
    label[which(labels == name[i])] <- 1
    label[-which(labels == name[i])] <- 0
    pred <- prediction(probalility[, i], label)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    auc <- ROCR::performance(pred, "auc")
    auc <- unlist(slot(auc, "y.values"))
    Auc <- c(Auc, auc)
    x <- unlist(perf@x.values)
    y <- unlist(perf@y.values)
    plotdata <- data.frame(x = x, y = y)
    plotdata$type <- name[i]
    Plotdata <- rbind(Plotdata, plotdata)
  }
  names(Auc) <- name
  Plotdata <- Plotdata[-1, ]
  Auc <- round(Auc, 4)
  names(Plotdata) <- c("fpr", "tpr", "type")
  g <- ggplot(Plotdata, aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
    scale_color_manual(breaks = c("CPCGENE", "VPC"),
                       values = c("#1B9E77" ,"#D95F02")) +
    labs(x = "False positive rate FPR=FP/N", y = "Ture positive rate TPR=TP/P") +
    theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()
  
  
  g <- g + annotate("text", x = 0.35, y = 0.85, label = paste("AUC", Auc[which(names(Auc) == "CPCGENE")]), color ="#1B9E77", size = 6) +
    annotate("text", x = 0.35, y = 0.75, label = paste("AUC", Auc[which(names(Auc) == "VPC")]), color = "#D95F02", size = 6) +p_theme+
    theme(legend.position = "top") + guides(color = guide_legend(title = NULL))
  #print(g)
  return(list(AUC = Auc, ROC = Plotdata))
  
}






data.barrier <- data.vd[,grep("barrier",colnames(data.vd))]













#test.result <- list(randomforestPre = pred, accuracy = accuracy, confusion = result,
#                    evaluate = temp, ROC = ROC)

test.result <- list(trainSet=trainSet,testSet=testSet,randomforestPre = pred, accuracy = accuracy, confusion = result,
                    evaluate = temp, ROC = ROC)








#saveRDS(test.result,file="FR_result_PE.Rdata")



plotROC <- function(probalility, labels) {
  # ROC <- plotROC(proba, testSet$class)
  name <- colnames(probalility)
  Plotdata <- data.frame(x = 1, y = 2, type = "aa")
  Auc <- c()
  for (i in 1:length(name)) {
    label <- as.character(labels)
    label[which(labels == name[i])] <- 1
    label[-which(labels == name[i])] <- 0
    pred <- prediction(probalility[, i], label)
    perf <- ROCR::performance(pred, "tpr", "fpr")
    auc <- ROCR::performance(pred, "auc")
    auc <- unlist(slot(auc, "y.values"))
    Auc <- c(Auc, auc)
    x <- unlist(perf@x.values)
    y <- unlist(perf@y.values)
    plotdata <- data.frame(x = x, y = y)
    plotdata$type <- name[i]
    Plotdata <- rbind(Plotdata, plotdata)
  }
  names(Auc) <- name
  Plotdata <- Plotdata[-1, ]
  Auc <- round(Auc, 4)
  names(Plotdata) <- c("fpr", "tpr", "type")
  g <- ggplot(Plotdata, aes(x = fpr, y = tpr, colour = type)) + geom_path(size = 1) +
    scale_color_manual(breaks = c("CPCGENE", "VPC"),
                       values = c("#1B9E77" ,"#D95F02")) +
    labs(x = "False positive rate FPR=FP/N", y = "Ture positive rate TPR=TP/P") +
    theme(plot.title = element_text(face = "bold", size = 15)) + theme_bw()
  
  
  g <- g + annotate("text", x = 0.35, y = 0.85, label = paste("AUC", Auc[which(names(Auc) == "CPCGENE")]), color ="#1B9E77", size = 6) +
    annotate("text", x = 0.35, y = 0.75, label = paste("AUC", Auc[which(names(Auc) == "VPC")]), color = "#D95F02", size = 6) +p_theme+
    theme(legend.position = "top") + guides(color = guide_legend(title = NULL))
  #print(g)
  return(list(AUC = Auc, ROC = Plotdata))
  
}
#colors <- brewer.pal(8,"Dark2")





























pe_hyper<- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/20200518_Figure_Control_Age/Prostate_PE_Total_24374_Hyprt_VPC_CPC_97.Rdata")
pe_hypo <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/20200518_Figure_Control_Age/Prostate_PE_Total_7662_Hypo_VPC_CPC_97.Rdata")
dim(pe_hypo)
mean.table <- data.frame(NULL)
for (i in 1:3){
  index <- match(data.combine[[i]],rownames(pe_hyper))
  print(which(is.na(index)))
  pe.hyper.temp <- pe_hyper[index,]
  mean.rpkm <- apply(pe.hyper.temp ,2,mean)
  mean.temp <- data.frame(mean.rpkm=mean.rpkm,sample=as.character(sapply(strsplit(colnames(pe.hyper.temp),"\\_"),"[[",2)) ,
                          type=names(data.combine)[i])
  mean.table <- rbind(mean.table,mean.temp)
}
mean.table$status <- "hyper"
mean.table$group <- "CPCGENE"
mean.table$group[grep("V",mean.table$sample)] <- "VPC"
mean.table.vpc <- subset(mean.table,group=="VPC")

mean.table.hypo <- data.frame(NULL)
for (i in 1:3){
  index <- match(data.combine.hypo[[i]],rownames(pe_hypo))
  pe.hypo.temp <- pe_hypo[index,]
  mean.rpkm <- apply(pe.hypo.temp ,2,mean)
  mean.temp <- data.frame(mean.rpkm=mean.rpkm,sample=as.character(sapply(strsplit(colnames(pe.hypo.temp),"\\_"),"[[",2)) ,
                          type=names(data.combine.hypo)[i])
  mean.table.hypo <- rbind(mean.table.hypo,mean.temp)
}
mean.table.hypo$status <- "hypo"
mean.table.hypo$group <- "CPCGENE"
mean.table.hypo$group[grep("V",mean.table.hypo$sample)] <- "VPC"
mean.table.hypo.vpc <- subset(mean.table.hypo,group=="VPC")

head(mean.table)
head(mean.table.hypo)
identical(mean.table$sample,mean.table.hypo$sample)

colnames(mean.table)[1] <- "hyper.mean"
mean.table$hypo.mean <- mean.table.hypo$mean.rpkm
mean.table$status <- NULL
length(unique(rownames(mean.table)));
length(unique(mean.table$sample))


vpc.clinical <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/Age_VPC_67cases_54patient.Rdata")
index <- match(mean.table$sample,vpc.clinical$sample.ID)
mean.table <- cbind(mean.table,vpc.clinical[index,])

CPCGENE.clinical <- read.csv("/Users/wenbinye/Dropbox/Jessica/WenbinYE/Data_information/Single_End_5mc/Clinical_data_30CPCGENE_update.csv",
                             row.names = 1)

index <- match(mean.table$sample,CPCGENE.clinical$MeDIP.sample)
mean.table$cpc_age <- CPCGENE.clinical$Age_at_Treatment[index]
table(mean.table$type)
mean.table$type <-factor(mean.table$type,levels = c("raw.67vpc.30cpc" ,"age.67vpc.30cpc","age.16vpc.16cpc"))

mean.table.vpc <- subset(mean.table,group=="VPC")
mean.table.cpc <- subset(mean.table,group=="CPCGENE")

table(mean.table.vpc$Timepoint)
mean.table.vpc$Timepoint <- factor(mean.table.vpc$Timepoint,levels=c("Baseline","1st progression","2nd progression"))
mean.table.vpc <- mean.table.vpc[order(mean.table.vpc$type,mean.table.vpc$`patient number`,mean.table.vpc$Timepoint),]
mean.table.vpc$label <- paste0(mean.table.vpc$type,mean.table.vpc$`patient number`)
mean.table.vpc$duplicated <- 0
mean.table.vpc$duplicated[duplicated(mean.table.vpc$label)]<-1
identical(mean.table.vpc$ctDNA,mean.table.vpc$`ctDNA%`)
mean.table.vpc$ctDNA <- as.numeric(mean.table.vpc$ctDNA)
#=================================================
#' @title Age and hyper mean
ggplot(subset(mean.table.vpc,duplicated ==0), 
       aes(x=Published_Age,y=hyper.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age",y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

#=================================================
#' @title Hyper and %ctDNA
colors <- brewer.pal(8,"Dark2")
ggplot(subset(mean.table.vpc,type=="age.16vpc.16cpc"), 
       aes(x=ctDNA,y=hyper.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))+geom_text_repel(
            data = subset(mean.table.vpc,type=="age.16vpc.16cpc"&hyper.mean>=4),
            aes(label=sample),
            size=4,
            box.padding=unit(0.35,"lines"),
            point.padding=unit(0.3,"lines")
          )
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



ggplot(mean.table.cpc, 
       aes(x=cpc_age,y=hyper.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age (CPCGENE)",y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

ggplot(mean.table.cpc, 
       aes(x=cpc_age,y=hypo.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age (CPCGENE)",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

#=================================================
#' @title Age and hypo mean
ggplot(subset(mean.table.vpc,duplicated ==0), 
       aes(x=Published_Age,y=hypo.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

#=================================================
#' @title %ctDNA and hypo mean
ggplot(subset(mean.table.vpc,duplicated ==0), 
       aes(x=ctDNA,y=hypo.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))

ggplot(subset(mean.table.vpc,type=="age.16vpc.16cpc"), 
       aes(x=ctDNA,y=hypo.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))+geom_text_repel(
            data = subset(mean.table.vpc,type=="age.16vpc.16cpc"&hypo.mean>=6),
            aes(label=sample),
            size=4,
            box.padding=unit(0.35,"lines"),
            point.padding=unit(0.3,"lines")
          )
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

#=================================================
#' @title hypo and hyper mean
ggplot(subset(mean.table.vpc,duplicated ==0), 
       aes(x=hyper.mean,y=hypo.mean,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Hyper-DMRs\nMean RPKM counts across sample",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

# ggplot(mean.table.vpc, 
#        aes(x=hyper.mean,y=hypo.mean,color=type))+geom_point()+theme_classic()+
#   p_theme+scale_color_manual(values=colors[c(1:3)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
#   stat_cor(method = "pearson")+
#   labs(x="Hyper-DMRs\nMean RPKM counts across sample",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=guide_legend(title=NULL))+
#   theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
#           legend.title=element_text(color="black",size=14,family="Helvetiva"),
#           axis.text.x = element_text(angle = 0,vjust=0.6),
#           strip.text.x = element_text(size=12),
#           legend.background = element_blank(),
#           legend.position = c(0.8,0.9))
# setwd(out.path)
# graph2ppt(file=out.figure, append=TRUE)

#=================================================
#' @title Age and %ctDNA
mean.table.vpc$`ctDNA%` <- as.numeric(mean.table.vpc$`ctDNA%`)
mean.table.vpc$ctDNA <- mean.table.vpc$`ctDNA%`
mean.table.temp <- subset(mean.table.vpc,duplicated ==0)
mean.table.temp <- mean.table.temp[!duplicated(mean.table.temp$sample),]
ggplot(mean.table.temp, 
       aes(x=Published_Age,y=ctDNA,color=type))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(5)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age",y="% ctDNA")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


#=================================================
mean.table.vpc.dup <- mean.table.vpc[which(mean.table.vpc$label %in% mean.table.vpc$label[which(mean.table.vpc$duplicated==1)]),]
mean.table.vpc.dup$patientID <- mean.table.vpc.dup$`patient number`
ggplot(mean.table.vpc.dup, 
       aes(x=Timepoint,y=hyper.mean,color=type,group=patientID))+geom_point()+theme_classic()+
  geom_line()+
  facet_wrap(~type)+p_theme+scale_color_manual(values=colors[c(1:3)])+
  labs(x=NULL,y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 45,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggplot(mean.table.vpc.dup, 
       aes(x=Timepoint,y=hypo.mean,color=type,group=patientID))+geom_point()+theme_classic()+
  geom_line()+
  facet_wrap(~type)+p_theme+scale_color_manual(values=colors[c(1:3)])+
  labs(x=NULL,y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 45,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

ggplot(mean.table.vpc.dup, 
       aes(x=Timepoint,y=hypo.mean,color=type,group=patientID))+geom_point()+theme_classic()+
  geom_line()+
  facet_wrap(~type)+p_theme+scale_color_manual(values=colors[c(1:3)])+
  labs(x=NULL,y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 45,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))+ylim(0,20)

setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

ggplot(mean.table.vpc.dup, 
       aes(x=Timepoint,y=ctDNA,color=type,group=patientID))+geom_point()+theme_classic()+
  geom_line()+
  facet_wrap(~type)+p_theme+scale_color_manual(values=colors[c(5:7)])+
  labs(x=NULL,y="% ctDNA")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 45,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.8,0.9))

setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



#============================================================
mean.rpkm <- subset(mean.table.vpc,duplicated ==0)
#mean.rpkm <- mean.table.vpc
mean.rpkm$OS.time <- as.numeric(mean.rpkm$`Published_Days to death or last followup`)
mean.rpkm$OS <- mean.rpkm$Published_Censored...14
mean.rpkm$OS[which(mean.rpkm$OS=="Yes")] <- 1
mean.rpkm$OS[which(mean.rpkm$OS=="No")] <- 2
mean.rpkm$OS <- as.numeric(mean.rpkm$OS)
mean.rpkm$RFS.time <- as.numeric(mean.rpkm$`Published_Days to progression or last followup`)
mean.rpkm$RFS <- mean.rpkm$Published_Censored...12
table(mean.rpkm$RFS )
mean.rpkm$RFS[which(mean.rpkm$RFS=="Yes")] <- 1
mean.rpkm$RFS[which(mean.rpkm$RFS=="No")] <- 2
mean.rpkm$RFS<-as.numeric(mean.rpkm$RFS)

mean.rpkm$hyper.cluster <- "1"
mean.rpkm$hypo.cluster <- "1"
mean.rpkm$ratio.cluster <- "1"
mean.rpkm$ratio <- as.numeric(mean.rpkm$hyper.mean)/as.numeric(mean.rpkm$hypo.mean)
type <- names(table(mean.rpkm$type))

mean.rpkm.renew <- data.frame(NULL)
for(i in type){
  mean.rpkm.temp <- subset(mean.rpkm,type==i)
  mean.rpkm.temp$hyper.cluster[which(mean.rpkm.temp$hyper.mean>=median(mean.rpkm.temp$hyper.mean))] <-2
  mean.rpkm.temp$hypo.cluster[which(mean.rpkm.temp$hypo.mean>=median(mean.rpkm.temp$hypo.mean))] <-2
  mean.rpkm.temp$ratio.cluster[which(mean.rpkm.temp$ratio>=median(mean.rpkm.temp$ratio))] <-2
  mean.rpkm.renew <- rbind(mean.rpkm.renew,mean.rpkm.temp)
}
table(mean.rpkm.renew$hyper.cluster)
table(mean.rpkm.renew$hypo.cluster)
table(mean.rpkm.renew$ratio.cluster)
#ref:https://github.com/kassambara/survminer/issues/64


mean.rpkm.renew.hyper <- mean.rpkm.renew
mean.rpkm.renew.hyper$cluter <- mean.rpkm.renew$hyper.cluster
mean.rpkm.renew.hyper$group <- "Hyper"

mean.rpkm.renew.hypo <- mean.rpkm.renew
mean.rpkm.renew.hypo$cluter <- mean.rpkm.renew$hypo.cluster
mean.rpkm.renew.hypo$group <- "Hypo"

mean.rpkm.renew.ratio <- mean.rpkm.renew
mean.rpkm.renew.ratio$cluter <- mean.rpkm.renew$ratio.cluster
mean.rpkm.renew.ratio$group <- "Ratio"
mean.rpkm.combine <- rbind(mean.rpkm.renew.hyper,mean.rpkm.renew.hypo,mean.rpkm.renew.ratio)



# OS.fit <- survfit(Surv(OS.time, OS) ~ hyper.cluster, data = subset(mean.rpkm.renew.ratio,type=="age.16vpc.16cpc"))
# ggsurvplot(OS.fit,size = 1,# change line size
#            palette = c("#E7B800", "#2E9FDF"), # custom color palette
#            conf.int = TRUE, # Add confidence interval
#            pval = TRUE, # Add p-value
#            #risk.table = TRUE, # Add risk table
#            risk.table.col = "Strata", # Risk table color by groups
#            ylab="Overall survival probability",
#            legend.labs=c("Low risk","High risk"),
#            ggtheme = p_theme# Change ggplot2 theme
# )

#RFS.fit <- survfit(Surv(RFS.time, RFS) ~ hyper.cluster+type, data = mean.rpkm)
#ggsurvplot(OS.fit)

mean.rpkm.combine.c <- subset(mean.rpkm.combine,type=="age.16vpc.16cpc")
library(survival)
library(survminer)
#https://rpkgs.datanovia.com/survminer/reference/ggsurvplot_facet.html
#https://rpkgs.datanovia.com/survminer/reference/ggsurvplot_facet.html
OS.fit <- survfit(Surv(OS.time, OS) ~ cluter, data = subset(mean.rpkm.combine.c,group=="Ratio"))
ggsurvplot(OS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("Ratio low","Ratio high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluter, data = subset(mean.rpkm.combine.c,group=="Ratio"))
ggsurvplot(RFS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progression-free survival probability",
           legend.labs=c("Ratio low","Ratio high"),
           ggtheme = p_theme# Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)






OS.fit <- survfit(Surv(OS.time, OS) ~ cluter, data =mean.rpkm.combine)
ggsurvplot(OS.fit, mean.rpkm.combine.c,size = 1, facet.by = c("type","group"),# change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           #risk.table = TRUE, # Add risk table
           risk.table.col = "Strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("Low methylation","High methylation"),
           ggtheme = p_theme# Change ggplot2 theme
)+theme(strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.position = "top")
setwd(out.path)

table(mean.rpkm.combine$type)
table(mean.rpkm.combine$group)

mean.rpkm.temp <- subset(mean.rpkm.combine, type=="age.67vpc.30cpc"&group=="Ratio")
mean.rpkm.temp$cluter <- "1"
mean.rpkm.temp$cluter[which(mean.rpkm.temp$ctDNA>median(mean.rpkm.temp$ctDNA))] <- "2"
OS.fit <- survfit(Surv(OS.time, OS) ~ cluter, data = mean.rpkm.temp)
ggsurvplot(OS.fit,size = 1,# change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("%ctDNA low","%ctDNA high"),
           ggtheme = p_theme# Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluter, data = mean.rpkm.temp)
ggsurvplot(RFS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progression-free survival probability",
           legend.labs=c("%ctDNA low","%ctDNA high"),
           ggtheme = p_theme# Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)




RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluter, data = mean.rpkm.combine)
ggsurvplot(RFS.fit, mean.rpkm.combine, size = 1,facet.by = c("type","group"), # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           # risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progression-free survival probability",
           legend.labs=c("Low methylation","High methylation"),
           ggtheme = p_theme# Change ggplot2 theme
)+theme(strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        legend.position = "top")
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



library(gridExtra)
#grid.arrange(p1,p2,p3,nrow=3)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
#p$plot+facet_grid(~type)+p_theme+
#  theme(strip.text.x = element_text(size=12),legend.position = "top")
library(ggquickeda)
#install.packages("ggquickeda")
OS.fit <- survfit(Surv(OS.time, OS) ~ ratio.cluster, data = mean.rpkm.renew)
p<- ggsurvplot(OS.fit, mean.rpkm.renew,facet.by = "type", conf.int = TRUE, pval = TRUE,
               palette = c("#E7B800", "#2E9FDF") )

p$plot+p_theme+theme(strip.text.x = element_text(size=12))
#+facet_grid(~type)
#fun = "cumhaz"
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
















library(VennDiagram)
data.combine <- list(raw.67vpc.30cpc = rownames(subset(pe.vpc,status=="hyper methylated") ),
                     age.67vpc.30cpc = rownames(subset(pe.vpc.age,status=="hyper methylated") ),
                     age.16vpc.16cpc = rownames(subset(pe.vpc.age.match,status=="hyper methylated") ))

lengths(data.combine)
cols=brewer.pal(9, "Set1")
venn.plot <- venn.diagram(list(A = data.combine[[1]],
                               B = data.combine[[2]],
                               C = data.combine[[3]] ),
                          col = "transparent",
                          fill=rev(cols[1:3]),
                          alpha=0.5,
                          label.col = c( "darkblue", "white", "darkblue",
                                         "white", "orange", "white", "darkblue"),
                          category.names = names(data.combine),
                          cex=1.4,
                          cat.fontface=3, 
                          filename=NULL,
                          cat.cex = 1.5,main=NULL,main.cex=3,
                          cat.fontfamily = "serif")
dev.off()
grid.draw(venn.plot)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

data.combine.hypo <- list(raw.67vpc.30cpc = rownames(subset(pe.vpc,status=="hypo methylated") ),
                          age.67vpc.30cpc = rownames(subset(pe.vpc.age,status=="hypo methylated") ),
                          age.16vpc.16cpc = rownames(subset(pe.vpc.age.match,status=="hypo methylated") ))

lengths(data.combine.hypo)
cols=brewer.pal(9, "Set1")
venn.plot <- venn.diagram(list(A = data.combine.hypo[[1]],
                               B = data.combine.hypo[[2]],
                               C = data.combine.hypo[[3]] ),
                          col = "transparent",
                          fill=rev(cols[1:3]),
                          alpha=0.5,
                          label.col = c( "darkblue", "white", "darkblue",
                                         "white", "orange", "white", "darkblue"),
                          category.names = names(data.combine.hypo),
                          cex=1.4,
                          cat.fontface=3, 
                          filename=NULL,
                          cat.cex = 1.5,main=NULL,main.cex=3,
                          cat.fontfamily = "serif")
dev.off()
grid.draw(venn.plot)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

1047/lengths(data.combine.hypo)[1]
47/lengths(data.combine.hypo)[2]
4113/lengths(data.combine.hypo)[3]








#=======================================================================
#===================================================
#Annotation hyperDMRs
library(annotatr)
pe.vpc.hyper <- subset(pe.vpc,status=="hyper methylated")
pe.vpc.hypo <- subset(pe.vpc,status=="hypo methylated")

#hypo.index <- as.numeric(rownames(pe.vpc.hypo))

load("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/annotation_Granges_genome.Rdata")
rm(annotation)
annotations <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/annotation.Rdata")

head(GRanges.genome)

deseq_annotated = annotate_regions(
  regions = GRanges.genome[as.numeric(rownames(pe.vpc.hyper))],
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)



annotated <- as.data.frame(deseq_annotated)


annotated$annot.type <- factor(annotated$annot.type,
                               levels=c(  "hg19_genes_promoters"  , "hg19_cpg_islands","hg19_cpg_shores" , "hg19_genes_1to5kb",
                                          "hg19_genes_3UTRs", "hg19_genes_5UTRs","hg19_cpg_shelves" , "hg19_genes_exons" ,"hg19_genes_introns",
                                          "hg19_genes_intronexonboundaries","hg19_genes_intergenic",  "hg19_cpg_inter"),
                               labels=c("Promoter","CpG island","CpG shore","Gene 1to5kb","Gene body","Gene body",
                                        "CpG shelves","Gene body","Gene body","Gene body","Intergenic","Intergenic"))
annotated <- annotated[order(annotated$index,annotated$annot.type),]
annotated <- annotated[!duplicated(annotated$index),]


annotated.need <- subset( annotated,annot.type %in% c("Promoter","CpG island","CpG shore"))

index.hyper.need <- as.numeric(annotated.need$index)
index.hyper.need <- index.hyper.need[order(index.hyper.need)]
length(unique(index.hyper.need ))
#saveRDS(index.hyper.need, file="/Users/wenbinye/Documents/study/methyaltion/result/AGE_matched_16VPC_CPC_with_age_control/PE_Hyper_13531_CpCisland_Shores_Promoter_index_16VPC_16CPC.Rdata")





table(annotated.need$annot.type)/nrow(pe.vpc.hyper)

table(annotated$annot.type)/nrow(pe.vpc.hyper)

annotate.frame <- as.data.frame(table(annotated$annot.type))
colnames(annotate.frame) <- c("group","value")
annotate.frame
sum(annotate.frame$value)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(export)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
annotate.frame$group <- as.character(annotate.frame$group)
annotate.frame$prop <- annotate.frame$value/sum(annotate.frame$value)

count.data <- annotate.frame%>%
  arrange(desc(group)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
count.data

p<- ggplot(count.data,aes(x=2,y=prop,fill=group))+
  geom_bar(stat="identity",color="white")+coord_polar("y",start=0)+
  geom_text(aes(y = lab.ypos, label = paste0(round(prop*100,2),"%")), 
            color = "white",size=4)+scale_fill_brewer(palette="Dark2")+
  theme_void()+xlim(0.5,2.5)
p+p_theme+
  guides(fill=guide_legend(title=NULL))+
  theme_void()+  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),#Helvetica
                         legend.title=element_text(color="black",size=14,family="Helvetica"),legend.position  ="top")
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE,width=5,height=5)


#============================================================================
#' @title DO Enrichment 
#' Promoter
library(org.Hs.eg.db)
annotation.promoter <- subset(annotated.need,annot.type=="Promoter")
#annotation.promoter <- subset(annotated,annot.type=="Promoter")
#annotation.promoter <- subset(annotated,annot.type=="Gene body")
length(unique(annotation.promoter$annot.id))
length(unique(annotation.promoter$annot.gene_id))


gene.df <- bitr(unique(annotation.promoter$annot.gene_id), fromType = "ENTREZID",
                toType = c("ENSEMBL", "ENTREZID","SYMBOL"),
                OrgDb = org.Hs.eg.db)
#write.csv(gene.df,file="/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/Gene_body_hypo_67VPC_30CPCG.csv")


kegg.resu<- enrichKEGG(gene         = unique(annotation.promoter$annot.gene_id),
                       organism     = 'hsa',
                       keyType       = "ncbi-geneid",  pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1)
dotplot(kegg.resu,showCategory=10)+ggtitle("Promoter")+p_theme
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
temp <- as.data.frame(kegg.resu)
write.csv(temp,file=  "PE_16946_hyperDMR_VPC_CPCGENE_Promoter_KEGG.csv")
y <- setReadable(kegg.resu, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
y<-as.data.frame(y)
write.csv(y,file=  "PE_16946_hyperDMR_VPC_CPCGENE_Promoter_KEGG_symbol.csv")

#==============================================
#CPCisland
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#txdb
gene_region <-  genes(txdb,columns="gene_id")
gene_region$type <- "gene"

cpg.file <- read.table("/Users/wenbinye/Documents/study/methyaltion/data/hg10_bed_cogIsland.txt",header=FALSE,
                       stringsAsFactors = FALSE)
#str(cpg.file)
cpgiland_region <- GRanges(
  seqnames=as.character(cpg.file$V1),
  ranges=IRanges(start=as.integer(cpg.file$V2),end=as.integer(cpg.file$V3)),
  strand="*")
cpgiland_region$type <- "cpg_island"


annotation.island <- subset(annotated.need,annot.type=="CpG island")
#annotation.island <- subset(annotated,annot.type=="CpG island")
dmr.region <- GRanges.genome[as.numeric(annotation.island$index)]

ov.cpgisland_region <- findOverlaps(dmr.region,cpgiland_region)
cpgiland_region_overlap <- cpgiland_region[unique(ov.cpgisland_region@to)]
ov.island.gene <- distanceToNearest(cpgiland_region_overlap,gene_region)
str((ov.island.gene@elementMetadata$distance))
index <- which(ov.island.gene@elementMetadata$distance<=5000)
ov.island.gene <- ov.island.gene[index]
gene_region_overlap <- gene_region[unique(ov.island.gene@to)]
length(unique(gene_region_overlap$gene_id))


gene.df <- bitr(unique(gene_region_overlap$gene_id), fromType = "ENTREZID",
                toType = c("ENSEMBL", "ENTREZID","SYMBOL"),
                OrgDb = org.Hs.eg.db)
dim(gene.df)
#write.csv(gene.df,file="/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/CpGisland_hypo_Gene_ID_67VPC_30CPC.csv")

kegg.resu<- enrichKEGG(gene         = unique(gene_region_overlap$gene_id),
                       organism     = 'hsa',
                       keyType       = "ncbi-geneid",  pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.1)
dotplot(kegg.resu,showCategory=10)+ggtitle("CpG island")+p_theme
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
temp <- as.data.frame(kegg.resu)
write.csv(temp,file=  "PE_16946_hyperDMR_VPC_CPCGENE_CpGisland_KEGG.csv")
y <- setReadable(kegg.resu, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
y<-as.data.frame(y)
write.csv(y,file=  "PE_16946_hyperDMR_VPC_CPCGENE_Gisland_KEGG_symbol.csv")

#===============================================================
#' @title validate by SE data
#=======================
#' #PE
rpkm.table <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/Prostate_PE_16946Hyper_VPC_CPC_153.Rdata")
dim(rpkm.table)
rpkm.table <-log2(rpkm.table+1)
sample.names <- colnames(rpkm.table)
group <- gsub("\\d+","",sample.names)
group[grep("Norm|BL",sample.names)] <- "Control"
group[grep("sorted_V",sample.names)] <- "VPC"
group[grep("sorted_J|G",sample.names)] <- "CPCGENE"
group[grep("sorted_C",sample.names)] <- "Benign"
group[grep("sorted_R",sample.names)] <- "Barrier"
group[grep("bam_Counts.wig.txt",group)] <- "Tony"

index <- which(group=="Tony")
gsub("sorted_|_L001_qcs.bam_Counts.wig.txt|_S\\d+","", colnames(rpkm.table)[index])


index <- which(group=="Benign")
group <- group[-index]

rpkm.table<- rpkm.table[,-index]




#index <- which(group=="VPC")
#group <- group[index]
#rpkm.table<- rpkm.table[,index]

renew.names <- as.character(sapply(strsplit(colnames(rpkm.table),"\\_"),"[[",2))
colnames(rpkm.table)<- renew.names

# sample.v.combine<- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/DMR_Fre21_2020/VPC_ctDNA_information_PE.Rdata")
# index <- match(colnames(rpkm.table),sample.v.combine$names)
# identical(colnames(rpkm.table),sample.v.combine$names[index])
# ctDNA <- sample.v.combine$ctDNA[index]

table(group)
group <- factor(group,levels=c("Control","CPCGENE","Tony","Barrier","VPC"))


mean.rpkm <-apply(rpkm.table,2,mean)
mean.rpkm <- data.frame(mean.rpkm=mean.rpkm,ctDNA=ctDNA,group=group,sample=colnames(rpkm.table))
index <- match(mean.rpkm$sample,names(cluster))
mean.rpkm$cluster <- cluster[index]
mean.rpkm <- mean.rpkm[order(mean.rpkm$ctDNA),]

annotation_col <- data.frame(
  class=mean.rpkm$cluster,
  ctDNA=mean.rpkm$ctDNA
)
#order(data$group2)
rownames(annotation_col)<- mean.rpkm$sample
library(pheatmap)
data <-as.matrix( t(mean.rpkm[,1]) )
dim(data)
colnames(data) <- mean.rpkm$sample
pheatmap(as.matrix(data),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE,show_colnames =FALSE,
         # color=colorRampPalette(c("navy","white","firebrick3"))(50),
         annotation_col = annotation_col,
         fontsize=12,height = 0.5)

mean.rpkm$cluster <- as.character(mean.rpkm$cluster)
ggplot(mean.rpkm , 
       aes(x=ctDNA,y=mean.rpkm.mean.rpkm,color=group))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Mean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.9,0.9)) +geom_text_repel(
            data = subset(mean.rpkm, (mean.rpkm.mean.rpkm>4.8) ),
            aes(label=sample),
            size=4,
            box.padding=unit(0.35,"lines"),
            point.padding=unit(0.3,"lines")
          )#+stat_n_text(size=4)



pdf("Heatmap_V2_hypo.pdf",width=7,height=7)

dev.off()

#=======================
#SE
rpkm.table<-readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/SE_RPKM_table_169946_hyper_index_67VPC_30CPC.Rdata.Rdata")
rpkm.table <-log2(rpkm.table+1)
group <- gsub("\\d+","",colnames(rpkm.table))
table(group)
group <- factor(group,levels=c("H","J","JP","R","V"),
                labels=c("Control","CPCGENE","6-CPCGENE","Barrier","VPC"))
table(group)

library(RColorBrewer)
library(ggpubr)
library(EnvStats)
library(Rtsne)
library(ggrepel)
set.seed(42)
tsne <- Rtsne(as.matrix(t(rpkm.table)),
              dims = 3,check_duplicates = FALSE,
              perplexity = 30)  #for VPC cohort ,perplexity=22

tsne_projection <- as.data.frame(tsne$Y)
tsne_projection$sample <- colnames(rpkm.table)
tsne_projection$group <- group
#tsne_projection$type <- type
tsne_projection$shape <- 1
index <- which(tsne_projection$sample  %in% names(cluster)[which(cluster==2)])
tsne_projection$shape[index] <- 2


sample.v.combine<- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/DMR_Fre21_2020/VPC_ctDNA_information_PE.Rdata")
index <- match(tsne_projection$sample,sample.v.combine$names)
identical(tsne_projection$sample[-which(is.na(index))],sample.v.combine$names[index[-which(is.na(index))]])
tsne_projection$ctDNA <- sample.v.combine$ctDNA[index]
#tsne_projection$ctDNA <- ctDNA
colors <- brewer.pal(8,"Dark2")

library(wesanderson)
pal <- wes_palette("Zissou1",100,type="continuous")

tsne_projection$shape <- as.character(tsne_projection$shape)
tsne_projection$class <- tsne_projection$shape 
ggplot(tsne_projection,aes(V1,V2))+
  geom_point(size=5,aes(color=ctDNA,shape=class)) + labs(title=NULL,x="Dim1",y="Dim2")+#scale_color_gradientn(colours=pal)
  guides(color=guide_legend(title="%ctDNA"))+p_theme+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetica"),
          legend.position  ="top")+scale_color_gradientn(colours=pal)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggplot(tsne_projection,aes(V1,V2,color=ctDNA))+
  geom_point(size=4) + labs(title=NULL,x="Dim1",y="Dim2")+
  guides(color=guide_legend(title="%ctDNA"))+p_theme+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetica"),
          legend.position  ="top")+scale_color_gradientn(colours=pal)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
ggplot(tsne_projection,aes(V1,V3,color=ctDNA))+
  geom_point(size=4) + labs(title=NULL,x="Dim1",y="Dim3")+
  guides(color=guide_legend(title="%ctDNA"))+p_theme+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetica"),
          legend.position  ="top")+scale_color_gradientn(colours=pal)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
ggplot(tsne_projection,aes(V2,V3,color=ctDNA))+
  geom_point(size=4,shape=tsne_projection$shape) + labs(title=NULL,x="Dim2",y="Dim3")+
  guides(color=guide_legend(title="%ctDNA"))+p_theme+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetica"),
          legend.position  ="top")+scale_color_gradientn(colours=pal)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

tsne_projection$pcolor <-"black"
tsne_projection$pcolor[tsne_projection$ctDNA<=0]<- "#3399BB"
tsne_projection$pcolor[tsne_projection$ctDNA>0 & tsne_projection$ctDNA<=20] <- "#77BBCC"
tsne_projection$pcolor[tsne_projection$ctDNA>20 & tsne_projection$ctDNA<=40] <- "#DDCC44"
tsne_projection$pcolor[tsne_projection$ctDNA>40 & tsne_projection$ctDNA<=60] <- "#EEBB00"
tsne_projection$pcolor[tsne_projection$ctDNA>60] <-"#EE5500"

library(scatterplot3d)
scatterplot3d(x=tsne_projection$V1,y=tsne_projection$V2,z=tsne_projection$V3,
              color=tsne_projection$pcolor ,
              pch=19,
              xlab="Dim1",
              ylab="Dim2",
              zlab="Dim3",grid=TRUE)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



ggplot(tsne_projection,aes(V2,V3,color=group))+
  geom_point(size=4) + labs(title=NULL,x="Dim2",y="Dim3")+
  guides(color=guide_legend(title=NULL))+
  #scale_color_brewer(palette = "Dark2")+
  scale_color_manual(values=c(colors[c(1:3,6,7)]))+p_theme+
  #scale_color_manual(values=c(colors[1],colors[c(2:3,6:8)],"#999999","#3288BD" ,"#C6CDF7"))+p_theme+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetica"),
          legend.position  ="top")+
  geom_text_repel(
    data = subset(tsne_projection, ctDNA <=10),
    #aes(label=as.character(sapply(strsplit(sample,"\\_"),"[[",2)) ),
    aes(label=paste0(sample,": ",ctDNA)),
    size=4,
    box.padding=unit(0.35,"lines"),
    point.padding=unit(0.3,"lines")
  )
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

library(plotly)
p<-plot_ly(tsne_projection,
           x=~V1,
           y=~V2,
           z=~V3,
           color=~group,
           colors=c(colors[c(1:3,6,7)]),
           #colors=c(colors[1],"#0F644B",colors[c(2:3,6:8)],"#999999","#3288BD" ,"#C6CDF7"),
           #alpha = 1,
           type="scatter3d",
           mode="markers",
           text= paste(
             "</br> ID: ", tsne_projection$sample,
             "</br> group: ", tsne_projection$group,
             "</br> %ctDNA", tsne_projection$ctDNA
           ))
p<-plot_ly(tsne_projection,
           x=~V1,
           y=~V2,
           z=~V3,
           color=~ctDNA,
           #color=~group,
           #colors=pal,
           #colors=c(colors[1],"#0F644B",colors[c(2:3,6:8)],"#999999","#3288BD" ,"#C6CDF7"),
           #alpha = 1,
           type="scatter3d",
           mode="markers",
           text= paste(
             "</br> ID: ", tsne_projection$sample,
             "</br> group: ", tsne_projection$group,
             "</br> %ctDNA", tsne_projection$ctDNA
           ))


t <- list(family="Helvetica",
          size=14,
          color="black")
p1 <- p %>% layout(font=t,
                   scene=list(
                     xaxis=list(title="Dim1",showticklabels=FALSE),
                     yaxis=list(title="Dim2",showticklabels=FALSE),
                     zaxis=list(title="Dim3",showticklabels=FALSE)
                   )) %>% 
  add_trace(  x=tsne_projection$V1,
              y=tsne_projection$V2,
              z=tsne_projection$V3,marker=list(size=1),showlegend=F)
print(p1)
setwd(out.path)
htmlwidgets::saveWidget(p1,"PE_TSNE_based_on_16946_PE_67VPC_30CPCGENE_hyperDMR_total135cases.html")
htmlwidgets::saveWidget(p1,"SE_TSNE_based_on_16946_PE_67VPC_30CPCGENE_hyperDMR_total114cases.html")


#======================================================
#==========================================================
#optional  STAR
library(factoextra)
res.hc <- eclust(t(rpkm.table),"hclust") 
fviz_dend(res.hc,rect=TRUE)
dd <- dist(t(rpkm.table),method="euclidean")

hc <- hclust(dd,method = "ward.D2")
hcd <- as.dendrogram(hc)
plot(hc,hang=01,cex=0.6)
plot(hcd,type="rectangle",ylab="Height")
cluster <- cutree(hc,k=2)
fviz_cluster(list(data=t(rpkm.table),cluster=cluster))+p_theme
library(ggdendro)
library(dendextend)
ggdendrogram(hc)
# Build dendrogram object from hclust results
dend <- as.dendrogram(hc)
# Extract the data (for rectangular lines)
# Type can be "rectangle" or "triangle"
dend_data <- dendro_data(dend, type = "rectangle")
# What contains dend_data
names(dend_data)
head(dend_data$segments)
head(dend_data$labels)

p <- ggplot(dend_data$segments) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = label),
            hjust = 1, angle = 90, size = 3)#+
#ylim(-3, 15)
print(p+p_theme)

dend %>% plot 

labels(dend)
colnames(rpkm.table)[index]

index <-match(labels(dend),colnames(rpkm.table))

pdf("VPC_HC_ward2D_based_16946_hyper_Mes_loc_ctDNA.pdf")
dend %>% set("labels",round(ctDNA[index],0) ) %>% 
  set("labels_col", value = c("#1B9E77", "#D95F02"), k=2) %>%
  set("labels_cex", 0.7) %>%
  set("branches_k_color",value=c("#1B9E77", "#D95F02"), k = 2) %>% plot(horiz = TRUE)
dev.off()

pdf("VPC_HC_ward2D_based_16946_hyper_Mes_loc_sampleID.pdf")
dend %>% 
  set("labels_col", value = c("#1B9E77", "#D95F02"), k=2) %>%
  set("labels_cex", 0.7) %>%
  set("branches_k_color",value=c("#1B9E77", "#D95F02"), k = 2) %>% plot(horiz = TRUE)
dev.off()

pdf("VPC_HC_ward2D_based_16946_hyper_Mes_loc_sampleID_ctDNA.pdf")
dend %>% set("labels",paste0(labels(dend),": ",round(ctDNA[index],0) )) %>% 
  set("labels_col", value = c("#1B9E77", "#D95F02"), k=2) %>%
  set("labels_cex", 0.6) %>%
  set("branches_k_color",value=c("#1B9E77", "#D95F02"), k = 2) %>% plot(horiz = FALSE)
dev.off()
ggd1 <- as.ggdend(dend)
ggplot(ggd1,horiz = TRUE, theme = NULL) 
#optional  END
#==========================================================

#Survival
options(stringsAsFactors =FALSE )
library(readxl)
sample.v <- read_excel("/Users/wenbinye/Dropbox/Jessica/WenbinYE/Data_information/Single_End_5mc/Batches_and_ALL_QCs_MeDIP.xlsx",
                       sheet="VPC Cohort")
sample.v <- sample.v[grep("V",sample.v $`Sequencing ID`),]
sample.v$`Sequencing ID` <- gsub("\\s+\\(33\\)","-33",sample.v$`Sequencing ID`)
index <-match(mean.rpkm$sample,sample.v$`Sequencing ID`)

mean.rpkm$sample[c(which(is.na(index)))]
mean.rpkm$patient <- sample.v$`patient number`[index]
mean.rpkm$tube <- sample.v$`Tube #`[index]
mean.rpkm$ctDNA_cp <-sample.v$`ctDNA%`[index]

#index <-match(sample.v$`Sequencing ID`,mean.rpkm$sample)
#sample.v$`Sequencing ID`[c(which(is.na(index)))]
#match(mean.rpkm$sample,sample.v$`Sequencing ID`)
v.clinical <- read_excel("/Users/wenbinye/Dropbox/Jessica/WenbinYE/Data_information/Single_End_5mc/VPC\ Cohort\ Clinical\ Data.xlsx",
                         sheet = "Table S1 - Cohort",skip=1)
index<- match(gsub("^0+","",mean.rpkm$patient), gsub("^0+","",v.clinical$Patient))

mean.rpkm <-cbind(mean.rpkm,v.clinical[index,])
#write.csv(mean.rpkm,file="/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/VPC_cohort_complete_clinical.csv")


#identical(mean.table.vpc$sample.ID,mean.table.hypo.vpc$sample.ID)
mean.rpkm <- mean.table.vpc
mean.rpkm$mean.rpkm <- as.numeric(mean.table.vpc$mean.rpkm)/as.numeric(mean.table.hypo.vpc$mean.rpkm)
#mean.rpkm$mean.rpkm <- as.numeric(mean.table.vpc$mean.rpkm)
#mean.rpkm.raw <- mean.rpkm
table(mean.rpkm$type)
mean.rpkm <- subset(mean.rpkm,type=="age.67vpc.30cpc") 
table(mean.table.vpc$type)
#mean.rpkm <- subset(mean.table.hypo.vpc,type=="age.67vpc.30cpc") 

colnames(mean.rpkm) <- gsub("Published_","",colnames(mean.rpkm))  
table(mean.rpkm$Timepoint)
mean.rpkm$Timepoint <- factor(mean.rpkm$Timepoint,levels=c("Baseline","1st progression","2nd progression"))
mean.rpkm <- mean.rpkm[order(mean.rpkm$Timepoint,mean.rpkm$`patient number`),]
mean.rpkm  <- mean.rpkm[!duplicated(mean.rpkm$`patient number`),]

mean.rpkm$OS.time <- as.numeric(mean.rpkm$`Days to death or last followup`)
mean.rpkm$OS <- mean.rpkm$Censored...14
mean.rpkm$OS[which(mean.rpkm$OS=="Yes")] <- 1
mean.rpkm$OS[which(mean.rpkm$OS=="No")] <- 2
mean.rpkm$OS <- as.numeric(mean.rpkm$OS)
mean.rpkm$RFS.time <- as.numeric(mean.rpkm$`Days to progression or last followup`)
mean.rpkm$RFS <- mean.rpkm$Censored...12
table(mean.rpkm$RFS )
mean.rpkm$RFS[which(mean.rpkm$RFS=="Yes")] <- 1
mean.rpkm$RFS[which(mean.rpkm$RFS=="No")] <- 2
mean.rpkm$RFS<-as.numeric(mean.rpkm$RFS)

mean.rpkm$cluster <- "1"
mean.rpkm$cluster[which(mean.rpkm$mean.rpkm>=median(mean.rpkm$mean.rpkm))] <- "2"
table(mean.rpkm$cluster)

OS.fit <- survfit(Surv(OS.time, OS) ~ cluster, data = mean.rpkm)
RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluster, data = mean.rpkm)

#ggsurvplot(OS.fit)
library(survminer)
ggsurvplot(OS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("High risk","Low risk"),
           ggtheme = p_theme# Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggsurvplot(RFS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progression-free survival probability",
           legend.labs=c("High risk","Low risk"),
           ggtheme = p_theme# Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



#temp <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/VPC_mean_RPKM_based_on_16946_hyperDMR_VPC_vs_CPC.Rdata")

#=======================================================================
#' @title adding 1930 hypo sites
#'  

pe_hypo <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/Cancer_Type/PE_Hypo_1930/Prostate_PE_1930Hypo_VPC_CPC_153.Rdata")
pe_hyper <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/Cancer_Type/PE_Hyper_16946/Prostate_PE_16946Hyper_VPC_CPC_153.Rdata")
sample.names <- colnames(pe_hypo)
group <- gsub("\\d+","",sample.names)
group[grep("Norm|BL",sample.names)] <- "Control"
group[grep("sorted_V",sample.names)] <- "VPC"
group[grep("sorted_J|G",sample.names)] <- "CPCGENE"
group[grep("sorted_C",sample.names)] <- "Benign"
group[grep("sorted_R",sample.names)] <- "Barrier"
group[grep("bam_Counts.wig.txt",group)] <- "Tony"
table(group)
hypo.mean <- apply(pe_hypo,2,mean)
hypo.table <- data.frame(mean.hypo=as.numeric(hypo.mean),group=group,sample=sample.names)
hypo.table$sampleID <-as.character(sapply(strsplit(hypo.table$sample,"\\_"),"[[",2))
summarise(group_by(hypo.table,group),mean(mean.hypo))
index <- match(mean.rpkm$sample,hypo.table$sampleID)
identical(mean.rpkm$sample,hypo.table$sampleID[index])


mean.rpkm$mean.hypo.rpkm <- hypo.table$mean.hypo[index]
#mean.rpkm$mean.hyper.rpkm.check <- hypo.table$mean.hypo[index]
identical(mean.rpkm$mean.hyper.rpkm.check,mean.rpkm$mean.hyper.rpkm)

colnames(mean.rpkm)[1:3] <- c("mean.hyper.rpkm","JessctDNA","Group")
#saveRDS(mean.rpkm,file="/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/67VPC_clinical_add_16946hyper_1930hypo.Rdata")
mean.rpkm <- readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/67VPC_clinical_add_16946hyper_1930hypo.Rdata")

library(readxl)
v.publised <- read_excel("/Users/wenbinye/Dropbox/Wenbin/Paired_End/Supplement.xlsx",
                         sheet = "Table S1 - Cohort",skip=1)
table(v.publised$`Genomic markers`)

index <- match(mean.rpkm$Patient,v.publised$Patient)
identical(mean.rpkm$Patient,v.publised$Patient[index])
identical(mean.rpkm$`LDH / ULN`,v.publised$`LDH / ULN`[index])
mean.rpkm$Genomic.markers <- v.publised$`Genomic markers`[index]
#saveRDS(mean.rpkm,file="/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/67VPC_clinical_add_16946hyper_1930hypo.Rdata")
temp <-readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/67VPC_clinical_add_16946hyper_1930hypo.Rdata")


str(mean.rpkm$mean.hyper.rpkm)
str(mean.rpkm$mean.hypo.rpkm)
cor.test(as.numeric(mean.rpkm$mean.hyper.rpkm),as.numeric(mean.rpkm$mean.hypo.rpkm),method=c("spearman"))
ggplot(mean.rpkm , 
       aes(x=ctDNA,y=mean.hypo.rpkm,color=group))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Hypo-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.9,0.9)) # +geom_text_repel(
#data = subset(mean.rpkm, (mean.hypo.rpkm>30) ),
#aes(label=sample),
#size=4,
#box.padding=unit(0.35,"lines"),
#point.padding=unit(0.3,"lines")
#)#+stat_n_text(size=4)
ggplot(mean.rpkm , 
       aes(x=ctDNA,y=mean.hyper.rpkm,color=group))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="% ctDNA",y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.9,0.9))

out.path <-"/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure"
out.figure <- "20200513_PE_DMRs.pptx"
ggplot(mean.rpkm , 
       aes(x=mean.hypo.rpkm,y=mean.hyper.rpkm,color=group))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Hypo-DMRs\nMean RPKM counts across sample",y="Hyper-DMRs\nMean RPKM counts across sample")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.9,0.9))+ylim(0,8)

ggplot(mean.rpkm , 
       aes(x=Age,y=ctDNA,color=group))+geom_point()+theme_classic()+
  p_theme+scale_color_manual(values=colors[c(7)])+geom_smooth(method=lm,se=TRUE,fullrange=FALSE)+
  stat_cor(method = "pearson")+
  labs(x="Age",y="% ctDNA")+guides(color=FALSE)+
  theme(  legend.text =element_text(color="black",size=12,family="Helvetica"),
          legend.title=element_text(color="black",size=14,family="Helvetiva"),
          axis.text.x = element_text(angle = 0,vjust=0.6),
          strip.text.x = element_text(size=12),
          legend.background = element_blank(),
          legend.position = c(0.9,0.9))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

cor.test(as.numeric(mean.rpkm$mean.hypo.rpkm),as.numeric(mean.rpkm$mean.hyper.rpkm),method=c("spearman"))


#=============================================
#Age
summary(mean.rpkm$ctDNA)

#raw.67vpc.30cpc age.67vpc.30cpc age.16vpc.16cpc
mean.rpkm <- subset(mean.rpkm.renew,type=="age.16vpc.16cpc")
table(mean.rpkm$type)
mean.rpkm$ctDNA <- as.numeric(mean.rpkm$ctDNA)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$Published_Age),method=c("spearman"))
summary(as.numeric(mean.rpkm$Published_Age))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$Published_Age<=70)],
            mean.rpkm$hyper.mean[which(mean.rpkm$Published_Age>70)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$Published_Age),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$Published_Age<=70)],
            mean.rpkm$hypo.mean[which(mean.rpkm$Published_Age>70)],
            alternative = c("greater"))

cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$Published_Age),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$Published_Age<=70)],
            mean.rpkm$ratio[which(mean.rpkm$Published_Age>70)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$Published_Age),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$Published_Age<=70)],
            mean.rpkm$ctDNA[which(mean.rpkm$Published_Age>70)],
            alternative = c("less"))

#=========================================
#PSA
summary(mean.rpkm$PSA)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$Published_PSA),method=c("pearson"))
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$Published_PSA),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$Published_PSA<=20)],
            mean.rpkm$hyper.mean[which(mean.rpkm$Published_PSA>20)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$Published_PSA),method=c("pearson"))
cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$Published_PSA),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$Published_PSA<=20)],
            mean.rpkm$hypo.mean[which(mean.rpkm$Published_PSA>20)],
            alternative = c("greater"))


cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$Published_PSA),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$Published_PSA<=20)],
            mean.rpkm$ratio[which(mean.rpkm$Published_PSA>20)],
            alternative = c("less"))
#summary(mean.rpkm$Published_PSA)

cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$Published_PSA),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$Published_PSA<=20)],
            mean.rpkm$ctDNA[which(mean.rpkm$Published_PSA>20)],
            alternative = c("less"))



#====================================
#LDH
summary(mean.rpkm$`Published_LDH / ULN`)
str(mean.rpkm$`Published_LDH / ULN`)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$`Published_LDH / ULN`),method=c("pearson"))
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$`Published_LDH / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_LDH / ULN`<=1)],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_LDH / ULN`>1)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$`Published_LDH / ULN`)                                      ,method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_LDH / ULN`<=1)],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_LDH / ULN`>1)],
            alternative = c("greater"))

cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$`Published_LDH / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_LDH / ULN`<=1)],
            mean.rpkm$ratio[which(mean.rpkm$`Published_LDH / ULN`>1)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$`Published_LDH / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_LDH / ULN`<=1)],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_LDH / ULN`>1)],
            alternative = c("less"))
#====================================================
#ALP
summary(mean.rpkm$`Published_ALP / ULN`)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$`Published_ALP / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ALP / ULN`<=1)],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ALP / ULN`>1)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$`Published_ALP / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_ALP / ULN`<=1)],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_ALP / ULN`>1)],
            alternative = c("greater"))

cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$`Published_ALP / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_ALP / ULN`<=1)],
            mean.rpkm$ratio[which(mean.rpkm$`Published_ALP / ULN`>1)],
            alternative = c("less"))


cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$`Published_ALP / ULN`),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_ALP / ULN`<=1)],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_ALP / ULN`>1)],
            alternative = c("less"))



#=========================================
#cfDNA
summary(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`<=45)],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`>45)],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`<=45)],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`>45)],
            alternative = c("greater"))

cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`<=45)],
            mean.rpkm$ratio[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`>45)],
            alternative = c("less"))


cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`<=45)],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_cfDNA yield (ng / mL plasma)`>45)],
            alternative = c("less"))

#========================================
#ctDNA, mean hyper, mean hypo
summary(mean.rpkm$ctDNA)
mean.rpkm$`ctDNA%` <- as.numeric(mean.rpkm$`ctDNA%`)
summary(mean.rpkm$`ctDNA%`)
mean.rpkm$`Published_ctDNA%` <- as.numeric(mean.rpkm$`Published_ctDNA%`)
summary(mean.rpkm$`Published_ctDNA%`)

cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$ctDNA),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$ctDNA<=38)],
            mean.rpkm$hyper.mean[which(mean.rpkm$ctDNA>38)],
            alternative = c("less"))
#wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$ctDNA<=43)],
#            mean.rpkm$hyper.mean[which(mean.rpkm$ctDNA>43)],
#            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$ctDNA),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$ctDNA<=38)],
            mean.rpkm$hypo.mean[which(mean.rpkm$ctDNA>38)],
            alternative = c("greater"))


cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$ctDNA),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$ctDNA<=38)],
            mean.rpkm$ratio[which(mean.rpkm$ctDNA>38)],
            alternative = c("less"))


cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$hypo.mean),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$hypo.mean>median(mean.rpkm$hypo.mean))],
            mean.rpkm$hyper.mean[which(mean.rpkm$hypo.mean<=median(mean.rpkm$hypo.mean))],
            alternative = c("less"))

median(mean.rpkm$ratio)
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$ratio),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$ratio<=median(mean.rpkm$ratio))],
            mean.rpkm$hyper.mean[which(mean.rpkm$ratio>median(mean.rpkm$ratio))],
            alternative = c("less"))


cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$hyper.mean),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$hyper.mean<median(mean.rpkm$hyper.mean))],
            mean.rpkm$hypo.mean[which(mean.rpkm$hyper.mean>=median(mean.rpkm$hyper.mean))],
            alternative = c("greater"))

median(mean.rpkm$ratio)
cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$ratio),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$ratio<=median(mean.rpkm$ratio))],
            mean.rpkm$hypo.mean[which(mean.rpkm$ratio>median(mean.rpkm$ratio))],
            alternative = c("greater"))


cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$hyper.mean),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$hyper.mean<=median(mean.rpkm$hyper.mean))],
            mean.rpkm$ratio[which(mean.rpkm$hyper.mean>median(mean.rpkm$hyper.mean))],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$hypo.mean),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$hypo.mean>median(mean.rpkm$hypo.mean))],
            mean.rpkm$ratio[which(mean.rpkm$hypo.mean<=median(mean.rpkm$hypo.mean))],
            alternative = c("less"))


cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$hyper.mean),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$hyper.mean<=median(mean.rpkm$hyper.mean))],
            mean.rpkm$ctDNA[which(mean.rpkm$hyper.mean>median(mean.rpkm$hyper.mean))],
            alternative = c("less"))
cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$ratio),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$ratio<=median(mean.rpkm$ratio))],
            mean.rpkm$ctDNA[which(mean.rpkm$ratio>median(mean.rpkm$ratio))],
            alternative = c("less"))
cor.test(as.numeric(mean.rpkm$ctDNA),as.numeric(mean.rpkm$hypo.mean),method=c("spearman"))
wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$hypo.mean>median(mean.rpkm$hypo.mean))],
            mean.rpkm$ctDNA[which(mean.rpkm$hypo.mean<=median(mean.rpkm$hypo.mean))],
            alternative = c("less"))

#-----------------------------------
#published  %ctDNA
cor.test(as.numeric(mean.rpkm$hyper.mean),as.numeric(mean.rpkm$`Published_ctDNA%`),method=c("spearman"))
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ctDNA%`<=38)],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ctDNA%`>38)],
            alternative = c("less"))

wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ctDNA%`<=median(mean.rpkm$`Published_ctDNA%`))],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_ctDNA%`>median(mean.rpkm$`Published_ctDNA%`))],
            alternative = c("less"))

cor.test(as.numeric(mean.rpkm$hypo.mean),as.numeric(mean.rpkm$`Published_ctDNA%`),method=c("spearman"))
wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_ctDNA%`<=median(mean.rpkm$`Published_ctDNA%`))],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_ctDNA%`>median(mean.rpkm$`Published_ctDNA%`))],
            alternative = c("greater"))


cor.test(as.numeric(mean.rpkm$ratio),as.numeric(mean.rpkm$`Published_ctDNA%`),method=c("spearman"))
wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_ctDNA%`<=median(mean.rpkm$`Published_ctDNA%`))],
            mean.rpkm$ratio[which(mean.rpkm$`Published_ctDNA%`>median(mean.rpkm$`Published_ctDNA%`))],
            alternative = c("less"))




#==================================================
#Bone
table(mean.rpkm$`Published_Bone metastases`)
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Bone metastases`=="No")],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Bone metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_Bone metastases`=="No")],
            mean.rpkm$ratio[which(mean.rpkm$`Published_Bone metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_Bone metastases`=="No")],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_Bone metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Bone metastases`=="No")],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Bone metastases`=="Yes")],
            alternative = c("greater"))

#===================================
table(mean.rpkm$`Published_Lung metastases`)
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Lung metastases`=="No")],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Lung metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_Lung metastases`=="No")],
            mean.rpkm$ratio[which(mean.rpkm$`Published_Lung metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_Lung metastases`=="No")],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_Lung metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Lung metastases`=="No")],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Lung metastases`=="Yes")],
            alternative = c("greater"))

#=============================================
table(mean.rpkm$`Published_Liver metastases`)
wilcox.test(mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Liver metastases`=="No")],
            mean.rpkm$hyper.mean[which(mean.rpkm$`Published_Liver metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ratio[which(mean.rpkm$`Published_Liver metastases`=="No")],
            mean.rpkm$ratio[which(mean.rpkm$`Published_Liver metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$ctDNA[which(mean.rpkm$`Published_Liver metastases`=="No")],
            mean.rpkm$ctDNA[which(mean.rpkm$`Published_Liver metastases`=="Yes")],
            alternative = c("less"))

wilcox.test(mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Liver metastases`=="No")],
            mean.rpkm$hypo.mean[which(mean.rpkm$`Published_Liver metastases`=="Yes")],
            alternative = c("greater"))

#===================================================================================
#Adding Genomic markers
temp <-readRDS("/Users/wenbinye/Documents/study/methyaltion/result/Rdata/PE_DMR/67VPC_clinical_add_16946hyper_1930hypo.Rdata")
index <- match(mean.rpkm$sample.ID,temp$sample)
mean.rpkm$Genomic.markers <- temp$Genomic.markers[index]
table(mean.rpkm$Genomic.markers )

mean.rpkm.sub <- subset(mean.rpkm,Genomic.markers=="1")
dim(mean.rpkm.sub)
mean.rpkm.sub$`Published_Germline HRR defect`[grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)]


wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_Germline HRR defect`)],
            alternative = c("greater"))
#==========================

mean.rpkm.sub$`Published_Somatic truncating HRR defect`[grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)]


wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_Somatic truncating HRR defect`)],
            alternative = c("greater"))

#==========================
mean.rpkm.sub$`Published_TP53 defect`[grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)]
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_TP53 defect`)],
            alternative = c("greater"))


#==========================
table(mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`[grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)])
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_BRCA2 / ATM deletion (mono- or biallelic)`)],
            alternative = c("greater"))


#==========================
table(mean.rpkm.sub$`Published_AR rearrangement`[grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)])
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_AR rearrangement`)],
            alternative = c("greater"))



#==========================
table(mean.rpkm.sub$`Published_AR mutation(s)`[grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)])
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_AR mutation(s)`)],
            alternative = c("greater"))


#==========================
mean.rpkm.sub$`Published_PI3K pathway defect(s)`[grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)]
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            alternative = c("greater"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_PI3K pathway defect(s)`)],
            alternative = c("less"))
#=============================================
mean.rpkm.sub$`Published_WNT pathway defect(s)`[grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)]
wilcox.test(mean.rpkm.sub$hyper.mean[-grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            mean.rpkm.sub$hyper.mean[grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$ratio[-grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            mean.rpkm.sub$ratio[grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            alternative = c("less"))
wilcox.test(mean.rpkm.sub$hypo.mean[-grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            mean.rpkm.sub$hypo.mean[grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            alternative = c("greater"))
wilcox.test(mean.rpkm.sub$ctDNA[-grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            mean.rpkm.sub$ctDNA[grep("\\S+",mean.rpkm.sub$`Published_WNT pathway defect(s)`)],
            alternative = c("less"))


#================================================================




feature <- read.table("/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/67VPC_DMR_with_Clinical.txt",
                      sep="\t",stringsAsFactors = FALSE,header = TRUE)

feature <- read.table("/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/60_out_67VPC_genomic.txt",
                      sep="\t",stringsAsFactors = FALSE,header = TRUE)

feature <- read.table("/Users/wenbinye/Documents/study/methyaltion/result/20200503_Figure/55_VPC_Control_AGE_Comparision.txt",
                      sep="\t",stringsAsFactors = FALSE,header = TRUE)
library(reshape2)
feature.table <- reshape2::melt(feature)
feature.table$variable <- gsub("X.ctDNA","%ctDNA",feature.table$variable)
feature.table$value <- -log10(feature.table$value)
table(feature.table$Feature)
feature.table$Feature <- factor(feature.table$Feature,levels =feature$Feature )
table(feature.table$variable)
feature.table$group <- as.character(sapply(strsplit(feature.table$variable,"\\."),"[[",1))
feature.table$variable1 <-gsub("A.|B.|C.|.1","",feature.table$variable)
table(feature.table$group)
feature.table$group <- factor(feature.table$group,levels = c("%ctDNA","A","B","C"))

table(feature.table$variable1)
feature.table$variable1 <- factor(feature.table$variable1,levels = c("%ctDNA"  ,"Hyper","Hypo",  "Ratio"))

unique(feature.table$Feature)

colors <- brewer.pal(8,"Dark2")
feature.table.new <- subset(feature.table,group %in%  c("%ctDNA","B","C"))
ggplot(subset(feature.table.new,Feature %in%  as.character(unique(feature.table.new$Feature)[25:29]) ),
       aes(x=group,y=value,fill=variable1))+scale_fill_manual(values=colors[c(1:3,6)])+
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(value,2)),vjust=1.6,
            position = position_dodge(0.9),size=4)+p_theme+
  geom_hline(yintercept=-log10(0.05),linetype="dotted",color="black")+
  facet_wrap(~Feature, scales="free",ncol=3)+labs(x=NULL,y="-Log10 p-value")+
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(size=12),legend.position ="top")+
  guides(fill=guide_legend(title=NULL))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)







ggplot(subset(feature.table,Feature %in%  as.character(unique(feature.table$Feature)[28:29]) ),
       aes(x=group,y=value,fill=variable1))+scale_fill_manual(values=colors[c(1:3,6)])+
  geom_bar(stat="identity",position=position_dodge())+
  geom_text(aes(label=round(value,2)),vjust=1.6,
            position = position_dodge(0.9),size=4)+p_theme+
  geom_hline(yintercept=-log10(0.05),linetype="dotted",color="black")+
  facet_wrap(~Feature, scales="free")+labs(x=NULL,y="-Log10 p-value")+
  theme(#axis.text.x = element_blank(),
    strip.text = element_text(size=16),legend.position ="top")+
  guides(fill=guide_legend(title=NULL))
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

#==========================================================================
#plotting survival curves
table(mean.rpkm$Censored...12)
table(mean.rpkm$RFS)
mean.rpkm$RFS[mean.rpkm$Censored...12=="Yes"]<-1
mean.rpkm$RFS[mean.rpkm$Censored...12=="No"]<-2

table(mean.rpkm$OS)
mean.rpkm$OS[mean.rpkm$Censored...14=="Yes"]<-1
mean.rpkm$OS[mean.rpkm$Censored...14=="No"]<-2


OS.fit <- survfit(Surv(OS.time, OS) ~ cluster, data = mean.rpkm)
RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluster, data = mean.rpkm)

summary(mean.rpkm$mean.hyper.rpkm)
mean.rpkm$hyper.status <- 1
mean.rpkm$hyper.status[which(mean.rpkm$mean.hyper.rpkm>=1)] <- 2
OS.fit.hyper <- survfit(Surv(OS.time, OS) ~ hyper.status, data = mean.rpkm)


summary(mean.rpkm$mean.hypo.rpkm)
mean.rpkm$hypo.status <- 1
mean.rpkm$hypo.status[which(mean.rpkm$mean.hypo.rpkm<12.53745)] <- 2
table(mean.rpkm$hypo.status)
OS.fit.hypo<- survfit(Surv(OS.time, OS) ~ hypo.status, data = mean.rpkm)

summary(mean.rpkm$ctDNA)
mean.rpkm$ctDNA.status <-1
mean.rpkm$ctDNA.status[which(mean.rpkm$ctDNA>=median(mean.rpkm$ctDNA))]<-2
table(mean.rpkm$ctDNA.status)
OS.fit.ctDNA<- survfit(Surv(OS.time, OS) ~ ctDNA.status, data = mean.rpkm)


RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluster, data = mean.rpkm)
OS.fit.hyper <- survfit(Surv(OS.time, OS) ~ hyper.status, data = mean.rpkm)
OS.fit.hypo<- survfit(Surv(OS.time, OS) ~ hypo.status, data = mean.rpkm)
OS.fit.ctDNA<- survfit(Surv(OS.time, OS) ~ ctDNA.status, data = mean.rpkm)
#ggsurvplot(OS.fit)
library(survminer)
ggsurvplot(OS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("Cluster1","Cluster2"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

ggsurvplot(OS.fit.hyper, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("Hyper low","Hyper high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggsurvplot(OS.fit.hypo, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("Hypo low","Hypo high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)

ggsurvplot(OS.fit.ctDNA, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Overall survival probability",
           legend.labs=c("%ctDNA low","%ctDNA high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)





RFS.fit <- survfit(Surv(RFS.time, RFS) ~ cluster, data = mean.rpkm)
RFS.fit.hyper <- survfit(Surv(RFS.time, RFS) ~ hyper.status, data = mean.rpkm)
RFS.fit.hypo<- survfit(Surv(RFS.time, RFS) ~ hypo.status, data = mean.rpkm)
RFS.fit.ctDNA<- survfit(Surv(RFS.time, RFS) ~ ctDNA.status, data = mean.rpkm)

ggsurvplot(RFS.fit, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progresstion-free survival probability",
           legend.labs=c("Cluster1","Cluster2"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)



ggsurvplot(RFS.fit.hyper, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progresstion-free survival probability",
           legend.labs=c("Hyper low","Hyper high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggsurvplot(RFS.fit.hypo, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progresstion-free survival probability",
           legend.labs=c("Hypo low","Hypo high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)


ggsurvplot(RFS.fit.ctDNA, size = 1, # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palette
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           ylab="Progresstion-free survival probability",
           legend.labs=c("%ctDNA low","%ctDNA high"),
           ggtheme = p_theme # Change ggplot2 theme
)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)






mean.rpkm$cluster <- as.character(mean.rpkm$cluster)

mean.rpkm <- mean.rpkm[order(mean.rpkm$mean.rpkm),]

annotation_col <- data.frame(
  Lung_met=mean.rpkm$`Lung metastases`,
  Bone_met=mean.rpkm$`Bone metastases`,
  Liver_met =mean.rpkm$`Bone metastases`,
  cfDNA_yield=mean.rpkm$`cfDNA yield (ng / mL plasma)`,
  ALP=mean.rpkm$`ALP / ULN`,
  LDH = mean.rpkm$`LDH / ULN`,
  ctDNA=mean.rpkm$ctDNA,
  meanRPKM=mean.rpkm$mean.rpkm,
  class=mean.rpkm$cluster,
  age =mean.rpkm$Age
  
  
)
annotation_col <- data.frame(
  Lung_met=mean.rpkm$`Lung metastases`,
  Bone_met=mean.rpkm$`Bone metastases`,
  Liver_met =mean.rpkm$`Liver metastases`,
  cfDNA_yield=mean.rpkm$`cfDNA yield (ng / mL plasma)`,
  ALP=mean.rpkm$`ALP / ULN`
  #LDH = mean.rpkm$`LDH / ULN`,
  #ctDNA=mean.rpkm$ctDNA,
  #meanRPKM=mean.rpkm$mean.rpkm,
  #class=mean.rpkm$cluster,
  #age =mean.rpkm$Age
)
#order(data$group2)
rownames(annotation_col)<- mean.rpkm$sample
library(pheatmap)
data <-as.matrix( t(mean.rpkm[,1]) )
dim(data)
colnames(data) <- mean.rpkm$sample
pheatmap(as.matrix(data),cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE,show_colnames =FALSE,
         # color=colorRampPalette(c("navy","white","firebrick3"))(50),
         annotation_col = annotation_col,
         fontsize=12,height = 0.5)
setwd(out.path)
graph2ppt(file=out.figure, append=TRUE)
```
