library("ggcorrplot")
library("xlsx")
library("ggplot2")
library("ggpubr")
library("gridExtra")
library("RColorBrewer")
library("patchwork")
library("pROC")

setwd("SHP2/ROC")
RC_mat <- read.xlsx("SHP2_data.xlsx",1)
head(RC_mat)
#   NUM POS SEQ protein_change cDNA_change      Type Domain Score        Si         MI       ΔΔG      RASA
# 1   1   2   T          p.T2I      c.5C>T RASopathy   <NA>     4 1.7203559 0.19447303  0.178354 95.007500
# 2   2  10   N         p.N10D     c.28A>G RASopathy  SH2_1     6 1.4419147 0.31591289  4.353530 13.753700
# 3   3  12   T         p.T12A     c.34A>G RASopathy  SH2_1     7 1.1196461 0.15424828  0.475647 36.936600
# 4   4  18   N         p.N18S     c.53A>G RASopathy  SH2_1     3 2.1819636 0.34615723 -0.176573 54.408800
# 5   5  35   K         p.K35E    c.103A>G RASopathy  SH2_1     5 1.3564694 0.26279907 -1.142700 51.944500
# 6   6  43   L         p.L43F    c.127C>T RASopathy  SH2_1     8 0.1726519 0.05485857  2.684760  0.141452
#       Degree Betweenness  Closeness       MSF Effectiveness Sensitivity
# 1 -0.0111473  0.05047757  0.6266855  9.543178     32.610301   0.2850602
# 2 -0.0111473  0.05047757 -0.7424476 15.328404     21.482988   0.5373204
# 3 -0.0111473  0.04966134 -0.7442245 21.566318     12.074674   0.8247089
# 4 -0.0111473  0.03602627 -0.3815926 29.140534      8.647312   1.0685560
# 5 -0.0111473  0.15392553 -1.5556249 35.624986      6.976254   1.2440535
# 6 -0.0111473 -0.22747229 -1.4904131 28.130999      9.677493   0.9818018


#相关性
M <- cor(RC_mat[,8:18],method = "spearman")

p1 <- ggcorrplot(M,
           method = "square", #设置相关性图展示类型
           show.legend = T, #设置是否展示图例
           legend.title = "Corr", #设置图例的标题
           colors = c("#6D9EC1", "white", "#E46726"), #设置相关性图的颜色
           ggtheme = ggplot2::theme_gray, #设置背景
           lab = T, #设置是否显示显关系数
           hc.order = T #设置排序
           )+
        labs(title = "Correlation of the Parameters",tags = "A")

#建模
set.seed(9123)

train_row <- sample(nrow(RC_mat), 7/10*nrow(RC_mat))
train_data <- RC_mat[train_row,]
test_data <- RC_mat[-train_row,]
RC_mat$Type <- as.factor(RC_mat$Type)

#分层抽样
library(sampling)
table(RC_mat)
# CANCER RASopathy 
#    181       127 
#按3:7分为测试集和训练集
train <- strata(RC_mat, stratanames = "Type", size = c(89,127), method = "srswor")
data_train = RC_mat[train$ID_unit,]
data_test = RC_mat[-train$ID_unit,]

library(caret)
rfeControls <- rfeControl(functions = rfFuncs, method = 'repeatedcv', repeats =50)
fs_nb <- rfe(x = data_train[,8:18], y = data_train[,"Type"], sizes = 1:11, rfeControl = rfeControls)
pdf("accuracy.pdf", height = 4, width = 6)
plot(fs_nb, type = "b", cex = 2, pt.cex = 1.2, col = "black")
dev.off()
fs_nb$optVariables
# [1] "Betweenness" "Degree"      "Closeness"   "Si" 

library(randomForest)
err <- as.numeric()
for (i in 1:4) {
  mtry_test <- randomForest(Type ~ Betweenness+Closeness+Degree+Si, data = data_train, mtry = i)
  err <- append(err,mean(mtry_test$err.rate))
}
mtry <- which.min(err)
ntree_fit <- randomForest(Type ~ Betweenness+Closeness+Degree+Si, data = data_train, mtry = mtry, ntree = 10000)
plot(ntree_fit)
rf <- randomForest(Type ~ Betweenness+Closeness+Degree+Si, data = data_train, mtry = mtry, ntree = 5000,importance = TRUE)
rf$importance
#                 CANCER  RASopathy MeanDecreaseAccuracy MeanDecreaseGini
# Betweenness 0.09944920 0.18203887           0.13313135         17.02616
# Closeness   0.06616207 0.03583261           0.05356388         13.39230
# Degree      0.20799078 0.24543959           0.22216629         45.96373
# Si          0.13619011 0.08901264           0.11662843         27.82527
varImpPlot(rf, main = "variable importance")
ipt[,5]<-rownames(ipt)
colnames(ipt)[5]<-"Features"
ipt2<-ipt[order(ipt$MeanDecreaseAccuracy),]
ipt$Features<-factor(ipt$Features, levels=ipt2$Features)
require(grid)
a<-ggplot(data = ipt,aes(x = Features, y =MeanDecreaseAccuracy))+ 
  theme_set(theme_bw())+ 
  geom_col(position="dodge",color="#000066",fill="#000066")+
  coord_flip()+
  theme(panel.grid.major=element_line(colour=NA))+
  xlab("Features")+
  ylab("Mean Decrease Accuracy")+
  theme(axis.text.x =element_text(size = 20),
        axis.text.y =element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)
  )
b<-ggplot(data = ipt,aes(x = Features, y =MeanDecreaseGini))+ 
  theme_set(theme_bw())+ 
  geom_col(position="dodge",color="#000066",fill="#000066")+
  coord_flip()+
  theme(panel.grid.major=element_line(colour=NA))+
  xlab("")+
  ylab("Mean Decrease Gini")+
  theme(axis.text.x =element_text(size = 20),
        axis.text.y =element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)
  )

pred <- predict(rf, newdata = data_test[,-6])
collect_test_data <- data.frame(prob = pred, obs = data_test$Type)
table(data_test$Type, pred, dnn = c("symbol","predict"))

library(pROC)
mroc <- roc(data_test$Type, as.numeric(pred), percent = T, levels = c("RASopathy","CANCER"))
plot(mroc, print.auc=T, auc.polygon=T, 
    grid=c(0.1,0.2), auc.polygon.col = "skyblue", main = "ROC Curve for the Final RF Model")

#建模
err <- as.numeric()
for (i in 1:11) {
  mtry_test <- randomForest(Type ~ Score+Si+MI+ΔΔG+RASA+Degree+Betweenness+Closeness+MSF+Effectiveness+Sensitivity, 
    data = data_train, mtry = i)
  err <- append(err,mean(mtry_test$err.rate))
}
mtry <- which.min(err)
ntree_fit <- randomForest(Type ~ Score+Si+MI+ΔΔG+RASA+Degree+Betweenness+Closeness+MSF+Effectiveness+Sensitivity, 
    data = data_train, mtry = mtry, ntree = 10000)
plot(ntree_fit)
ntree_fit <- randomForest(Type ~ Score+Si+MI+ΔΔG+RASA+Degree+Betweenness+Closeness+MSF+Effectiveness+Sensitivity, 
    data = data_train, mtry = mtry, ntree = 1000)
plot(ntree_fit)
rf <- randomForest(Type ~ Score+Si+MI+ΔΔG+RASA+Degree+Betweenness+Closeness+MSF+Effectiveness+Sensitivity, 
    data = data_train, mtry = mtry, ntree = 500,importance = TRUE)
plot(rf)
rf$importance
#                    CANCER     RASopathy MeanDecreaseAccuracy MeanDecreaseGini
# Score         0.009062134 -9.717459e-04          0.004935991        0.8019283
# Si            0.083448204  4.976393e-02          0.069526845       15.5781351
# MI            0.065586593  4.133321e-02          0.055580502       11.4929907
# ΔΔG           0.008407840 -2.462735e-03          0.003976341        2.8790637
# RASA          0.016693470  3.071922e-03          0.011067088        2.3772777
# Degree        0.172021859  2.117388e-01          0.187504608       39.0584225
# Betweenness   0.088319858  1.696983e-01          0.121557058       15.3753467
# Closeness     0.055620561  3.830277e-02          0.048589217       11.4630305
# MSF           0.006922767  4.537669e-03          0.005959592        1.8632387
# Effectiveness 0.005228204  3.476906e-03          0.004490475        1.6731040
# Sensitivity   0.006903774  6.213327e-05          0.004170094        1.6835547
varImpPlot(rf, main = "variable importance")
ipt<-as.data.frame(rf$importance)
ipt[,5]<-rownames(ipt)
colnames(ipt)[5]<-"Features"
ipt2<-ipt[order(ipt$MeanDecreaseAccuracy),]
ipt$Features<-factor(ipt$Features, levels=ipt2$Features)
require(grid)
a<-ggplot(data = ipt,aes(x = Features, y =MeanDecreaseAccuracy))+ 
  theme_set(theme_bw())+ 
  geom_col(position="dodge",color="#000066",fill="#000066")+
  coord_flip()+
  theme(panel.grid.major=element_line(colour=NA))+
  xlab("Features")+
  ylab("Mean Decrease Accuracy")+
  theme(axis.text.x =element_text(size = 20),
        axis.text.y =element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)
  )+
  labs(tags = "C")
b<-ggplot(data = ipt,aes(x = Features, y =MeanDecreaseGini))+ 
  theme_set(theme_bw())+ 
  geom_col(position="dodge",color="#000066",fill="#000066")+
  coord_flip()+
  theme(panel.grid.major=element_line(colour=NA))+
  xlab("Features")+
  ylab("Mean Decrease Gini")+
  theme(axis.text.x =element_text(size = 20),
        axis.text.y =element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20)
  )+
  labs(tags = "D")

pred <- predict(rf, newdata = data_test[,-6])
collect_test_data <- data.frame(prob = pred, obs = data_test$Type)
table(data_test$Type, pred, dnn = c("symbol","predict"))

library(pROC)
mroc <- roc(data_test$Type, as.numeric(pred), percent = T, levels = c("RASopathy","CANCER"))

#建模后预测的ROC
plot(mroc, print.auc=T, auc.polygon=T, 
    grid=c(0.1,0.2), auc.polygon.col = "skyblue", main = "ROC Curve for the Final RF Model")

pa <- ggroc(mroc)+
        geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100),color="darkgrey", linetype="dashed")+
        labs(title="ROC Curve for the Final RF Model",subtitle = "AUC = 93.27%",tags = "E")+
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()
            )+
        xlab("Specificity(%)")+
        ylab("Sensitivity(%)")


roc.list <- roc(Type ~ Score+Si+MI+ΔΔG+RASA+Degree+Betweenness+Closeness+MSF+Effectiveness+Sensitivity, 
    data = RC_mat,
    percent = T,
    levels = c("RASopathy","CANCER"))
#读取每条曲线的auc
auc_list <- c()
for(i in roc.list){
    auc_list <- c(auc_list,i$auc)
}
auc_list <- round(auc_list,2)
para_list <- colnames(RC_mat[,8:18])
n_auc_list <- as.numeric(auc_list)
tab <- cbind(para_list,n_auc_list)
tab <- tab[order(tab[,2],decreasing = TRUE),]
tab[,2] <- as.numeric(tab[,2])


lab_list <- c()
for(i in 1:length(auc_list)){
    lab_text <- paste(tab[,1][i],"(AUC=",tab[,2][i],")",sep="") 
    lab_list <- c(lab_list,lab_text)
}
lab_list

g.list <- ggroc(roc.list)+
    geom_segment(aes(x = 100, xend = 0, y = 0, yend = 100), color="darkgrey", linetype="dashed")+
    scale_colour_discrete(breaks = para_list,labels = lab_list)+
    xlab("Specificity(%)")+
    ylab("Sensitivity(%)")+
    theme(
        legend.title=element_blank(),
        legend.position = c(.95, .0),
        legend.justification = c("right", "bottom"),
        legend.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()
        )+
    labs(title = "ROC Curve for the Parameters",tags = "B")

g.list

#拼图
(p1|a/b)/(g.list|pa)


