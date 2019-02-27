##################################################################
#          Predicting Long-term Outcome of Metabolic Syndrome on 
#          Cognitive Impairment using Machine Learning
##################################################################

# 1.load packages and combine data ------------------------------------------------
library(dplyr)
library(foreign)
library(stringr)
library(epiDisplay)
library(VIM)
library(mice)
library(naniar)
library(caret)
biom11<-read.spss("biom12.sav",to.data.frame =T)
cross11<-read.spss("lds1114.sav",to.data.frame =T)
biom14<-read.spss("biom14.sav",to.data.frame =T)
cross14<-read.spss("cross_sectional_14.sav",to.data.frame =T)
# combine dataset
met<-cross11%>%filter(trueage>=65)%>%
  right_join(biom11,by="id")%>%
  mutate(g102c=as.numeric(g102c),g511=as.numeric(g511),g512=as.numeric(g512),
          g521=as.numeric(g521),g522=as.numeric(g522),hdlc=as.numeric(hdlc),
          glu=as.numeric(glu),tg=as.numeric(tg))
met14<-biom14%>%left_join(cross14,by="id")

# 2. data manipulation ----------------------------------------------------
#2.1 MMSE 
varnames<-c('c11','c12','c13','c14','c15','c21a','c21b','c21c','c31a','c31b',
            'c31c','c31d','c31e','c32','c41a','c41b','c41c','c51a','c51b','c52',
            'c53a','c53b','c53c')
nvarnames<-c('nc11','nc12','nc13','nc14','nc15','nc21a','nc21b','nc21c','nc31a','nc31b',
             'nc31c','nc31d','nc31e','nc32','nc41a','nc41b','nc41c','nc51a','nc51b','nc52',
             'nc53a','nc53b','nc53c')

met14$c16[met14$c16=='not able to answer']<-0
met14$c16[met14$c16=='missing']<-0
met14$nnc16<-as.numeric(met14$c16)
met14$nc16[met14$nnc16==1]<-1
met14$nc16[met14$nnc16==2]<-2
met14$nc16[met14$nnc16==3]<-4
met14$nc16[met14$nnc16==4]<-4
met14$nc16[met14$nnc16==5]<-5
met14$nc16[met14$nnc16==6]<-6
met14$nc16[met14$nnc16>6]<-7

met14$nc11<-NA
met14$nc11[met14$c11=='correct']<-1
met14$nc11[met14$c11=='wrong']<-0
met14$nc11[met14$c11=='not able to answer']<-0
met14$nc11[met14$c11=='missing']<-NA

met14$nc12<-NA
met14$nc12[met14$c12=='correct']<-1
met14$nc12[met14$c12=='wrong']<-0
met14$nc12[met14$c12=='not able to answer']<-0
met14$nc12[met14$c12=='missing']<-NA

met14$nc13<-NA
met14$nc13[met14$c13=='correct']<-1
met14$nc13[met14$c13=='wrong']<-0
met14$nc13[met14$c13=='not able to answer']<-0
met14$nc13[met14$c13=='missing']<-NA

met14$nc14<-NA
met14$nc14[met14$c14=='correct']<-1
met14$nc14[met14$c14=='wrong']<-0
met14$nc14[met14$c14=='not able to answer']<-0
met14$nc14[met14$c14=='missing']<-NA

met14$nc15<-NA
met14$nc15[met14$c15=='correct']<-1
met14$nc15[met14$c15=='wrong']<-0
met14$nc15[met14$c15=='not able to answer']<-0
met14$nc15[met14$c15=='missing']<-NA

met14$nc21a<-NA
met14$nc21a[met14$c21a=='correct']<-1
met14$nc21a[met14$c21a=='wrong']<-0
met14$nc21a[met14$c21a=='not able to answer']<-0
met14$nc21a[met14$c21a=='missing']<-NA

met14$nc21b<-NA
met14$nc21b[met14$c21b=='correct']<-1
met14$nc21b[met14$c21b=='wrong']<-0
met14$nc21b[met14$c21b=='not able to answer']<-0
met14$nc21b[met14$c21b=='missing']<-NA

met14$nc21c<-NA
met14$nc21c[met14$c21c=='correct']<-1
met14$nc21c[met14$c21c=='wrong']<-0
met14$nc21c[met14$c21c=='not able to answer']<-0
met14$nc21c[met14$c21c=='missing']<-NA

met14$nc22<-0
met14$nc22[met14$c22=='0']<-1
met14$nc22[met14$c22=='missing']<-NA

met14$nc31a<-NA
met14$nc31a[met14$c31a=='correct']<-1
met14$nc31a[met14$c31a=='wrong']<-0
met14$nc31a[met14$c31a=='not able to answer']<-0
met14$nc31a[met14$c31a=='missing']<-NA

met14$nc31b<-NA
met14$nc31b[met14$c31b=='correct']<-1
met14$nc31b[met14$c31b=='wrong']<-0
met14$nc31b[met14$c31b=='not able to answer']<-0
met14$nc31b[met14$c31b=='missing']<-NA

met14$nc31c<-NA
met14$nc31c[met14$c31c=='correct']<-1
met14$nc31c[met14$c31c=='wrong']<-0
met14$nc31c[met14$c31c=='not able to answer']<-0
met14$nc31c[met14$c31c=='missing']<-NA

met14$nc31d<-NA
met14$nc31d[met14$c31d=='correct']<-1
met14$nc31d[met14$c31d=='wrong']<-0
met14$nc31d[met14$c31d=='unable to do']<-0
met14$nc31d[met14$c31d=='missing']<-0

met14$nc31e<-NA
met14$nc31e[met14$c31e=='correct']<-1
met14$nc31e[met14$c31e=='wrong']<-0
met14$nc31e[met14$c31e=='not able to answer']<-0
met14$nc31e[met14$c31e=='missing']<-NA

met14$nc32<-NA
met14$nc32[met14$c32=='correct']<-1
met14$nc32[met14$c32=='wrong']<-0
met14$nc32[met14$c32=="can't use pen to draw the figure"]<-0
met14$nc32[met14$c32=="not able to do this (disabled)"]<-0
met14$nc32[met14$c32=='missing']<-NA

met14$nc41a<-NA
met14$nc41a[met14$c41a=='correct']<-1
met14$nc41a[met14$c41a=='wrong']<-0
met14$nc41a[met14$c41a=='not able to answer']<-0
met14$nc41a[met14$c41a=='missing']<-NA

met14$nc41b<-NA
met14$nc41b[met14$c41b=='correct']<-1
met14$nc41b[met14$c41b=='wrong']<-0
met14$nc41b[met14$c41b=='not able to answer']<-0
met14$nc41b[met14$c41b=='missing']<-NA

met14$nc41c<-NA
met14$nc41c[met14$c41c=='correct']<-1
met14$nc41c[met14$c41c=='wrong']<-0
met14$nc41c[met14$c41c=='not able to answer']<-0
met14$nc41c[met14$c41c=='missing']<-NA

met14$nc51a<-NA
met14$nc51a[met14$c51a=='correct']<-1
met14$nc51a[met14$c51a=='wrong']<-0
met14$nc51a[met14$c51a=='not able to answer']<-0
met14$nc51a[met14$c51a=='missing']<-NA

met14$nc51b<-NA
met14$nc51b[met14$c51b=='correct']<-1
met14$nc51b[met14$c51b=='wrong']<-0
met14$nc51b[met14$c51b=='not able to answer']<-0
met14$nc51b[met14$c51b=='missing']<-NA

met14$nc52<-NA
met14$nc52[met14$c52=='correct']<-3
met14$nc52[met14$c52=='wrong']<-0
met14$nc52[met14$c52=='not able to answer']<-0
met14$nc52[met14$c52=='missing']<-NA

met14$nc53a<-NA
met14$nc53a[met14$c53a=='correct']<-1
met14$nc53a[met14$c53a=='wrong']<-0
met14$nc53a[met14$c53a=='not able to do']<-0
met14$nc53a[met14$c53a=='missing']<-NA

met14$nc53b<-NA
met14$nc53b[met14$c53b=='correct']<-1
met14$nc53b[met14$c53b=='wrong']<-0
met14$nc53b[met14$c53b=='not able to do']<-0
met14$nc53b[met14$c53b=='missing']<-NA

met14$nc53c<-NA
met14$nc53c[met14$c53c=='correct']<-1
met14$nc53c[met14$c53c=='wrong']<-0
met14$nc53c[met14$c53c=='not able to do']<-0
met14$nc53c[met14$c53c=='missing']<-NA

met14$mmsescore<-rowSums(met14[,nvarnames])
summary(factor(met14$mmsescore))

met14$mmse<-NA
met14$mmse[met14$mmsescore<18]<-2
met14$mmse[met14$mmsescore<24 & met14$mmsescore>=18]<-1
met14$mmse[met14$mmsescore>=24]<-0

metswant14<-met14[,c("id","mmsescore","mmse")]
#2.2 met
#waist-obesity
met<-met%>%mutate(waist = ifelse(a1.y == 1 & g102c>=90, 1, ifelse(a1.y == 2 & g102c>=85,1, 0)))
# 2.2空腹血糖≥6.1mmol/L 或糖负荷后2 小时血糖≥7.8mmol/L 和(或)已确诊为糖尿病并治疗者。 -----------------

met$gluind<-0
met$gluind[met$glu>=6.1 | met$g15b2=='yes']<-1
met$gluind[is.na(met$glu) & is.na(met$g15b2)]<-NA

# 2.3高血压:血压≥130/85mmHg 及(或)已确认为高血压并治疗者。  ------------------------------------

met$sbp<-(met$g511+met$g521)/2
met$sbp[which(is.na(met$sbp))]<-met$g511[which(is.na(met$sbp))]
met$sbp[which(is.na(met$sbp))]<-met$g521[which(is.na(met$sbp))]
met$dbp<-(met$g512+met$g522)/2
met$dbp[which(is.na(met$sbp))]<-met$g512[which(is.na(met$sbp))]
met$dbp[which(is.na(met$sbp))]<-met$g522[which(is.na(met$sbp))]
met$bp<-0
met$bp[met$sbp>=130 | met$dbp>=85 | met$g15a2=='yes']<-1
met$bp[is.na(met$sbp) & is.na(met$dbp) & is.na(met$g15a2)]<-NA

# 2.4空腹TG≥1.70mmol/L。  ----------------------------------------------------
met$tgind<-NA
met$tgind[met$tg>=1.7]<-1
met$tgind[met$tg<1.7]<-0

# 2.5空腹HDL-C<1.04mmol/L。  -------------------------------------------------
met$hdlcind<-NA
met$hdlcind[met$hdlc<1.04]<-1
met$hdlcind[met$hdlc>=1.04]<-0

# 2.6met -----------------------------------------------------------------
met$index5<-met$waist+met$gluind+met$tgind+met$bp+met$hdlcind
met$ms<-NA
met$ms[met$index5<3]<-0
met$ms[met$index5>=3]<-1
# 3.control variable --------------------------------------------------------
##ADL
met$adl1<-NA
met$adl1[met$e1=='without assistance']<-0
met$adl1[met$e1=='one part assistance']<-1
met$adl1[met$e1=='more than one part assistance']<-2

met$adl2<-NA
met$adl2[met$e2=='without assistance']<-0
met$adl2[met$e2=='need assistance for trying shoes']<-1
met$adl2[met$e2==' assistance in getting clothes and getting dressed']<-2

met$adl3<-NA
met$adl3[met$e3=='without assistance']<-0
met$adl3[met$e3=='assistance in cleaning or arranging clothes']<-1
met$adl3[met$e3=="don't use toilet"]<-2

met$adl4<-NA
met$adl4[met$e4=='without assistance']<-0
met$adl4[met$e4=='with assistance']<-1
met$adl4[met$e4=='bedridden']<-2

met$adl5<-NA
met$adl5[met$e5=='without assistance']<-0
met$adl5[met$e5==' occasional accidents']<-1
met$adl5[met$e5=='incontinent']<-2

met$adl6<-NA
met$adl6[met$e6=='without assistance']<-0
met$adl6[met$e6=='with some help']<-1
met$adl6[met$e6=='need feeding']<-2

met$adl<-met$adl1+met$adl2+met$adl3+met$adl4+met$adl5+met$adl6
met$adl_2<-met$adl
met$adl_2[met$adl_2>0]<-1

# marriage status
met$f41<-str_trim(met$f41,side ="both")
met$marriage<-NA
met$marriage[met$f41=='currently married and living with spouse']<-1
met$marriage[met$f41=='married but not living with spouse']<-1
met$marriage[met$f41=='divorced']<-0
met$marriage[met$f41=='never married']<-0
met$marriage[met$f41=='widowed']<-0

#nationality
met$a2<-str_trim(met$a2,side ="both")
met$nationality<-0
met$nationality[met$a2=='han']<-1
met$nationality[is.na(met$a2)]<-NA
met$nationality[met$a2=='missing']<-NA

#sex
met$a1.y<-str_trim(met$a1.y,side ="both")
met$sex<-NA
met$sex[met$a1.y==1]<-1
met$sex[met$a1.y==2]<-0

#smoke
met$d71<-str_trim(met$d71,side ="both")
met$smoke<-NA
met$smoke[met$d71=='yes']<-1
met$smoke[met$d71=='no']<-0
#drinking
met$d81<-str_trim(met$d81,side ="both")
met$acohol<-NA
met$acohol[met$d81=='yes']<-1
met$acohol[met$d81=='no']<-0
#exercise
met$d71<-str_trim(met$d71)
met$exercise<-NA
met$exercise[met$d91=='yes']<-1
met$exercise[met$d91=='no']<-0

#education
met$f2<-str_trim(met$f2,side ="both")
met$education<-1
met$education[met$f1=='0']<-0
met$education[met$f1=="don't know"]<-NA
met$education[is.na(met$f1)]<-NA

#resident
met$a51<-str_trim(met$a51,side ="both")
met$reside<-NA
met$reside[met$a51=='with household member(s)']<-1
met$reside[met$a51=='alone']<-0
met$reside[met$a51=='in an institution']<-0
#income
met$f31<-str_trim(met$f31,side ="both")
met$income<-0
met$income[met$f31=='retirement wages']<-1
met$income[met$f31=='missing']<-NA
met$income[is.na(met$f31=='missing')]<-NA


# 4.final dataset -------------------------------------------------------
varwant<-c("id","trueage.x","waist","gluind","glu","sbp","dbp","bp","tgind","tg","hdlcind",
           "hdlc","adl","adl_2","marriage","sex","smoke","acohol","exercise","education",
           "reside","income")
metswant11<-met[,varwant]
data<-metswant11%>%inner_join(metswant14,by="id")

# 5.exploring analysis --------------------------------------------------
data$outcome<-NA
data$outcome[data$mmsescore<24 ]<-1
data$outcome[data$mmsescore>=24]<-0
#计算各样本缺失状况
miss_var_summary(data)
vis_miss(data)
vis_miss(data,cluster = T)
gg_miss_var(data)
gg_miss_case(data)
gg_miss_upset(data)
aggr(data, prop=TRUE, numbers=TRUE)
x <- as.data.frame(abs(is.na(data)))
y <- x[which(apply(x,2,sum)>0)]
cor(y)
cor(data, y, use="pairwise.complete.obs")
#单变量描述
sapply(data, function(x) {summ.factor(x)})

# 6.multiple imputation ------------------------------------------------------------
data1<-as.data.frame(sapply(data, function(x) {as.factor(x)}))
data1$trueage.x<-as.numeric(data1$trueage.x)
data1$sbp<-as.numeric(data1$sbp)
data1$dbp<-as.numeric(data1$dbp)
data1$adl<-as.numeric(data1$adl)
data1$mmsescore<-as.numeric(data1$mmsescore)
data1$tg<-as.numeric(data1$tg)
data1$hdlc<-as.numeric(data1$hdlc)
data1$glu<-as.numeric(data1$glu)
summary(data1)

# Perform mice imputation, excluding certain less-than-useful variables:
set.seed(129)

mice_mod <- mice(data1[, !names(data1) %in%c('id')], method='rf') 

summary(mice_mod)
mice_output <- complete(mice_mod)
summary(mice_output)


# 7.Prediction ------------------------------------------------------------

#random forest

#Split the data back into a train set and a test set
set.seed(300)

trainIndex <- createDataPartition(mice_output$outcome,p = .8, 
                                  list = FALSE, 
                                  times = 1)
Train <- mice_output[ trainIndex,]
Test  <- mice_output[-trainIndex,]
#Building the model
ctrl <- trainControl(method = "repeatedcv",number = 10, repeats = 10)
grid_rf <- expand.grid(.mtry = c(1, 2, 4,6))
set.seed(300)

m_rf <- train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                acohol+exercise+education+reside+income,
              data=Train,method = "rf",metric = "Kappa", 
              trControl = ctrl,tuneGrid = grid_rf)
m_rf
rf_imp <- varImp(m_rf, scale = FALSE)
rf_imp <- rf_imp$importance
rf_gini <- data.frame(Variables = row.names(rf_imp), MeanDecreaseGini = rf_imp$Overall)
library(ggplot2)
ggplot(rf_gini, aes(x=reorder(Variables, MeanDecreaseGini), y=MeanDecreaseGini, fill=MeanDecreaseGini)) +
  geom_bar(stat='identity') + coord_flip() + theme(legend.position="none") + labs(x="") +
  ggtitle('Variable Importance of Random Forest') + theme(plot.title = element_text(hjust = 0.5))
plot(m_rf)
#decision tree
set.seed(300)
m_rpart <- train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                acohol+exercise+education+reside+income, 
                 data=Train, method = "rpart",metric = "Kappa", 
                 tuneLength = 30,trControl = ctrl)
m_rpart
plot(m_rpart)
#KNN
m_KNN <- train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                 acohol+exercise+education+reside+income, 
               data=Train, method = "knn",metric = "Kappa", 
               trControl=trainControl(method = "repeatedcv", number = 10, repeats = 10),
               tuneLength = 10)
m_KNN
plot(m_KNN)
#SVM
set.seed(300)
m_SVM <- train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                 acohol+exercise+education+reside+income, 
               data=Train, method = "svmLinear",metric = "Kappa", 
               trControl=trainControl(method = "repeatedcv", number = 10, repeats = 10))
m_SVM
plot(m_SVM)
#Naïve Bayes
set.seed(300)
m_nb <- train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                acohol+exercise+education+reside+income, 
              data=Train, method = "nb",metric = "Kappa", 
              trControl = ctrl)
m_nb
plot(m_nb)
#nnet
set.seed(300)
m_nnet<-train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
                acohol+exercise+education+reside+income, 
              data=Train, method = "nnet", trace = FALSE,metric = "Kappa", 
              trControl = ctrl)
m_nnet
plot(m_nnet)
#glm
set.seed(300)
m_glm<-train(outcome~trueage.x+waist+glu+sbp+dbp+tg+hdlc+adl_2+marriage+sex+smoke+
               acohol+exercise+education+reside+income, 
              data=Train, method = "glm", trace = FALSE,metric = "Kappa", 
              trControl = ctrl)
m_glm
plot(m_glm)
###Accuracy
#Lets take a look at how well each prediction does with itself, and the test set we partitioned. We partitioned the data so we could see how much over fitting was resulting from our model choice. We can now estimate what the accuracy would be on a real data set where we didn't know the outcome. 

rf.mat <- confusionMatrix(Train$outcome, predict(m_rf, Train))
rpart.mat <- confusionMatrix(Train$outcome, predict(m_rpart, Train))
knn.mat <- confusionMatrix(Train$outcome, predict(m_KNN, Train))
nnet.mat <- confusionMatrix(Train$outcome, predict(m_nnet, Train))
svm.mat <- confusionMatrix(Train$outcome, predict(m_SVM, Train))
glm.mat <- confusionMatrix(Train$outcome, predict(m_glm, Train))
train.accuracies <- c(rf.mat$overall[1], rpart.mat$overall[1], knn.mat$overall[1],
                       nnet.mat$overall[1], svm.mat$overall[1], 
                      glm.mat$overall[1])

rf.mat <- confusionMatrix(Test$outcome, predict(m_rf, Test))
rpart.mat <- confusionMatrix(Test$outcome, predict(m_rpart, Test))
knn.mat <- confusionMatrix(Test$outcome, predict(m_KNN, Test))
nnet.mat <- confusionMatrix(Test$outcome, predict(m_nnet, Test))
svm.mat <- confusionMatrix(Test$outcome, predict(m_SVM, Test))
glm.mat <- confusionMatrix(Test$outcome, predict(m_glm, Test))
test.accuracies <- c(rf.mat$overall[1], rpart.mat$overall[1], knn.mat$overall[1],
                       nnet.mat$overall[1], svm.mat$overall[1], 
                      glm.mat$overall[1])
model<-c("rf","rpart","knn","nnet","svm","glm")
model.metas <- data.frame(model,train.accuracies, test.accuracies)
model.metas








