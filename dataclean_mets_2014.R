library(foreign)
library(utils)
library(haven)
library(epiDisplay)
library(survival)
library(gmodels)
# 1. --------------------------------------------------------------
#1.1导入数据
biom14<-read.spss("D:/master/clhls/biom/biom14.sav",to.data.frame =T)
cross14<-read.spss("D:/master/clhls/biom/cross_sectional_14.sav",to.data.frame =T)
summary(cross14$trueage)
subcross14<-cross14[cross14$trueage>=65,]
#1.2连接数据
mets<-merge.data.frame(biom14,subcross14,by.x = 'id',by.y ='id')
summ.factor(mets$a1.x)

# 2 METS诊断 ------------------------------------------------------------------
varnames<-c('g102c','g511','g512','g521','g522','hdlc12','glu12','tg12')
for (i in varnames) {
  mets[,i]<-as.numeric(mets[,i])
}
mets_f<-mets[which(mets$a1.y=='female'),]
mets_m<-mets[which(mets$a1.y=='male'),]


# 2.1腹型肥胖(即中心型肥胖):腰围男性≥90cm, 女性≥85cm -----------------------------------------
summary(mets$g102c)
mets_f$waist<-NA
mets_f$waist[mets_f$g102c>=85]<-1
mets_f$waist[mets_f$g102c<85]<-0
summ.factor(mets_f$waist)

mets_m$waist<-NA
mets_m$waist[mets_m$g102c>=90]<-1
mets_m$waist[mets_m$g102c<90]<-0
summ.factor(mets_m$waist)

mets<-rbind(mets_f,mets_m)

# 2.2空腹血糖≥6.1mmol/L 或糖负荷后2 小时血糖≥7.8mmol/L 和(或)已确诊为糖尿病并治疗者。 -----------------
summ.factor(mets$g15b2)
mets$gluind<-0
mets$gluind[mets$glu>=6.1 | mets$g15b2=='yes']<-1
mets$gluind[is.na(mets$glu) & is.na(mets$g15b2)]<-NA
summ.factor(mets$gluind)

# 2.3高血压:血压≥130/85mmHg 及(或)已确认为高血压并治疗者。  ------------------------------------
summ.factor(mets$g15a2)
mets$sbp<-(mets$g511+mets$g521)/2
mets$sbp[which(is.na(mets$sbp))]<-mets$g511[which(is.na(mets$sbp))]
mets$sbp[which(is.na(mets$sbp))]<-mets$g521[which(is.na(mets$sbp))]
mets$dbp<-(mets$g512+mets$g522)/2
mets$dbp[which(is.na(mets$sbp))]<-mets$g512[which(is.na(mets$sbp))]
mets$dbp[which(is.na(mets$sbp))]<-mets$g522[which(is.na(mets$sbp))]
mets$bp<-0
mets$bp[mets$sbp>=130 | mets$dbp>=85 | mets$g15a2=='yes']<-1
mets$bp[is.na(mets$sbp) & is.na(mets$dbp) & is.na(mets$g15a2)]<-NA
summ.factor(mets$bp)

# 2.4空腹TG≥1.70mmol/L。  ----------------------------------------------------
mets$tgind<-NA
mets$tgind[mets$tg>=1.7]<-1
mets$tgind[mets$tg<1.7]<-0
summ.factor(mets$tgind)



# 2.5空腹HDL-C<1.04mmol/L。  -------------------------------------------------
mets$hdlcind<-NA
mets$hdlcind[mets$hdlc<1.04]<-1
mets$hdlcind[mets$hdlc>=1.04]<-0
summ.factor(mets$hdlcind)


# 2.6mets -----------------------------------------------------------------
mets$index5<-mets$waist+mets$gluind+mets$tgind+mets$bp+mets$hdlcind
summ.factor(mets$index5)
mets$ms<-NA
mets$ms[mets$index5<3]<-0
mets$ms[mets$index5>=3]<-1
summ.factor(mets$ms)


# 3.control variable --------------------------------------------------------

summ.factor(mets$f41)
mets$f41<-str_trim(mets$f41,side ="both")
mets$marriage<-NA
mets$marriage[mets$f41=='currently married and living with spouse']<-1
mets$marriage[mets$f41=='married but not living with spouse']<-1
mets$marriage[mets$f41=='divorced']<-0
mets$marriage[mets$f41=='never married']<-0
mets$marriage[mets$f41=='widowed']<-0
summ.factor(mets$marriage)

summ.factor(mets$a2)
mets$a2<-str_trim(mets$a2,side ="both")
mets$nationality<-0
mets$nationality[mets$a2=='han']<-1
mets$nationality[is.na(mets$a2)]<-NA
mets$nationality[mets$a2=='missing']<-NA
summ.factor(mets$nationality)

summ.factor(mets$a1.y)
mets$a1.y<-str_trim(mets$a1.y,side ="both")
mets$sex<-NA
mets$sex[mets$a1.y=='male']<-1
mets$sex[mets$a1.y=='female']<-0
summ.factor(mets$sex)

summ.factor(mets$d71)
mets$d71<-str_trim(mets$d71,side ="both")
mets$smoke<-NA
mets$smoke[mets$d71=='yes']<-1
mets$smoke[mets$d71=='no']<-0
summ.factor(mets$smoke)

summ.factor(mets$d81)
mets$d81<-str_trim(mets$d81,side ="both")
mets$acohol<-NA
mets$acohol[mets$d81=='yes']<-1
mets$acohol[mets$d81=='no']<-0
summ.factor(mets$acohol)

summ.factor(mets$d91)
mets$d71<-str_trim(mets$d71)
mets$exercise<-NA
mets$exercise[mets$d91=='yes']<-1
mets$exercise[mets$d91=='no']<-0
summ.factor(mets$exercise)

?str_trim()
summary(factor(mets$f2))
mets$f2<-str_trim(mets$f2,side ="both")
mets$occupation<-0
mets$occupation[mets$f2=='agriculture, forestry, animal husbandry or fishery worker']<-1
mets$occupation[is.na(mets$f2)]<-NA
summ.factor(mets$occupation)

summary(factor(mets$f1))
mets$f2<-str_trim(mets$f2,side ="both")
mets$education<-1
mets$education[mets$f1=='0']<-0
mets$education[mets$f1=="don't know"]<-NA
mets$education[is.na(mets$f1)]<-NA
summ.factor(mets$education)

summ.factor(mets$a43)
mets$a43<-str_trim(mets$a43,side ="both")
mets$urban<-NA
mets$urban[mets$a43=='urban']<-1
mets$urban[mets$a43=='rural']<-0
summ.factor(mets$urban)

summ.factor(mets$a51)
mets$a51<-str_trim(mets$a51,side ="both")
mets$reside<-NA
mets$reside[mets$a51=='with household member(s)']<-1
mets$reside[mets$a51=='alone']<-0
mets$reside[mets$a51=='in an institution']<-0
summ.factor(mets$reside)

summ.factor(mets$f31)
mets$f31<-str_trim(mets$f31,side ="both")
mets$income<-0
mets$income[mets$f31=='retirement wages']<-1
mets$income[mets$f31=='missing']<-NA
mets$income[is.na(mets$f31=='missing')]<-NA
summ.factor(mets$income)

# 最终使用变量 ------------------------------------------------------------------

varwant<-c("id","waist","gluind","sbp","dbp","bp","tgind","hdlcind","index5","ms", 
           "marriage","nationality","sex","smoke","acohol","exercise")
metswant14<-mets[,varwant]
metswant14$year<-2014
write.csv(metswant14, file="D:/master/clhls/biom/MetS/MetS/mets14only.csv")

# 合并11、14年数据库 -------------------------------------------------------------
mets_11_14<-rbind(metswant11,metswant14)
write.csv(metstemp, file="D:/master/clhls/biom/MetS/MetS/mets_11_14.csv")