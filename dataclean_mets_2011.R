library(foreign)
library(utils)
library(haven)
library(epiDisplay)
library(survival)
library(gmodels)
library(stringr)
# 1 --------------------------------------------------------------
#1.1导入数据
biom11<-read.spss("D:/master/clhls/biom/biom12.sav",to.data.frame =T)
cross11<-read.spss("D:/master/clhls/biom/lds1114.sav",to.data.frame =T)
summary(cross11$trueage)
subcross11<-cross11[cross11$trueage>=65,]
#1.2连接数据
mets<-merge.data.frame(biom11,subcross11,by.x = 'id',by.y ='id')
summ.factor(mets$a1.x)
# 2 METS诊断 ------------------------------------------------------------------
varnames<-c('g102c','g511','g512','g521','g522','hdlc12','glu12','tg12')
temp<-tapply(mets[,varnames],as.numeric())

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

# 3.MMSE --------------------------------------------------------------------
##mmse
varnames<-c('c11','c12','c13','c14','c15','c21a','c21b','c21c','c31a','c31b',
            'c31c','c31d','c31e','c32','c41a','c41b','c41c','c51a','c51b','c52',
            'c53a','c53b','c53c')
nvarnames<-c('nc11','nc12','nc13','nc14','nc15','nc21a','nc21b','nc21c','nc31a','nc31b',
             'nc31c','nc31d','nc31e','nc32','nc41a','nc41b','nc41c','nc51a','nc51b','nc52',
             'nc53a','nc53b','nc53c')

summ.factor(factor(mets$c16))
mets$c16[mets$c16=='not able to answer']<-0
mets$c16[mets$c16=='missing']<-0
mets$nnc16<-as.numeric(mets$c16)
mets$nc16[mets$nnc16==1]<-1
mets$nc16[mets$nnc16==2]<-2
mets$nc16[mets$nnc16==3]<-4
mets$nc16[mets$nnc16==4]<-4
mets$nc16[mets$nnc16==5]<-5
mets$nc16[mets$nnc16==6]<-6
mets$nc16[mets$nnc16>6]<-7
summ.factor(factor(mets$nc16))

summ.factor(mets$c11)
mets$nc11<-NA
mets$nc11[mets$c11=='correct']<-1
mets$nc11[mets$c11=='wrong']<-0
mets$nc11[mets$c11=='not able to answer']<-0
mets$nc11[mets$c11=='missing']<-NA
summ.factor(mets$nc11)


summ.factor(mets$c12)
mets$nc12<-NA
mets$nc12[mets$c12=='correct']<-1
mets$nc12[mets$c12=='wrong']<-0
mets$nc12[mets$c12=='not able to answer']<-0
mets$nc12[mets$c12=='missing']<-NA
summ.factor(mets$nc12)

summ.factor(mets$c13)
mets$nc13<-NA
mets$nc13[mets$c13=='correct']<-1
mets$nc13[mets$c13=='wrong']<-0
mets$nc13[mets$c13=='not able to answer']<-0
mets$nc13[mets$c13=='missing']<-NA
summ.factor(mets$nc13)

summ.factor(mets$c14)
mets$nc14<-NA
mets$nc14[mets$c14=='correct']<-1
mets$nc14[mets$c14=='wrong']<-0
mets$nc14[mets$c14=='not able to answer']<-0
mets$nc14[mets$c14=='missing']<-NA
summ.factor(mets$nc14)

summ.factor(mets$c15)
mets$nc15<-NA
mets$nc15[mets$c15=='correct']<-1
mets$nc15[mets$c15=='wrong']<-0
mets$nc15[mets$c15=='not able to answer']<-0
mets$nc15[mets$c15=='missing']<-NA
summ.factor(mets$nc15)

summ.factor(mets$c21a)
mets$nc21a<-NA
mets$nc21a[mets$c21a=='correct']<-1
mets$nc21a[mets$c21a=='wrong']<-0
mets$nc21a[mets$c21a=='not able to answer']<-0
mets$nc21a[mets$c21a=='missing']<-NA
summ.factor(mets$nc21a)

summ.factor(mets$c21b)
mets$nc21b<-NA
mets$nc21b[mets$c21b=='correct']<-1
mets$nc21b[mets$c21b=='wrong']<-0
mets$nc21b[mets$c21b=='not able to answer']<-0
mets$nc21b[mets$c21b=='missing']<-NA
summ.factor(mets$nc21b)

summ.factor(mets$c21c)
mets$nc21c<-NA
mets$nc21c[mets$c21c=='correct']<-1
mets$nc21c[mets$c21c=='wrong']<-0
mets$nc21c[mets$c21c=='2']<-0
mets$nc21c[mets$c21c=='not able to answer']<-0
mets$nc21c[mets$c21c=='missing']<-NA
summ.factor(mets$nc21c)

summary(factor(mets$c22))
mets$nc22<-0
mets$nc22[mets$c22=='0']<-1
mets$nc22[mets$c22=='missing']<-NA
summ.factor(mets$nc22)

summ.factor(mets$c31a)
mets$nc31a<-NA
mets$nc31a[mets$c31a=='correct']<-1
mets$nc31a[mets$c31a=='wrong']<-0
mets$nc31a[mets$c31a=='not able to answer']<-0
mets$nc31a[mets$c31a=='missing']<-NA
summ.factor(mets$nc31a)

summ.factor(mets$c31b)
mets$nc31b<-NA
mets$nc31b[mets$c31b=='correct']<-1
mets$nc31b[mets$c31b=='wrong']<-0
mets$nc31b[mets$c31b=='not able to answer']<-0
mets$nc31b[mets$c31b=='missing']<-NA
summ.factor(mets$nc31b)

summ.factor(mets$c31c)
mets$nc31c<-NA
mets$nc31c[mets$c31c=='correct']<-1
mets$nc31c[mets$c31c=='wrong']<-0
mets$nc31c[mets$c31c=='not able to answer']<-0
mets$nc31c[mets$c31c=='missing']<-NA
summ.factor(mets$nc31c)

summ.factor(mets$c31d)
mets$nc31d<-NA
mets$nc31d[mets$c31d=='correct']<-1
mets$nc31d[mets$c31d=='wrong']<-0
mets$nc31d[mets$c31d=='unable to do']<-0
mets$nc31d[mets$c31d=='missing']<-0
summ.factor(mets$nc31d)

summ.factor(mets$c31e)
mets$nc31e<-NA
mets$nc31e[mets$c31e=='correct']<-1
mets$nc31e[mets$c31e=='wrong']<-0
mets$nc31e[mets$c31e=='not able to answer']<-0
mets$nc31e[mets$c31e=='missing']<-NA
summ.factor(mets$nc31e)

summ.factor(mets$cc32)
mets$cc32<-as.numeric(mets$c32)
tail(mets$c32,n=20)
mets$nc32<-NA
mets$nc32[mets$c32=='correct']<-1
mets$nc32[mets$c32=='wrong']<-0
mets$nc32[mets$c32=="not able to answer"]<-0
mets$nc32[mets$c32=="not able to do this (disabled)"]<-0
mets$nc32[mets$c32=='missing']<-0
summ.factor(mets$nc32)

summ.factor(mets$c41a)
mets$nc41a<-NA
mets$nc41a[mets$c41a=='correct']<-1
mets$nc41a[mets$c41a=='wrong']<-0
mets$nc41a[mets$c41a=='not able to answer']<-0
mets$nc41a[mets$c41a=='missing']<-NA
summ.factor(mets$nc41a)

summ.factor(mets$c41b)
mets$nc41b<-NA
mets$nc41b[mets$c41b=='correct']<-1
mets$nc41b[mets$c41b=='wrong']<-0
mets$nc41b[mets$c41b=='not able to answer']<-0
mets$nc41b[mets$c41b=='missing']<-NA
summ.factor(mets$nc41b)

summ.factor(mets$c41c)
mets$nc41c<-NA
mets$nc41c[mets$c41c=='correct']<-1
mets$nc41c[mets$c41c=='wrong']<-0
mets$nc41c[mets$c41c=='not able to answer']<-0
mets$nc41c[mets$c41c=='missing']<-NA
summ.factor(mets$nc41c)

summ.factor(mets$c51a)
mets$nc51a<-NA
mets$nc51a[mets$c51a=='correct']<-1
mets$nc51a[mets$c51a=='wrong']<-0
mets$nc51a[mets$c51a=='not able to answer']<-0
mets$nc51a[mets$c51a=='missing']<-NA
summ.factor(mets$nc51a)

summ.factor(mets$c51b)
mets$nc51b<-NA
mets$nc51b[mets$c51b=='correct']<-1
mets$nc51b[mets$c51b=='wrong']<-0
mets$nc51b[mets$c51b=='not able to answer']<-0
mets$nc51b[mets$c51b=='missing']<-NA
summ.factor(mets$nc51b)

summ.factor(mets$c52)
mets$nc52<-NA
mets$nc52[mets$c52=='correct']<-3
mets$nc52[mets$c52=='wrong']<-0
mets$nc52[mets$c52=='not able to answer']<-0
mets$nc52[mets$c52=='missing']<-NA
summ.factor(mets$nc52)

summ.factor(mets$c53a)
mets$nc53a<-NA
mets$nc53a[mets$c53a=='correct']<-1
mets$nc53a[mets$c53a=='wrong']<-0
mets$nc53a[mets$c53a=='not able to do']<-0
mets$nc53a[mets$c53a=='missing']<-NA
summ.factor(mets$nc53a)

summ.factor(mets$c53b)
mets$nc53b<-NA
mets$nc53b[mets$c53b=='correct']<-1
mets$nc53b[mets$c53b=='wrong']<-0
mets$nc53b[mets$c53b=='not able to do']<-0
mets$nc53b[mets$c53b=='missing']<-NA
summ.factor(mets$nc53b)

summ.factor(mets$c53c)
mets$nc53c<-NA
mets$nc53c[mets$c53c=='correct']<-1
mets$nc53c[mets$c53c=='wrong']<-0
mets$nc53c[mets$c53c=='not able to do']<-0
mets$nc53c[mets$c53c=='missing']<-NA
summ.factor(mets$nc53c)

mets$mmsescore<-rowSums(mets[,nvarnames])
summ.factor(mets$mmsescore)
summary(factor(mets$mmsescore))

# 24无障碍 23~18轻度<=17中重度
mets$mmse<-NA
mets$mmse[mets$mmsescore<18]<-2
mets$mmse[mets$mmsescore<24 & mets$mmsescore>=18]<-1
mets$mmse[mets$mmsescore>=24]<-0
summ.factor(mets$mmse)

# 4.control variable --------------------------------------------------------
##ADL
summ.factor(mets$e1)
mets$adl1<-NA
mets$adl1[mets$e1=='without assistance']<-0
mets$adl1[mets$e1=='one part assistance']<-1
mets$adl1[mets$e1=='more than one part assistance']<-2
summ.factor(mets$adl1)

summ.factor(mets$e2)
mets$adl2<-NA
mets$adl2[mets$e2=='without assistance']<-0
mets$adl2[mets$e2=='need assistance for trying shoes']<-1
mets$adl2[mets$e2==' assistance in getting clothes and getting dressed']<-2
summ.factor(mets$adl2)

summ.factor(mets$e3)
mets$adl3<-NA
mets$adl3[mets$e3=='without assistance']<-0
mets$adl3[mets$e3=='assistance in cleaning or arranging clothes']<-1
mets$adl3[mets$e3=="don't use toilet"]<-2
summ.factor(mets$adl3)

summ.factor(mets$e4)
mets$adl4<-NA
mets$adl4[mets$e4=='without assistance']<-0
mets$adl4[mets$e4=='with assistance']<-1
mets$adl4[mets$e4=='bedridden']<-2
summ.factor(mets$adl4)

summ.factor(mets$e5)
mets$adl5<-NA
mets$adl5[mets$e5=='without assistance']<-0
mets$adl5[mets$e5==' occasional accidents']<-1
mets$adl5[mets$e5=='incontinent']<-2
summ.factor(mets$adl5)

summ.factor(mets$e6)
mets$adl6<-NA
mets$adl6[mets$e6=='without assistance']<-0
mets$adl6[mets$e6=='with some help']<-1
mets$adl6[mets$e6=='need feeding']<-2
summ.factor(mets$adl6)

mets$adl<-mets$adl1+mets$adl2+mets$adl3+mets$adl4+mets$adl5+mets$adl6
summ.factor(mets$adl)
mets$adl_2<-mets$adl
mets$adl_2[mets$adl_2>0]<-1
summ.factor(mets$adl_2)


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

varwant<-c("id","trueage.x","waist","gluind","sbp","dbp","bp","tgind","hdlcind",
           "index5","ms","mmsescore","mmse","adl","adl_2",
           "marriage","nationality","sex","smoke","acohol","exercise",
           "occupation", "education","urban","reside","income")
metswant11<-mets[,varwant]
metswant11$year<-2011
write.csv(metswant11, file="D:/master/clhls/biom/MetS/MetS/mets11only.csv")


