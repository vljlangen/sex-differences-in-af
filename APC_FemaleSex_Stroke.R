library(dplyr)
library(gridExtra)
library(Epi)
library(haven)
library(ggplot2)
library(RColorBrewer)
library(showtext)
library(splines)

APC.analysis <- data.frame(read_sav("APC_analysis.sav"))
head(APC.analysis)

dput(names(APC.analysis))

APC.analysis$SexBin <- 1 - APC.analysis$SexBin

tmp.dt<-Lexis(entry=list(age=cal.yr(CohortEntryDate,format = "%Y-%m-%d")-BirthYear,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorEndOrDeath,format = "%Y-%m-%d")-cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = 0,
              exit.status = ifelse(ISaftercohortall==1,2,0),
              id=SID,
              data=APC.analysis[,c("SID", "deathDateSPSSdate", "purchasedate", "HypertensionBOAC", "CongestiveHeartFailureBOAC", 
                                   "HyperlipidemiaBOAC", "DiabetesBOAC", "IschemicStrokeOrTIABOAC", 
                                   "AnyVascularDiseaseBOAC", "CohortEntryYear", "SexBin", 
                                   "Age", "ISaftercohortall", "DateISorEndOrDeath", "LastAKdateplus120days",
                                   "Incometertiles", "CohortEntryDate", "BirthYear")])

aux.cut<-subset(tmp.dt,!is.na(purchasedate))[,c("SID","purchasedate")]
aux.cut$lex.id<-aux.cut$SID
aux.cut$cut<-cal.yr(aux.cut$purchasedate,format = "%Y-%m-%d")
aux.cut$new.state<-1
tmp.dt$age <- tmp.dt$age + tmp.dt$lex.dur / 2

tmp.dt1<-cutLexis(data = tmp.dt,
                  cut=aux.cut[,c("lex.id","cut","new.state")],
                  new.scale="ak",
                  timescale = "per")
summary(tmp.dt1)

tmp.dt2<-splitLexis(lex = tmp.dt1,
                    breaks=c(2006:2018),
                    time.scale = "per")

APC.analysis.Lx<-tmp.dt2
save(APC.analysis.Lx,file="APCanalysisLx.RData")
table(APC.analysis.Lx$lex.Cst)


range(APC.analysis.Lx$BirthYear)

tmp.m0<-glm(cbind(lex.Xst==2,lex.dur) ~
              DiabetesBOAC+HypertensionBOAC+HyperlipidemiaBOAC+CongestiveHeartFailureBOAC+
              IschemicStrokeOrTIABOAC+
              AnyVascularDiseaseBOAC+Incometertiles+
              SexBin+
              ns(BirthYear,df=10)+ 
              ns(per,df=5)+
              ns(age,df=10),
            data=APC.analysis.Lx,family=poisreg())

tmp.m1<-update(tmp.m0,~.+SexBin:ns(BirthYear,df=10))
anova(tmp.m0,tmp.m1,test="Chisq")
# ci.exp(tmp.m1,subset="Sex")

range(APC.analysis.Lx$age)
tmp.df<-with(APC.analysis.Lx,
             expand.grid( lex.Cst=sort(unique(lex.Cst))[1],
                          DiabetesBOAC=sort(unique(DiabetesBOAC))[1],
                          HypertensionBOAC=sort(unique(HypertensionBOAC))[1],
                          HyperlipidemiaBOAC=sort(unique(HyperlipidemiaBOAC))[1],
                          CongestiveHeartFailureBOAC=sort(unique(CongestiveHeartFailureBOAC))[1],
                          IschemicStrokeOrTIABOAC=sort(unique(IschemicStrokeOrTIABOAC))[1],
                          AnyVascularDiseaseBOAC=sort(unique(AnyVascularDiseaseBOAC))[1],
                          Incometertiles=sort(unique(Incometertiles))[1],
                          SexBin=sort(unique(SexBin)),
                          BirthYear=1902:1998,
                          age=seq(from=20,to=110,by=10),
                          per=2006:2019,
                          lex.Xst=0,
                          lex.dur=1
             ))

tmp.mtr<-model.matrix(tmp.m1$formula,data=tmp.df)

aux.all<-as.data.frame(ci.exp(tmp.m1,
                              ctr.mat=(tmp.mtr[tmp.df$SexBin==1,]-
                                         tmp.mtr[tmp.df$SexBin!=1,])))

aux.in<-with(tmp.df[tmp.df$SexBin==1,],(per==2014)&(age==70) )
aux.all.1<-aux.all[aux.in,]

aux.all.2<-subset(tmp.df,(per==2014)&(age==70)&(SexBin==1))
aux.all.2$IRR<-aux.all.1[,1]
aux.all.2$IRR.lo<-aux.all.1[,2]
aux.all.2$IRR.hi<-aux.all.1[,3]

head(aux.all.2)

library(ggplot2)

ggplot(aux.all.2,aes(x=BirthYear,y=IRR))+geom_line()+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(aes(ymin=IRR.lo,ymax=IRR.hi),alpha=0.2,data=aux.all.2)+
  coord_cartesian(ylim=c(0,2),xlim=c(1910,1980))+
  labs(x="Birth year",y="IRR, females vs. males")

font_add_google("Rosario", "rosario")
showtext_auto()

pastel_colors <- brewer.pal(n = 3, name = "Set1")

adjusted_fig<- ggplot(aux.all.2, aes(x = BirthYear, y = IRR)) +
  geom_ribbon(aes(ymin = IRR.lo, ymax = IRR.hi),
              fill = pastel_colors[2], alpha = 0.3) +
  geom_line(color = pastel_colors[1], size = 1.2) +
  scale_x_continuous(breaks = seq(1910, 1980, by = 10)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", size = 0.6) +
  annotate("text", x = 1910, y = 2, label = "B", hjust = 0, vjust = 1.2, size = 6, fontface = "bold") + 
  coord_cartesian(ylim = c(0, 2), xlim = c(1910, 1980)) +
  labs(
    x = "Birth year",
    y = "Incidence rate ratio (women vs. men)",
  ) +
  theme_minimal(base_family = "rosario") +
  theme(
    text = element_text(color = "gray10"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.caption = element_text(size = 10, color = "gray40", hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

adjusted_fig

######unadjsted analysis###

tmp.dt<-Lexis(entry=list(age=cal.yr(CohortEntryDate,format = "%Y-%m-%d")-BirthYear,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorEndOrDeath,format = "%Y-%m-%d")-cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = 0,
              exit.status = ifelse(ISaftercohortall==1,2,0),
              id=SID,
              data=APC.analysis[,c("SID", "deathDateSPSSdate", "purchasedate", "HypertensionBOAC", "CongestiveHeartFailureBOAC", 
                                   "HyperlipidemiaBOAC", "DiabetesBOAC", "IschemicStrokeOrTIABOAC", 
                                   "AnyVascularDiseaseBOAC", "CohortEntryYear", "SexBin", 
                                   "Age", "ISaftercohortall", "DateISorEndOrDeath", "LastAKdateplus120days",
                                   "Incometertiles", "CohortEntryDate", "BirthYear")])

tmp.dt$age <- tmp.dt$age + tmp.dt$lex.dur / 2

tmp.dt2<-splitLexis(lex = tmp.dt,
                    breaks=c(2006:2018),
                    time.scale = "per")

APC.analysis.Lx<-tmp.dt2
table(APC.analysis.Lx$lex.Cst)

library(splines)
range(APC.analysis.Lx$BirthYear)

tmp.m0<-glm(cbind(lex.Xst==2,lex.dur) ~
              SexBin+
              ns(BirthYear,df=10)+ 
              ns(per,df=5)+
              ns(age,df=10),
            data=APC.analysis.Lx,family=poisreg())

tmp.m1<-update(tmp.m0,~.+SexBin:ns(BirthYear,df=10))
anova(tmp.m0,tmp.m1,test="Chisq")

range(APC.analysis.Lx$age)
tmp.df<-with(APC.analysis.Lx,
             expand.grid( lex.Cst=sort(unique(lex.Cst))[1],
                          SexBin=sort(unique(SexBin)),
                          BirthYear=1902:1998,
                          age=seq(from=20,to=110,by=10),
                          per=2006:2019,
                          lex.Xst=0,
                          lex.dur=1
             ))

tmp.mtr<-model.matrix(tmp.m1$formula,data=tmp.df)

aux.all<-as.data.frame(ci.exp(tmp.m1,
                              ctr.mat=(tmp.mtr[tmp.df$SexBin==1,]-
                                         tmp.mtr[tmp.df$SexBin!=1,])))

aux.in<-with(tmp.df[tmp.df$SexBin==1,],(per==2014)&(age==70) )
aux.all.1<-aux.all[aux.in,]

aux.all.2<-subset(tmp.df,(per==2014)&(age==70)&(SexBin==1))
aux.all.2$IRR<-aux.all.1[,1]
aux.all.2$IRR.lo<-aux.all.1[,2]
aux.all.2$IRR.hi<-aux.all.1[,3]

head(aux.all.2)

library(ggplot2)

ggplot(aux.all.2,aes(x=BirthYear,y=IRR))+geom_line()+
  geom_hline(yintercept = 1,linetype=2)+
  geom_ribbon(aes(ymin=IRR.lo,ymax=IRR.hi),alpha=0.2,data=aux.all.2)+
  coord_cartesian(ylim=c(0,2),xlim=c(1910,1980))+
  labs(x="Birth year",y="IRR, females vs. males")

font_add_google("Rosario", "rosario")
showtext_auto()

pastel_colors <- brewer.pal(n = 3, name = "Set1")

# Plot
unadjusted_fig<- ggplot(aux.all.2, aes(x = BirthYear, y = IRR)) +
  geom_ribbon(aes(ymin = IRR.lo, ymax = IRR.hi),
              fill = pastel_colors[2], alpha = 0.3) +
  geom_line(color = pastel_colors[1], size = 1.2) +
  scale_x_continuous(breaks = seq(1910, 1980, by = 10)) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", size = 0.6) +
  annotate("text", x = 1910, y = 2, label = "A", hjust = 0, vjust = 1.2, size = 6, fontface = "bold") + 
  coord_cartesian(ylim = c(0, 2), xlim = c(1910, 1980)) +
  labs(
    x = "Birth year",
    y = "Incidence rate ratio (women vs. men)",
  ) +
  theme_minimal(base_family = "rosario") +
  theme(
    text = element_text(color = "gray10"),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 13),
    plot.caption = element_text(size = 10, color = "gray40", hjust = 1),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

unadjusted_fig
adjusted_fig

combined_fig <- grid.arrange(unadjusted_fig, adjusted_fig, ncol = 2)
