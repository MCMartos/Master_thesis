# ANOVA test
# Maria del Carmen Martos Contreras: mariadelcarmen.martos@e-campus.uab.cat
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(ggpubr)
library(rstatix)
library(WRS2)
library(tidyverse)
library(vegan)
library(dplyr)

AR_T<-read.csv('C:/Users/mcarm/Desktop/2019_2020/Master_bioinformatics/06_Module/02_metagenomics/Results/P_a_reads.csv', header = T, sep = ';')

##Import data
# AR<-read.csv('C:/Users/mcarm/Desktop/2019_2020/Master_bioinformatics/06_Module/02_metagenomics/Results/ART/Informative_reads.csv', header = T, sep = ';')
# AR_1<-read.csv('C:/Users/mcarm/Desktop/2019_2020/Master_bioinformatics/06_Module/02_metagenomics/Results/ISS/Informative_reads.csv', header = T, sep = ';')
# AR_2<-read.csv('C:/Users/mcarm/Desktop/2019_2020/Master_bioinformatics/06_Module/02_metagenomics/Results/REAL/Informative_reads.csv', header = T, sep = ';')
# Simulator<-rep("ART", 132)
# Simulator<- as.data.frame(Simulator)
# AR<-cbind(AR, Simulator)
# Simulator<-rep("ISS", 132)
# Simulator<- as.data.frame(Simulator)
# AR_1<-cbind(AR_1, Simulator)
# Simulator<-rep("REAL", 132)
# Simulator<- as.data.frame(Simulator)
# AR_2<-cbind(AR_2, Simulator)
# AR_T<-full_join(AR, AR_1)
# AR_T<-full_join(AR_T, AR_2)

##Subset of data
G99D98 <- AR_T[ which(AR_T$G.D=="99_98"), ]
G99D97 <- AR_T[ which(AR_T$G.D=="99_97"), ]
G99D96 <- AR_T[ which(AR_T$G.D=="99_96"), ]
G98D97 <- AR_T[ which(AR_T$G.D=="98_97"), ]
G98D96 <- AR_T[ which(AR_T$G.D=="98_96"), ]
G97D96 <- AR_T[ which(AR_T$G.D=="97_96"), ]

print(G99D98 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))

print(G99D97 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))
print(G99D96 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))
print(G98D97 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))

print(G98D96 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))

print(G97D96 %>%
        group_by(as.factor(Simulator)) %>%
        summarise(Min = min(P_reads_sp),
                  Max = max(P_reads_sp),
                  Median = median(P_reads_sp),
                  Mean = mean(P_reads_sp),
                  IQRange = IQR(P_reads_sp)))


##Correlation test
correlation <- function(data){
  p<-c()
  x<-cor.test(data[ which(data$Simulator=="ISS"), ]$P_reads_sp, data[ which(data$Simulator=="ART"), ]$P_reads_sp , method = "pearson")
  p<-c(p,x$p.value)
  y<-cor.test(data[ which(data$Simulator=="ISS"), ]$P_reads_sp, data[ which(data$Simulator=="REAL"), ]$P_reads_sp , method = "pearson")
  p<-c(p,y$p.value)
  z<-cor.test(data[ which(data$Simulator=="ART"), ]$P_reads_sp, data[ which(data$Simulator=="REAL"), ]$P_reads_sp , method = "pearson")
  p<-c(p,z$p.value)
  print(x)
  print(y)
  print(z)
  return(p)
}

 p<-c()
 G99D98_cor <- correlation(G99D98)
 p<-c(p,G99D98_cor)
 G99D97_cor <- correlation(G99D97)
 p<-c(p,G99D97_cor)
 G99D96_cor <- correlation(G99D96)
 p<-c(p,G99D96_cor)
 G98D97_cor <- correlation(G98D97)
 p<-c(p,G98D97_cor)
 G98D96_cor <- correlation(G98D96)
 p<-c(p,G98D96_cor)
 G97D96_cor <- correlation(G97D96)
 p<-c(p,G97D96_cor)
 
cor_yes=0
for(i in p){
  if(i >=0.05){
    cor_yes=cor_yes+1
  }else{
    cor_yes=cor_yes+0
  }
}
if (cor_yes==length(p)){
  cor=T
}else{
  cor=F
}

##Normality test Shapiro-Wilk test

normal<- function(data){
  p<-c()
  x<-shapiro.test(data[ which(data$Simulator=="ISS"), ]$P_reads_sp)
  p<-c(p,x$p.value)
  y<-shapiro.test(data[ which(data$Simulator=="ART"), ]$P_reads_sp)
  p<-c(p,y$p.value)
  z<-shapiro.test(data[ which(data$Simulator=="REAL"), ]$P_reads_sp)
  p<-c(p,z$p.value)
  print(x)
  print(y)
  print(z)
  return(p)
}

p<-c()
G99D98_nor <- normal(G99D98)
p<-c(p,G99D98_nor)
G99D97_nor <- normal(G99D97)
p<-c(p,G99D97_nor)
G99D96_nor <- normal(G99D96)
p<-c(p,G99D96_nor)
G98D97_nor <- normal(G98D97)
p<-c(p,G98D97_nor)
G98D96_nor <- normal(G98D96)
p<-c(p,G98D96_nor)
G97D96_nor <- normal(G97D96)
p<-c(p,G97D96_nor)

nor_yes=0
for(i in p){
  if(i >=0.05 ){
    nor_yes=nor_yes+1
  }else{
    nor_yes=nor_yes+0
  }
}
if (nor_yes==length(p)){
  nor=T
}else{
  nor=F
}

#Normality plot
plote<- function(datos){
  ggplot(data = datos, P_reads_spping = aes(x = P_reads_sp, colour = Simulator)) +
    geom_histogram() +
    theme_bw() +
    facet_grid(. ~ Simulator) +
    theme(legend.position = "none")}

# plote(G99D98)
plote(AR_T)
##Homocedasticity test
# p<-c()
x<-fligner.test(AR_T$P_reads_sp ~ AR_T$Simulator)
print(x)
p<-x$p.value
x<-fligner.test(G99D98$P_reads_sp ~ G99D98$Simulator)
print(x)
p<-c(p,x$p.value)
x<-fligner.test(G99D97$P_reads_sp ~ G99D97$Simulator)
print(x)
p<-c(p,x$p.value)
x<-fligner.test(G99D96$P_reads_sp ~ G99D96$Simulator)
print(x)
p<-c(p,x$p.value)
x<-fligner.test(G98D97$P_reads_sp ~ G98D97$Simulator)
print(x)
p<-c(p,x$p.value)
x<-fligner.test(G98D96$P_reads_sp ~ G98D96$Simulator)
print(x)
p<-c(p,x$p.value)
x<-fligner.test(G97D96$P_reads_sp ~ G97D96$Simulator)
print(x)
p<-c(p,x$p.value)

hom_yes=0
for(i in p){
  if(i >=0.05 ){
    hom_yes=hom_yes+1
  }else{
    hom_yes=hom_yes+0
  }
}
if (hom_yes==length(p)){
  hom=T
}else{
  hom=F
}

##ANOVA test choice
test<-function(data){
  if (cor==F & nor==F & hom==F){
    stat.test <- data %>%
      group_by(G.D) %>%
      games_howell_test(P_reads_sp ~ Simulator) %>%
      add_significance("p.adj")
    print(stat.test)
    
     jpeg(filename="C:/Users/mcarm/Desktop/2019_2020/Master_bioinformatics/06_Module/02_metagenomics/Results/Assigned_reads_mapped.jpeg" ,res=300, width=2000, height=1500)
    bxp <- ggboxplot(
      AR_T, x = "G.D", y = "P_reads_sp", 
      color = "Simulator", palette = c("#00AFBB", "#E7B800", "#4e0041")
    )
    bxp <- bxp + labs(x = "", y = "Assigned mapped reads")
    # Add p-values onto the box plots
    stat.test <- stat.test %>%
      add_xy_position(x = "G.D", dodge = 0.8)
    bxp <- bxp + stat_pvalue_manual(
      stat.test, step.group.by = "G.D", label = "p.adj",
      tip.length = 0, step.increase = 0.1,bracket.nudge.y = 0.001
    )
    print(bxp)
     dev.off()
    
  }else if(cor==F & nor==F & hom==T){
    stat_test <- data %>%
      group_by(G.D) %>%
      kruskal_test(P_reads_sp ~ Simulator)
    print(stat_test)
    stat.test1 <- data %>%
      group_by(G.D) %>%
      dunn_test(P_reads_sp ~ Simulator)
    print(stat.test1)
    jpeg(filename="Assigned_reads_mapped.jpeg" ,res=300, width=2000, height=1500)
    bxp <- ggboxplot(
      data, x = "G.D", y = "P_reads_sp",
      color = "Simulator", palette = c("#00AFBB", "#E7B800", "#4e0041")
    )
    bxp <- bxp + labs(x = "Assigned mapped reads", y = "")
    # Add p-values onto the box plots
    stat.test1 <- stat.test1 %>%
      add_xy_position(x = "G.D", dodge = 0.8)
    bxp <- bxp + stat_pvalue_manual(
      stat.test1, step.group.by = "G.D", #label = "p.adj",
      tip.length = 0, step.increase = 0.1,bracket.nudge.y = 0.001
    )
    print(bxp)
    dev.off()
    
  }else if(cor==T){
    print("Samples aren't independents")
  }else if(cor==T & nor==T & hom==T){
    print("Try an anova test")
  }else{
    print("ERROR")
  }
}
test(AR_T)