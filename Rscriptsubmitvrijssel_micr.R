rm(list=ls())

#install.packages("devtools")
#library(devtools)
#install_github("vmikk/metagMisc")
#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#biocLite("genefilter")
library(data.table)
library("ape")
library("car")
library(plyr)
library("dplyr")
library("scales")
library("grid")
library(ade4)
library(tidyr)
library(jsonlite)   
library(survival)           
library(stringi)
library(phyloseq)

library(vegan)

library(corrplot)
#library(ggpubr)
#library(ggtern)

library(RColorBrewer)
library(gridExtra)
library("permute")
library("lattice")
library("cowplot")
library("reshape2")
library("VennDiagram")
library(lmerTest)
library(geosphere)
library(gridExtra)
#library(plotrix)
library(ggplot2)

set.seed(123)

#################### WORKING DRIVE ###################################################################################
#setwd("C:\\Users\\sophi\\Documents\\PhD") #Sophiethuis

#setwd("N:\\Dep.TE\\TE-Share\\Vital Soils\\Sophie\\Chapter 1 Sophie Bacterial and fungal communities\\analysis R and data")
setwd("C:\\Users\\sophi\\Documents\\PhD\\Chapter 1 Sophie Bacterial and fungal communities\\analysis R and data")
data_phylum4<-read.table("dataframe.txt")

#################### FUNCTIONS #######################################################################################

merge_samples_mean <- function(physeq, group){
  #Calculate the number of samples in each group
  group_sums <- as.matrix(table(sample_data(physeq)[ ,group]))[,1]
  #Merge samples by summing
  merged <- merge_samples(physeq, group)
  #Divide summed OTU counts by number of samples in each group to get mean
  #Calculation is done while taxa are columns, but then transposed at the end
  x <- as.matrix(otu_table(merged)) 
  if(taxa_are_rows(merged)){ x<-t(x) }
  out <- t(x/group_sums)
  #Return new phyloseq object with taxa as rows
  out <- otu_table(out, taxa_are_rows = TRUE)
  otu_table(merged) <- out
  return(merged)
}
######################### P H Y L O S E Q for 16S ####################################################################
Chronodata3<-read.csv("data.csv",header=T,row.names=1,sep=";",text=text,dec=",",stringsAsFactors=FALSE,check.names = FALSE)

Chronodata3<-Chronodata3[rownames(Chronodata3)!="172",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="173",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="174",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="31",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="32",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="33",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="49",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="50",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="51",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="175",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="176",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="177",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="34",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="35",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="36",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="52",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="53",]
Chronodata3<-Chronodata3[rownames(Chronodata3)!="54",]


Chronodata3$percSOM<-Chronodata3$SOM*100

Chronodata3$pH.std<-scale(Chronodata3$pH)
Chronodata3$SOM.std<-scale(Chronodata3$SOM)
Chronodata3$Soil_Origin<-factor(Chronodata3$Soil_Origin,levels=c("Clay Middle","Clay SOuth West","Sand North","Sand Middle","Sand South"))

OTU.16Stable <- read.csv("16SOTUtable.csv",header=TRUE,row.names=1,sep=";",text=text,dec=".",stringsAsFactors=FALSE,check.names = FALSE)
TAX.16Stable <- read.csv("16STAXtable.csv",header=TRUE,row.names=1,sep=";",text=text,stringsAsFactors=FALSE,check.names = FALSE) 
OTU.16S <- as.matrix(OTU.16Stable) 
TAX.16S <- as.matrix(TAX.16Stable) 

OTU16SSamples <- Chronodata3

OTU.16S = otu_table(OTU.16S, taxa_are_rows = TRUE)
TAX.16S = tax_table(TAX.16S)
ENV.16S = sample_data(OTU16SSamples)


physeq16S = phyloseq(OTU.16S, TAX.16S, ENV.16S)
physeq16S
physeq16S= subset_samples(physeq16S,  sample_names(physeq16S)!= "81")
physeq16S= subset_samples(physeq16S,  sample_names(physeq16S)!= "193")
physeq16S= subset_samples(physeq16S,  sample_names(physeq16S)!= "194")
physeq16S= subset_samples(physeq16S,  sample_names(physeq16S)!= "195")
physeq16S= subset_samples(physeq16S,  sample_names(physeq16S)!= "189")

physeq16S = subset_taxa(physeq16S, Kingdom == "k__Archaea" | Kingdom == "k__Bacteria") 
physeq16S
physeq16S = prune_samples(sample_sums(physeq16S)>=5000, physeq16S)
physeq16S
physeq16S = prune_taxa(taxa_sums(physeq16S) > 3, physeq16S)
physeq16S
physeq16Srar = rarefy_even_depth(physeq16S)
physeq16Srar

OTU.ITStable <- read.csv("ITSOTUtable.csv",header=TRUE,row.names=1,sep=";",text=text,dec=".",stringsAsFactors=FALSE,check.names = FALSE)
TAX.ITStable <- read.csv("ITSTAXtablefunguild.csv",header=TRUE,row.names=1,sep=";",text=text,stringsAsFactors=FALSE,check.names = FALSE) 

OTU.ITS <- as.matrix(OTU.ITStable) 
TAX.ITS <- as.matrix(TAX.ITStable) 

#listofOTUs16S<-row.names(OTU.16Stable)
samplesITS<-colnames(OTU.ITStable)
OTUITSSamples <- Chronodata3

OTU.ITS = otu_table(OTU.ITS, taxa_are_rows = TRUE)
TAX.ITS = tax_table(TAX.ITS)

ENV.ITS = sample_data(OTUITSSamples)

physeqITS = phyloseq(OTU.ITS, TAX.ITS, ENV.ITS)
physeqITS

physeqITS= subset_samples(physeqITS,  sample_names(physeqITS)!= "193")
physeqITS= subset_samples(physeqITS,  sample_names(physeqITS)!= "194")
physeqITS= subset_samples(physeqITS,  sample_names(physeqITS)!= "195")

physeqITS2 = subset_taxa(physeqITS, Kingdom == "k__Fungi")
physeqITS2 = prune_samples(sample_sums(physeqITS2)>=1500, physeqITS2)
physeqITS2 = prune_taxa(taxa_sums(physeqITS2) > 3, physeqITS2)


physeqITSrar = rarefy_even_depth(physeqITS2)
physeqITSrar



#####Group physeqs#####
physeqITS2
physeq16Sbac = subset_taxa(physeq16S, Kingdom == "k__Bacteria") 
physeq16Sarch = subset_taxa(physeq16S, Kingdom == "k__Archaea" ) 

Fungi.reads = sample_sums(physeqITS2)
Bacteria.reads = sample_sums(physeq16Sbac)
Archaea.reads = sample_sums(physeq16Sarch)

PhyseqBacteria<-subset_taxa(physeq16Srar,Kingdom == "k__Bacteria")
PhyseqArchaea<-subset_taxa(physeq16Srar,Kingdom == "k__Archaea")
PhyseqFungi<-physeqITSrar
sampledata.Bacteria<-data.frame(sample_data(PhyseqBacteria))
sampledata.Archaea<-data.frame(sample_data(PhyseqArchaea))
sampledata.Fungi<-data.frame(sample_data(PhyseqFungi))


Bacteria.div<-estimate_richness(PhyseqBacteria,measures=c("Observed","Shannon"))
sampledata.Bacteria<-cbind(sampledata.Bacteria,Bacteria.div)
sampledata.Bacteriareads<-cbind(sampledata.Bacteria,Bacteria.div,Bacteria.reads)
Archaea.div<-estimate_richness(PhyseqArchaea,measures=c("Observed","Shannon"))
sampledata.Archaea<-cbind(sampledata.Archaea,Archaea.div)
sampledata.Archaeareads<-cbind(sampledata.Archaea,Archaea.div,Archaea.reads)
Fungi.div<-estimate_richness(PhyseqFungi,measures=c("Observed","Shannon"))
sampledata.Fungi<-cbind(sampledata.Fungi,Fungi.div)
sampledata.Fungireads<-cbind(sampledata.Fungi,Fungi.div,Fungi.reads)


#####no of reads#####
aggregate(Bacteria.reads~Soil,data=sampledata.Bacteriareads,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Archaea.reads~Soil,data=sampledata.Archaeareads,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Fungi.reads~Soil,data=sampledata.Fungireads,function(x) c(mean = mean(x), sd = sd(x)))

aggregate(Bacteria.reads~Soil*ConOrg,data=sampledata.Bacteriareads,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Archaea.reads~Soil*ConOrg,data=sampledata.Archaeareads,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Fungi.reads~Soil*ConOrg,data=sampledata.Fungireads,function(x) c(mean = mean(x), sd = sd(x)))

anova(lmer(Bacteria.reads~Soil*ConOrg+(1|Fieldnr),data=sampledata.Bacteriareads))
anova(lmer(Archaea.reads~Soil*ConOrg+(1|Fieldnr),data=sampledata.Archaeareads))
anova(lmer(Fungi.reads~Soil*ConOrg+(1|Fieldnr),data=sampledata.Fungireads))


#####pH, SOM, Organic inputs#####
anova(lmer(SOM~Soil*ConOrg+(1|Fieldnr),data=Chronodata3))
anova(lmer(pH~Soil*ConOrg+(1|Fieldnr),data=Chronodata3))
anova(lmer(external.om.inputs~Soil*ConOrg+(1|Fieldnr),data=Chronodata3))
summary(aov(SOM~Soil*ConOrg+Error(Fieldnr),data=Chronodata3))
summary(aov(pH~Soil*ConOrg+Error(Fieldnr),data=Chronodata3))
summary(aov(external.om.inputs~Soil*ConOrg+Error(Fieldnr),data=Chronodata3))

aggregate(SOM~Soil,data=Chronodata3,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(pH~Soil,data=Chronodata3,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(SOM~Soil*ConOrg,data=Chronodata3,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(pH~Soil*ConOrg,data=Chronodata3,function(x) c(mean = mean(x), sd = sd(x)))


ggplot(Chronodata3,aes(x=Soil,y=SOM,fill=ConOrg))+
  geom_boxplot()+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y='Soil organic matter (%)') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggplot(Chronodata3,aes(x=Soil,y=pH,fill=ConOrg))+
  geom_boxplot()+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y='pH (%)') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())






#####diversity clay vs sand and conventional vs organic#####
aggregate(Observed~Soil,data=sampledata.Bacteria,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Observed~Soil,data=sampledata.Archaea,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Observed~Soil,data=sampledata.Fungi,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil,data=sampledata.Bacteria,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil,data=sampledata.Archaea,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil,data=sampledata.Fungi,function(x) c(mean = mean(x), sd = sd(x)))

aggregate(Observed~Soil*ConOrg,data=sampledata.Bacteria,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Observed~Soil*ConOrg,data=sampledata.Archaea,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Observed~Soil*ConOrg,data=sampledata.Fungi,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil*ConOrg,data=sampledata.Bacteria,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil*ConOrg,data=sampledata.Archaea,function(x) c(mean = mean(x), sd = sd(x)))
aggregate(Shannon~Soil*ConOrg,data=sampledata.Fungi,function(x) c(mean = mean(x), sd = sd(x)))

anova(lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Fungi))

Mb1<-lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Bacteria)
Mb2<-lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Archaea)
Mb3<-lmer(Observed~Soil*ConOrg+(1|Fieldnr),data=sampledata.Fungi)
Mb4<-lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Bacteria)
Mb5<-lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Archaea)
Mb6<-lmer(Shannon~Soil*ConOrg+(1|Fieldnr),data=sampledata.Fungi)
#testmetpaarerin
anova(lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Bacteria))
anova(lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Archaea))
anova(lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Fungi))
Mbp1<-lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Bacteria)
Mbp2<-lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Archaea)
Mbp3<-lmer(Observed~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Fungi)
anova(lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Bacteria))
anova(lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Archaea))
anova(lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Fungi))
Mbp4<-lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Bacteria)
Mbp5<-lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Archaea)
Mbp6<-lmer(Shannon~Soil*ConOrg+(1|Pair/Fieldnr),data=sampledata.Fungi)
anova(Mb1,Mbp1)
anova(Mb2,Mbp2)
anova(Mb3,Mbp3)
anova(Mb4,Mbp4)
anova(Mb5,Mbp5)
anova(Mb6,Mbp6)

anova(lmer(Observed~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Observed~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Observed~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(Shannon~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Shannon~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Shannon~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Fungi))

anova(lmer(Observed~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Observed~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Observed~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(Shannon~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Shannon~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Shannon~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Fungi))




#####general NMDS & diversity#####
Physeq.ordBacteria <- ordinate(PhyseqBacteria, "NMDS", "bray", k = 2,trymax=100)
Physeq.ordBacteria
Physeq.ordArchaea <- ordinate(PhyseqArchaea, "NMDS", "bray", k = 3,trymax=100)
Physeq.ordArchaea
Physeq.ordFungi <- ordinate(PhyseqFungi, "NMDS", "bray", k = 3,trymax=100)
Physeq.ordFungi

sampledata.Bacteria$NMDS1<-Physeq.ordBacteria$points[,1]
sampledata.Bacteria$NMDS2<-Physeq.ordBacteria$points[,2]
sampledata.Bacteria$Group<-"Bacteria"
sampledata.Archaea$NMDS1<-Physeq.ordArchaea$points[,1]
sampledata.Archaea$NMDS2<-Physeq.ordArchaea$points[,2]
sampledata.Archaea$Group<-"Archaea"
sampledata.Fungi$NMDS1<-Physeq.ordFungi$points[,1]
sampledata.Fungi$NMDS2<-Physeq.ordFungi$points[,2]
sampledata.Fungi$Group<-"Fungi"
#names(sampledata.Bacteria)[names(sampledata.Bacteria)=="Bacteria.reads"]<-"reads"
#names(sampledata.Archaea)[names(sampledata.Archaea)=="Archaea.reads"]<-"reads"
#names(sampledata.Fungi)[names(sampledata.Fungi)=="Fungi.reads"]<-"reads"

datafig1<-rbind(sampledata.Bacteria,sampledata.Archaea,sampledata.Fungi)
datafig1$Group<-factor(datafig1$Group,levels=c("Bacteria","Archaea","Fungi"))

#####compositional lmers NMDS1&2#####
anova(lmer(NMDS1~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS1~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS1~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(NMDS2~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS2~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS2~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=sampledata.Fungi))

anova(lmer(NMDS1~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS1~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS1~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(NMDS2~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS2~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS2~Soil*ConOrg*(pH.std+SOM.std)+(1|Fieldnr),data=sampledata.Fungi))

anova(lmer(Observed~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Observed~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Observed~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(Shannon~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(Shannon~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(Shannon~Soil*ConOrg*NMDS2+(1|Fieldnr),data=sampledata.Fungi))

anova(lmer(NMDS1~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS1~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS1~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Fungi))
anova(lmer(NMDS2~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Bacteria))
anova(lmer(NMDS2~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Archaea))
anova(lmer(NMDS2~Soil*ConOrg*(lat.num+lon.num)+(1|Fieldnr),data=sampledata.Fungi))


#####perMANOVAs#####

sampledata.Bacteriasub<-data.frame(sample_data(PhyseqBacteriasub))
sampledata.Archaeasub<-data.frame(sample_data(PhyseqArchaeasub))
sampledata.Fungisub<-data.frame(sample_data(PhyseqFungisub))
OTU.Bacteriasub<-as.data.frame(t(otu_table(PhyseqBacteriasub)))
OTU.Archaeasub<-as.data.frame(t(otu_table(PhyseqArchaeasub)))
OTU.Fungisub<-as.data.frame(t(otu_table(PhyseqFungisub)))

#betadispersion in sand vs clay
betadisper(vegdist(OTU.Bacteriasub, method="bray"),sampledata.Bacteriasub$Soil)
betadisper(vegdist(OTU.Archaeasub, method="bray"),sampledata.Archaeasub$Soil)
betadisper(vegdist(OTU.Fungisub, method="bray"),sampledata.Fungisub$Soil)

#permanova soil origin
adonis((vegdist(OTU.Bacteriasub, method="bray"))~sampledata.Bacteriasub$Fieldnr+sampledata.Bacteriasub$Soil_Origin)
adonis((vegdist(OTU.Archaeasub, method="bray"))~sampledata.Archaeasub$Fieldnr+sampledata.Archaeasub$Soil_Origin)
adonis((vegdist(OTU.Fungisub, method="bray"))~sampledata.Fungisub$Fieldnr+sampledata.Fungisub$Soil_Origin)


#pairwise permanova; whether effects of management are stronger/less strong in diff. soil types
adonis((vegdist(OTU.Bacteriasub, method="bray"))~sampledata.Bacteriasub$Fieldnr+sampledata.Bacteriasub$Soil_Origin)
adonis((vegdist(OTU.Archaeasub, method="bray"))~sampledata.Archaeasub$Fieldnr+sampledata.Archaeasub$Soil_Origin)
adonis((vegdist(OTU.Fungisub, method="bray"))~sampledata.Fungisub$Fieldnr+sampledata.Fungisub$Soil_Origin)
pairwise.adonis((vegdist(OTU.Bacteriasub, method="bray")), fact, test = c("Pillai", "Wilks",
                                                                          "Hotelling-Lawley", "Roy", "Spherical"), nperm = 999, 
                progress = TRUE, p.method = "fdr", F = FALSE, R2 = FALSE)

#permanova soil&management
adonis((vegdist(OTU.Bacteriasub, method="bray"))~sampledata.Bacteriasub$Fieldnr+sampledata.Bacteriasub$Soil*sampledata.Bacteriasub$ConOrg)
adonis((vegdist(OTU.Archaeasub, method="bray"))~sampledata.Archaeasub$Fieldnr+sampledata.Archaeasub$Soil*sampledata.Archaeasub$ConOrg)
adonis((vegdist(OTU.Fungisub, method="bray"))~sampledata.Fungisub$Fieldnr+sampledata.Fungisub$Soil*sampledata.Fungisub$ConOrg)

#permanova time
adonis((vegdist(OTU.Bacteriasub, method="bray"))~sampledata.Bacteriasub$Fieldnr+sampledata.Bacteriasub$Soil*sampledata.Bacteriasub$ConOrg*sampledata.Bacteriasub$Yearssince.as.if)
adonis((vegdist(OTU.Archaeasub, method="bray"))~sampledata.Archaeasub$Fieldnr+sampledata.Archaeasub$Soil*sampledata.Archaeasub$ConOrg*sampledata.Archaeasub$Yearssince.as.if)
adonis((vegdist(OTU.Fungisub, method="bray"))~sampledata.Fungisub$Fieldnr+sampledata.Fungisub$Soil*sampledata.Fungisub$ConOrg*sampledata.Fungisub$Yearssince.as.if)

#permanova other factors
adonis((vegdist(OTU.Bacteriasub, method="bray"))~sampledata.Bacteriasub$Fieldnr+sampledata.Bacteriasub$Soil*sampledata.Bacteriasub$ConOrg*(sampledata.Bacteriasub$pH.std+sampledata.Bacteriasub$SOM.std))
adonis((vegdist(OTU.Archaeasub, method="bray"))~sampledata.Archaeasub$Fieldnr+sampledata.Archaeasub$Soil*sampledata.Archaeasub$ConOrg*(sampledata.Archaeasub$pH.std+sampledata.Archaeasub$SOM.std))
adonis((vegdist(OTU.Fungisub, method="bray"))~sampledata.Fungisub$Fieldnr+sampledata.Fungisub$Soil*sampledata.Fungisub$ConOrg*(sampledata.Fungisub$pH.std+sampledata.Fungisub$SOM.std))



#####conorg braycurtis dissimilarities & time#####
#1 create matrix of braycurtis
#2 melt to long format
#3 select paired data
#4 get time data
#5 plot

bcdisbac<-distance(PhyseqBacteria,method="bray",type="samples",diag = T,upper = T)
bcdisarch<-distance(PhyseqArchaea,method="bray",type="samples",diag = T,upper = T)
bcdisfungi<-distance(PhyseqFungi,method="bray",type="samples",diag = T,upper = T)



m1 <- as.matrix(bcdisbac)
m3 <- as.matrix(bcdisarch)
m5 <- as.matrix(bcdisfungi)


m2.1 <- melt(m1)[melt(upper.tri(m1))$value,]
names(m2.1) <- c("c1", "c2", "distance")

m2.3 <- melt(m3)[melt(upper.tri(m3))$value,]
names(m2.3) <- c("c1", "c2", "distance")

m2.5 <- melt(m5)[melt(upper.tri(m5))$value,]
names(m2.5) <- c("c1", "c2", "distance")

m2.1$Group<-rep("Bacteria")

m2.3$Group<-rep("Archaea")

m2.5$Group<-rep("Fungi")

braycurtislongformat<-rbind(m2.1,m2.3,m2.5)
#write.csv(braycurtislongformat,"braycurtislongformatnew.csv")
#edit in excel: select paired fields only
bcpairs<-read.csv("braycurtispairsnew.csv",header=T,sep=";",text=text,dec=",",stringsAsFactors=FALSE,check.names = FALSE)

bcpairs.av<-aggregate(bcpairs$distance,by=list(bcpairs$Pair,bcpairs$Group,bcpairs$Soil,bcpairs$Yearssince),na.rm=TRUE,FUN="mean")
names(bcpairs.av) <- c("Pair", "Group", "Soil", "Yearssince","distance")

bcpairs.av$Group<-factor(bcpairs.av$Group,levels=c("Bacteria","Archaea","Fungi"))

bcpairsarchclay.av<-subset(bcpairs.av,Group=="Archaea" & Soil=="Clay")
bcpairsarchsand.av<-subset(bcpairs.av,Group=="Archaea" & Soil=="Sand")
bcpairsbacclay.av<-subset(bcpairs.av,Group=="Bacteria" & Soil=="Clay")
bcpairsbacsand.av<-subset(bcpairs.av,Group=="Bacteria" & Soil=="Sand")
bcpairsfungiclay.av<-subset(bcpairs.av,Group=="Fungi" & Soil=="Clay")
bcpairsfungisand.av<-subset(bcpairs.av,Group=="Fungi" & Soil=="Sand")

summary(lm(distance~Yearssince,data=bcpairsarchclay.av))
summary(lm(distance~Yearssince,data=bcpairsarchsand.av))
summary(lm(distance~Yearssince,data=bcpairsbacclay.av))
summary(lm(distance~Yearssince,data=bcpairsbacsand.av))
summary(lm(distance~Yearssince,data=bcpairsfungiclay.av))
summary(lm(distance~Yearssince,data=bcpairsfungisand.av))

ggplot(bcpairs.av,aes(x=Yearssince,y=distance,color=Soil))+
  geom_point()+
  xlim(0,25)+
  #geom_smooth(method="lm")+
  facet_wrap(vars(Soil,Group),scales="free",ncol=3,dir="h")+
  scale_color_manual(values = c("#336666","#CC9966"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15),
        strip.text=element_blank()
  )
ggsave("bctime.pdf",width = 12.4, height = 9.3)
#####correlations year#####
#these dont work since the model ran out of degrees of freedom (due to confounded factors?)
anova(lmer(Yearssince.as.if~Soil*ConOrg*SOM+(1|Fieldnr),data=Chronodata3))
anova(lmer(Yearssince.as.if~Soil*ConOrg*pH+(1|Fieldnr),data=Chronodata3))
anova(lmer(Yearssince.as.if~Soil*ConOrg*lat.num+(1|Fieldnr),data=Chronodata3))
anova(lmer(Yearssince.as.if~Soil*ConOrg*lon.num+(1|Fieldnr),data=Chronodata3))

#but shouldnt these be the other way around? yearssince as effect, not as response...? That doesnt give problems seems..
anova(lmer(SOM~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(pH~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(lat.num~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(lon.num~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(external.om.inputs.std~Soil*ConOrg*Yearssince.as.if+(1|Fieldnr),data=Chronodata3)) #so this one should be done per field, not per sample, right?!


anova(lmer(SOM~Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(pH~Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(lat.num~Yearssince.as.if+(1|Fieldnr),data=Chronodata3))
anova(lmer(lon.num~Yearssince.as.if+(1|Fieldnr),data=Chronodata3))


ggplot(datafig1,aes(x=Yearssince.as.if,y=SOM,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~.,scales="free",ncol=1)+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Soil organic matter (%)') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("yearsom.pdf",width = 5, height = 9.3)
ggplot(datafig1,aes(x=Yearssince.as.if,y=pH,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~.,scales="free",ncol=1)+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='pH') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("yearph.pdf",width = 5, height = 9.3)



ggplot(datafig1,aes(x=Yearssince.as.if,y=lat.num,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~.,scales="free",ncol=1)+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Latitude') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("yearlat.pdf",width = 5, height = 9.3)
ggplot(datafig1,aes(x=Yearssince.as.if,y=lon.num,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~.,scales="free",ncol=1)+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Longitude') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("yearlon.pdf",width = 5, height = 9.3)




#####new plots black&white#####
ggplot(datafig1,aes(x=NMDS1,y=NMDS2,shape=ConOrg,color=Soil))+
  geom_point(size=2)+
  #xlim(-1.5,1.5)+
  facet_wrap(.~Group,scales="free")+
  scale_color_manual(values = c("#336666","#CC9966"))+
  scale_shape_manual(values = c(21,16))+
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1.pdf",width = 12.4, height = 4.65)
ggplot(datafig1,aes(x=NMDS1,y=NMDS2,color=Soil_Origin))+
  geom_point(size=2)+
  facet_wrap(.~Group,scales="free")+
  #scale_color_manual(values = c("#99CCCC","#336666","#FFCC99","#CC9966"))+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmdsorigin.pdf",width = 12.4, height = 4.65)


ggplot(datafig1,aes(x=ConOrg,y=Observed,fill=Soil_Management))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  facet_wrap(Soil~Group,scales="free")+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y='Observed OTU richness') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("richnessbox.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=ConOrg,y=Shannon,fill=Soil_Management))+
  geom_boxplot()+
  geom_jitter(position=position_jitter(width=.1, height=0))+
  facet_wrap(Soil~Group,scales="free")+
  scale_fill_manual(values = c("white","#336666","white","#CC9966"))+
  labs(x="Management",y="Shannon's diversity index") + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("shannonbox.pdf",width = 12.4, height = 9.3)


ggplot(datafig1,aes(x=Yearssince.as.if,y=Observed,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Observed OTU richness') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())

#new try
ggplot(datafig1,aes(x=Yearssince.as.if,y=Observed,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Observed OTU richness') + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15),
        strip.text=element_blank()
  )

ggplot(datafig1,aes(x=Yearssince.as.if,y=Observed,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='Observed OTU richness') + 
  theme_bw()+
  theme(legend.position = "none",
        axis.title=element_text(size=20,face="bold"),
        axis.text=element_text(size=15),
        strip.text=element_blank()
  )

ggsave("richnesstime.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=Yearssince.as.if,y=Shannon,shape=ConOrg,linetype=ConOrg,color="black"))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y="Shannon's diversity index") + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("shannontime.pdf",width = 12.4, height = 9.3)


ggplot(datafig1,aes(x=NMDS2,y=Shannon,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  facet_wrap(Soil~Group,scales="free")+
  labs(x="NMDS2",y="Shannon's diversity index") + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("shannonnmds2.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=NMDS2,y=Observed,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  facet_wrap(Soil~Group,scales="free")+
  labs(x="NMDS2",y='Observed OTU richness') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("richnessnmds2.pdf",width = 12.4, height = 9.3)




ggplot(datafig1,aes(x=Yearssince.as.if,y=NMDS1,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  xlim(0,25)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='NMDS1') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1time.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=Yearssince.as.if,y=NMDS2,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  xlim(0,25)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Years of organic management",y='NMDS2') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds2time.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=percSOM,y=NMDS1,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Soil Organic Matter (%)",y='NMDS1') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1som.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=percSOM,y=NMDS2,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="Soil Organic Matter (%)",y='NMDS2') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds2som.pdf",width = 12.4, height = 9.3)

ggplot(datafig1,aes(x=pH,y=NMDS1,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="pH",y='NMDS1') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1ph.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=pH,y=NMDS2,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="pH",y='NMDS2') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds2ph.pdf",width = 12.4, height = 9.3)

ggplot(datafig1,aes(x=lat.num,y=NMDS1,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  #xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="latitude",y='NMDS1') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1lat.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=lat.num,y=NMDS2,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  #xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="latitude",y='NMDS2') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds2lat.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=lon.num,y=NMDS1,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  #xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="longitude",y='NMDS1') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds1lon.pdf",width = 12.4, height = 9.3)
ggplot(datafig1,aes(x=lon.num,y=NMDS2,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  #xlim(0,8)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="longitude",y='NMDS2') + 
  theme_bw()+
  theme(legend.position = "none",axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("nmds2lon.pdf",width = 12.4, height = 9.3)





#####plots diversity ph,som,inputs#####

ggplot(datafig1,aes(x=pH,y=Observed,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="pH",y='Observed') + 
  theme_bw()

ggplot(datafig1,aes(x=SOM,y=Observed,color=Soil,shape=ConOrg,linetype=ConOrg))+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  facet_wrap(Soil~Group,scales="free")+
  scale_color_manual(values = c("black","black"))+
  scale_shape_manual(values = c(21,16))+
  scale_linetype_manual(values=c("longdash","solid"))+
  labs(x="SOM",y='Observed') + 
  theme_bw()


#####composition groups#####

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861")

#bacteria phyla level
data_phylum0<-transform_sample_counts(PhyseqBacteria, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-subset_taxa(data_phylum0,Phylum_selection=="Keep")
data_phylum2<-tax_glom(data_phylum1,taxrank="Phylum",NArm=FALSE)

melt1<-psmelt(data_phylum2)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
modelsHSD <- list()
phylbac<-as.character(unique(melt1_clay$Phylum))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_bac_clay <- significance 


significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_bac_sand <- significance 

data_phylum2<-merge_samples_mean(data_phylum2,"Soil_Management")

data_phylum4 <- psmelt(data_phylum2) # create dataframe from phyloseq object
#data_phylum4 <- subset(data_phylum4, select=c(""))
data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character

# group dataframe by Phylum, calculate median rel. abundance
maximums <- ddply(data_phylum4, ~Phylum, function(x) c(max=max(x$Abundance)))
maximums# find Phyla whose rel. abund. is less than 1%
remainder <- maximums[maximums$max <= 0.01,]$Phylum
# change their name to "Remainder"
data_phylum4[data_phylum4$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacphylum.pdf",width = 15, height = 8)
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"), 
        legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("compbacphylumwhlegend.pdf",width = 6.5, height = 8)



#bacteria genus level
data_phylum0<-transform_sample_counts(PhyseqBacteria, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-subset_taxa(data_phylum0,Genus_selection=="Keep")
data_phylum2<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])

melt1<-psmelt(data_phylum2)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand") 
phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")

#forallgenera
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_bac_clay_all <- significance 

for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_bac_sand_all <- significance 


#fortop10
gentop10  = prune_taxa(phylumtop10,  data_phylum2)
melt1<-psmelt(gentop10)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_bac_clay_top10 <- significance 


significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Genus)), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Genus))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_bac_sand_top10 <- significance 


data_phylum3<-merge_samples_mean(data_phylum2,"Soil_Management")
phylumtop10  = prune_taxa(phylumtop10,  data_phylum3)

data_phylum4 <- psmelt(phylumtop10) # create dataframe from phyloseq object
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacgenus.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacgenuswhlegend.pdf",width = 6.5, height = 8)

#Archaea phyla level
data_phylum0<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-tax_glom(data_phylum0,taxrank="Genus",NArm=FALSE)

melt1<-psmelt(data_phylum1)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_clay$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_clay$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_arch_clay_all <- significance 

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_arch_sand_all <- significance 

data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Management")
data_phylum4 <- psmelt(data_phylum2) # create dataframe from phyloseq object

data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character

# group dataframe by Phylum, calculate median rel. abundance
maximums <- ddply(data_phylum4, ~Phylum, function(x) c(max=max(x$Abundance)))
maximums# find Phyla whose rel. abund. is less than 1%
remainder <- maximums[maximums$max <= 0.01,]$Phylum
# change their name to "Remainder"
data_phylum4[data_phylum4$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("comparchphylum.pdf")
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("comparchphylumwhlegend.pdf",width = 6.5, height = 8)



#Archaea genus level - identified genera
PhyseqArchaearel<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU))
data_phylum1<-subset_taxa(PhyseqArchaearel,Genus_selection=="Keep")
data_phylum1<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)

melt1<-psmelt(data_phylum1)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
modelsHSD <- list()
str(melt1_clay$Abundance)
phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")


for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_arch_clay_all <- significance 

phylbac<-as.character(unique(melt1_sand$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")


for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_arch_sand_all <- significance 

data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Management") #takes the mean of the samples
data_phylum4<-psmelt(data_phylum2)

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchgenusidentified.pdf")
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchgenusidentifiedwhlegend.pdf",width=6.5,height = 8)

#Archaea OTU level - identified species
PhyseqArchaearel<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU))

data_phylum2<-merge_samples_mean(PhyseqArchaearel,"Soil_Management") #takes the mean of the samples
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])
phylumtop10  = prune_taxa(phylumtop10,  data_phylum2)
plot_bar(phylumtop10,fill="OTU")+
  geom_bar(aes(color=OTU, fill=OTU), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchOTUunidentified.pdf")
plot_bar(phylumtop10,fill="OTU")+
  geom_bar(aes(color=OTU, fill=OTU), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchOTUunidentifiedwhlegend.pdf",width = 6.5, height = 8)


#fungi phyla level
data_phylum0<-transform_sample_counts(PhyseqFungi, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-tax_glom(data_phylum0,taxrank="Phylum",NArm=FALSE)
data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Management")
data_phylum4 <- psmelt(data_phylum2) # create dataframe from phyloseq object
data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character

melt1<-psmelt(data_phylum1)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_clay$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_clay$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_fungi_clay <- significance 

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_phylum_fungi_sand <- significance 

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungiphylum.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Management",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungiphylumwhlegend.pdf",width = 6.5, height = 8)

#fungi genus level
data_phylum0<-transform_sample_counts(PhyseqFungi, function(X) X/sum(X))
data_phylum1<-subset_taxa(data_phylum0,Genus_selection=="keep")

data_phylum2<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])

#OTUtop25<-as.data.frame(t(otu_table(top25data)))
#Sampletop25<-as.data.frame(sample_data(top25data))

gentop10  = prune_taxa(phylumtop10,  data_phylum2)
melt1<-psmelt(gentop10)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_fungi_clay_top10 <- significance 


significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Genus)), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Genus))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~ConOrg,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_genus_fungi_sand_top10 <- significance 

phylumtop10  = prune_taxa(phylumtop10,  data_phylum1)
data_phylum2<-merge_samples_mean(phylumtop10,"Soil_Management") #takes the sum of the samples
data_phylum4<-psmelt(data_phylum2)

plot_bar(data_phylum2,fill="Genus")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)


ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungigenus.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Management",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungigenuswhlegend.pdf",width = 6.5,height=8)



#####composition flevoland vs zeeland#####
#bacteria phyla level
data_phylum0<-transform_sample_counts(PhyseqBacteria, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-subset_taxa(data_phylum0,Phylum_selection=="Keep")
data_phylum2<-tax_glom(data_phylum1,taxrank="Phylum",NArm=FALSE)

melt1<-psmelt(data_phylum2)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
modelsHSD <- list()
phylbac<-as.character(unique(melt1_clay$Phylum))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_bac_clay <- significance 

phylbac<-as.character(unique(melt1_sand$Phylum))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_bac_sand <- significance 
#data_phylum2<-subset_samples(data_phylum2,Soil=="Clay")
data_phylum3<-merge_samples_mean(data_phylum2,"Soil_Origin")

data_phylum4 <- psmelt(data_phylum3) # create dataframe from phyloseq object
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))
#data_phylum4 <- subset(data_phylum4, select=c(""))
data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character

# group dataframe by Phylum, calculate median rel. abundance
maximums <- ddply(data_phylum4, ~Phylum, function(x) c(max=max(x$Abundance)))
#maximums# find Phyla whose rel. abund. is less than 1%
remainder <- maximums[maximums$max <= 0.01,]$Phylum
# change their name to "Remainder"
data_phylum4[data_phylum4$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacphylum_soilorigin.pdf",width = 15, height = 8)
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampiling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"), 
        legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("compbacphylumwhlegend_soilorigin.pdf",width = 6.5, height = 8)



#bacteria genus level
data_phylum0<-transform_sample_counts(PhyseqBacteria, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-subset_taxa(data_phylum0,Genus_selection=="Keep")
data_phylum2<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])

melt1<-psmelt(data_phylum2)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand") 
phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")

#forallgenera
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_bac_clay_all <- significance 

for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~SOil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_bac_sand_all <- significance 


#fortop10
gentop10  = prune_taxa(phylumtop10,  data_phylum2)
melt1<-psmelt(gentop10)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_bac_clay_top10 <- significance 


significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Genus)), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Genus))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_bac_sand_top10 <- significance 


data_phylum3<-merge_samples_mean(data_phylum2,"Soil_Origin")
phylumtop10  = prune_taxa(phylumtop10,  data_phylum3)

data_phylum4 <- psmelt(phylumtop10) # create dataframe from phyloseq object
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacgenus_soilorigin.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compbacgenuswhlegend_soilorigin.pdf",width = 6.5, height = 8)

#Archaea phyla level
data_phylum0<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-tax_glom(data_phylum0,taxrank="Genus",NArm=FALSE)

melt1<-psmelt(data_phylum1)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_clay$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_clay$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_arch_clay_all <- significance 

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_arch_sand_all <- significance 

data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Origin")
data_phylum4 <- psmelt(data_phylum2) # create dataframe from phyloseq object
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))
data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character

# group dataframe by Phylum, calculate median rel. abundance
maximums <- ddply(data_phylum4, ~Phylum, function(x) c(max=max(x$Abundance)))
maximums# find Phyla whose rel. abund. is less than 1%
remainder <- maximums[maximums$max <= 0.01,]$Phylum
# change their name to "Remainder"
data_phylum4[data_phylum4$Phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("comparchphylum_soilorigin.pdf")
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("comparchphylumwhlegend_soilorigin.pdf",width = 6.5, height = 8)



#Archaea genus level - identified genera
PhyseqArchaearel<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU))
data_phylum1<-subset_taxa(PhyseqArchaearel,Genus_selection=="Keep")
data_phylum1<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)

melt1<-psmelt(data_phylum1)

melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")
modelsaov <- list()
modelsHSD <- list()
str(melt1_clay$Abundance)
phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")


for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_arch_clay_all <- significance 

phylbac<-as.character(unique(melt1_sand$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")


for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_arch_sand_all <- significance 

data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Origin") #takes the mean of the samples
data_phylum4<-psmelt(data_phylum2)
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchgenusidentified_soilorigin.pdf")
ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchgenusidentifiedwhlegend_soilorigin.pdf",width=6.5,height = 8)

#Archaea OTU level - identified species
PhyseqArchaearel<-transform_sample_counts(PhyseqArchaea, function(OTU) OTU/sum(OTU))

data_phylum2<-merge_samples_mean(PhyseqArchaearel,"Soil_Origin") #takes the mean of the samples
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])
phylumtop10  = prune_taxa(phylumtop10,  data_phylum2)

plot_bar(phylumtop10,fill="OTU")+
  geom_bar(aes(color=OTU, fill=OTU), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchOTUunidentified_soilorigin.pdf")
plot_bar(phylumtop10,fill="OTU")+
  geom_bar(aes(color=OTU, fill=OTU), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())
ggsave("comparchOTUunidentifiedwhlegend_soilorigin.pdf",width = 6.5, height = 8)


#fungi phyla level
data_phylum0<-transform_sample_counts(PhyseqFungi, function(OTU) OTU/sum(OTU)) # normalization
data_phylum1<-tax_glom(data_phylum0,taxrank="Phylum",NArm=FALSE)
data_phylum2<-merge_samples_mean(data_phylum1,"Soil_Origin")
data_phylum4 <- psmelt(data_phylum2) # create dataframe from phyloseq object
data_phylum4$Phylum <- as.character(data_phylum4$Phylum) #convert to character
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))

melt1<-psmelt(data_phylum1)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_clay$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_clay$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_fungi_clay <- significance 

significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Phylum)), ncol = 4))
colnames(significance) <- c("Phylum","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Phylum))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Phylum == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_phylum_fungi_sand <- significance 

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungiphylum_soilorigin.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Phylum),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per phylum\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungiphylumwhlegend_soilorigin.pdf",width = 6.5, height = 8)

#fungi genus level
data_phylum0<-transform_sample_counts(PhyseqFungi, function(X) X/sum(X))
data_phylum1<-subset_taxa(data_phylum0,Genus_selection=="keep")

data_phylum2<-tax_glom(data_phylum1,taxrank="Genus",NArm=FALSE)
phylumtop10 = names(sort(taxa_sums(data_phylum2), TRUE)[1:10])

#OTUtop25<-as.data.frame(t(otu_table(top25data)))
#Sampletop25<-as.data.frame(sample_data(top25data))

gentop10  = prune_taxa(phylumtop10,  data_phylum2)
melt1<-psmelt(gentop10)
melt1_clay <- subset(melt1,Soil =="Clay")
melt1_sand <- subset(melt1,Soil =="Sand")

phylbac<-as.character(unique(melt1_clay$Genus))
significance <- as.data.frame(matrix(0,nrow = length(phylbac), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_clay%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_fungi_clay_top10 <- significance 


significance <- as.data.frame(matrix(0,nrow = length(unique(melt1_sand$Genus)), ncol = 4))
colnames(significance) <- c("Genus","F","P","P.adj")
phylbac<-as.character(unique(melt1_sand$Genus))
for (i in 1:length(phylbac)){
  modelsaov [[i]] <- melt1_sand%>%filter(Genus == phylbac[i])%>%
    aov(Abundance~Soil_Origin,data=.)
  
  significance[i,1] <- phylbac[i]
  significance[i,2] <- summary(modelsaov[[i]])[[1]][["F value"]][1]
  significance[i,3] <- summary(modelsaov[[i]])[[1]][["Pr(>F)"]][1]
  significance[i,4] <- p.adjust(significance[i,3], method = "bonferroni", n = length(phylbac))
}
sign_soilorigin_genus_fungi_sand_top10 <- significance 

phylumtop10  = prune_taxa(phylumtop10,  data_phylum1)
data_phylum2<-merge_samples_mean(phylumtop10,"Soil_Origin") #takes the sum of the samples
data_phylum4<-psmelt(data_phylum2)
data_phylum4$Sample<-factor(data_phylum4$Sample,levels=c("Flevoland","Zeeland","North","Middle","South"))
plot_bar(data_phylum2,fill="Genus")+
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)


ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.title=element_text(size=15), 
        legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungigenus_soilorigin.pdf")

ggplot(data=data_phylum4)+
  geom_bar(aes(x=Sample, y=Abundance, fill=Genus),stat="identity") +
  labs(x="Sampling Region",y="relative abundance per genus\n (average of samples)") + 
  scale_fill_manual(values = phylum_colors)+
  scale_color_manual(values = phylum_colors)+
  theme_bw()+
  theme(axis.title=element_text(size=20,face="bold"),legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),strip.text=element_blank())

ggsave("compfungigenuswhlegend_soilorigin.pdf",width = 6.5,height=8)

