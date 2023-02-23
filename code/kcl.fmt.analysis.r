require(momr)
require(matrixStats)
require(gplots)
require(vegan)
require(ape)
require(ggplot2)
require(ggsci)

gutMsp1992Data = "C:\\Data/comparative.analysis.healthy.sweden/hs_10.4_1992_MSP_freeze2_20180905.RData"
load(gutMsp1992Data)

mgsFile = "merged.final.mgs.med.vec.10M.RData"
load(mgsFile)
mgs_med_vec_10m = mgs_med_vec_10m[,colnames(mgs_med_vec_10m) != "FMT102.1"]
mgs_richness = colSums( mgs_med_vec_10m > 0)

corMat = cor(mgs_med_vec_10m)
heatmap.2(corMat, trace="none")

sampleFile = "P15952_20_03_sample_info_stool.txt"
sampleTab = read.delim(sampleFile, sep="\t")
sampleTab$fmtId = gsub("P15952_", "FMT", sampleTab$NGI.ID)
sampleTab$User.ID = gsub("Donor ","Donor", sampleTab$User.ID)
sampleTab$User.ID = gsub("Stool","", sampleTab$User.ID)

sampleTab$User.ID = gsub("Day 7","Day7", sampleTab$User.ID)
sampleTab$User.ID = gsub("Day 30","Day30", sampleTab$User.ID)
sampleTab$User.ID = gsub("Day 90","Day90", sampleTab$User.ID)

sampleTab$User.ID = gsub(" \\+ "," ", sampleTab$User.ID)
sampleTab$subjects = sapply(sampleTab$User.ID, function(x) strsplit(x, split = " ")[[1]][1])
sampleTab$condition = sapply(sampleTab$User.ID, function(x) strsplit(x, split = " ")[[1]][2])
sampleTab$richness = mgs_richness[match(sampleTab$fmtId, names(mgs_richness))]
sampleTab$category = rep("recipient",90)
sampleTab$category[grepl("Donor",sampleTab$subjects)] = "donor"
sampleTab_donor = sampleTab[sampleTab$category == "donor",]
richness_donor_types = split(sampleTab_donor$richness, sampleTab_donor$condition)
boxplot(richness_donor_types)
wilcox.test(richness_donor_types$glycerol, richness_donor_types$raw)
fullSubjects = names(table(sampleTab$subjects)[table(sampleTab$subjects) == 4])
sampleTab_full = sampleTab[sampleTab$subjects %in% fullSubjects, ]
sampleTab_full = sampleTab_full[order(sampleTab_full$condition),]
sampleTab_full = sampleTab_full[order(sampleTab_full$subjects),]

richnessByTimes = split(sampleTab_full$richness, sampleTab_full$condition)
boxplot(richnessByTimes[c("Baseline", "Day7","Day30","Day90")])

wilcox.test(richnessByTimes$Baseline, richnessByTimes$Day7, paired = T)
wilcox.test(richnessByTimes$Baseline, richnessByTimes$Day30, paired = T )
wilcox.test(richnessByTimes$Baseline, richnessByTimes$Day90, paired = T )


richness_donor_recipient = split(sampleTab$richness, sampleTab$category)
boxplot(richness_donor_recipient)



ind_baseline = grepl("Baseline", sampleTab$User.ID)

tail(sampleTab)


### richness analysis ###

### pcoa analysis ###

### all dataset ###
mgsDist=vegdist(t(mgs_med_vec_10m), method="bray")
pcoaOut=pcoa(mgsDist)
pcoaMat = pcoaOut$vectors 
pcoaTab = data.frame(pcoaMat[,1:2],
                     group = sampleTab$condition[match(rownames(pcoaMat), sampleTab$fmtId)],
                     stringsAsFactors = F)
ggplot(pcoaTab, aes(x=Axis.1, y=Axis.2, colour=group)) + geom_point(size=3)+ theme_bw() + scale_colour_jco()

### subjects with full visits ###
mgs_med_vec_10m_full = mgs_med_vec_10m[,colnames(mgs_med_vec_10m) %in% sampleTab_full$fmtId]
mgsDistFull=vegdist(t(mgs_med_vec_10m_full), method="bray")
pcoaOutFull=pcoa(mgsDistFull)
pcoaMatFull = pcoaOutFull$vectors 
pcoaTabFull = data.frame(pcoaMatFull[,1:2],
                         subj = sampleTab_full$subjects[match(rownames(pcoaMatFull), sampleTab_full$fmtId)],
                         group = sampleTab_full$condition[match(rownames(pcoaMatFull), sampleTab_full$fmtId)],
                         stringsAsFactors = F)
pcoaTabFull$time = 0
pcoaTabFull$time[pcoaTabFull$group == "Day7"] = 1
pcoaTabFull$time[pcoaTabFull$group == "Day30"] = 2
pcoaTabFull$time[pcoaTabFull$group == "Day90"] = 3
pcoaTabFull$time = factor(pcoaTabFull$time)

#pcoaTabFull$group = factor(pcoaTabFull$group, levels = c("Baseline", "Day7", "Day30", "Day90"))
ggplot(pcoaTabFull, aes(x=Axis.1, y=Axis.2, colour=time))+ geom_path(aes(group=subj)) + geom_point(size=3)+ theme_bw() + scale_colour_npg()


testRelations()

### longitudinal analysis ###
samplesByTimes = split(sampleTab_full$fmtId, sampleTab_full$condition)
pCut= 0.05

targetMat = cbind(mgs_med_vec_10m[,samplesByTimes$Baseline], mgs_med_vec_10m[,samplesByTimes$Day7])
traits = rep(c("baseline", "day7"), each=16)
stats.baseline.day7 = testRelations(targetMat,traits, type = "wilcoxon",paired = T )
stats.baseline.day7$species = taxo$species[match(rownames(stats.baseline.day7), rownames(taxo))]
stats.baseline.day7.sig = stats.baseline.day7[!is.na(stats.baseline.day7$p) & stats.baseline.day7$p < pCut, ]

targetMat = cbind(mgs_med_vec_10m[,samplesByTimes$Baseline], mgs_med_vec_10m[,samplesByTimes$Day30])
traits = rep(c("baseline", "day30"), each=16)
stats.baseline.day30 = testRelations(targetMat,traits, type = "wilcoxon",paired = T )
stats.baseline.day30$species = taxo$species[match(rownames(stats.baseline.day30), rownames(taxo))]
stats.baseline.day30.sig = stats.baseline.day30[!is.na(stats.baseline.day30$p) & stats.baseline.day30$p < pCut, ]

targetMat = cbind(mgs_med_vec_10m[,samplesByTimes$Baseline], mgs_med_vec_10m[,samplesByTimes$Day90])
traits = rep(c("baseline", "day90"), each=16)
stats.baseline.day90 = testRelations(targetMat,traits, type = "wilcoxon",paired = T )
stats.baseline.day90$species = taxo$species[match(rownames(stats.baseline.day90), rownames(taxo))]
stats.baseline.day90.sig = stats.baseline.day90[!is.na(stats.baseline.day90$p) & stats.baseline.day90$p < pCut, ]
