library("ALDEx2")

metadata <- read.table("data/clean_metadata_gg.txt",header=TRUE,sep="\t",quote="")

# healthy is represented by 0, SS/NASH is represented by 1
binaryGroups <- metadata$SSvsNASH
binaryGroups[which(!is.na(metadata$SSvsNASH))] <- 1
binaryGroups[which(is.na(metadata$SSvsNASH))] <- 0

# NASH is 1, healthy is 0, SS is NA
nashHealthyGroups <- metadata$SSvsNASH
nashHealthyGroups[which(metadata$SSvsNASH==0)] <- NA
nashHealthyGroups[which(is.na(metadata$SSvsNASH))] <- 0

# SS is 1, healthy is 0, NASH is NA
ssHealthyGroups <- metadata$SSvsNASH
ssHealthyGroups[which(metadata$SSvsNASH==1)] <- NA
ssHealthyGroups[which(metadata$SSvsNASH==0)] <- 1
ssHealthyGroups[which(is.na(metadata$SSvsNASH))] <- 0

metagenomicGroups <- metadata$SSvsNASH
metagenomicGroups[] <- NA
metagenomicNASH <- c("CL-166-BL", "CL-169-BL", "CL-139-BL-2", "CL-173-2", "CL-177", "CL-160", "CL-165", "CL-119")
metagenomicGroups[which(metadata$X%in%metagenomicNASH)] <- 1
metagenomicHealthy <- c("HLD-100", "HLD-102", "HLD-111-2", "HLD-80", "HLD-85", "HLD-28", "HLD-112", "HLD-23")
metagenomicGroups[which(metadata$X%in%metagenomicHealthy)] <- 0

data <- read.table("data/td_OTU_tag_mapped_lineage.txt", header=TRUE,sep="\t",quote="",row.names=1)

taxonomy <- data.taxonomy[,ncol(data)]
taxonomy <- as.character(taxonomy)

data <- data[,c(1:(ncol(data)-1))]

taxonomy <- strsplit(taxonomy,";")

genus <- lapply(taxonomy,function(x) { return(x[6]) })

genus <- unlist(genus)

grouplist <- list()
grouplist[[1]] <- genus

data.genus <- aggregate(data, by=grouplist, sum,simplify=TRUE)

rownames(data.genus) <- data.genus$Group.1

data.genus <- data.genus[,c(2:ncol(data.genus))]

healthy.ssnash <- aldex(data.genus, binaryGroups, mc.samples=128, test="t", effect=TRUE,
    include.sample.summary=FALSE, verbose=FALSE)

healthy.nash <- aldex(data.genus, nashHealthyGroups, mc.samples=128, test="t", effect=TRUE,
    include.sample.summary=FALSE, verbose=FALSE)

healthy.ss <- aldex(data.genus, ssHealthyGroups, mc.samples=128, test="t", effect=TRUE,
    include.sample.summary=FALSE, verbose=FALSE)

metagenomic.nash.data.genus <- aldex(data, metagenomicGroups, mc.samples=128, test="t", effect=TRUE,
    include.sample.summary=FALSE, verbose=FALSE)

# healthy.ssnash <- aldex(data, binaryGroups, mc.samples=128, test="t", effect=TRUE,
#     include.sample.summary=FALSE, verbose=FALSE)
# 
# healthy.nash <- aldex(data, nashHealthyGroups, mc.samples=128, test="t", effect=TRUE,
#     include.sample.summary=FALSE, verbose=FALSE)
# 
# healthy.ss <- aldex(data, ssHealthyGroups, mc.samples=128, test="t", effect=TRUE,
#     include.sample.summary=FALSE, verbose=FALSE)
# 
# metagenomic.nash.healthy <- aldex(data, metagenomicGroups, mc.samples=128, test="t", effect=TRUE,
#     include.sample.summary=FALSE, verbose=FALSE)

#sort
healthy.ssnash <- healthy.ssnash[order(-abs(healthy.ssnash$effect)),]
healthy.nash <- healthy.nash[order(-abs(healthy.nash$effect)),]
healthy.ss <- healthy.ss[order(-abs(healthy.ss$effect)),]
metagenomic.nash.healthy <- metagenomic.nash.healthy[order(-abs(metagenomic.nash.healthy$effect)),]

healthy.ssnash.20 <- rownames(healthy.ssnash)[1:20]
healthy.nash.20 <- rownames(healthy.nash)[1:20]
healthy.ss.20 <- rownames(healthy.ss)[1:20]
metagenomic.nash.healthy.20 <- metagenomic.nash.healthy[1:20,]

ssnash.nash <- intersect(healthy.ssnash.20, healthy.nash.20)
ssnash.ss <- intersect(healthy.ssnash.20, healthy.ss.20)
ss.nash <- intersect(healthy.nash.20, healthy.ss.20)
ss.nash.ssnash <- intersect(ssnash.nash, healthy.ss.20)

nash.metagenomic <- intersect(healthy.nash.20,rownames(metagenomic.nash.healthy.20))

x_axis <- c(1:(length(ss.nash.ssnash)*3))

effect <- c(1:(length(ss.nash.ssnash)*3))
effect[c(TRUE, FALSE, FALSE)] <- healthy.ssnash$effect[which(rownames(healthy.ssnash)%in%ss.nash.ssnash)]
effect[c(FALSE, TRUE, FALSE)] <- healthy.nash$effect[which(rownames(healthy.nash)%in%ss.nash.ssnash)]
effect[c(FALSE, FALSE, TRUE)] <- healthy.ss$effect[which(rownames(healthy.ss)%in%ss.nash.ssnash)]

effect.975 <- c(1:(length(ss.nash.ssnash)*3))
effect.975[c(TRUE, FALSE, FALSE)] <- healthy.ssnash$effect.975[which(rownames(healthy.ssnash)%in%ss.nash.ssnash)]
effect.975[c(FALSE, TRUE, FALSE)] <- healthy.nash$effect.975[which(rownames(healthy.nash)%in%ss.nash.ssnash)]
effect.975[c(FALSE, FALSE, TRUE)] <- healthy.ss$effect.975[which(rownames(healthy.ss)%in%ss.nash.ssnash)]

effect.025 <- c(1:(length(ss.nash.ssnash)*3))
effect.025[c(TRUE, FALSE, FALSE)] <- healthy.ssnash$effect.025[which(rownames(healthy.ssnash)%in%ss.nash.ssnash)]
effect.025[c(FALSE, TRUE, FALSE)] <- healthy.nash$effect.025[which(rownames(healthy.nash)%in%ss.nash.ssnash)]
effect.025[c(FALSE, FALSE, TRUE)] <- healthy.ss$effect.025[which(rownames(healthy.ss)%in%ss.nash.ssnash)]

# get taxa names
taxa_table <- read.table("data/td_OTU_tag_mapped_lineage.txt",header=TRUE,sep="\t",quote="",row.names=1)

taxa_labels <- c(1:(length(ss.nash.ssnash)))
taxa_labels <- as.character(taxa_table$taxonomy[which(rownames(taxa_table)%in%ss.nash.ssnash)])
taxa_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;","",taxa_labels)

pdf("common_otus_in_hc_vs_nash_ss_effect_size_with_confidence_interval.pdf")
plot(NA, xlim=c(0,length(ss.nash.ssnash)*3+1),ylim=c(-10,10),main="Confidence Intervals for OTUs with High Effect Size",xlab="OTU",ylab="Effect size",xaxt='n')
points(x_axis, effect, col=c(1,2,3),ann=FALSE,pch=19)
arrows(x_axis, effect.025, x_axis, effect.975, length=0.05, angle=90, code=3,col=c(1,2,3))
text(x =x_axis[c(TRUE,FALSE,FALSE)], y = par("usr")[3] - 0.5, labels = taxa_labels,srt = 45, pos = 1, xpd = TRUE, cex=0.5)
legend(x=1, y=9, legend=c("Healthy vs. SS/NASH", "Healthy vs. NASH", "Healthy vs. SS"),col=c(1,2,3),pch=19)
dev.off()

#plot metagenomic effect sizes
x_axis <- c(1:20)
effect <- metagenomic.nash.healthy.20$effect
effect.975 <- metagenomic.nash.healthy.20$effect.975
effect.025 <- metagenomic.nash.healthy.20$effect.025

taxa_labels <- c(1:20)
taxa_labels <- as.character(taxa_table$taxonomy[which(rownames(taxa_table)%in%rownames(metagenomic.nash.healthy.20))])
taxa_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;","",taxa_labels)


pdf("metagenomic_samples_otu_effect_size_with_confidence_intervals.pdf")
plot(NA, xlim=c(0,21),ylim=c(-15,15),main="Confidence Intervals for OTUs with High Effect Size",xlab="OTU",ylab="Effect size",xaxt='n')
points(x_axis, effect,ann=FALSE,pch=19)
arrows(x_axis, effect.025, x_axis, effect.975, length=0.05, angle=90, code=3)
text(x =x_axis, y = par("usr")[3] - 1.75, labels = taxa_labels,srt = 45, pos = 1, xpd = TRUE, cex=0.5)
dev.off()

#plot effect sizes for top OTUs in metagenomic and healthy vs. nash
x_axis <- c(1:(length(nash.metagenomic)*2))

effect <- c(1:(length(nash.metagenomic)*2))
effect[c(TRUE, FALSE)] <- healthy.nash$effect[which(rownames(healthy.nash)%in%nash.metagenomic)]
effect[c(FALSE, TRUE)] <- metagenomic.nash.healthy.20$effect[which(rownames(metagenomic.nash.healthy.20)%in%nash.metagenomic)]

effect.975 <- c(1:(length(nash.metagenomic)*2))
effect.975[c(TRUE, FALSE)] <- healthy.nash$effect.975[which(rownames(healthy.nash)%in%nash.metagenomic)]
effect.975[c(FALSE, TRUE)] <- metagenomic.nash.healthy.20$effect.975[which(rownames(metagenomic.nash.healthy.20)%in%nash.metagenomic)]

effect.025 <- c(1:(length(nash.metagenomic)*2))
effect.025[c(TRUE, FALSE)] <- healthy.nash$effect.025[which(rownames(healthy.nash)%in%nash.metagenomic)]
effect.025[c(FALSE, TRUE)] <- metagenomic.nash.healthy.20$effect.025[which(rownames(metagenomic.nash.healthy.20)%in%nash.metagenomic)]

taxa_labels <- as.character(taxa_table$taxonomy[which(rownames(taxa_table)%in%nash.metagenomic)])
taxa_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;","",taxa_labels)

pdf("metagenomic_and_healthy_vs_nash_samples_otu_effect_size_with_confidence_intervals.pdf")
plot(NA, xlim=c(0,((length(nash.metagenomic)*2)+1)),ylim=c(-15,15),main="Confidence Intervals for OTUs with High Effect Size",xlab="OTU",ylab="Effect size",xaxt='n')
points(x_axis, effect,ann=FALSE,pch=19,col=c(1,2))
arrows(x_axis, effect.025, x_axis, effect.975, length=0.05, angle=90, code=3,col=c(1,2))
text(x =x_axis[c(TRUE,FALSE)], y = par("usr")[3] - 1.25, labels = taxa_labels,srt = 45, pos = 1, xpd = TRUE, cex=0.5)
legend(x=1, y=14, legend=c("Healthy vs. NASH", "Metagenomic study samples only"),col=c(1,2),pch=19)
dev.off()

#plot effect sizes for top OTUs in metagenomic
x_axis <- c(1:20)

effect <- c(1:20)
effect <- metagenomic.nash.healthy.20$effect

effect.975 <- c(1:20)
effect.975 <- metagenomic.nash.healthy.20$effect.975

effect.025 <- c(1:20)
effect.025 <- metagenomic.nash.healthy.20$effect.025

taxa_labels <- as.character(taxa_table$taxonomy[which(rownames(taxa_table)%in%rownames(metagenomic.nash.healthy.20))])
taxa_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;[^;]*;","",taxa_labels)

pdf("metagenomic_samples_top_20_otu_effect_size_with_confidence_intervals.pdf")
plot(NA, xlim=c(0,21),ylim=c(-20,20),main="Confidence Intervals for OTUs with High Effect Size",xlab="OTU",ylab="Effect size",xaxt='n')
points(x_axis, effect,ann=FALSE,pch=19)
arrows(x_axis, effect.025, x_axis, effect.975, length=0.05, angle=90, code=3)
text(x =x_axis, y = par("usr")[3] - 2.25, labels = taxa_labels,srt = 45, pos = 1, xpd = TRUE, cex=0.5)
dev.off()


write.table(metagenomic.nash.healthy,file="metagenomic_samples_aldex_results.txt",sep="\t",quote=FALSE)

pdf("metagenomic_and_healthy_vs_nash_effect_size_correlation.pdf")
metagenomic.effect <- metagenomic.nash.healthy$effect[order(rownames(metagenomic.nash.healthy))]
nash.effect <- healthy.nash$effect[order(rownames(healthy.nash))]
plot(metagenomic.effect, nash.effect, pch=19, col=rgb(1,0,1,0.3), main="Metagenomic samples' effect size against all healthy vs. nash effect size", xlab="metagenomic samples effect size", ylab="all samples effect size")
aldex.plot(metagenomic.nash.healthy, type="MA", test="welch")
aldex.plot(metagenomic.nash.healthy, type="MW", test="welch")

dev.off()

magenta <- c("2", "19", "21", "31", "40", "54", "99")
cyan <- c("42","48","200","1367")
green <- c("27","76","82","89","92","110","114","140","154","178","194","226","362","1095","1243")
magenta_labels <- taxa_table$taxonomy[which(rownames(taxa_table)%in%magenta)]
magenta_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;","",magenta_labels)
cyan_labels <- taxa_table$taxonomy[which(rownames(taxa_table)%in%cyan)]
cyan_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;","",cyan_labels)
green_labels <- taxa_table$taxonomy[which(rownames(taxa_table)%in%green)]
green_labels <- gsub("^[^;]*;[^;]*;[^;]*;[^;]*;","",green_labels)

pdf("Healthy_vs_NASH_metagenomic_samples.pdf")
aldex.plot(metagenomic.nash.healthy, type="MW", test="welch")
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%magenta)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%magenta)],col="magenta",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%cyan)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%cyan)],col="cyan",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%green)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%green)],col="green",pch=20)
aldex.plot(metagenomic.nash.healthy, type="MW", test="welch")
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%magenta)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%magenta)],col="magenta",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%cyan)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%cyan)],col="cyan",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%green)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%green)],col="green",pch=20)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%magenta)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%magenta)],labels=magenta_labels[which(magenta%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%magenta)]))],pos=4,cex=0.3)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%cyan)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%cyan)],labels=cyan_labels[which(cyan%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%cyan)]))],pos=4,cex=0.3)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%green)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%green)],labels=green_labels[which(green%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%green)]))],pos=4,cex=0.3)
aldex.plot(metagenomic.nash.healthy, type="MW", test="welch")
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%magenta)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%magenta)],col="magenta",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%cyan)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%cyan)],col="cyan",pch=20)
points(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%green)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%green)],col="green",pch=20)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%magenta)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%magenta)],labels=magenta_labels[which(magenta%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%magenta)]))],pos=4,cex=0.5)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%cyan)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%cyan)],labels=cyan_labels[which(cyan%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%cyan)]))],pos=4,cex=0.5)
text(metagenomic.nash.healthy$diff.win[which(rownames(metagenomic.nash.healthy)%in%green)],metagenomic.nash.healthy$diff.btw[which(rownames(metagenomic.nash.healthy)%in%green)],labels=green_labels[which(green%in%(rownames(metagenomic.nash.healthy)[which(rownames(metagenomic.nash.healthy)%in%green)]))],pos=4,cex=0.5)
dev.off()

pdf("Healthy_vs_NASH_all_samples.pdf")
aldex.plot(healthy.nash, type="MW", test="welch")
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%magenta)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%magenta)],col="magenta",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%cyan)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%cyan)],col="cyan",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%green)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%green)],col="green",pch=20)
aldex.plot(healthy.nash, type="MW", test="welch")
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%magenta)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%magenta)],col="magenta",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%cyan)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%cyan)],col="cyan",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%green)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%green)],col="green",pch=20)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%magenta)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%magenta)],labels=magenta_labels[which(magenta%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%magenta)]))],pos=4,cex=0.3)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%cyan)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%cyan)],labels=cyan_labels[which(cyan%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%cyan)]))],pos=4,cex=0.3)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%green)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%green)],labels=green_labels[which(green%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%green)]))],pos=4,cex=0.3)
aldex.plot(healthy.nash, type="MW", test="welch")
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%magenta)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%magenta)],col="magenta",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%cyan)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%cyan)],col="cyan",pch=20)
points(healthy.nash$diff.win[which(rownames(healthy.nash)%in%green)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%green)],col="green",pch=20)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%magenta)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%magenta)],labels=magenta_labels[which(magenta%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%magenta)]))],pos=4,cex=0.5)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%cyan)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%cyan)],labels=cyan_labels[which(cyan%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%cyan)]))],pos=4,cex=0.5)
text(healthy.nash$diff.win[which(rownames(healthy.nash)%in%green)],healthy.nash$diff.btw[which(rownames(healthy.nash)%in%green)],labels=green_labels[which(green%in%(rownames(healthy.nash)[which(rownames(healthy.nash)%in%green)]))],pos=4,cex=0.5)
dev.off()

#
# ssnash.nash
#  [1] "30"   "40"   "54"   "2"    "99"   "1367" "24"   "15"   "41"   "140"
# [11] "67"   "89"   "31"   "362"
#
#
# ssnash.ss
#  [1] "40"   "54"   "2"    "1367" "21"   "15"   "41"   "1334" "19"   "89"
# [11] "31"
#
#
# ss.nash
# [1] "54"   "40"   "1367" "89"   "41"   "2"    "15"   "31"
#
# ss.nash.ssnash
# [1] "40"   "54"   "2"    "1367" "15"   "41"   "89"   "31"
