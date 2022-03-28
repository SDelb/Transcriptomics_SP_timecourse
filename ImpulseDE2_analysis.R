library(ImpulseDE2)

d = read.delim("data/design_updated.txt", header = T, as.is = T)
dat = read.delim("data/raw_corrected_counts_filtered_updated.txt", header = T, as.is = T)

# Comparisons of interest:
CvsV = dat[,c(11:30, 41:60, 71:87)]
SvsV = dat[,c(1:20, 31:50, 61:79)]
SvsC = dat[,c(1:10, 21:40, 51:70, 80:87)]

#######################################
#                                     #
#         C vs V                      #
#                                     #
#######################################

# annotation info:
Sample = rownames(d)[c(11:30, 41:60, 71:87)]
Condition = c(rep(rep(c("control", "case"), each = 10), 2), rep("control", 9), rep("case", 8))
Time = d[c(11:30, 41:60, 71:87), 3]
design = cbind(Sample, Condition, Time)
colnames(design) = c("Sample", "Condition", "Time")
write.table(design, file="d.txt", row.names = F, quote=F, sep="\t")
design = read.delim("d.txt", header=T, as.is=T)

ImpulseDE2_Results = runImpulseDE2(matCountData = as.matrix(CvsV),
                                   dfAnnotation = design,
                                   boolCaseCtrl = T,
                                   scaQThres = 0.05,
                                   boolIdentifyTransients = T,
                                   scaNProc = 6)
save(ImpulseDE2_Results, file= "CvsV.RData")

#######################################
#                                     #
#         S vs V                      #
#                                     #
#######################################

# annotation info:
Sample = rownames(d)[c(1:20, 31:50, 61:79)]
Condition = c(rep(rep(c("case", "control"), each = 10), 2), rep("case", 10), rep("control", 9))
Time = d[c(1:20, 31:50, 61:79), 3]
design = cbind(Sample, Condition, Time)
colnames(design) = c("Sample", "Condition", "Time")
write.table(design, file="d.txt", row.names = F, quote=F, sep="\t")
design = read.delim("d.txt", header=T, as.is=T)

ImpulseDE2_Results = runImpulseDE2(matCountData = as.matrix(SvsV),
                                   dfAnnotation = design,
                                   boolCaseCtrl = T,
                                   scaQThres = 0.05,
                                   boolIdentifyTransients = T,
                                   scaNProc = 6)
save(ImpulseDE2_Results, file= "SvsV.RData")

#######################################
#                                     #
#         S vs C                      #
#                                     #
#######################################

# annotation info:
Sample = rownames(d)[c(1:10, 21:40, 51:70, 80:87)]
Condition = c(rep(rep(c("case", "control"), each = 10), 2), rep("case", 10), rep("control", 8))
Time = d[c(1:10, 21:40, 51:70, 80:87), 3]
design = cbind(Sample, Condition, Time)
colnames(design) = c("Sample", "Condition", "Time")
write.table(design, file="d.txt", row.names = F, quote=F, sep="\t")
design = read.delim("d.txt", header=T, as.is=T)

ImpulseDE2_Results = runImpulseDE2(matCountData = as.matrix(SvsC),
                                   dfAnnotation = design,
                                   boolCaseCtrl = T,
                                   scaQThres = 0.05,
                                   boolIdentifyTransients = T,
                                   scaNProc = 6)
save(ImpulseDE2_Results, file= "SvsC.RData")


