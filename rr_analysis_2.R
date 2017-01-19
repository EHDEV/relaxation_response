
# Loading GE and Disease Data

setwd("~/columbia/fall_2015/binf/")
rm(list=ls())
# GDS3416.soft - Relaxation response practice effect on blood

x = readLines("./data/r_data/GDS3416.soft")
read.GEO <- function(filename){
  # x <- readLines(filename) 
  start.line<-grep("!dataset_table_begin",x)+1
  end.line<-grep("!dataset_table_end",x)-1
  return (read.table(textConnection(x[start.line:end.line]),header = T, stringsAsFactors=F,sep="\t"))
}
data = read.GEO("./data/GDS3416.soft")


diabetes = read.csv("./data/diabetes.csv", header=T, stringsAsFactors = F)
disease2 = read.csv("./data/obesity.tsv", header=T, stringsAsFactors = F, sep="\t")
disease3 = read.csv("./data/alzhymers.tsv", header=T, stringsAsFactors = F, sep="\t")
disease4 = read.csv("./data/t2diabetes.tsv", header=T, stringsAsFactors = F, sep="\t")

d = readline("./data/diabetes.tsv")

# 
lb = which(data$IDENTIFIER %in% diseases$Symbol)
length(lb)
rownames(data) = data$ID_REF
data$ID_REF = NULL
data$IDENTIFIER = NULL

data.t <- na.omit(data.frame(t(data[lb,])))
colnames(data.t) = rownames(data[lb, ])

length(which(is.na(data.t)))
# premdata = data[c(neg.sample, pos.sample), ]

##### Main Data Filtered  -- some munging ######

for (i in 1:dim(data.t)[2])
  if (class(data.t[,i]) != "numeric")
    data.t[,i] = as.numeric(as.character(data.t[,i]))

### Extracting Annotation of Samples

first.line<-grep('!subset_description', x)+1  
number_of_comparators_data = length(first.line)
# subset the sample ids that are control/disease case
group_ids <- x[first.line]
group_ids <- gsub('!subset_sample_id = ','', group_ids)
group_ids <-strsplit(group_ids, ",")
# Keep only the disease cases and discard the control and other cases
ctl_ids = unlist(group_ids[[1]]) # control
st_ids = unlist(group_ids[[2]]) # Short ter RR training
lt_ids = unlist(group_ids[[3]]) # Long term RR training

ctl.idx = which(rownames(data.t) %in%ctl_ids) # control indices
st.idx = which(rownames(data.t) %in%st_ids)
lt.idx = which(rownames(data.t) %in%lt_ids)
################

# Combined class labels for each of the 72 samples
sample_class = c(rep(0, length(ctl.idx)), rep(1, length(st.idx)), rep(2, length(lt.idx)))

# Normalizing the data and performing t-test (N0 vs. N2, N0 vs. N1, N1 vs. N2)
# norm_data = scale(data.t)

ct.data = data.t[which(sample_class == 0), ]
st.data = data.t[which(sample_class == 1), ]
lt.data = data.t[which(sample_class == 2), ]

tt.st.lt = Map(t.test, x=lt.data, y=st.data)

### Start ### Control vs. Short-term 

tt.ct.st = Map(t.test, x=st.data, y=ct.data)
p.ct.st = rep(NA, length(tt.ct.st))
for(i in 1:length(p.ct.st)){
  p.ct.st[i] = tt.ct.st[[colnames(ct.data)[i]]][[3]]
}
length(which(p.ct.st <= 0.01))
##### End #####

### Start ### Control vs. Long-term 
tt.ct.lt = Map(t.test, x=ct.data, y=lt.data)
p.ct.lt = rep(NA, length(tt.ct.lt))
for(i in 1:length(p.ct.lt)){
  p.ct.lt[i] = tt.ct.lt[[colnames(ct.data)[i]]][[3]]
}
length(which(p.ct.lt <= 0.03))

diff_exp_genes_lt_idx = which(p.ct.lt <= 0.03)
diff_exp_genes_names_lt = colnames(data.t)[diff_exp_genes_lt_idx]

write(diff_exp_genes_names_lt, "signif_exp_nominal_genes.txt")

##### End #####

### Start ### Short-term vs. Long-term 

tt.st.lt = Map(t.test, x=lt.data, y=st.data)
p.st.lt = rep(NA, length(tt.st.lt))
for(i in 1:length(p.st.lt)){
  p.st.lt[i] = tt.st.lt[[colnames(st.data)[i]]][[3]]
}
length(which(p.st.lt <= 0.05))
##### End #####

tt.genes.ltst = colnames(lt.data)[which(p.st.lt <= 0.01)]
tt.genes.ltct = colnames(lt.data)[which(p.ct.lt <= 0.01)]
tt.genes.stct = colnames(lt.data)[which(p.ct.st <= 0.01)]


################# qFDR #############

library(qvalue)

pfdrtest.st.lt <- qvalue(p.st.lt) 
pfdrtest.ct.st <- qvalue(p.ct.st)
pfdrtest.ct.lt <- qvalue(p.ct.lt)

summary(pfdrtest.ct.lt)
hist(pfdrtest.ct.lt)
hist(pfdrtest.ltct)




# ct.data.log = log(data.t[ctl.idx, ])
# st.data.log = log(data.t[st.idx, ])
# lt.data.log = log(data.t[lt.idx, ])

# ttres = Map(t.test, x = st.data.log, y = ct.data.log)

# ttest = function(x) {
#   tt = t.test(as.matrix(x[1:23, ]), as.matrix(x[24:47, ]))
#   return (tt$p.value)
# }
# mm.pvals<-apply(mm,2, 
#                 function(x) {out <- t.test(x[1:15],x[16:25]);out$p.value}) 

# fun <- function(d){return(t.test(d[st.idx],d[ctl.idx])$p.value)}

# pval.ct.st = NULL
# for (c in 1:dim(norm_data)[2])
#   pval.ct.st = c(pval.ct.st, ttest(norm_data[, c]))
#            
# 
# pval.ct.lt = NULL
# for (c in 1:dim(norm_data)[2])
#   pval.ct.lt = c(pval.ct.lt, ttest(norm_data[, c]))
#            
# 
# pval.st.lt = NULL
# for (c in 1:dim(norm_data)[2])
#   pval.st.lt = c(pval.st.lt, ttest(norm_data[, c]))
#            

pval2 = lapply(norm_data, ttest)           
           
ct.data = data.t[ctl.idx, ]
st.data = data.t[st.idx, ]
lt.data = data.t[lt.idx, ]


hist(as.matrix(ct.data.log[1,]))

ttall= Map(t.test, x = lt.data.log, y = ct.data.log)

ct.data$group = rep(1, length(ctl.idx))
lt.data$group = rep(2, length(lt.idx))
ctl_m = merge.data.frame(ct.data, lt.data, by=0, all=T)

k = unlist(ttall)
p = rep(NA, length(ttall))
for(i in 1:length(ttall)){
  p[i] = ttall[[colnames(lt.data.log)[i]]][[3]]
}
?Map
length(which(p <= 0.01))

hist(p)
hist(pfdrtest)



### FDR ###

psort <- sort(p.ct.lt)
fdrtest <- NULL
for (i in 1:dim(data.t)[2])
  fdrtest <- c(fdrtest, p[i] > match(p[i],psort) * .05/1000)

length(sig)
sig = which(fdrtest)

pp <- p.adjust(p, method="BH")
table(pp)
library(qvalue)

pfdrtest <- qvalue(p.ct.lt) 
summary(pfdrtest)
hist(p.ltct)
hist(pfdrtest)


gene_list = colnames(mdata)
write(gene_list, "./t2d_gene_list.txt")

########################################-c(1295, 3731)
library(penalized)
data.t$group = sample_class
data.t.log = log(data.t[,which(names(data.t) != 'group')])
data.t.log$group = data.t$group

# Lasso - testing cv options

fit <- penalized(data.t.log$group[1:60], data.t.log[1:60,-c(which(names(data.t.log) %in% c("group", "207145_at")))], data=data.t, lambda1=1.5)
fit_cvl = cvl(data.t.log$group, data.t.log[,-c(which(names(data.t.log) %in% c("group", "207145_at")))], lambda1=1.5, fold=10)
fit_opt_lt <- optL1(data.t$group[c(ctl.idx, lt.idx)], penalized = data.t[c(ctl.idx, lt.idx),-c(which(names(data.t) %in% c("group", "207145_at")))], minlambda1 = 0.5, maxlambda1 = 25,fold = 10)
fit_cvl_lt <- cvl(data.t.log$group[c(ctl.idx, lt.idx)], penalized = data.t.log[c(ctl.idx, lt.idx),-c(which(names(data.t.log) %in% c("group", "207145_at")))], fold = 10)


fit_cvl_lt <- cvl(data.t.log$group[c(ctl.idx, lt.idx)], penalized = data.t.log[c(ctl.idx, lt.idx),-c(which(names(data.t.log) %in% c("group", "207145_at")))], fold = 10)


length(which(abs(coef(fit)) > .01))


##### NMF - Experimenting

res <- nmf(t(data.t[,-c(which(names(data.t.log) %in% c("group", "207145_at")))]), 3)
V.hat <- fitted(res)
dim(V.hat)

summary(res, class = data.t$group)
W = basis(res)
H = coef(res)
s = featureScore(res)
summary(s)
s = extractFeatures(res)
str(s)

# When the seeding method is stochastic, multiple runs are usually required to achieve stability or a
# resonable result. This can be done by setting argument nrun to the desired value. For performance
# reason we use nrun=5 here, but a typical choice would lies between 100 and 200:

res.multirun <- nmf(t(data.t[,-c(which(names(data.t.log) %in% c("group", "207145_at")))]), 10, nrun = 100)

# Estimating the rank factorization (r)
estim.r <- nmf(t(data.t[,-c(which(names(data.t.log) %in% c("group", "207145_at")))]), 5, nrun = 50)

W = basis(estim.r)
H = coef(estim.r)
s = featureScore(estim.r)
summary(s)
s = extractFeatures(estim.r)
str(s)

plot(estim.r)
consensusmap(estim.r, annCol = data.t[,-c(which(names(data.t.log) %in% c("group", "207145_at")))], labCol = NA, labRow = NA)


# Getting list of genes for further analysis

s1_genes = colnames(data.t)[s[[1]]]
write(s1_genes, "s1_genes.txt")

