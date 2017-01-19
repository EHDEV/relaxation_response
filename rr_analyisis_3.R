

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
obesity = read.csv("./data/obesity.tsv", header=T, stringsAsFactors = F, sep="\t")
alzhymers = read.csv("./data/alzhymers.tsv", header=T, stringsAsFactors = F, sep="\t")
t2diabetes = read.csv("./data/t2diabetes.tsv", header=T, stringsAsFactors = F, sep="\t")

disease_genes_all = c(diabetes$Symbol, obesity$Symbol, alzhymers$Symbol, t2diabetes$Symbol)

lb = which(data$IDENTIFIER %in% disease_genes_all)

length(disease_genes_all)


rownames(data) = data$ID_REF
data$ID_REF = NULL
data$IDENTIFIER = NULL

data.t <- na.omit(data.frame(t(data[lb,])))
colnames(data.t) = rownames(data[lb, ])

length(which(is.na(data.t)))

for (i in 1:dim(data.t)[2])
  if (class(data.t[,i]) != "numeric")
    data.t[,i] = as.numeric(as.character(data.t[,i]))
data.t = na.omit(data.t)

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

ct.data = data.t[which(sample_class == 0), ]
st.data = data.t[which(sample_class == 1), ]
lt.data = data.t[which(sample_class == 2), ]

#### Begin
# Testing null hypothesis that variables have no relationship with response
# source("http://bioconductor.org/biocLite.R")
# biocLite("globaltest")
library(globaltest)
ass.test = gt(response = sample_class, data.t)
# p-value of .0307 < .05 thus null hypothesis is rejected for the alternative.
#### End 

####Plot
par(mfrow=c(1,3)) 
boxplot(log(data.t[sample_class==0,21]), main="Boxplot for a Diff Expressed Gene - Control", ylim=c(3,5), col=8)
boxplot(log(data.t[sample_class==1,21]), main="Boxplot for a Diff Expressed Gene - Short Term", ylim=c(3,5), col=5)
boxplot(log(data.t[sample_class==2,21]), main="Boxplot for a Diff Expressed Gene - Long Term",  ylim=c(3,5), col=3)
#####

set.seed(1000)
lt_ct_ids = which(sample_class %in% c(0,2))
st_ct_ids = which(sample_class %in% c(0,1))
lt_st_ids = which(sample_class %in% c(1,2))

#binary responses
lt_ct_sample = as.factor(as.character(sample_class[lt_ct_ids]))
st_ct_sample = as.factor(as.character(sample_class[st_ct_ids]))
lt_st_sample = as.factor(as.character(sample_class[lt_st_ids]))


# take a random test data
fn_test_data = function(x){
  set.seed(1001)
  return (sample(1:length(x), round(.2 * length(x))))
}

testids_lt.ct = fn_test_data(lt_ct_sample)
testids_st.ct = fn_test_data(st_ct_sample)
testids_lt.st = fn_test_data(lt_st_sample)

# rand_subjects = sample(1:length(ctl_ids)+length(lt_ids)-1, length(ctl_ids)+length(lt_ids))


fit_cvl = cvl(data.t.log$group, data.t.log[,-c(which(names(data.t.log) %in% c("group", "207145_at")))], lambda1=1.5, fold=10)
fit_opt_lt <- optL1(sample_class[rand_subjects], penalized = data.t[rand_subjects,lt.idx],fold = 10)

fit_opt_lt <- optL1(lt_ct_sample[-test_ids], penalized = data.t[lt_ct_ids[-test_ids],], minlambda1 = .1, maxlambda1 = 1000,fold = 10)
fit_opt_lt <- optL2(lt_ct_sample[-test_ids], penalized = data.t[lt_ct_ids[-test_ids],], minlambda2 = .1, maxlambda2 = 1000,fold = 10)

fit_prof = profL1(lt_ct_sample[-test_ids], penalized = data.t[lt_ct_ids[-test_ids],], minlambda1 = .1, maxlambda1 = 1000,fold = 10)
fit_prof2 = profL2(lt_ct_sample, data.t[lt_ct_ids,], fold=10, plot=TRUE, minlambda2 = 0.01, maxlambda2 = 1000)

fit <- cvl(lt_ct_sample[-test_ids], penalized = data.t[lt_ct_ids[-test_ids],], lambda1 = 1000, lambda2 = 50,fold = 10)

pen_fit <- penalized(lt_ct_sample[-test_ids], penalized = data.t[lt_ct_ids[-test_ids],], lambda1 = 1000, lambda2 = 5)

pred2 = predict(pen_fit, data.t[testids_lt.ct,])
table(pred2 >= .5)
table(pred >= .5)
table(lt_ct_sample[test_ids] == 2)

pred = predict(fit_prof, data.t[test_ids,])
sample_class = as.factor(as.character(sample_class))
par(mfrow=c(1,1))

plotpath(fit_prof$fullfit, log="x")
fit_opt_lt$predictions

table(fit_opt_lt$predictions >= .5)

table(lt_ct_sample)

#### ElasticNet Train, CV, and test for Short-term vs. Control ####

fit_opt_lt <- optL1(st_ct_sample[-testids_st.ct], penalized = data.t[st_ct_ids[-testids_st.ct],], minlambda1 = .1, maxlambda1 = 1000,fold = 10)
fit_opt_lt <- optL2(st_ct_sample[-testids_st.ct], penalized = data.t[st_ct_ids[-testids_st.ct],], minlambda2 = .1, maxlambda2 = 1000,fold = 10)

pen_fit.st <- penalized(st_ct_sample[-testids_st.ct], penalized = data.t[st_ct_ids[-testids_st.ct],], lambda1 = 1000, lambda2 = 1000)

pred_st = predict(pen_fit.st, data.t[testids_st.ct,])
table(pred_st >= .5)
table(st_ct_sample[testids_st.ct] == 1)

rc.st = roc(st_ct_sample[testids_st.ct], pred_st)
pr.st = prediction(pred_st, st_ct_sample[testids_st.ct] == 1)
pf.st = performance(pr.st, "tpr", "fpr")
plot(pf.st)
plot(rc.st, main="Short term vs. Control", col=3)

par(mfrow=c(1,3))

#### ElasticNet Train, CV, and test for Short-term vs. Long-term ####

fit_opt_lt.st1 <- optL1(lt_st_sample[-testids_lt.st], penalized = data.t[lt_st_ids[-testids_lt.st],], minlambda1 = .1, maxlambda1 = 1000,fold = 10)
fit_opt_lt.st2 <- optL2(lt_st_sample[-testids_lt.st], penalized = data.t[lt_st_ids[-testids_lt.st],], minlambda2 = .1, maxlambda2 = 1000,fold = 10)

pen_fit.st.lt <- penalized(lt_st_sample[-testids_lt.st], penalized = data.t[lt_st_ids[-testids_lt.st],], lambda1 = 1000, lambda2 = 1000)

pred_st.lt = predict(pen_fit.st.lt, data.t[testids_lt.st,])
table(pred_st.lt >= .5)
table(pred2 >= .5)
table(lt_st_sample[testids_lt.st] == 2)

rc.lt.st = roc(lt_st_sample[testids_lt.st], pred_st.lt)
pr.st.lt = prediction(pred_st.lt, lt_st_sample[testids_lt.st] == 1)
pf.lt.st = performance(pr.st.lt, "tpr", "fpr")

plot(pf.st, colorize = TRUE)
plot(pf.lt.st, add = TRUE, colorize = TRUE)

plot(rc.lt.st, col=2, main="Long term vs. Short term")

#### ElasticNet Train, CV, and test for Long-term vs. Control ####

fit_opt_lt <- optL1(lt_ct_sample[-testids_lt.ct], penalized = data.t[lt_ct_ids[-testids_lt.ct],], minlambda1 = .1, maxlambda1 = 1000,fold = 10)
fit_opt_lt <- optL2(lt_ct_sample[-testids_lt.ct], penalized = data.t[lt_ct_ids[-testids_lt.ct],], minlambda2 = .1, maxlambda2 = 1000,fold = 10)

pen_fit.lt <- penalized(lt_ct_sample[-testids_lt.ct], penalized = data.t[lt_ct_ids[-testids_lt.ct],], lambda1 = 1000, lambda2 = 1000)

pred_lt = predict(pen_fit.lt, data.t[testids_lt.ct,])

table(pred_lt >= .69)
table(pred_lt >= .5)
table(lt_ct_sample[testids_lt.ct] == 2)

rc.lt = roc(as.numeric(lt_ct_sample[testids_lt.ct] == 2), as.numeric(pred_lt>0.5))
plot(rc.lt)


pr.lt = prediction(pred_lt, lt_ct_sample[testids_lt.ct] == 2)
pf.lt = performance(pr.lt, "tpr", "fpr")

plot(pf.st, colorize = TRUE)
plot(pf.lt.st, colorize = TRUE)
plot(pf.lt, colorize = TRUE)

plot(rc.lt, col=4, main="Long term vs. Control")


##############



library(pROC)
library(ROCR)
?roc

rc = roc(lt_st_sample[testids_lt.st], pred_st.lt)
plot(rc)

detach("package:ROCR", unload=TRUE)

######### Post-Elastic Net ##########
# capture the list of genes

lt.st.diff.genes = sort(coef(pen_fit.st.lt), decreasing = T)
lt.ct.diff.genes = sort(coef(pen_fit.lt), decreasing = T)
st.ct.diff.genes = sort(coef(pen_fit.st), decreasing = T)

all.diff.genes = unique(c(names(lt.st.diff.genes)[-1], names(lt.ct.diff.genes)[-1], names(st.ct.diff.genes)[-1]))

write(names(lt.st.diff.genes)[-1], "./out/lt.st.diff.genes.txt")
write(names(lt.ct.diff.genes)[-1], "./out/lt.ct.diff.genes.txt")
write(names(st.ct.diff.genes)[-1], "./out/st.ct.diff.genes.txt")
write(all.diff.genes, "./out/all.diff.genes.txt")

#clustering and some plots
?kmeans

kx = kmeans(data.t[,which(colnames(data.t) %in%all.diff.genes)], centers=4, iter.max = 1000)
d = dist(data.t[,which(colnames(data.t) %in%all.diff.genes)], method='euclidean', na.action(na.omit))
h = hclust(d, "complete")
h.dend = as.dendrogram(h)
plot(h, hang=-1)

hc.cut = cutree(h, 3)

heatmap(as.matrix(data.t[,which(colnames(data.t) %in%all.diff.genes)]), Rowv=h.dend, Colv=sample_class, scale="column", main="Heatmap of Differentially Expressed Genes from Control, 8 week and Lifelong Subjects", ylab="Samples", xlab="Genes")
plot(kx, data=data.t[, which(colnames(data.t) %in%all.diff.genes)], title = "PCA Result with K-Means Clustering (k=3)")

################ Post-Elastic Net #############
tt.genes.ltst = colnames(lt.data)[which(p.st.lt <= 0.01)]
tt.genes.ltct = colnames(lt.data)[which(p.ct.lt <= 0.01)]
tt.genes.stct = colnames(lt.data)[which(p.ct.st <= 0.01)]


length(which(names(st.ct.diff.genes)[-1] %in% tt.genes.stct))
length(which(names(lt.ct.diff.genes)[-1] %in% tt.genes.ltct))
length(which(names(st.ct.diff.genes)[-1] %in% tt.genes.stct))



mean_ct_values = apply(ct.data[,which(colnames(ct.data) %in% names(st.ct.diff.genes))], 2, mean)
mean_st_values = apply(st.data[,which(colnames(st.data) %in% names(st.ct.diff.genes))], 2, mean)

mean_lt_values = apply(lt.data[,which(colnames(lt.data) %in% names(lt.ct.diff.genes))], 2, mean)
mean_ctlt_values = apply(ct.data[,which(colnames(ct.data) %in% names(lt.ct.diff.genes))], 2, mean)
length(which(mean_lt_values > mean_ctlt_values))
length(which(!(mean_lt_values > mean_ctlt_values)))


median_ct_values = apply(ct.data[,which(colnames(ct.data) %in% names(st.ct.diff.genes))], 2, median)
median_st_values = apply(st.data[,which(colnames(st.data) %in% names(st.ct.diff.genes))], 2, median)

