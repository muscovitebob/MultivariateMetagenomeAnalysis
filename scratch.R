library(zCompositions)
library(tidyverse)
library(ggfortify)
set.seed(100)
data = read_tsv("data.r")

functional_space = filter(data, FeatureType %in% c("GMM", "SEED"))
functional_space

united_functional_space = functional_space %>%
  unite(SampleInfo, Sample, Dataset, Status, sep="~") %>% 
  unite(FeatureInfo, FeatureType, Feature, sep="~")

X = united_functional_space %>% group_by(FeatureInfo) %>% spread(FeatureInfo, Abundance)

X = drop_na(X)
separated_X = X %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")

# train and test set code if i need it later

train = X %>% sample_frac(0.8)
test = anti_join(X, train, by="SampleInfo")

# log transform the counts

X.df.transposed = X %>% column_to_rownames("SampleInfo") %>% as.data.frame()
X.df = t(X.df.transposed)
X.df.czm = cmultRepl(t(X.df),  label=0, method="CZM")
X.df.clr = t(apply(X.df.czm, 1, function(x){log(x) - mean(log(x))}))
X.clr = X.df.clr %>% as_tibble(rownames = "SampleInfo") %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")

# do a little PCA with the log transformed results
X.df.clr.PCA = prcomp(X.df.clr)
PCA.rotation.transpose = t(X.df.clr.PCA$rotation)
plot(X.df.clr.PCA)

# without normalising, can we detect any status patterns separately?
danish = filter(separated_X, Dataset %in% c("MHD"))
autoplot(prcomp(danish[,-c(1,2,3)]), data=danish, colour="Status", x=1, y=2)

swedish = filter(separated_X, Dataset %in% c("SWE"))
autoplot(prcomp(swedish[,-c(1,2,3)]), data=swedish, colour="Status", x=1, y=2)

chinese = filter(separated_X, Dataset %in% c("CHN"))
autoplot(prcomp(chinese[,-c(1,2,3)]), data=chinese, colour="Status", x=1, y=2)

# treatment status does not align with the variance captured by the principal components at all

# what if we use transformed results?

danish.clr = filter(X.clr, Dataset %in% c("MHD"))
autoplot(prcomp(danish.clr[,-c(1,2,3)]), data=danish.clr, colour="Status", x=1, y=2)

swedish.clr = filter(X.clr, Dataset %in% c("SWE"))
autoplot(prcomp(swedish.clr[,-c(1,2,3)]), data=swedish.clr, colour="Status", x=1, y=2)

# here in swedish we see an interesting pattern not previously observed - the metformin- is mostly not in PC1, instead being explained by PC2 mostly, and in PC2 it aligns with control subjects. 

chinese.clr = filter(X.clr, Dataset %in% c("CHN"))
autoplot(prcomp(chinese.clr[,-c(1,2,3)]), data=chinese.clr, colour="Status", x=1, y=2)

# its hard to see a separation with danish and chinese, but swedish does have some

# what if we attempt normalisation using the sum of each row?

X.df.transposed[,1:ncol(X.df.transposed)] = X.df.transposed[,1:ncol(X.df.transposed)]/rowSums(X.df.transposed[,1:ncol(X.df.transposed)])
X.sumrownormal =  X.df.transposed %>% as_tibble(rownames = "SampleInfo") %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")

danish.norm = filter(X.sumrownormal, Dataset %in% c("MHD"))
autoplot(prcomp(danish.norm[,-c(1,2,3)]), data=danish.norm, colour="Status", x=1, y=2)

swedish.norm = filter(X.sumrownormal, Dataset %in% c("SWE"))
autoplot(prcomp(swedish.norm[,-c(1,2,3)]), data=swedish.norm, colour="Status", x=1, y=2)

chinese.norm = filter(X.sumrownormal, Dataset %in% c("CHN"))
autoplot(prcomp(chinese.norm[,-c(1,2,3)]), data=chinese.norm, colour="Status", x=1, y=2)

# we succeed in making PC1 capture large propotions of the variance, but again treatment status does not seem to play a large role

# what if we log transform the normalised data?

X.sumrownormal_united = X.sumrownormal %>% unite(SampleInfo, Sample, Dataset, Status, sep="~") %>% unite(FeatureInfo, FeatureType, Feature, sep="~")

X.sumrownormal_united= X %>% column_to_rownames("SampleInfo") %>% as.data.frame()
X.sumrownormal.df = t(X.sumrownormal_united)
X.sumrownormal.df.czm = cmultRepl(t(X.sumrownormal.df),  label=0, method="CZM")
X.sumrownormal.df.clr = t(apply(X.sumrownormal.df.czm, 1, function(x){log(x) - mean(log(x))}))
X.sumrownormal.clr = X.sumrownormal.df.clr %>% as_tibble(rownames = "SampleInfo") %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")

danish.sumrownormal.clr = filter(X.sumrownormal.clr, Dataset %in% c("MHD"))
autoplot(prcomp(danish.sumrownormal.clr[,-c(1,2,3)]), data=danish.sumrownormal.clr, colour="Status", x=1, y=2)

swedish.sumrownormal.clr = filter(X.sumrownormal.clr, Dataset %in% c("SWE"))
autoplot(prcomp(swedish.sumrownormal.clr[,-c(1,2,3)]), data=swedish.sumrownormal.clr, colour="Status", x=1, y=2)

chinese.sumrownormal.clr = filter(X.sumrownormal.clr, Dataset %in% c("CHN"))
autoplot(prcomp(chinese.sumrownormal.clr[,-c(1,2,3)]), data=chinese.sumrownormal.clr, colour="Status", x=1, y=2)

# similar results as if we just log transformed the data, perhaps predictably

# it seems that log transforming the data is the most interesting way forward
# also the data has already been normalised before so the normalisation steps were spurious

# OK, we carry log transformed forward for ordination analysis via dbRDA

library(vegan)
constrainedModel3 = capscale(X.clr[, -c(1,2,3)] ~ Status + Condition(Dataset), data=X.clr, method="canberra")
constrainedModel3Summary = summary(constrainedModel1)
ggpairs(as.tibble(constrainedModel3Summary$species), progress = FALSE)

