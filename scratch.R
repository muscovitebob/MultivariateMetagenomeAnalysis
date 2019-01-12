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
X.clr.abundances = X.clr[,-c(1,2,3)]
annotated_X.clr = bind_cols(separated_X[1:3], as.tibble(X.df.clr))

# do a little PCA with the log transformed results
X.df.clr.PCA = prcomp(X.df.clr)
autoplot(X.df.clr.PCA, data=annotated_X.clr, colour="Status", x=1, y=2)

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
constrainedModel1 = capscale(separated_X[, -c(1,2,3)] ~ Status + Condition(Dataset), data=separated_X, method="canberra")
constrainedModel1Summary = summary(constrainedModel1)
CAP1_results = bind_cols(separated_X[,1:3],as.tibble(constrainedModel1$Ybar))
ggplot(CAP1_results, aes(x=Dim1, y=Dim2,colour=Status)) + geom_point() 

constrainedModel3 = capscale(X.clr.abundances ~ Status + Condition(Dataset), data=X.clr, distance="canberra")
constrainedModel3Summary = summary(constrainedModel3)
CAP3_results = bind_cols(X.clr[,1:3],as.tibble(constrainedModel3$Ybar))
ggplot(CAP3_results, aes(x=Dim1, y=Dim2,colour=Status)) + geom_point() 

# what if we try the origin dataset format instead of the abundance matrix?

constrainedModel4 = capscale(Abundance ~ Status + Condition(Dataset), data=functional_space)
# this is too memory intensive to be viable

constrainedModel5 = capscale(separated_X[,-c(1,2,3)] ~ Status + Dataset, data=separated_X, distance="canberra")
CAP5_results = bind_cols(separated_X[,1:3],as.tibble(constrainedModel5$Ybar))
ggplot(CAP5_results, aes(x=Dim1, y=Dim2,colour=Status)) + geom_point() 

constrainedModel6 = capscale(X.clr.abundances ~ Status + Dataset, data=X.clr, distance="canberra")
CAP6_results = bind_cols(X.clr[,1:3],as.tibble(constrainedModel6$Ybar))
ggplot(CAP6_results, aes(x=Dim1, y=Dim2,colour=Status)) + geom_point() 

constrainedModel7 = capscale(X.clr.abundances ~ Condition(Status), data=X.clr, distance="canberra")
constrainedModel7Summary = summary(constrainedModel3)
CAP7_results = bind_cols(X.clr[,1:3],as.tibble(constrainedModel7$Ybar))
ggplot(CAP7_results, aes(x=Dim1, y=Dim2,colour=Status)) + geom_point()

# what about k means clustering

library(cluster) 
library(factoextra)
fourMeansLog = kmeans(X.clr.abundances, centers=4)
# this plotter does PCA and colours points in PCA by cluster
fviz_cluster(fourMeansLog, data=X.clr.abundances)
# these clusters are not very convincing
fviz_nbclust(X.clr.abundances, kmeans, method="silhouette")
# apparently two clusters are the best but this doesnt really correspond to what we want - which is discrimination by clinical status, of which there are four levels
# if we look back to the PCA plot coloured by Status or Dataset we can clearly see that kmeans actually does a pretty solid job of identifying MHD and SWE as separate Datasets. It's even pretty good at seeing CHN as a separate cloud. Nothing to do with clinical status, though.

fourMeans = kmeans(X[-1], centers=4)
fviz_cluster(fourMeans, data=X[-1])
# produces complete garbage
