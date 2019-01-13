library(zCompositions)
library(tidyverse)
library(ggfortify)
library(PPtreeViz)
set.seed(100)
data = read_tsv("data.r")

functional_space = filter(data, FeatureType %in% c("GMM", "SEED"))

united_functional_space = functional_space %>%
  unite(SampleInfo, Sample, Dataset, Status, sep="~") %>% 
  unite(FeatureInfo, FeatureType, Feature, sep="~")

X = united_functional_space %>% group_by(FeatureInfo) %>% spread(FeatureInfo, Abundance)

X = drop_na(X)
separated_X = X %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")
X.df = X %>% column_to_rownames("SampleInfo") %>% as.data.frame()
X.df.czm = cmultRepl(X.df,  label=0, method="CZM")
X.df.clr = t(apply(X.df.czm, 1, function(x) {log(x) - mean(log(x))} ))
annotated_X.clr = bind_cols(separated_X[1:3], as.tibble(X.df.clr))

proj1 = LDAopt(annotated_X.clr$Status, X.df.clr, q=2)
proj2 = LDAopt(annotated_X.clr$Dataset, X.df.clr, q=2)
PPoptViz(proj1)
PPoptViz(proj2)
# amazing clusters with projection pursuit

# what if we try doing this on the taxonomic space?

taxonomic_space = filter(data, FeatureType %in% c("Genus"))
united_taxonomic_space = taxonomic_space %>%
  unite(SampleInfo, Sample, Dataset, Status, sep="~") %>% 
  unite(FeatureInfo, FeatureType, Feature, sep="~")

X.tax = united_taxonomic_space %>% group_by(FeatureInfo) %>% spread(FeatureInfo, Abundance)

X.tax = drop_na(X.tax)
separated_X.tax = X.tax %>% separate(SampleInfo, into=c("Sample", "Dataset", "Status"), sep="~")
X.tax.df = X.tax %>% column_to_rownames("SampleInfo") %>% as.data.frame()
X.tax.df.czm = cmultRepl(X.tax.df,  label=0, method="CZM")
X.tax.df.clr = t(apply(X.tax.df.czm, 1, function(x) {log(x) - mean(log(x))} ))
annotated_X.tax.clr = bind_cols(separated_X.tax[1:3], as.tibble(X.tax.df.clr))

proj3 = LDAopt(annotated_X.tax.clr$Status, X.tax.df.clr, q=2)
proj4 = LDAopt(annotated_X.tax.clr$Dataset, X.tax.df.clr, q=2)
PPoptViz(proj3)
PPoptViz(proj4)

# drop in some Genus level capscale analysis for fun too
library(vegan)
constrainedModel8 = capscale(X.tax[,-1] ~ Status + Condition(Dataset), data=separated_X.tax, distance="bray")
CAP8_results=bind_cols(separated_X.tax[,1:3],as.tibble(scores(constrainedModel8)$sites))
ggplot(CAP8_results, aes(x=CAP1, y=CAP2,colour=Status)) + geom_point()

