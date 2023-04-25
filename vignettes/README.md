Place .Rmd files of analyses here.
## R Markdown

```{r, message=FALSE, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::version()
BiocManager::install("remotes", dependencies = TRUE)
BiocManager::install("waldronlab/bugSigSimple")
BiocManager::install("curatedMetagenomicData")
```
```{r}
suppressPackageStartupMessages({
  library(bugSigSimple)
  library(BugSigDBStats)
  library(bugsigdbr)
  library(tidyverse)
  library(stringr)
  library(kableExtra)
  library(dplyr)
  library(ComplexHeatmap)
})
```
## Import data from Bugsigdb database

```{r}

bsdb <- importBugSigDB(version = 'devel', cache = FALSE)
dim(bsdb)
names(bsdb)
```


# tax.levels <- c("mixed", "genus", "species")
# id.types <- c("ncbi", "metaphlan", "taxname")
tax.levels <- c("mixed")
id.types <- c("metaphlan")
exact.tax.levels <- c(TRUE, FALSE)


## Subset Signatures by Curator and Conditions

```{r}
bsdb.sub <- subsetByCurator(bsdb, "Sharmilac") %>%
  dplyr::filter(`Body site` == "Feces") %>%
  mutate(`Body site` = tolower(`Body site`)) %>%
  mutate(`Condition` = tolower(`Condition`))
table(bsdb.sub[,"Condition"])
```
# 2 of the 10 papers are missing when the above code is used. 

## Import data from Bugsigdb database (this produces all 10 papers)
```{r, messages=FALSE}
dat <-
  bsdb[which(
    bsdb$PMID == "29302014" |
    
      bsdb$PMID == "31026576" |
      bsdb$PMID == "28368458" |
      bsdb$PMID == "31990790" |
      bsdb$PMID == "32010563" |
     bsdb$PMID == "32733788" |
      bsdb$PMID == "29097493" |
      bsdb$PMID == "31337439" |
      bsdb$PMID == "28923537" ), ]

dim(dat)
names(dat)
```
#Recoding
```{r}

library(dplyr)
bsdb.subnew <- bsdb.sub %>% 
  mutate(Condition = recode_factor(Condition, `metastatic melanoma` = "response to immunochemotherapy", `response to immune checkpoint inhibitor,response to immunochemotherapy` ="response to immunochemotherapy", `response to immunochemotherapy,response to immune checkpoint inhibitor` = "response to immunochemotherapy"))



library(dplyr)
dat.recoded <- dat %>% 
  mutate(Condition = recode_factor(Condition, `metastatic melanoma` = "response to immunochemotherapy", `Response to immune checkpoint inhibitor,Response to immunochemotherapy` ="response to immunochemotherapy", `Response to immunochemotherapy,Response to immune checkpoint inhibitor` = "response to immunochemotherapy", `Response to immunochemotherapy` = "response to immunochemotherapy"))

select(dat.recoded, all_of(c("PMID","16S variable region","Location of subjects", "Source", "Group 0 name", "Group 1 name", "Condition",  "Abundance in Group 1", "Year", "Sequencing type", "Significance threshold" )))

```

## Table of studies
```{r}
bugSigSimple::createStudyTable(dat.recoded) %>%
  kbl() %>%
  kable_styling()
```

## Create a Taxon Frequency Table by response to immunochemotherapy
```{r}
dat.recoded_im <- filter(dat.recoded, dat.recoded$Condition == "response to immunochemotherapy") 
dim(dat.recoded_im)
bugSigSimple::createTaxonTable(dat.recoded_im)  %>%
  kbl() %>%
  kable_styling()
```
#The above table displays a total of 10 taxon names

## Cluster Analysis for signatures of increased abundance and decreased abundance in group 0 (non-responders) calculated using Jaccard Distance to create a distance matrix.

#did I specify g0? how?

```{r}
allsigs <- bugsigdbr::getSignatures(dat.recoded_im , tax.id.type = "taxname")
allsigs <- allsigs[sapply(allsigs, length) > 1] #require length > 1
length(allsigs)

mydists <- BugSigDBStats::calcPairwiseOverlaps(allsigs)
dim(mydists)
```

## Visualize the distribution of the signature lengths

```{r}
library(ggplot2)
siglengths <- sapply(allsigs, length)
siglengths.df <- data.frame(siglengths = siglengths)

ggplot(siglengths.df, aes(x=siglengths)) +
  geom_bar()
table(siglengths)
```

## Create a matrix of Jaccard similarities (0 for no overlap, 1 for 100% overlap)
```{r}
jmat <- BugSigDBStats::calcJaccardSimilarity(allsigs)
jmat1 <- 1-jmat
```
##Create a Clustered heatmap
```{r, fig.width=8.5, fig.height=11}
library(ComplexHeatmap)
ha <- HeatmapAnnotation(`Signature Length` = anno_barplot(siglengths))
hr <- rowAnnotation(
  `Signature Length` = anno_barplot(siglengths)
  )
hm <- Heatmap(
  jmat1,
# top_annotation = ha,
#  left_annotation = hr,
#  column_names_max_height = unit(23, "cm"),
  column_names_rot = 45,
#  row_names_max_width = unit(15, "cm"),
#  get rid of study labels
  row_labels = sub("bsdb:", "", sub("_.+", "", rownames(jmat)), fixed = TRUE),  
  column_labels = sub("bsdb:", "", sub("_.+", "", colnames(jmat)), fixed = TRUE)
)
hm
```
#This similarity is a result of having 3 different tests and multiple signatures within a single study. The same taxon were found to be increased abundantly between different tests. 

## Create a wide format dataframe for Regression Analysis
```{r}
dat_withsigs <- filter(dat.recoded_im , !is.na(dat.recoded_im$`NCBI Taxonomy IDs`))
sigs <- bugsigdbr::getSignatures(dat_withsigs, tax.id.type = "taxname")
cmat <- t(safe::getCmatrix(sigs, as.matrix = TRUE, min.size = 0, prune = FALSE))
cdf <- data.frame(cmat, stringsAsFactors = FALSE, check.names = FALSE)
cdf <- cbind(dat_withsigs, cdf)
dim(cdf)
```
An arbitrary example of meta-regression:

```{r}
fit <-
  glm(
    Prevotella ~ `Location of subjects` + `Sequencing type` + `Abundance in Group 1`,
    family = binomial(link = "logit"),
    data = cdf
  )
summary(fit)
```

## Create another heatmap on correlations of presence/absence of taxa.
```{r}
sigcors <- cor(t(cmat))
#view(sigcors)
siglengths <- sapply(sigs, length)
ha <- HeatmapAnnotation(`Signature Length` = anno_barplot(siglengths))
hr <- rowAnnotation(`Signature Length` = anno_barplot(siglengths))
hm <- Heatmap(
  sigcors,
  top_annotation = ha, left_annotation = hr,
  row_names_max_width = unit(.05, "cm"),
  column_names_max_height = unit(.1, "cm"),
 # row_labels = sub(".+:", "", rownames(sigcors)), ##removing study just to make signature names legible
  column_labels = sub(".+:", "", colnames(sigcors))
)
hm
```
