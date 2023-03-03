library(wCorr)
library(ggplot2)
library(reshape2)

# load the Python-processed files (still have negative values)
wd <- file.path("/Users", "alison", "Documents", "FishEtho", "data")
potential <- read.csv(file.path(wd, "potential.csv"))
likelihood <- read.csv(file.path(wd, "likelihood.csv"))
certainty <- read.csv(file.path(wd, "certainty.csv"))

# clean up negative values
likelihood[likelihood < 0] <- NA
potential[potential < 0] <- NA
certainty[certainty < 0] <- NA

criteria <- c("X17", "X18", "X29", "X30", "X31", "X32", "X33", "X34", "X35", "X36")
criteria.labels <- c("Home range", "Depth range", "Migration", "Reproduction",
                     "Aggregation", "Aggression", "Substrate", "Stress",
                     "Malformation", "Slaughter")

C <- length(criteria)
N <- nrow(likelihood)
index <- rownames(likelihood)

wcorrs <- array(, c(C, C))
present.counts <- array(, c(C, C))

for (i in 1:C){
  for (j in 1:C){
    # what to calculate correlation for
    column1 <- criteria[i]
    column2 <- criteria[j]
    
    # get vars (removing NULLS)
    x <- likelihood[, column1]
    y <- likelihood[, column2]
    w <- certainty[, column1] + certainty[, column1] + 1
    idx <- which(!is.na(x))
    idy <- which(!is.na(y))
    ids <- intersect(idx, idy)
    x <- x[ids]
    y <- y[ids]
    w <- w[ids]
    
    # set up nan-correcting
    p = length(ids)
    present.counts[i, j] <- p / N
    present.counts[j, i] <- p / N
    
    if (p > 0){
      # calculate weighted correlation (NaN for all zeros)
      wcorr <- weightedCorr(x, y, method="Spearman", weights=w)
      wcorrs[i, j] <- wcorr
      wcorrs[j, i] <- wcorr
    }
    
  }
}

wcorrs
write.csv(wcorrs, "/Users/alison/Documents/FishEtho/data/spearman/likelihood/total_corrs.csv")
write.csv(present.counts, "/Users/alison/Documents/FishEtho/data/spearman/likelihood/percentage_present.csv")

# plotting corrs for full dataset
melted.cormat <- melt(wcorrs, na.rm=TRUE, value.name="correlation", varnames=c("Criterion1", "Criterion"))
mean.heatmap <- ggplot(data=melted.cormat, aes(x=Criterion1, y=Criterion, fill=correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1,1), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
mean.heatmap

# plotting %present for full dataset
melted.cormat <- melt(present.counts, na.rm=TRUE, value.name="correlation", varnames=c("Criterion1", "Criterion"))
mean.heatmap <- ggplot(data=melted.cormat, aes(x=Criterion1, y=Criterion, fill=correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0.25, limit=c(0, 1), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  ggtitle("% present") + 
  coord_fixed()
mean.heatmap

# multiplying by %present
wcorr.present <- wcorrs * present.counts
melted.cormat <- melt(wcorr.present, na.rm=TRUE, value.name="correlation", varnames=c("Criterion1", "Criterion"))
mean.heatmap <- ggplot(data=melted.cormat, aes(x=Criterion1, y=Criterion, fill=correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1, 1), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
mean.heatmap


# set up bootstrap
#########
B <- 100  # number of bootstrap samples
n <- 30  # size of bootstrap samples
wcorrs <- array(, c(B, C, C))

# start bootstrap sampling
set.seed(1)
for (b in 1:B){
  # sample dataset
  sample.b <- sample(index, size=n, replace=TRUE)
  for (i in 1:C){
    for (j in 1:C){
      # what to calculate correlation for
      column1 <- criteria[i]
      column2 <- criteria[j]
      
      # get vars (removing NULLS)
      x <- likelihood[sample.b, column1]
      y <- likelihood[sample.b, column2]
      w <- certainty[sample.b, column1] + certainty[sample.b, column1] + 1
      idx <- which(!is.na(x))
      idy <- which(!is.na(y))
      ids <- intersect(idx, idy)
      x <- x[ids]
      y <- y[ids]
      w <- w[ids]
      
      # set up nan-correcting
      p = length(ids)
      t = length(sample.b)
      
      if (p > 0){
        # calculate weighted correlation (NaN for all zeros)
        wcorr <- weightedCorr(x, y, method="Spearman", weights=w)
        wcorrs[b, i, j] <- wcorr
      }
    }
  }
}

means <- apply(wcorrs, c(2, 3), mean, na.rm=TRUE)
stds <- apply(wcorrs, c(2, 3), sd, na.rm=TRUE)

# plotting means
melted.cormat <- melt(means, na.rm=TRUE, value.name="correlation", varnames=c("Criterion1", "Criterion"))
mean.heatmap <- ggplot(data=melted.cormat, aes(x=Criterion1, y=Criterion, fill=correlation)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1,1), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_bw() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
mean.heatmap


# plotting standard deviations
melted.cormat <- melt(stds, na.rm = TRUE, varnames=c("Criterion1", "Criterion"))
stds.heatmap <- ggplot(data=melted.cormat, aes(x=Criterion1, y=Criterion, fill=value)) +
  geom_tile(color="white") +
  # geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(0.0,0.5), space="Lab", 
                       name="std. dev") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
stds.heatmap

