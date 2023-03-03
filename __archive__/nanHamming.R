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
index <- rownames(likelihood)

# set up bootstrap
B <- 100  # number of bootstrap samples
n <- 30  # size of bootstrap samples
nanhams <- array(, c(B, C, C))

# start bootstrap sampling
set.seed(1)
for (b in 1:B){
  for (i in 1:C){
    for (j in 1:C){
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
      
      nanhams[b, i, j] <- (t / p) * sum(abs(x - y))
    }   
  }
}

means <- apply(nanhams, c(2, 3), mean, na.rm=TRUE)
stds <- apply(nanhams, c(2, 3), sd, na.rm=TRUE)

# plotting
melted.cormat <- melt(means, na.rm = TRUE)
mean.heatmap <- ggplot(data=melted.cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=5, limit=c(0,10), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
mean.heatmap


# plotting standard deviations
melted.cormat <- melt(stds, na.rm = TRUE)
stds.heatmap <- ggplot(data=melted.cormat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="white") +
  # geom_text(aes(label = round(value, 2))) +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0.5, limit=c(0, 1), space="Lab", 
                       name="std. dev") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
stds.heatmap


