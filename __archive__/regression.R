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
B <- 100 # number of bootstrap samples
n <- 30  # size of bootstrap samples
betas = array(, c(B, C, C))
p.vals = array(, c(B, C, C))

# start bootstrap sampling
set.seed(1)
for (b in 1:B){
  # sample dataset
  sample.b <- sample(index, size=n, replace=TRUE)
  for (i in 1:C){
    for (j in 1:C){
  
      column1 <- criteria[i]
      column2 <- criteria[j]
      
      # get vars (removing NULLS)
      x <- likelihood[sample.b, column1]
      y <- likelihood[sample.b, column2]
      w <- certainty[sample.b, column1] + certainty[sample.b, column1] + 1
      
      # get pairwise complete samples
      idx <- which(!is.na(x))
      idy <- which(!is.na(y))
      ids <- intersect(idx, idy)
      x <- x[ids]
      y <- y[ids]
      w <- w[ids]
      
      if ((sum(x) > 0) && sum(y) >0){
        lm.no.cert <- lm(y ~ x)
        summary(lm.no.cert)
        
        betas[b, i, j] <- summary(lm.no.cert)$coefficients[2,1]
        p.vals[b, i, j] <- summary(lm.no.cert)$coefficients[2,4]
      }
    }
  }
}

means <- apply(betas, c(2, 3), mean, na.rm=TRUE)
stds <- apply(betas, c(2, 3), sd, na.rm=TRUE)

means.p <- apply(p.vals, c(2, 3), mean, na.rm=TRUE)
stds.p <- apply(p.vals, c(2, 3), sd, na.rm=TRUE)

# plotting
melted.bmat <- melt(means, na.rm = TRUE)
beta.heatmap <- ggplot(data=melted.bmat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="blue", high="red", mid="white",
                       midpoint=0, limit=c(-1,1.01), space='Lab', 
                       name="Weighted\nSpearman\nCorrelation") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
beta.heatmap

melted.pmat <- melt(means.p, na.rm = TRUE)
p.heatmap <- ggplot(data=melted.pmat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="red", high="blue", mid="white",
                       midpoint=0.05, limit=c(0, 0.1), space='Lab', 
                       name="p-value") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
p.heatmap

melted.pmat <- melt(stds.p, na.rm = TRUE)
p.heatmap <- ggplot(data=melted.pmat, aes(x=Var1, y=Var2, fill=value)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="red", high="blue", mid="white",
                       midpoint=0.05, limit=c(0, 0.1), space='Lab', 
                       name="p-value") + 
  scale_y_reverse(breaks=c(1:10), labels=criteria.labels) + 
  scale_x_continuous(breaks=c(1:10), labels=criteria.labels) + 
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=45, vjust=1, size=12, hjust=1)) +
  theme(axis.text.y=element_text(angle=0, vjust=1, size=12, hjust=1)) +
  coord_fixed()
p.heatmap
