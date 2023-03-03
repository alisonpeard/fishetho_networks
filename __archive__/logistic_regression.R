library(wCorr)
library(ggplot2)
library(reshape2)
library(MASS)

# glm giving huge coefficients, might be good to revisit

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

likelihood[criteria] <- likelihood[criteria]
C <- length(criteria)
index <- rownames(likelihood)

# set up bootstrap
n <- 30  # size of bootstrap samples
betas = array(, c(C, C))
p.vals = array(, c(C, C))

for (i in 1:C){
  for (j in 1:C){

    column1 <- criteria[i]
    column2 <- criteria[j]
    
    # get vars (removing NULLS)
    x <- likelihood[, column1] / 2
    y <- likelihood[, column2] / 2
    w <- certainty[, column1] + certainty[, column1] + 1
    
    # get pairwise complete samples
    idx <- which(!is.na(x))
    idy <- which(!is.na(y))
    ids <- intersect(idx, idy)
    x <- x[ids]
    y <- y[ids]
    w <- w[ids]
    
    if ((sum(x) > 0) && sum(y) > 0){
      x <- as.factor(x)
      y <- as.factor(y)
      # lm.no.cert <- polr(y ~ x)  # for potential
      lm.no.cert <- glm(y ~ x, family='binomial', weights=w)
      summary(lm.no.cert)
      
      betas[i, j] <- summary(lm.no.cert)$coefficients[2,1]
      p.vals[i, j] <- summary(lm.no.cert)$coefficients[2,4]
    }
  }
}


# means <- apply(betas, c(2, 3), mean, na.rm=TRUE)
#mstds <- apply(betas, c(2, 3), sd, na.rm=TRUE)

# means.p <- apply(p.vals, c(2, 3), mean, na.rm=TRUE)
# stds.p <- apply(p.vals, c(2, 3), sd, na.rm=TRUE)

# plotting
melted.bmat <- melt(betas, na.rm = TRUE)
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
