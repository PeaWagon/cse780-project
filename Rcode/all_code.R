# to run: Rscript all_code.R

#####################################################
#              hierarchical clustering              #
#####################################################

# to compare the clustering results
library(fpc)

files <- list(
    "../data_processing/fp2_dissimilarity.csv",
    "../data_processing/fp3_dissimilarity.csv",
    "../data_processing/fp4_dissimilarity.csv",
    "../data_processing/maccs_dissimilarity.csv"
)

fingerprints <- list("fp2", "fp3", "fp4", "maccs")

linkages <- list("single", "complete", "average", "ward.D", "ward.D2")

fpnum <- 1 # keep track of fingerprint

for (fname in files) {

    # read in data
    d <- read.csv(fname)
    # get rid of first column (cid)
    d <- d[,-1]
    # make the matrix a distance/dissimilarity matrix
    d <- as.dist(d)

    # cannot compare the results for the different
    # fingerprints since they use different distance
    # matrices - can compare results for one
    # fingerprint
    # let's choose FP3 since that performed well for
    # kmedoids analysis
    # do a comparison of ward.D, complete, avg
    if (fpnum == 2) {
        # number of cuts to make: 50, 100, 150
        for (cut in 1:3) {
            hc1 <- hclust(d, "ward.D")
            hc2 <- hclust(d, "complete")
            hc3 <- hclust(d, "average")
            fit1 <- cutree(hc1, cut*50)
            fit2 <- cutree(hc2, cut*50)
            fit3 <- cutree(hc3, cut*50)
            c1 <- cluster.stats(d, fit1, fit2)
            oname <- paste("fp3_wardD_complete_",
                           cut*50, "-cuts.rds", sep="")
            saveRDS(c1, oname)
            c2 <- cluster.stats(d, fit1, fit3)
            oname <- paste("fp3_wardD_avg_",
                           cut*50, "-cuts.rds", sep="")
            saveRDS(c2, oname)
            c3 <- cluster.stats(d, fit2, fit3)
            oname <- paste("fp3_complete_avg_",
                           cut*50, "-cuts.rds", sep="")
            saveRDS(c3, oname)
        }
    }

    # make some hierarchical clustering plots
    for (linkage in linkages) {
        hc <- hclust(d, linkage)
        outfile <- paste(fingerprints[fpnum],
                 "_hclust_", linkage, ".png",
                 sep="")
        png(outfile)
        plot(hc, sub="", xlab=paste(linkage,
             fingerprints[fpnum]), ylab="",
             cex.axis=1.5, cex.lab=3, cex.main=3)
        dev.off()
    }

    # increment fingerprint name
    fpnum <- fpnum + 1
}

# example of how to get the clustering results
data <- readRDS("fp3_complete_avg_50-cuts.rds")
data$vi             # 1.313243
data$corrected.rand # 0.5259344

#####################################################
#              kmedoids clustering                  #
#####################################################

# for using silhouette function
library(cluster)

# example silhouette plot
d <- read.csv(files[[2]])
d <- d[,-1]
# need to put as data.matrix otherwise the title
# will read dist=d not dmatrix=d (and using a data
# frame or distance matrix - as.dist - will result
# in an error from the cluster package)
d <- data.matrix(d)
kmed <- pam(d, 50, diss=TRUE)
sil1 <- silhouette(kmed$clustering, dmatrix=d)
# png/jpeg don't work (don't show silhouettes)
pdf("silhouette_kmedoids_k50_fp3.pdf")
plot(sil1)
dev.off()
#kmed$silinfo$avg.width

# make a matrix of zeros with 4 rows and 10 columns
silwidths <- matrix(0, 4, 10)
# access with silwidths[row, column]

# make silhouette plots comparing the average
# silhouette width for cluster sizes from 25-250 in
# increments of 25
numClusters <- c(1:10*25)

findex <- 1 # keep track of current file index
for (fname in files) {
    # read in data
    d <- read.csv(fname)
    # get rid of first column (cid)
    d <- d[,-1]
    # convert from data frame to matrix
    d <- data.matrix(d)
    for (k in numClusters) {
        # kmedoids clustering with dissimilarity matrix
        kmed <- pam(d, k, diss=TRUE)
        # get silhouette information
        avsil <- kmed$silinfo$avg.width
        silwidths[findex, k/25] <- avsil
    }
    findex <- findex + 1
}

saveRDS(silwidths, "silwidths.rds")

# basic way to plot the average silhouette width (y) for values of k (x)
png("silwidths_plot.png")
matplot(t(silwidths), type="l")
dev.off()

#####################################################
#                 arules analysis                   #
#####################################################

library(arules)      # for apriori function
source("std_lift.R") # professor's code to calculate
                     # standardised lift

infile <- "../data_processing/drug_data_trimmed_text_len.csv"

d <- read.csv(infile)

# we can consider:
# drugName (1), conditon (2), rating (4), revLenDesc(9)
x <- d[, c(1,2,4,9)]

# put the rating in two categories: >= 7 and < 7
# sets rating variable to 0 for values below 7
x$rating[x$rating<7] <- 0

# sets rating variable to 1 for values above or equal to 7
x$rating[x$rating>0] <- 1

# make the option 0 or 1
x$rating <- as.factor(x$rating)

# minlen is number of parameters to consider
# maxlen is how many parameters there are for the rules
# if no rules are generated, reduce support and/or
# confidence

# set the right-hand side or the left-hand side
# to see what rules there are (a)
params <- list(support=0.01, confidence=0.8, minlen=2,
               maxlen=4)

app <- list(rhs=c("rating=0", "rating=1"),
            default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# see where drugs are poorly rated (b)
params <- list(support=0.005, confidence=0.6, minlen=2,
               maxlen=4)

app <- list(rhs=c("rating=0"), default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# see if length of text is a rule for anything (c)
params <- list(support=0.0001, confidence=0.7,
               minlen=2, maxlen=4)

app <- list(rhs=c("revLenDesc=v.long",
                  "revLenDesc=long",
                  "revLenDesc=medium",
                  "revLenDesc=short",
                  "revLenDesc=v.short"), default="lhs")

fit <- apriori(x, parameter=params, appearance=app)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# remove rhs restrictions (d)
params <- list(support=0.01, confidence=0.7, minlen=2,
               maxlen=4)

fit <- apriori(x, parameter=params)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

# remove rhs restrictions again (e)
params <- list(support=0.0075, confidence=0.75,
               minlen=2, maxlen=4)

fit <- apriori(x, parameter=params)
qual <- quality(fit)
inspect(sort(fit, by = "lift"))

fit2 <- fit
quality(fit2) <- std_lift(fit2, x)
inspect(sort(fit2, by="slift"))

