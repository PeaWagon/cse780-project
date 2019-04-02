
# to run: Rscript hclust_analysis.R

# to compare the clustering results
library(fpc)

# clustering of fingerprint data

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
    # uncomment to produce the plots
    #for (linkage in linkages) {
        #hc <- hclust(d, linkage)
        #outfile <- paste(fingerprints[fpnum],
        #         "_hclust_", linkage, ".png",
        #         sep="")
        #png(outfile)
        #plot(hc, sub="", xlab=paste(linkage,
        #     fingerprints[fpnum]), ylab="",
        #     cex.axis=1.5, cex.lab=3, cex.main=3)
        #dev.off()
        # now make some cuts to compare the sizes
    #}

    # increment fingerprint name
    fpnum <- fpnum + 1
}

