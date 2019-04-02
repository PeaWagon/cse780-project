
# to run: Rscript kmedoids_analysis.R

library(cluster) # for using silhouette function

files <- list(
    "../data_processing/fp2_dissimilarity.csv",
    "../data_processing/fp3_dissimilarity.csv",
    "../data_processing/fp4_dissimilarity.csv",
    "../data_processing/maccs_dissimilarity.csv"
)

fingerprints <- list("fp2", "fp3", "fp4", "maccs")

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

