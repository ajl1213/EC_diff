

library(covNormRpkg)



args <- commandArgs(TRUE)
inFile=args[1]
sampleID=args[2]
chrID=args[3]
resolution=args[4]
maxDist=as.integer(args[5])


raw_data <- read.table(gzfile(inFile),head=TRUE)
print("1: Data Loaded.")

raw_data_filter <- covNormRpkg::filterInputDF(raw_data)
print("2: Data Filtered.")

cov_result <- covNormRpkg::normCoverage(raw_data_filter)
cov_result$coeff_cov1
cov_result$coeff_cov2
cov_df <- cov_result$result_df

outFile=paste0('CovNorm/',sampleID, '.', chrID, '.',resolution, '.CovNormFeat.gz')
write.table(cov_df, file=gzfile(outFile), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE) # Coverage normalization result
print("3: Coverage normalized.")

covNormRpkg::checkFreqCovPCC(cov_df, outpdfname=paste0('Figs/',sampleID, '.', chrID, '.', resolution,'.QCplot_coverage_PCC.pdf'))
covNormRpkg::plotCovNormRes( cov_df, outpdfname=paste0('Figs/',sampleID, '.', chrID, '.', resolution,'.QCplot_coverage_heatmap.pdf'))
print("4: Plot coverage normalization results.")

dist_result <- covNormRpkg::normDistance(cov_df, max_dist=maxDist)
dist_df <- dist_result$result_df
print("5: Distance normalized.")

covNormRpkg::checkFreqDistPCC(dist_df, outpdfname=paste0('Figs/',sampleID, '.', chrID, '.', resolution,'.QCplot_dist_PCC.pdf'))
covNormRpkg::plotDistNormRes( dist_df, outpdfname=paste0('Figs/',sampleID, '.', chrID, '.', resolution,'.QCplot_dist_hexmap.pdf'))
print("6: Plot distance normalization results.")

final_df <- covNormRpkg::contactPval(dist_df, paste0('Figs/',sampleID, '.', chrID, '.', resolution,'.fit.pdf'))
print("7: Significant interactions called.")

#Uncomment 'saveEachChr' to split-save file for each chromosome.
#covNormRpkg::saveEachChr(final_df, "./outputFolder", "outputSampleName") 
outFile=paste0('DistNorm/',sampleID, '.', chrID, '.',resolution, '.DistNormFeat.gz')
write.table(final_df, file=gzfile(outFile), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE) # Distance normalization & significant interactions







