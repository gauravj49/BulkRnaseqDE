# ****************************************************************
# USAGE: Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/csf_all_samples_diagnosis_LS100K_age40_pirnaome_normalized.txt output/test/pairwise_heatmap/R_csf_all_pirnaome.png 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

#Rscript scripts/R_pairwise_correlation_heatmap.R input/modelingIN/01_csf_filter/01_normalized_expFiltered/csf_AD_samples_diagnosis_LS100K_age40_pirnaome_normalized.txt output/test/pairwise_heatmap/R_csf_AD_pirnaome.png 11 2 initials p "initals,age,gender,detailed_diagnosis,library_size,diagnosis"

# ****************************************************************

################ load libraries ###############
suppressPackageStartupMessages(library("argparse"))

## Main logic
main <- function(){    

	##### Parse command line arguments ############
	args        <- check_options()
	inputfile   <- args$inputfile
	outputfile  <- args$outputfile
	xc          <- args$xc
	yc          <- args$yc
	sample_pc   <- args$sample_pc
	read_cutoff <- args$read_cutoff
	nf          <- args$nf
	nclust      <- args$nclust
	vmin        <- args$vmin
	vmax        <- args$vmax

	# Create output dir
	system(paste("mkdir -p ", dirname(outputfile), sep=''))

	if(args$scale == 0){
		scale <- "none"
	}else if(args$scale == 1){
		scale <- "row"
	}else if(args$scale == 2){
		scale <- "column"
	}
	if(length(args$annCols)){
		annCols <- as.vector(unlist(strsplit(args$annCols, ','))) # "initals,age,gender,detailed_diagnosis,library_size,diagnosis"
	}else{
		annCols <- "NONE"
	}
	if(length(args$vmin)){
		vmin    <- args$vmin
	}else{
		vmin    <- NULL
	}
	if(length(args$vmax)){
		vmax    <- args$vmax
	}else{
		vmax    <- NULL
	}

	################ load libraries ###############
	load_libraries()

	## Get the input data
	cat("- Reading input file ...\n")
	allDataOrig     <- data.frame(read.table(inputfile, header=TRUE, sep="\t", na.strings=c(NA, NaN, Inf), row.names=1))

	####### ignore NA rows
	allDataOrig[is.na(allDataOrig)] <- 0
	
	if (!nf){
		# Remove all the rows which 95% of it's entries are less than certain threshold
		allData <- allDataOrig[rowSums(allDataOrig >= read_cutoff) >= sample_pc*ncol(allDataOrig),] 
	} else{
		allData <- allDataOrig
	}

	if(annCols != "NONE"){
		anndf   <- allData[annCols]

		# # START
		# # Comment out this entire block from START to END if not planning to plot clusters
		# library(RColorBrewer)
		# aColors <- brewer.pal(n = length(unique(anndf$cluster)), name = "Set3")
		# # aColors        <- colorRampPalette(grDevices::rainbow(length(unique(anndf$cluster))))(length(unique(anndf$cluster)))
		# names(aColors) <- paste("cluster_",unique(anndf$cluster), sep='')
		# ofile    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_cluster_legend.pdf",sep='');
		# pdf(ofile)
		# pie(unique(anndf$cluster), col=aColors)
		# dev.off()
		# aColors        <- list(cluster = aColors)
		# # END

		# Delete annCols
		`%ni%`  <- Negate(`%in%`)
		allData <- subset(allData, select = names(allData) %ni% annCols)
	}else{
		anndf   <- data.frame()
		aColors <- NULL
	}
	# Get the number of colums
	ncols  <- length(names(allData))
	nrows  <- length(rownames(allData))
	cat("- Total samples : ", ncols,"\n- Total features: ", nrows, "\n")

	show_rownames = TRUE
	if(nrows > 300){
		show_rownames = FALSE
	}

	#	scaledData <- log2(allData+1)
	psize=10
	
	# Plot without any scaling
	ofile     <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_none_scaled.pdf" ,sep='');
	colorP    <-  colorRampPalette(c('white','gold','darkorange','red','darkred'))(256)
	if(annCols != "NONE"){
		ohnescRes <- pheatmap(allData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=file_path_sans_ext(basename(inputfile)), annotation_row = anndf,  annotation_colors = aColors, annotation_names_row = TRUE)
	}else{
		ohnescRes <- pheatmap(allData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=file_path_sans_ext(basename(inputfile)))
	}

	# Convert the pdf to png
	system(paste("convert -background white -alpha remove -density 300 ",ofile, " ", paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_expression.png" ,sep='')))

	ofile    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_values.pdf",sep='');
	colorP   <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(256)

	# If color scale is set, then plot the scaled data
	if(length(vmin) > 0 & length(vmax) > 0){
		# Sets the minimum (-2), the maximum (+2), and the increasing steps (+1) for the color scale
		# Note: if some of your genes are outside of this range, they will appear white on the heatmap
		breaksList = seq(vmin, vmax, by = 0.1)

		# Define the vector of colors for the legend (it has to be of the same lenght of breaksList)
		ofile    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_values_vmin",vmin,"_vmax",vmax,".pdf",sep='');
		pngof    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_values_vmin",vmin,"_vmax",vmax,".png",sep='');
		colorP   <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(length(breaksList))
		# colorP    <-  colorRampPalette(c('white','gold','darkorange','red','darkred'))(length(breaksList))
		if(annCols != "NONE"){
			# Sets the breaks of the color scale as in breaksList
			mitscRes <-pheatmap(allData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=paste(file_path_sans_ext(basename(inputfile))," (", scale, " scaled) ", sep=''), scale=scale, breaks = breaksList, annotation_row = anndf,  annotation_colors = aColors, annotation_names_row = TRUE) 
		}else{
			mitscRes <-pheatmap(allData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=paste(file_path_sans_ext(basename(inputfile))," (", scale, " scaled) ", sep=''), scale=scale, breaks = breaksList)
		}
	}else{
		# Automatic scale
		ofile    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_values.pdf",sep='');
		pngof    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_values.png",sep='');
		colorP   <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(256)
		if(annCols != "NONE"){
			cat("pheatmap(allclustlData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = \"grey50\", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=paste(file_path_sans_ext(basename(inputfile)),\" (\", scale, \" scaled) \", sep=''), scale=scale, annotation_row = anndf,  annotation_colors = aColors, annotation_names_row = TRUE)")

			mitscRes <-pheatmap(allclustlData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=paste(file_path_sans_ext(basename(inputfile))," (", scale, " scaled) ", sep=''), scale=scale, annotation_row = anndf,  annotation_colors = aColors, annotation_names_row = TRUE)
		}else{
			mitscRes <-pheatmap(allData, cluster_rows=xc, cluster_cols=yc, color = colorP, fontsize = psize/2, filename = ofile, border_color = "grey50", show_rownames=show_rownames, width=psize, height=psize, treeheight_row=psize*10, treeheight_col=psize*10, main=paste(file_path_sans_ext(basename(inputfile))," (", scale, " scaled) ", sep=''), scale=scale)
		}
	}

	# Convert the pdf to png
	system(paste("convert -background white -alpha remove -density 300 ",ofile, " ", pngof,sep=''))

	# Save the clusters from scaled data
	if (nclust > 0){
		ofile    <- paste(dirname(outputfile),"/",file_path_sans_ext(basename(outputfile)),"_", scale, "_scaled_clusters_n",nclust,".txt",sep='');
		allclust <- cbind(allData, cluster = cutree(mitscRes$tree_row, k = nclust))
		# # Sort the datafame on cluster values
		# allclust <- allclust[order(allclust$cluster),]
		write.table(allclust, file = ofile, col.names = NA, row.names = T, sep = '\t', quote = F)
	}

	# Turn off device driver (to flush output to PNG file)
	dev.off()

	# Print session info to the session's log file
	logdir <- paste(dirname(outputfile),"/logs", sep='')
	system(paste("mkdir -p ", logdir, sep=''))
	session_logfile <- paste(logdir,"/session_info_",file_path_sans_ext(basename(outputfile)),".log" ,sep='');
	print_session_info(session_logfile)
}

############ USER DEFINED FUNCTIONS ##########
# Print session info as log file formatted in tabular format
# Source: https://stackoverflow.com/questions/21967254/how-to-write-a-reader-friendly-sessioninfo-to-text-file
print_session_info <- function(session_logfile){
	suppressPackageStartupMessages(library("devtools"))
	suppressPackageStartupMessages(library("knitr"))

	# Get all the session info to the variable
	my_session_info <- devtools::session_info()

	# Print it in the tabular format using knitr
	writeLines(text = {
	    paste(sep = "\n", collapse = "",
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 32), collapse = ""),
	                "R environment",
	                paste0(rep("-", 33), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(data.frame(setting = names(my_session_info$platform),
	                                  value = as.character(my_session_info$platform))), collapse = "\n"),
	          paste0(rep("-", 80), collapse = ""),
	          paste(paste0(rep("-", 35), collapse = ""),
	                "packages",
	                paste0(rep("-", 35), collapse = "")),
	          paste0(rep("-", 80), collapse = ""),
	          paste(knitr::kable(my_session_info$packages), collapse = "\n")
	          )
	}, con = session_logfile)
}

# Load Libraries
load_libraries <- function(){
	# Load libraries at start
	suppressMessages( library(ggplot2))
	suppressMessages( library(gplots))
	suppressMessages( library(RColorBrewer))
	suppressMessages( library(pheatmap))
	suppressMessages( library(tools))        # for file path, basename_without_extension(file_path_sans_ext) and extensions(file_ext)
	suppressPackageStartupMessages(library("argparse"))
	#suppressPackageStartupMessages(library("BHC"))
	suppressPackageStartupMessages(library("dendsort"))
	suppressPackageStartupMessages(library("coop"))
	suppressPackageStartupMessages(library("effects"))
	suppressPackageStartupMessages(library("grid"))
	suppressPackageStartupMessages(library("scales"))
	suppressPackageStartupMessages(library("gtable"))
	suppressPackageStartupMessages(library("stats"))
	suppressPackageStartupMessages(library("grDevices"))
	suppressPackageStartupMessages(library("graphics"))
}

check_options <- function(){
	# Description of the script
	desc <- sprintf("
	----------------- SAMPLE USAGE ------------------
		- Rscript scripts/R_heatmap.R -if=output/filtered_data/deseq2/mirna/results_DEseq2/ACC18mwt_over_ACC4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/ACC18mwt_over_ACC4mwt_DE_heatmaps.txt -xc -yc
		- Rscript scripts/R_heatmap.R -if=output/filtered_data/deseq2/mirna/results_DEseq2/Blood18mwt_over_Blood4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/Blood18mwt_over_Blood4mwt_DE_heatmaps.txt
	-------------------------------------------------
	CONTACT: 
		Gaurav Jain
		gaurav.jain@dzne.de
	-------------------------------------------------\n
	")
	# create parser object
	parser <- ArgumentParser(description=cat(desc))

	# Add arguments 
	parser$add_argument("-if", "--inputfile"  , dest="inputfile"  , help="*Input matrix file", type="character", required=TRUE)
	parser$add_argument("-of", "--outputfile" , dest="outputfile" , help="*Output file name" , type="character", required=TRUE)
	parser$add_argument("-ac", "--annCols"    , dest="annCols"    , help=" Comma separated list of annotation columns to plot along with heatmap. Ex: the file columns may be - initals,age,gender,detailed_diagnosis,library_size,diagnosis etc. " , type="character")
	parser$add_argument('-xc', "--xc"         , dest="xc"         , help=" if set, cluster x-rows (Default=%(default)s)", action='store_true', default=FALSE)
	parser$add_argument('-yc', "--yc"         , dest="yc"         , help=" if set, cluster y-rows (Default=%(default)s)", action='store_true', default=FALSE)
	parser$add_argument('-sc', "--scale"      , dest="scale"      , help=" Number indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Please enter: 0 for NONE, 1 for ROW or 2 for COLUMN (Default=%(default)s)", default=0)
    parser$add_argument('-pc', "--samper"     , dest="sample_pc"  , help=" Percentage of sample for read cutoff (Default=%(default)s)", default=95)
    parser$add_argument('-ct', "--rdcutf"     , dest="read_cutoff", help=" Reads (normalized/unnormalized) below which to discard (Default=%(default)s)", default=5, type="integer")
	parser$add_argument('-nf', "--nf"         , dest="nf"         , help=" if set, do not filter reads (Default=%(default)s)", action='store_true', default=FALSE)
	parser$add_argument('-nc', "--nclust"     , dest="nclust"     , help=" Extract number of clusters.ex: 10 clusters (Default=%(default)s)", default=0, type="integer")
	parser$add_argument('-mn', "--vmin"       , dest="vmin"       , help=" Set the min value of the heatmap colorbar. Ex: -2", type="double")
	parser$add_argument('-mx', "--vmax"       , dest="vmax"       , help=" Set the max value of the heatmap colorbar. EX: +2", type="double")

	# Print the help message only if no arguments are supplied
	if(length(commandArgs(TRUE))==0){
		cat(desc)
		parser$print_help()
		quit()
	}

	# get command line options, if help option encountered print help and exit,
	# otherwise if options not found on command line then set defaults,
	args <- parser$parse_args()
	return(args)
}

## Call the main function in the end
main()

