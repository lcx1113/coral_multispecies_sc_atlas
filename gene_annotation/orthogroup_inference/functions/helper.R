# libraries
require("ComplexHeatmap")
require("circlize")
require("zoo")

# create a dictionary from two vectors
dic_from_vecs = function(names, terms, make_unique_on_names = TRUE) {
	dic = terms
	names(dic) = names
	if (make_unique_on_names) {
		dic = dic [ !duplicated(names(dic)) ]
	}
	return(dic)
}

# create a transcript to gene dictionary from a GTF annotation file
dictionary_t2g = function(gtf_fn, vector_to_fix, t2g = TRUE, transcript_field = "transcript", transcript_id = "transcript_id", gene_id = "gene_id", return_elements_not_in_gtf = TRUE) {
	
	# import gtf
	gene_txgtf = rtracklayer::import(gtf_fn)
	
	if (t2g) {
		dic = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,gene_id]
		names(dic) = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,transcript_id]
	} else {
		dic = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,transcript_id]
		names(dic) = as.data.frame(GenomicRanges::mcols(gene_txgtf[gene_txgtf$type == transcript_field]))[,gene_id]
	}
	
	# return object
	
	out = dic [ vector_to_fix ]
	
	# return elements not in GTF dictionary, unaltered
	if (return_elements_not_in_gtf) {
		ixs_to_keep = is.na(out)
		out[ixs_to_keep] = vector_to_fix[ixs_to_keep]
		names(out[ixs_to_keep]) = out[ixs_to_keep]
	}
	
	# return
	out
	
}

# Heatmaps
plot_complex_heatmap = function(
		mat,
		name = "heatmap",
		color_mat = c("white","#d6e72e","#6fb600","#003f4d"),
		color_min = 0,
		color_max = 1,
		fontsize = 10,
		hm_body_height = NULL,
		hm_body_width = NULL,
		categories_col = NULL,
		categories_row = NULL,
		separate_col = FALSE,
		separate_row = FALSE,
		colors_col = NULL,
		colors_row = NULL,
		title_row = NULL,
		title_col = NULL,
		name_row_show = TRUE,
		name_col_show = TRUE,
		size_row_annot = NULL,
		size_col_annot = NULL,
		cluster_row = TRUE,
		cluster_col = TRUE,
		use_raster = TRUE,
		raster_quality = 1,
		show_legend_row = FALSE,
		show_legend_col = FALSE,
		show_legend_hm = TRUE,
		both_sides_row = TRUE,
		both_sides_col = TRUE,
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = NA, lwd = 1, lty = 1),
		do_dotplot = FALSE,
		dot_size_mat = NULL,
		dot_size_min = NULL,
		dot_size_max = NULL,
		cex_dotplot = 0.02,
		col_dotplot_border = NA,
		label_mat = NULL,
		show_cell_label = FALSE,
		label_cex = 0.6,
		label_col = "black"
) {
	
	require("ComplexHeatmap")
	require("circlize")
	ht_opt$message = FALSE
	
	# color function
	col_fun = circlize::colorRamp2(seq(color_min, color_max, length.out = length(color_mat)), color_mat)
	
	
	# # vector of clusters (row/col annotation)
	# cluster_vector = as.character(mot_merge_d$cluster)
	# # categorical colors for clusters (get recycled)
	# catcol_vec = rep(catcol_lis, length.out = length(unique(cluster_vector)))
	# names(catcol_vec) = unique(cluster_vector)
	# cluster_colors = catcol_vec [ cluster_vector ]
	# names(cluster_colors) = cluster_colors
	# # quantitative colorscale for heatmap
	# catcol_fun = circlize::colorRamp2(1:length(catcol_vec), catcol_vec)
	# names(cluster_colors) = cluster_vector
	
	if (is.null(title_row)) {
		title_row = sprintf("n = %i rows", nrow(mat))
	}
	if (is.null(title_col)) {
		title_col = sprintf("n = %i columns", ncol(mat))
	}
	
	# left row annotations
	if (is.null(categories_row)) {
		categories_row = rownames(mat)
	}
	if (is.null(colors_row)) {
		colors_row = rep(NA, nrow(mat))
		right_annotation = NULL
	} else {
		if (is.null(names(colors_row))) { names(colors_row) = categories_row }
		if (is.null(size_row_annot)) { size_row_annot = ht_opt$simple_anno_size }
		right_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_row, name = title_row, col = list(c = colors_row), which = "row", show_legend = show_legend_row, simple_anno_size = size_row_annot, show_annotation_name = FALSE)
	}
	
	# top col annotations
	if (is.null(categories_col)) {
		categories_col = colnames(mat)
	}
	if (is.null(colors_col)) {
		colors_col = rep(NA, ncol(mat))
		top_annotation = NULL
	} else {
		if (is.null(names(colors_col))) { names(colors_col) = categories_col }
		if (is.null(size_col_annot)) { size_col_annot = ht_opt$simple_anno_size }
		top_annotation = ComplexHeatmap::HeatmapAnnotation(c = categories_col, name = title_col, col = list(c = colors_col), which = "column", show_legend = show_legend_col, simple_anno_size = size_col_annot, show_annotation_name = FALSE)
	}
	
	# add annotations to both sides of the rows or columns?
	if (both_sides_col) {
		bottom_annotation = top_annotation
	} else {
		bottom_annotation = NULL
	}
	if (both_sides_row) {
		left_annotation = right_annotation
	} else {
		left_annotation = NULL
	}
	
	# split columns and rows?
	if (separate_row) {
		split_row = categories_row
	} else {
		split_row = NULL
	}
	if (separate_col) {
		split_col = categories_col
	} else {
		split_col = NULL
	}
	
	cell_fun_use <- NULL
	# should this be a dot plot?
	if (do_dotplot || show_cell_label) {
		if (do_dotplot) {
			if (is.null(dot_size_mat)) {
				dot_size_mat = mat
			}
			if (!is.null(dot_size_max) & !is.null(dot_size_min)) {
				# dot_size_mat [ dot_size_mat > dot_size_max ] = dot_size_max
				# dot_size_min [ dot_size_mat < dot_size_min ] = dot_size_min
				rangesize = function(x) { (x - dot_size_min) / (dot_size_max - dot_size_min) }
			} else {
				rangesize = function(x) { (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) }
			}
		}
		# label_mat 默认由 mat 转字符（可配 label_format）
		if (show_cell_label) {
			if (is.null(label_mat)) {
			if (is.null(label_format)) {
				label_mat <- matrix(as.character(round(mat, 2)), nrow = nrow(mat), dimnames = dimnames(mat))
			} else {
				label_mat <- matrix(label_format(mat), nrow = nrow(mat), dimnames = dimnames(mat))
			}
			} else {
			  stopifnot(all(dim(label_mat) == dim(mat)))
			}
		}
		cell_fun_use <- function(j, i, x, y, width, height, fill) {
			# dotplot
			if (do_dotplot) {
			grid.circle(
				x = x, y = y,
				r = sqrt(rangesize(dot_size_mat)[i, j]) * cex_dotplot,
				gp = gpar(col = col_dotplot_border, fill = col_fun(mat[i, j]))
			)
			}

			# text label
			if (show_cell_label) {
			lab <- label_mat[i, j]
			if (!is.na(lab) && lab != "") {
				grid.text(
				label = lab, x = x, y = y,
				gp = gpar(col = label_col, fontsize = fontsize * label_cex)
				)
			}
			}
		}

		if (do_dotplot) cell_border = gpar(type = "none")
}
	
	# plot
	hm = ComplexHeatmap::Heatmap(
		mat,
		name = name,
		cell_fun = cell_fun_use,
		height = hm_body_height,
		width = hm_body_width,
		use_raster = use_raster,
		raster_quality = raster_quality,
		cluster_rows = cluster_row,
		cluster_columns = cluster_col,
		row_title = title_row,
		row_title_gp = gpar(fontsize = fontsize),
		column_title = title_col,
		column_title_gp = gpar(fontsize = fontsize),
		show_row_names = name_row_show,
		show_column_names = name_col_show,
		row_names_gp = gpar(fontsize = fontsize),
		column_names_gp = gpar(fontsize = fontsize),
		top_annotation = top_annotation,
		left_annotation = left_annotation,
		right_annotation = right_annotation,
		bottom_annotation = bottom_annotation,
		column_split = split_col,
		row_split = split_row,
		row_gap = unit(0.5, "mm"),
		column_gap = unit(0.5, "mm"),
		rect_gp = cell_border,
		border_gp = heatmap_border,
		show_heatmap_legend = show_legend_hm,
		col = col_fun)
	
	# return heatmap
	return(hm)
	
}




# inflection point of a knee plot
inflection_point <- function(counts_per_cell, min_value = 100) {
	
	counts_per_cell = counts_per_cell [!is.na(counts_per_cell)]
	# get log10 data
	dat = data.frame(
		umis = log10(sort(counts_per_cell, decreasing = TRUE)), 
		rank = log10(1:length(counts_per_cell)))
	
	# ignore cells below min value
	if (!is.null(min_value) & min_value != 0) {
		dat = dat[dat$umis > log10(min_value), ]
	}
	
	# rank difference
	dif_r = diff(dat$umis) / diff(dat$rank)
	min_ix = which.min(dif_r)
	
	# calculate inflection point
	inf = 10 ^ (dat$umis[min_ix])
	return(inf)
	
}


# valleys in a series of data
find_peaks = function (x, m = 10){
	shape <- diff(sign(diff(x, na.pad = FALSE)))
	pks <- sapply(which(shape < 0), FUN = function(i){
		z <- i - m + 1
		z <- ifelse(z > 0, z, 1)
		w <- i + m + 1
		w <- ifelse(w < length(x), w, length(x))
		if (all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
	})
	pks <- unlist(pks)
	pks
}

# wrapper to return peaks or valleys from an histogram
find_peaks_histogram = function(histogram, m = 10, valleys = FALSE) {
	
	x = histogram$counts
	
	if (valleys) {
		x = -x 
	}
	i = find_peaks(x, m = m)
	c = histogram$breaks [ i ]
	return(c)
	
}

# cell name cleaning function
# takes as input a vector of cellranger-style cell names and returns just the barcodes
cell_name_clean = function(v) {
	apply(stringr::str_split(v, "_", simplify = TRUE), 1, function(i) {
		i = c(i,"")
		ix = which(i == "") - 1
		c = i[ix][1]
		c = gsub("-\\d+$", "", c)
	}
	)
}


clean_og_pairs = function(og_pairs_fn, sp1, sp2, t2g = TRUE, t2g_sp1 = NULL, t2g_sp2 = NULL, t2g_sp1_field = "gene_id", t2g_sp2_field = "gene_id", header = FALSE) {
	
	# load
	if (is.character(og_pairs_fn)) {
		og_pairs = read.table(og_pairs_fn, header = header, sep = "\t", col.names = c("g1", "g2"), stringsAsFactors = FALSE)
	} else {
		og_pairs = og_pairs_fn
		colnames(og_pairs) = c("g1", "g2")
	}
	og_pairs = og_pairs [ (grepl(sprintf("^%s_",sp1),og_pairs[,1]) & grepl(sprintf("^%s_",sp2),og_pairs[,2])) | (grepl(sprintf("^%s_",sp1),og_pairs[,2]) & grepl(sprintf("^%s_",sp2),og_pairs[,1])) , ]
	og_pairs = data.frame(
		sp1 = c(og_pairs [ grepl(sprintf("^%s_",sp1),og_pairs[,1]), 1 ] , og_pairs [ grepl(sprintf("^%s_",sp1),og_pairs[,2]), 2 ] ),
		sp2 = c(og_pairs [ grepl(sprintf("^%s_",sp2),og_pairs[,2]), 2 ] , og_pairs [ grepl(sprintf("^%s_",sp2),og_pairs[,1]), 1 ] )
	)
	
	# load gene to transcript dictionary
	if (!is.null(t2g_sp1)) {
		og_pairs$sp1 = dictionary_t2g(gtf_fn = t2g_sp1, vector_to_fix = og_pairs$sp1, gene_id = t2g_sp1_field)
	}
	if (!is.null(t2g_sp2)) {
		og_pairs$sp2 = dictionary_t2g(gtf_fn = t2g_sp2, vector_to_fix = og_pairs$sp2, gene_id = t2g_sp2_field)
	}
	og_pairs = og_pairs [ !is.na(og_pairs$sp1) & !is.na(og_pairs$sp2) , ]
	
	return(og_pairs)
	
}


# subsample matrix
matrix_subsample = function(mat, subsample_rate = 0.1, subsample_rate_row = NULL, subsample_rate_col = NULL) {
	
	# subsample rates
	if (is.null(subsample_rate_row)) { subsample_rate_row = subsample_rate }
	if (is.null(subsample_rate_col)) { subsample_rate_col = subsample_rate }
	
	# num bins and widths
	nbin_r = nrow(mat) * subsample_rate_row
	nbin_c = ncol(mat) * subsample_rate_col
	wbin_r = round(nrow(mat) / nbin_r)
	wbin_c = round(ncol(mat) / nbin_c)
	
	# init subsampled matrix
	smat = matrix(nrow = nbin_r, ncol = nbin_c)
	
	# log
	message(sprintf("subsample matrix at %.2f x %.2f rate: %i x %i to %i x %i", subsample_rate_row, subsample_rate_col, nrow(mat), ncol(mat), nrow(smat), ncol(smat)))
	
	# loop
	for (rs in 0:(nbin_r - 1)) {
		ro_s = (rs * wbin_r) + 1
		ro_e = rs * wbin_r + wbin_r
		for (cs in 0:(nbin_c - 1)) {
			co_s = (cs * wbin_c) + 1
			co_e = cs * wbin_c + wbin_c
			smat_i = mat [ ro_s : ro_e , co_s : co_e ]
			smat[rs + 1, cs + 1]   = mean(smat_i)
		}
	}
	
	# out
	return(smat)
	
}

# Helper function to round values without changing total sum
round_smart <- function(x, digits = 0) {
	up <- 10 ^ digits
	x <- x * up
	y <- floor(x)
	indices <- tail(order(x-y), round(sum(x)) - sum(y))
	y[indices] <- y[indices] + 1
	y / up
}

# jaccard from vectors
jaccard_index = function(v1, v2) {
	
	l_int = length(intersect(v1, v2))
	l_uni = length(union(v1, v2))
	return(l_int / l_uni)
	
}


# function
table_to_matrix = function(table) {
  mat = matrix(table, nrow = nrow(table))
  rownames(mat) = rownames(table)
  colnames(mat) = colnames(table)
  return(mat)
}

# modified from RaceID
raceid_plotgraph_mod = function (object, showCells = FALSE, showMap = TRUE, tp = 0.5, scthr = 0, cex = 1) 
{
    if (length(object@cdata) <= 0) 
        stop("run comppvalue before plotgraph")
    if (!is.numeric(showCells) & !is.logical(showCells)) 
        stop("argument showCells has to be logical (TRUE/FALSE)")
    if (!is.numeric(scthr)) 
        stop("scthr has to be a non-negative number")
    else if (scthr < 0 | scthr > 1) 
        stop("scthr has to be a number between 0 and 1")
    ramp <- colorRamp(c("#d4d7b9","#d6e72e","#6fb600","#003f4d"))
    mcol <- rgb(ramp(seq(0, 1, length = 101)), maxColorValue = 255)
    co <- object@cdata
    fc <- (co$counts/(co$counts.br + 0.5)) * (co$pv.e < object@par$pthr)
    fc <- fc * (fc > t(fc)) + t(fc) * (t(fc) >= fc)
    fc <- log2(fc + (fc == 0))
    if (object@par$nmode | object@par$fast) {
        k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e)) < 
            object@par$pthr])) + 1e-05)
        if (length(k) == 1) 
            k <- c(k - k/100, k)
        mlpv <- -log10(co$pvn.e + 1e-05)
    }
    else {
        k <- -log10(sort(unique(as.vector(t(co$pvn.e))[as.vector(t(co$pv.e)) < 
            object@par$pthr])) + 1/object@par$pdishuf)
        if (length(k) == 1) 
            k <- c(k - k/100, k)
        mlpv <- -log10(co$pvn.e + 1/object@par$pdishuf)
    }
    diag(mlpv) <- min(mlpv, na.rm = TRUE)
    dcc <- t(apply(round(100 * (mlpv - min(k))/(max(k) - min(k)), 
        0) + 1, 1, function(x) {
        y <- c()
        for (n in x) y <- append(y, if (n < 1) 
            NA
        else mcol[n])
        y
    }))
    cx <- c()
    cy <- c()
    va <- c()
    m <- object@ldata$m
    for (i in 1:(length(m) - 1)) {
        for (j in (i + 1):length(m)) {
            if (min(co$pv.e[i, j], co$pv.e[j, i], na.rm = TRUE) < 
                object@par$pthr) {
                if (mlpv[i, j] > mlpv[j, i]) {
                  va <- append(va, dcc[i, j])
                }
                else {
                  va <- append(va, dcc[j, i])
                }
                cx <- append(cx, i)
                cy <- append(cy, j)
            }
        }
    }
    cnl <- object@ldata$cnl
    u <- object@ltcoord[, 1]
    v <- object@ltcoord[, 2]
    pardefault <- par()
    layout(cbind(c(1, 1), c(2, 3)), widths = c(5, 1, 1), heights = c(5, 
        5, 1))
    par(mar = c(10, 5, 1, 1))
    xlim <- c(min(u), max(u))
    ylim <- c(min(v), max(v))
    if (showMap) {
        part <- object@ldata$lp
        f <- object@sc@cpart %in% unique(part)
        if (object@par$fr) {
            d <- object@sc@fr[f, ]
        }
        else if (object@par$um) {
            d <- object@sc@umap[f, ]
        }
        else {
            d <- object@sc@tsne[f, ]
        }
        xlim <- c(min(u, d[, 1]), max(u, d[, 1]))
        ylim <- c(min(v, d[, 2]), max(v, d[, 2]))
    }
    plot(0, 0, cex = 0, xlim = xlim, ylim = ylim, xlab = "", 
        ylab = "", axes = FALSE)
    if (showMap) {
        for (i in 1:max(part)) {
            if (sum(part == i) > 0) 
                points(d[part == i, 1], d[part == i, 2], col = adjustcolor(object@sc@fcol[i], 
                  tp), pch = 20, cex = cex)
        }
    }
    if (showCells) 
        points(u, v, cex = 1.5, col = "grey", pch = 20)
    if (length(va) > 0) {
        f <- order(va, decreasing = TRUE)
        for (i in 1:length(va)) {
            if (object@cdata$linkscore[cx[i], cy[i]] > scthr) {
                if (showCells) {
                  lines(cnl[c(cx[i], cy[i]), 1], cnl[c(cx[i], 
                    cy[i]), 2], col = va[i], lwd = 2)
                }
                else {
                  lines(cnl[c(cx[i], cy[i]), 1], cnl[c(cx[i], 
                    cy[i]), 2], col = va[i], lwd = 5 * object@cdata$linkscore[cx[i], 
                    cy[i]])
                }
            }
        }
    }
    en <- aggregate(object@entropy, list(object@sc@cpart), median)
    en <- en[en$Group.1 %in% m, ]
    mi <- min(en[, 2], na.rm = TRUE)
    ma <- max(en[, 2], na.rm = TRUE)
    w <- round((en[, 2] - mi)/(ma - mi) * 99 + 1, 0)
    ramp <- colorRamp(c("#bbccd1", "#3ab2da","dodgerblue3","midnightblue"))
    ColorRamp <- rgb(ramp(seq(0, 1, length = 101)), maxColorValue = 255)
    ColorLevels <- seq(mi, ma, length = length(ColorRamp))
    if (mi == ma) {
        ColorLevels <- seq(0.99 * mi, 1.01 * ma, length = length(ColorRamp))
    }
    for (i in m) {
        f <- en[, 1] == m
        points(cnl[f, 1], cnl[f, 2], cex = 5, col = ColorRamp[w[f]], 
            pch = 20)
    }
    text(cnl[, 1], cnl[, 2], m, cex = 1.25, font = 2, col = "white")
    par(mar = c(5, 4, 1, 2))
    image(1, ColorLevels, matrix(data = ColorLevels, ncol = length(ColorLevels), nrow = 1), col = ColorRamp, xlab = "", ylab = "entropy", xaxt = "n")
    coll <- seq(min(k), max(k), length = length(mcol))
    image(1, coll, matrix(data = coll, ncol = length(mcol), nrow = 1), col = mcol, xlab = "", ylab = "", xaxt = "n")
    layout(1)
    par(mar = pardefault$mar)
}



# find PCA elbow function from findPC: Automatic selection of number of principal components
# Citation: https://doi.org/10.1093/bioinformatics/btac235
find_pca_elbow = function(sdev){

	x<-1:length(sdev)
	eb<-data.frame(x,sdev)

	# piecewise linear model
	sse<-NULL
	for (k in 1:length(sdev)) {
	D=ifelse(x<=k,0,1)
	x2<-(x-k)*D
	sse[k]=sum(lm(sdev~x+x2)$residuals^2)
	}
	dim_plm<-which.min(sse)
	names(dim_plm)<-'Piecewise linear model'

	# first derivative
	df1<-diff(sdev,differences=1)
	den<-density(df1)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_fid<-max(which(df1<mean(df1)))+1
	} else {
	cutoff<-sum(rlev[-((length(rlev)-1):length(rlev))])
	dim_fid<-max(which(df1<den$x[cutoff]))+1 }
	names(dim_fid)<-'First derivative'

	# second derivative
	df2<-diff(sdev,differences=2)
	df2p<-df2[df2>0]
	den<-density(df2p)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_sed<-which.max(df2)+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_sed<-max(which(df2>den$x[cutoff]))+1 }
	names(dim_sed)<-'Second derivative'

	# preceding residual
	fit<-NULL;res<-NULL
	for (i in 1:(length(sdev)-2)) {
	fit[[i]]<-lm(sdev~x,data=eb[(i+1):length(sdev),])
	res[i]<-sdev[i]-predict(fit[[i]],newdata=data.frame(x=i))
	}
	den<-density(res)
	rlev<-rle(diff(den$y)>0)$lengths
	if(length(rlev)<=2) {dim_pr<-max(which(res>mean(res)))+1
	} else {
	cutoff<-sum(rlev[1:2])
	dim_pr<-max(which(res>den$x[cutoff]))+1 }
	names(dim_pr)<-'Preceding residual'

	# perpendicular line
	A<-c(1,sdev[1]);B<-c(length(sdev),sdev[length(sdev)]);Dist<-NULL
	for (i in 2:(length(sdev)-1)) {
	C<-c(i,sdev[i]);D<-cbind(rbind(A,B,C),rep(1,3))
	S<-1/2*abs(det(D));Dist[i]<-2*S/dist(rbind(A,B))
	}
	dim_perl<-which.max(Dist)
	names(dim_perl)<-'Perpendicular line'

	# k-means clustering
	set.seed(2022)
	dim_clu<-min(kmeans(sdev,2)$size)+1
	names(dim_clu)<-'K-means clustering'

	dim_all<-c(dim_plm,dim_fid,dim_sed,dim_pr,dim_perl,dim_clu)
	return(dim_all)
}

select_top_markers = function(matrix, matrix_thr = 0, n_top_markers = 20, n_markers_rollmean = 2) {
	
	markers = unique(as.vector(unlist( apply(matrix, 2, function(c) { names( head(sort(-c[ c >= matrix_thr ]), n = n_top_markers ) ) } ))))
	markers_order = order(apply(matrix[markers,], 1, function(r) which.max(zoo::rollmean(r, n_markers_rollmean) )))
	markers_ordered = markers [ markers_order ]
	return(markers_ordered)
	
}

plot_complex_dotmap <- function(
  mat,
  dot_size_mat = NULL,
  status_mat = NULL,      # U/O/B
  pval_mat = NULL,
  p_cutoff = 0.05,

  name = "expr",
  color_mat = c("navy", "white", "firebrick3"),
  color_min = -2,
  color_max = 2,

  cex_dot = 0.45,
  dot_size_min = NULL,
  dot_size_max = NULL,

  cluster_row = FALSE,
  cluster_col = FALSE,
  show_row_names = TRUE,
  show_col_names = TRUE,
  fontsize = 8,
  na_col = "#F0F0F0",

  status_border_col = c(U = "black", O = "#2C7FB8", B = "#7A0177"),
  status_border_lwd = c(U = 0.8, O = 0.8, B = 1.8),
  status_border_lty = c(U = "solid", O = "dashed", B = "solid"),

  sig_ring_col = "#1a9850",
  sig_ring_lwd = 1.1,
  sig_ring_expand = 1.18,

  row_split = NULL,
  column_split = NULL,
  top_annotation = NULL,
  left_annotation = NULL,
  right_annotation = NULL,
  bottom_annotation = NULL
) {
  require(ComplexHeatmap)
  require(circlize)
  require(grid)

  stopifnot(is.matrix(mat), is.numeric(mat))
  nr <- nrow(mat); nc <- ncol(mat)

  if (is.null(dot_size_mat)) dot_size_mat <- abs(mat)
  stopifnot(is.matrix(dot_size_mat), all(dim(dot_size_mat) == c(nr, nc)))
  if (!is.null(status_mat)) stopifnot(is.matrix(status_mat), all(dim(status_mat) == c(nr, nc)))
  if (!is.null(pval_mat))   stopifnot(is.matrix(pval_mat),   all(dim(pval_mat) == c(nr, nc)))

  col_fun <- circlize::colorRamp2(
    seq(color_min, color_max, length.out = length(color_mat)),
    color_mat
  )
  draw_circle_border <- function(x, y, r, col = "black", lwd = 1, lty = "solid", n = 120) {
	theta <- seq(0, 2*pi, length.out = n)
	xx <- x + unit(r * cos(theta), "snpc")
	yy <- y + unit(r * sin(theta), "snpc")
	grid::grid.polyline(
		x = xx, y = yy,
		gp = gpar(col = col, lwd = lwd, lty = lty, fill = NA)
	)
  }

  # size scaling
  if (!is.null(dot_size_min) && !is.null(dot_size_max)) {
    if (dot_size_max <= dot_size_min) stop("dot_size_max must be > dot_size_min")
    size_scaled <- (dot_size_mat - dot_size_min) / (dot_size_max - dot_size_min)
  } else {
    rng <- range(dot_size_mat, na.rm = TRUE)
    if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
      size_scaled <- matrix(0.7, nrow = nr, ncol = nc, dimnames = dimnames(dot_size_mat))
    } else {
      size_scaled <- (dot_size_mat - rng[1]) / diff(rng)
    }
  }
  size_scaled[size_scaled < 0] <- 0
  size_scaled[size_scaled > 1] <- 1

cell_fun <- function(j, i, x, y, width, height, fill) {
  v <- mat[i, j]
  s <- size_scaled[i, j]
  if (is.na(v) || is.na(s)) return(NULL)

  r <- sqrt(s) * cex_dot
  st <- if (!is.null(status_mat)) as.character(status_mat[i, j]) else NA_character_

  bcol <- if (!is.na(st) && st %in% names(status_border_col)) status_border_col[[st]] else NA
  blwd <- if (!is.na(st) && st %in% names(status_border_lwd)) status_border_lwd[[st]] else 1
  blty <- if (!is.na(st) && st %in% names(status_border_lty)) status_border_lty[[st]] else "solid"

  # 先画填充圆（无边框）
  grid::grid.circle(
    x = x, y = y, r = unit(r, "snpc"),
    gp = gpar(fill = col_fun(v), col = NA)
  )

  # 再单独画边框（支持虚线）
  if (!is.na(bcol)) {
    draw_circle_border(x, y, r = r, col = bcol, lwd = blwd, lty = blty)
  }

  # 显著性外圈（保持你原来的逻辑）
  if (!is.null(pval_mat)) {
    pv <- pval_mat[i, j]
    if (!is.na(pv) && pv < p_cutoff) {
      grid::grid.circle(
        x = x, y = y, r = unit(r * sig_ring_expand, "snpc"),
        gp = gpar(fill = NA, col = sig_ring_col, lwd = sig_ring_lwd)
      )
    }
  }
}

  ht <- Heatmap(
    mat,
    name = name,
    col = col_fun,
    na_col = na_col,
    rect_gp = gpar(type = "none"),
    cell_fun = cell_fun,
    cluster_rows = cluster_row,
    cluster_columns = cluster_col,
    show_row_names = show_row_names,
    show_column_names = show_col_names,
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    row_split = row_split,
    column_split = column_split,
    top_annotation = top_annotation,
    left_annotation = left_annotation,
    right_annotation = right_annotation,
    bottom_annotation = bottom_annotation
  )

  # ---------- custom legends ----------
 # 1) Dot size legend: 与 cell_fun 完全一致
# 允许 dot_size_min/max 为空时自动使用实际范围
if (is.null(dot_size_min) || is.null(dot_size_max)) {
  rng_legend <- range(dot_size_mat, na.rm = TRUE)
  ds_min_lgd <- rng_legend[1]
  ds_max_lgd <- rng_legend[2]
} else {
  ds_min_lgd <- dot_size_min
  ds_max_lgd <- dot_size_max
}

# 你可按需要调整断点（保证在范围内）
size_break_vals <- c(0, 2, 4, 6)
size_break_vals <- size_break_vals[size_break_vals >= ds_min_lgd & size_break_vals <= ds_max_lgd]
if (length(size_break_vals) == 0) {
  size_break_vals <- signif(seq(ds_min_lgd, ds_max_lgd, length.out = 4), 2)
}

size_break_scaled <- (size_break_vals - ds_min_lgd) / (ds_max_lgd - ds_min_lgd)
size_break_scaled <- pmax(0, pmin(1, size_break_scaled))
size_break_r <- sqrt(size_break_scaled) * cex_dot * 8

lgd_size <- ComplexHeatmap::Legend(
  title = "|log2FC|",
  labels = sprintf("%.1f", size_break_vals),
  type = "points",
  pch = 16,
  size = unit(size_break_r * 2, "snpc"),  # 与图中半径公式对应
  legend_gp = gpar(col = "grey25", fill = "grey25")
)

# 2) Status legend: 放大点 + 显示虚线/实线
# 用 graphics 自绘，保证 lty 生效
uob_labels <- intersect(c("U", "O", "B"), names(status_border_col))
status_legend_dot_r <- cex_dot * 8  # 可再调大，例如 2.0

make_status_grob <- function(st) {
  function(x, y, w, h) {
    r_mm <- min(convertWidth(w, "mm", valueOnly = TRUE),
                convertHeight(h, "mm", valueOnly = TRUE)) * status_legend_dot_r * 0.9

    # 填充圆
    grid.circle(
      x, y, r = unit(r_mm, "mm"),
      gp = gpar(fill = "white", col = NA)
    )

    # 边框圆（用 polyline，虚线稳定）
    theta <- seq(0, 2 * pi, length.out = 160)
    xx <- x + unit(r_mm * cos(theta), "mm")
    yy <- y + unit(r_mm * sin(theta), "mm")

    grid.polyline(
      x = xx, y = yy,
      gp = gpar(
        col = status_border_col[[st]],
        lwd = status_border_lwd[[st]],
        lty = status_border_lty[[st]]
      )
    )
  }
}

lgd_uob <- ComplexHeatmap::Legend(
  title = "Status",
  labels = uob_labels,
  graphics = lapply(uob_labels, make_status_grob)
)

# 3) significance ring legend（保留）
lgd_sig <- ComplexHeatmap::Legend(
  title = "Significance",
  labels = paste0("p < ", p_cutoff),
  graphics = list(function(x, y, w, h) {
    r0 <- min(convertWidth(w, "mm", valueOnly = TRUE),
              convertHeight(h, "mm", valueOnly = TRUE)) * 0.18
    grid.circle(x, y, r = unit(r0, "mm"),
                gp = gpar(fill = "white", col = "white", lwd = 0.8))
    grid.circle(x, y, r = unit(r0 * sig_ring_expand, "mm"),
                gp = gpar(fill = NA, col = sig_ring_col, lwd = sig_ring_lwd))
  })
)

lgd_list <- ComplexHeatmap::packLegend(lgd_size, lgd_sig, direction = "vertical")

  list(
    heatmap = ht,
    extra_legends = lgd_list
  )
}

collapse_unique <- function(x, sep = "/") {
  x <- unique(na.omit(as.character(x)))
  x <- x[x != ""]
  if (length(x) == 0) return(NA_character_)
  paste(x, collapse = sep)
}

merge_match_type <- function(x) {
  x <- unique(na.omit(as.character(x)))
  if ("B" %in% x) return("B")
  if (all(c("U", "O") %in% x)) return("B")
  if ("U" %in% x) return("U")
  if ("O" %in% x) return("O")
  NA_character_
}
match_plot_target <- function(main_gene, target_vec) {
  if (is.na(main_gene) || main_gene == "") return(NA_character_)
  hit <- target_vec[sapply(target_vec, function(pt) {
    parts <- strsplit(pt, "/", fixed = TRUE)[[1]]
    main_gene %in% parts
  })]
  if (length(hit) == 0) NA_character_ else hit[1]  # 命中多个时取第一个
}
fix_gene_name <- function(gene_names, prefix) {
    gene_names<-gsub("-","_", gene_names)
    ifelse(grepl(paste0("^", prefix), gene_names), gene_names, paste0(prefix, gene_names))
}