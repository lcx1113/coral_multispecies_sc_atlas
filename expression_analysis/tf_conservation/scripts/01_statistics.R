#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(data.table);library(Matrix);library(jsonlite);library(digest)})
setDTthreads(4)
args <- commandArgs(TRUE); root <- normalizePath(args[1]); spi <- args[2]
cfg <- fromJSON(file.path(root,'config/analysis_config.json'))
out <- file.path(root,'results/by_species',spi); dir.create(out,recursive=TRUE,showWarnings=FALSE)
write_tsv <- function(x,path) fwrite(x,path,sep='\t',na='NA',quote=FALSE,compress=if(grepl('gz$',path))'gzip' else 'none')
started <- Sys.time()
objfile <- file.path(cfg$source_root,'intermediate/02_og_sum_objects',paste0(spi,'.OG_SUM.seurat.rds'))
sourcehash <- digest(file=objfile,algo='sha256',serialize=FALSE)
obj <- readRDS(objfile)
counts <- LayerData(obj,assay='OG_SUM',layer='counts'); dat <- LayerData(obj,assay='OG_SUM',layer='data')
stopifnot(identical(dimnames(counts),dimnames(dat)),all(Matrix::rowSums(counts)>0))
linear <- dat; linear@x <- expm1(linear@x)
labels <- as.character(obj$cell_type)
manifest <- fread(file.path(root,'config/celltype_species_manifest.tsv'))[object_species==spi]
targets <- cfg$celltype_order[cfg$celltype_order %in% unique(labels)]
stopifnot(setequal(targets,manifest[include_analysis==TRUE,cell_type]),setequal(unique(labels),manifest$cell_type))
annotation <- fread(file.path(cfg$annotation_root,'orthogroup_final_annotations.csv'))
tfogs <- annotation[has_TF==TRUE,orthogroup]
brocfile <- '../../gene_annotation/orthogroup_inference/reference_results/Anthozoa/orthogroup_conservation.csv'
broc <- fread(brocfile,select=c('orthogroup','species'))
orthosp <- if(spi=='Spin')'Spis' else spi
present <- unique(broc[species==orthosp,orthogroup])
fmta <- fread(file.path(cfg$source_root,'intermediate/02_og_sum_objects',paste0(spi,'.OG_SUM.feature_metadata.tsv.gz')))
mapped_positive <- fmta$orthogroup
mapping <- fread(cfg$gene_og_mapping,select=c('object_species','orthogroup','feature_total_umi'))[object_species==spi]
mapped_any <- unique(mapping$orthogroup)
allstats <- list(); tfdetails <- list(); summaries <- list()
for(ct in targets){
  ii <- which(!is.na(labels)&labels==ct); jj <- setdiff(seq_along(labels),ii)
  cat(spi,ct,length(ii),'vs',length(jj),'\n');flush.console()
  if(ct %in% cfg$reuse_comparisons){
    src <- file.path(cfg$source_root,'results/comparisons',ct,'markers/full',paste0(spi,'.OG_SUM.all_tested.tsv.gz'))
    fm <- fread(src); stopifnot(all(fm$target_labels==ct),all(fm$n_target_cells==length(ii)),all(fm$n_background_cells==length(jj)))
    stats_origin <- 'reused_formal_OG_SUM_comparison'
  } else {
    Idents(obj) <- factor(ifelse(labels==ct,'target','background'),levels=c('background','target'))
    set.seed(cfg$markers$seed)
    fm <- as.data.table(FindMarkers(obj,ident.1='target',ident.2='background',assay='OG_SUM',
      slot='data',fc.slot='data',min.pct=cfg$markers$min_pct,logfc.threshold=0,
      only.pos=FALSE,test.use='wilcox',random.seed=cfg$markers$seed,verbose=FALSE),keep.rownames='orthogroup')
    if('avg_logFC'%in%names(fm))setnames(fm,'avg_logFC','avg_log2FC')
    setnames(fm,c('pct.1','pct.2'),c('pct_target','pct_background'))
    src <- objfile;stats_origin <- 'calculated_from_OG_SUM'
  }
  keep <- c('orthogroup','p_val','avg_log2FC','pct_target','pct_background','p_val_adj')
  stopifnot(!anyDuplicated(fm$orthogroup),all(keep%in%names(fm)))
  fm <- fm[,..keep]
  fm[,`:=`(object_species=spi,cell_type=ct,n_target_cells=length(ii),n_background_cells=length(jj),
            statistics_origin=stats_origin,statistics_source=src)]
  fm[,is_up:=avg_log2FC>.25 & pct_target>.05 & p_val_adj<.05]
  fm[,is_down:=avg_log2FC<(-.25) & pct_background>.05 & p_val_adj<.05]
  stopifnot(!any(fm$is_up & fm$is_down),all(is.finite(fm$avg_log2FC)))
  allstats[[ct]] <- fm
  avail <- intersect(tfogs,rownames(dat))
  mt <- Matrix::rowMeans(linear[avail,ii,drop=FALSE]); mb <- Matrix::rowMeans(linear[avail,jj,drop=FALSE])
  desc <- data.table(orthogroup=avail,mean_normalized_target=as.numeric(mt),mean_normalized_background=as.numeric(mb),
    expression_pct_target=as.numeric(Matrix::rowMeans(dat[avail,ii,drop=FALSE]>0)),
    expression_pct_background=as.numeric(Matrix::rowMeans(dat[avail,jj,drop=FALSE]>0)),
    descriptive_log2FC=log2((as.numeric(mt)+1/length(ii))/(as.numeric(mb)+1/length(jj))))
  d <- merge(data.table(orthogroup=tfogs),fm[orthogroup%in%tfogs],by='orthogroup',all.x=TRUE,sort=FALSE)
  d <- merge(d,desc,by='orthogroup',all.x=TRUE,sort=FALSE)
  d[,`:=`(object_species=spi,orthology_species=orthosp,cell_type=ct,
    n_target_cells=length(ii),n_background_cells=length(jj),statistics_origin=stats_origin,statistics_source=src,
    present_in_broccoli=orthogroup%in%present,mapped_to_RNA=orthogroup%in%mapped_any,
    mapped_positive_total=orthogroup%in%mapped_positive,tested_as_marker=!is.na(p_val))]
  d[is.na(is_up),is_up:=FALSE];d[is.na(is_down),is_down:=FALSE]
  d[,coverage_status:=fcase(mapped_positive_total,'mapped_positive_total',mapped_to_RNA,'mapped_but_zero_species_total',
                           present_in_broccoli,'unmapped_to_RNA',default='not_in_species_broccoli')]
  d[,marker_state:=fcase(is_up,'up',is_down,'down',tested_as_marker,'tested_not_threshold',default='not_tested')]
  d[,plot_log2FC:=fcoalesce(avg_log2FC,descriptive_log2FC)]
  d[,plot_value_source:=fcase(!is.na(avg_log2FC),'original_marker_statistics',!is.na(descriptive_log2FC),'descriptive_from_same_OG_SUM_data',default='missing')]
  d[,column_id:=paste(orthosp,ct,sep='|')]
  d[,include_plot:=!column_id%in%cfg$filter_ct]
  tfdetails[[ct]]<-d
  summaries[[ct]]<-data.table(object_species=spi,cell_type=ct,n_target_cells=length(ii),n_background_cells=length(jj),
    statistics_origin=stats_origin,n_ogs_tested=nrow(fm),n_tf_tested=sum(d$tested_as_marker),
    n_tf_up=sum(d$is_up),n_tf_down=sum(d$is_down),source=src,
    max_descriptive_vs_test_difference=max(abs(d$avg_log2FC-d$descriptive_log2FC),na.rm=TRUE))
}
write_tsv(rbindlist(allstats),file.path(out,'all_celltypes.OG_SUM.all_tested.tsv.gz'))
write_tsv(rbindlist(tfdetails),file.path(out,'TF_OG.species_celltype_details.tsv.gz'))
write_tsv(rbindlist(summaries),file.path(out,'run_summary.tsv'))
write_tsv(data.table(object_species=spi,source_object=objfile,sha256=sourcehash,
  started_utc=format(started,tz='UTC',usetz=TRUE),
  finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE),Seurat=as.character(packageVersion('Seurat')),
  n_ogs=nrow(obj),n_cells=ncol(obj)),file.path(out,'run_record.tsv'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('COMPLETED',spi,'\n')
