#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(4);root<-normalizePath(commandArgs(TRUE)[1]);cfg<-fromJSON(file.path(root,'config/analysis_config.json'))
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',na='NA',quote=FALSE)}
raw<-fread(file.path(root,'results/TF_OG.original_species_celltype_details.tsv.gz'))
tf<-fread(file.path(root,'results/TF_OG.comparison_conservation.tsv.gz'))
cols<-fread(file.path(root,'results/plot_columns.tsv'));units<-fread(file.path(root,'config/original_combination_manifest.tsv'))
keys<-data.table(figure_key=paste0('majority_all_celltypes/',c('up','down')),matrix_name='plot_log2FC_matrix.tsv',overview=TRUE,direction=c('up','down'),comparison_id='all_celltypes')
for(id in cfg$comparison_order)for(direction in c('up','down')){
 cc<-copy(cols[column_id%in%units[comparison_id==id & include_plot==TRUE,column_id]])
 if(!nrow(cc))next
 if(id%in%names(cfg$lineage_groups))cc[,plot_column_group:=cell_type]
 rr<-copy(tf[comparison_id==id & get(paste0('is_',direction,'_selected'))==TRUE])
 rr[,`:=`(support=if(direction=='up')n_up else n_down,original_unit_support=if(direction=='up')n_up_units else n_down_units,primary_cell_type=id,effect=if(direction=='up')median_log2FC else -median_log2FC)]
 setorder(rr,-support,-effect,orthogroup);rr[,row_order:=seq_len(.N)]
 key<-paste('by_comparison',direction,id,sep='/');ddir<-file.path(root,'results/plot_data',key)
 put(rr,file.path(ddir,'row_annotations.tsv'));put(cc,file.path(ddir,'column_annotations.tsv'))
 keys<-rbind(keys,data.table(figure_key=key,matrix_name='matrix.tsv',overview=FALSE,direction=direction,comparison_id=id))
 if(!nrow(rr))next
 dd<-raw[orthogroup%in%rr$orthogroup & column_id%in%cc$column_id]
 for(field in c('plot_log2FC','marker_state','plot_value_source')){
  mm<-dcast(dd,orthogroup~column_id,value.var=field)[match(rr$orthogroup,orthogroup)]
  setcolorder(mm,c('orthogroup',cc$column_id))
  put(mm,file.path(ddir,if(field=='plot_log2FC')'matrix.tsv' else paste0(field,'_matrix.tsv')))
 }
}
for(k in seq_len(nrow(keys))){
 f<-keys[k];ddir<-file.path(root,'results/plot_data',f$figure_key);out<-file.path(root,'results/transposed_plot_data',f$figure_key)
 rr<-fread(file.path(ddir,'row_annotations.tsv'));cc<-fread(file.path(ddir,'column_annotations.tsv'))
 put(cc,file.path(out,'row_annotations.tsv'));put(rr,file.path(out,'column_annotations.tsv'))
 if(!nrow(rr))next
 for(name in c(f$matrix_name,'marker_state_matrix.tsv','plot_value_source_matrix.tsv')){
  mm<-fread(file.path(ddir,name));x<-as.matrix(mm[,-1]);rownames(x)<-mm$orthogroup
  tx<-t(x);put(data.table(column_id=rownames(tx),as.data.table(tx)),file.path(out,name))
 }
}
put(keys,file.path(root,'config/plotting_manifest.tsv'))
cat('PLOT TABLES COMPLETE',nrow(keys),'figures\n')
