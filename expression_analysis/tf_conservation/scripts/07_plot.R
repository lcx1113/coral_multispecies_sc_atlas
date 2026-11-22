#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite);library(ComplexHeatmap);library(circlize);library(grid)})
setDTthreads(4);root<-normalizePath(commandArgs(TRUE)[1]);layout<-fromJSON(file.path(root,'config/layout_config.json'))
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',na='NA',quote=FALSE)}
mapping<-fread(file.path(root,'config/original_celltype_to_lineage.tsv'));celltype_to_group<-setNames(mapping$row_group,mapping$original_cell_type)
palette<-fread(file.path(root,'config/celltype_group_palette.tsv'));group_colors<-setNames(palette$color,palette$celltype_group)
species_labels<-fread(file.path(root,'config/species_label_mapping.tsv'));display_map<-fread(file.path(root,'results/annotation_display_mapping.tsv'))
color_fun<-colorRamp2(layout$color_breaks,c(layout$negative_color,layout$zero_color,layout$positive_color))
put(data.table(log2FC=seq(-3,3,by=.25),color=color_fun(seq(-3,3,by=.25))),file.path(root,'results/layout_color_mapping.tsv'))
figures<-list();source(file.path(root,'scripts/plot_renderer.R'))
keys<-fread(file.path(root,'config/plotting_manifest.tsv'))
read_matrix<-function(path,first){d<-fread(path);stopifnot(names(d)[1]==first);x<-as.matrix(d[,-1]);rownames(x)<-d[[1]];x}
for(k in seq_len(nrow(keys))){
 f<-keys[k];figure_key<-f$figure_key;overview<-f$overview;direction<-f$direction
 olddir<-file.path(root,'results/plot_data',figure_key);newdir<-file.path(root,'results/transposed_plot_data',figure_key)
 rr<-fread(file.path(olddir,'row_annotations.tsv'));cc<-fread(file.path(olddir,'column_annotations.tsv'))
 prefix<-if(overview)file.path(root,'figures/majority_all_celltypes',paste0('TF_OG.',direction,'.all_celltypes')) else file.path(root,'figures/by_comparison',direction,f$comparison_id)
 title<-if(overview)paste0('TF-associated OGs (TF members >=80%): ',direction,'regulated; one vote per species') else paste0(f$comparison_id,': ',direction,'regulated; TF members >=80%')
 width<-if(overview){if(direction=='up')max(10,20.25+(nrow(rr)-82)*.125) else max(10,10+(nrow(rr)-21)*.125)} else max(8,nrow(rr)*.19+5)
 height<-if(overview)28.84+(nrow(cc)-157)*.12 else max(12,nrow(cc)*.42+8)
 if(!nrow(rr)){
  dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE);pdf(paste0(prefix,'.pdf'),width=width,height=height,useDingbats=FALSE)
  grid.text(title,y=.7);grid.text('No TF-related OG passes the selected threshold',y=.45);dev.off()
  figures[[length(figures)+1]]<-data.table(figure=paste0(prefix,'.pdf'),n_rows=nrow(cc),n_columns=0L,n_OGs=0L,n_species_celltype_combinations=nrow(cc),missing_cells=0L,marker_cells=0L,width_inches=width,height_inches=height,axes_transposed=TRUE)
  next
 }
 rr[,display_annotation:=display_map$display_annotation[match(orthogroup,display_map$orthogroup)]]
 rr[,plot_annotation:=display_map$plot_annotation[match(orthogroup,display_map$orthogroup)]]
 native<-read_matrix(file.path(olddir,f$matrix_name),'orthogroup');storage.mode(native)<-'numeric'
 native_states<-read_matrix(file.path(olddir,'marker_state_matrix.tsv'),'orthogroup')
 x<-read_matrix(file.path(newdir,f$matrix_name),'column_id');storage.mode(x)<-'numeric'
 st<-read_matrix(file.path(newdir,'marker_state_matrix.tsv'),'column_id')
 stopifnot(identical(x,t(native)),identical(st,t(native_states)))
 if(!overview)title<-paste0(title,'; >=',rr$selected_min_species[1],'/',rr$n_species[1],' species',if(rr$n_species[1]==1)' (single species)' else '')
 render_transposed(x,st,rr,cc,prefix,title,overview,width,height)
 cat('RENDERED',figure_key,nrow(x),'x',ncol(x),'\n');flush.console()
}
put(rbindlist(figures,fill=TRUE),file.path(root,'results/figure_manifest.tsv'))
cat('PLOTTING COMPLETE',length(figures),'PDFs\n')
