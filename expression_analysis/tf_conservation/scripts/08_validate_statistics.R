suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(4);root<-normalizePath(commandArgs(TRUE)[1]);cfg<-fromJSON(file.path(root,'config/analysis_config.json'))
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',quote=FALSE,na='NA')}
# Verify every retained vote using original native-celltype marker decisions.
d<-fread(file.path(root,'results/TF_OG.comparison_original_combination_details.tsv.gz'))
tf<-fread(file.path(root,'results/TF_OG.comparison_conservation.tsv.gz'))
v<-d[orthogroup%in%tf$orthogroup,.(up=any(is_up),down=any(is_down),tested=any(tested_as_marker)),by=.(orthogroup,comparison_id,object_species)]
s<-v[,.(up=sum(up),down=sum(down),tested=sum(tested),N=.N),by=.(orthogroup,comparison_id)]
s<-merge(s,tf,by=c('orthogroup','comparison_id'))
stopifnot(nrow(s)==nrow(tf),all(s$up==s$n_up_species),all(s$down==s$n_down_species),all(s$tested==s$n_tested_species),all(s$N==s$n_species),
 all(s$is_up_selected==(s$up>=ifelse(s$N<=2,s$N,ceiling(s$N/2)))),all(s$is_down_selected==(s$down>=ifelse(s$N<=2,s$N,ceiling(s$N/2)))))
# Frozen member annotation and complete, unmodified first_annotation.
a<-fread(file.path(cfg$annotation_root,'orthogroup_final_annotations.csv'))
stopifnot(all(tf$first_annotation==a$first_annotation[match(tf$orthogroup,a$orthogroup)]),all(5*tf$n_tf_genes>=4*tf$n_unique_species_gene))
mem<-fread(file.path(root,'results/tf_og_member_evidence.tsv'))
now<-fread(file.path(cfg$annotation_root,'orthogroup_conservation_final.csv'))[orthogroup%in%mem$orthogroup]
setorderv(mem,names(mem));setorderv(now,names(mem));stopifnot(identical(mem,now))
put(data.table(check=c('species_OR_votes_and_fixed_denominators','TF80_and_source_first_annotation','frozen_member_annotations'),n_rows=c(nrow(s),nrow(tf),nrow(mem)),status='PASS'),file.path(root,'results/statistical_validation.tsv'))
# Diagonal block indices are in inherited group order; only membership and within-group sorting are refreshed.
blocks<-list()
for(direction in c('up','down')){
 rr<-fread(file.path(root,'results/plot_data/majority_all_celltypes',direction,'row_annotations.tsv'))
 cc<-fread(file.path(root,'results/plot_columns.tsv'));mp<-fread(file.path(root,'config/original_celltype_to_lineage.tsv'))
 cc[,lineage:=mp$row_group[match(cell_type,mp$original_cell_type)]]
 stopifnot(!anyDuplicated(rr$orthogroup),!is.unsorted(match(rr$primary_cell_type,cfg$comparison_order)),!is.unsorted(match(cc$lineage,cfg$comparison_order)))
 for(g in unique(rr$primary_cell_type)){
  ii<-which(rr$primary_cell_type==g);jj<-which(cc$lineage==g)
  blocks[[length(blocks)+1]]<-data.table(direction=direction,group=g,OG_column_start=min(ii),OG_column_end=max(ii),native_row_start=min(jj),native_row_end=max(jj),N_species=rr$n_species[ii[1]])
 }
}
put(rbindlist(blocks),file.path(root,'results/diagonal_block_audit.tsv'))
cat('STATISTICAL VALIDATION PASS\n')
