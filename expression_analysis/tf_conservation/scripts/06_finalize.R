#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(4)
root<-normalizePath(commandArgs(TRUE)[1]);cfg<-fromJSON(file.path(root,'config/analysis_config.json'))
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',quote=FALSE,na='NA',compress=if(grepl('gz$',p))'gzip' else 'none')}
contrib<-rbindlist(lapply(cfg$species_order,function(s)fread(file.path(root,'results/member_contribution',s,'OG_TF_member_contribution.tsv.gz'))))
valid<-rbindlist(lapply(cfg$species_order,function(s)fread(file.path(root,'results/member_contribution',s,'validation.tsv'))))
stopifnot(nrow(valid)==11,all(valid$raw_UMI_reconstruction_max_error==0),all(valid$normalized_mean_reconstruction_max_error<1e-8))
put(contrib,file.path(root,'results/all_species.OG_TF_member_contribution.tsv.gz'))
put(valid,file.path(root,'results/all_species.member_contribution_validation.tsv'))
details<-fread(file.path(root,'results/TF_OG.comparison_original_combination_details.tsv.gz'))
tf<-fread(file.path(root,'results/TF_OG.comparison_conservation.tsv.gz'))
wanted<-unique(tf[is_up_selected|is_down_selected,orthogroup])
ck<-c('orthogroup','object_species','cell_type')
positive<-unique(details[orthogroup%in%wanted & mapped_positive_total,..ck])
stopifnot(nrow(fsetdiff(positive,contrib[,..ck]))==0)
native<-merge(details[orthogroup%in%wanted & (is_up|is_down)],contrib,by=ck,all.x=TRUE,sort=FALSE)
native[,direction:=ifelse(is_up,'up','down')]
put(native,file.path(root,'results/TF_OG.marker_support_with_native_member_evidence.tsv.gz'))
med<-function(x)if(all(is.na(x)))NA_real_ else median(x,na.rm=TRUE)
evidence<-native[,.(n_marker_support_units=.N,n_marker_support_species=uniqueN(object_species),
 n_target_TF_fraction_available=sum(!is.na(tf_marked_contribution_fraction)),
 n_target_TF_fraction_ge_half=sum(tf_marked_contribution_fraction>=.5,na.rm=TRUE),
 median_target_TF_fraction=med(tf_marked_contribution_fraction),
 n_dominant_target_gene_TF_marked=sum(dominant_gene_tf_marked%in%TRUE),
 dominant_target_genes=paste(paste(column_id,dominant_canonical_genes,sep=':'),collapse=';')),
 by=.(orthogroup,comparison_id,direction)]
selected<-rbindlist(list(copy(tf[is_up_selected==TRUE])[,direction:='up'],copy(tf[is_down_selected==TRUE])[,direction:='down']))
selected<-merge(selected,evidence,by=c('orthogroup','comparison_id','direction'),all.x=TRUE,sort=FALSE)
stopifnot(all(selected$n_marker_support_species==ifelse(selected$direction=='up',selected$n_up_species,selected$n_down_species)),
 all(selected$n_marker_support_units==ifelse(selected$direction=='up',selected$n_up_units,selected$n_down_units)))
put(selected,file.path(root,'results/TF_OG.selected_associations.with_member_support.tsv'))
species_support<-native[,.(n_supporting_original_subtypes=.N,source_cell_types=paste(cell_type,collapse=';'),
 median_target_TF_fraction_across_supporting_subtypes=med(tf_marked_contribution_fraction),
 n_subtypes_with_TF_fraction=sum(!is.na(tf_marked_contribution_fraction)),
 dominant_target_genes=paste(paste(cell_type,dominant_canonical_genes,sep=':'),collapse=';')),
 by=.(orthogroup,comparison_id,object_species,direction)]
put(species_support,file.path(root,'results/TF_OG.species_marker_support.with_member_evidence.tsv.gz'))
cat('MEMBER EVIDENCE FINALIZED\n');print(valid)
