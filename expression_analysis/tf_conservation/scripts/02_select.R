#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(jsonlite)})
setDTthreads(4)
root<-normalizePath(commandArgs(TRUE)[1]);cfg<-fromJSON(file.path(root,'config/analysis_config.json'))
broad<-root
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',quote=FALSE,na='NA',compress=if(grepl('gz$',p))'gzip' else 'none')}
comps<-fread(file.path(root,'config/comparison_thresholds.tsv'))
mapping<-fread(file.path(root,'config/original_celltype_to_lineage.tsv'))
lookup<-setNames(mapping$row_group,mapping$original_cell_type)
group_for<-function(x){z<-unname(lookup[x]);stopifnot(!anyNA(z));z}
raw<-rbindlist(lapply(cfg$species_order,function(s)fread(file.path(root,'results/by_species',s,'TF_OG.species_celltype_details.tsv.gz'))))
stopifnot(!any(raw$cell_type%in%cfg$excluded_analysis_targets))
stopifnot(!anyDuplicated(raw[,.(orthogroup,object_species,cell_type)]))
put(raw,file.path(root,'results/TF_OG.original_species_celltype_details.tsv.gz'))
details<-copy(raw);details[,comparison_id:=group_for(cell_type)]
stopifnot(!anyDuplicated(details[,.(orthogroup,comparison_id,object_species,cell_type)]))
put(details,file.path(root,'results/TF_OG.comparison_original_combination_details.tsv.gz'))
ordered_species<-function(x)paste(cfg$species_order[cfg$species_order%in%x],collapse=';')
med<-function(x)if(all(is.na(x)))NA_real_ else median(x,na.rm=TRUE)

# Reuse each native FindMarkers call, then apply logical OR within each species.
allvotes<-list()
for(s in cfg$species_order){
 x<-fread(file.path(broad,'results/by_species',s,'all_celltypes.OG_SUM.all_tested.tsv.gz'),
  select=c('orthogroup','object_species','cell_type','is_up','is_down'))
 x<-x[cell_type%in%mapping$original_cell_type]
 stopifnot(all(x$object_species==s),!anyDuplicated(x[,.(orthogroup,cell_type)]))
 x[,comparison_id:=group_for(cell_type)]
 v<-x[,.(tested_as_marker=TRUE,is_up=any(is_up),is_down=any(is_down),
  n_tested_subtypes=.N,n_up_subtypes=sum(is_up),n_down_subtypes=sum(is_down)),
  by=.(orthogroup,comparison_id,object_species)]
 put(v,file.path(root,'results/by_species',s,'all_OG.tested_species_votes.tsv.gz'))
 allvotes[[s]]<-v
}
av<-rbindlist(allvotes,use.names=TRUE)
counts<-av[,.(n_tested_species=.N,n_up_species=sum(is_up),n_down_species=sum(is_down),
 n_tested_units=sum(n_tested_subtypes),n_up_units=sum(n_up_subtypes),n_down_units=sum(n_down_subtypes),
 up_species=ordered_species(object_species[is_up]),down_species=ordered_species(object_species[is_down]),
 n_species_with_both_directions=sum(is_up&is_down)),by=.(orthogroup,comparison_id)]
pa<-fread(file.path(root,cfg$og_annotation_file))
af<-c('orthogroup','has_TF','first_annotation','first_support','second_annotation','second_support','final_annotation')
annotation<-unique(pa[,..af]);universe<-sort(annotation$orthogroup,method='radix')
stopifnot(length(universe)==21743)
grid<-CJ(orthogroup=universe,comparison_id=cfg$comparison_order,sorted=FALSE)
all<-merge(grid,counts,by=c('orthogroup','comparison_id'),all.x=TRUE,sort=FALSE)
for(f in c('n_tested_species','n_up_species','n_down_species','n_tested_units','n_up_units','n_down_units','n_species_with_both_directions'))set(all,which(is.na(all[[f]])),f,0L)
for(f in c('up_species','down_species'))set(all,which(is.na(all[[f]])),f,'')
all<-merge(all,comps,by='comparison_id',all.x=TRUE,sort=FALSE)
all<-merge(all,annotation,by='orthogroup',all.x=TRUE,sort=FALSE)
classify<-function(k,n,b,m)fcase(k>=n,'core',k>=b,'broad_only',k>=m,'majority_only',k>=2,'shared_only',k==1,'species_specific',default='not_positive')
all[,`:=`(conservation_class=classify(n_up_species,n_species,broad_min_species,majority_min_species),
 down_conservation_class=classify(n_down_species,n_species,broad_min_species,majority_min_species),
 is_up_selected=n_up_species>=selected_min_species,is_down_selected=n_down_species>=selected_min_species,
 up_fraction=n_up_species/n_species,down_fraction=n_down_species/n_species,
 direction_conflict=n_up_species>0&n_down_species>0)]
all[down_conservation_class=='not_positive',down_conservation_class:='not_downregulated']
stopifnot(all(all$n_tested_species<=all$n_species),all(all$n_up_species<=all$n_tested_species),
 all(all$n_down_species<=all$n_tested_species))
setorder(all,comparison_order,orthogroup)
put(all,file.path(root,'results/all_OG.comparison_conservation.tsv.gz'))
cat('ALL OG SPECIES VOTES COMPLETE\n');flush.console()

sv<-details[,.(n_original_units=.N,n_tested_subtypes=sum(tested_as_marker),n_up_subtypes=sum(is_up),n_down_subtypes=sum(is_down),
 tested_as_marker=any(tested_as_marker),is_up=any(is_up),is_down=any(is_down),
 present_in_broccoli=any(present_in_broccoli),mapped_to_RNA=any(mapped_to_RNA),mapped_positive_total=any(mapped_positive_total),
 coverage_status=first(coverage_status),
 source_cell_types=paste(cell_type,collapse=';'),up_source_cell_types=paste(cell_type[is_up],collapse=';'),
 down_source_cell_types=paste(cell_type[is_down],collapse=';'),
 up_original_combinations=paste(column_id[is_up],collapse=';'),down_original_combinations=paste(column_id[is_down],collapse=';')),
 by=.(orthogroup,comparison_id,object_species)]
sv[,species_vote_state:=fcase(is_up&is_down,'up_and_down',is_up,'up_only',is_down,'down_only',
 tested_as_marker,'tested_not_threshold',default='not_tested')]
stopifnot(!anyDuplicated(sv[,.(orthogroup,comparison_id,object_species)]))
put(sv,file.path(root,'results/TF_OG.comparison_species_votes.tsv.gz'))
template<-fread(file.path(root,cfg$tf_member_summary_file))
ef<-c('orthogroup','n_member_records','n_canonical_genes','n_tf_member_records','n_unflagged_records','n_unknown_records','n_species_with_tf','tf_support_species','tf_label')
ev<-unique(template[,..ef])
tf<-merge(all[has_TF%in%TRUE],ev,by='orthogroup',all.x=TRUE,sort=FALSE)
state<-sv[,.(n_species_votes=.N,n_tested_direct=sum(tested_as_marker),n_up_direct=sum(is_up),n_down_direct=sum(is_down),
 n_up_only_species=sum(species_vote_state=='up_only'),n_down_only_species=sum(species_vote_state=='down_only'),
 n_up_and_down_species=sum(species_vote_state=='up_and_down'),
 n_tested_not_threshold_species=sum(species_vote_state=='tested_not_threshold'),n_not_tested_species=sum(species_vote_state=='not_tested'),
 n_broccoli_present_species=sum(present_in_broccoli),n_mapped_RNA_species=sum(mapped_to_RNA),n_mapped_positive_species=sum(mapped_positive_total)),
 by=.(orthogroup,comparison_id)]
effects<-details[,.(median_log2FC=med(plot_log2FC),
 min_log2FC=if(all(is.na(plot_log2FC)))NA_real_ else min(plot_log2FC,na.rm=TRUE),
 max_log2FC=if(all(is.na(plot_log2FC)))NA_real_ else max(plot_log2FC,na.rm=TRUE),
 up_units=paste(column_id[is_up],collapse=';'),down_units=paste(column_id[is_down],collapse=';')),
 by=.(orthogroup,comparison_id)]
tf<-merge(tf,state,by=c('orthogroup','comparison_id'),all.x=TRUE,sort=FALSE)
tf<-merge(tf,effects,by=c('orthogroup','comparison_id'),all.x=TRUE,sort=FALSE)
stopifnot(all(tf$n_species_votes==tf$n_species),all(tf$n_tested_direct==tf$n_tested_species),
 all(tf$n_up_direct==tf$n_up_species),all(tf$n_down_direct==tf$n_down_species),
 all(tf$n_up_and_down_species==tf$n_species_with_both_directions),
 all(tf$n_up_only_species+tf$n_down_only_species+tf$n_up_and_down_species+tf$n_tested_not_threshold_species==tf$n_tested_species),
 all(tf$n_tested_species+tf$n_not_tested_species==tf$n_species))
tf[,`:=`(n_up=n_up_species,n_down=n_down_species,n_tested=n_tested_species,eligible_for_plot=n_display_original_units>0,
 is_up_broad=n_up_species>=broad_min_species,is_down_broad=n_down_species>=broad_min_species)]
# TF identity is member annotation, independent of expression conservation.
mem<-fread(file.path(root,'results/tf_og_member_evidence.tsv'))
genes<-mem[,.(tf=any(is_tf%in%TRUE),n_distinct_flags=uniqueN(is_tf)),by=.(orthogroup,species,gene)]
stopifnot(all(genes$n_distinct_flags==1))
fractions<-genes[,.(n_unique_species_gene=.N,n_tf_genes=sum(tf),tf_member_fraction=mean(tf),passes_TF80=5*sum(tf)>=4*.N),by=orthogroup]
put(fractions,file.path(root,'results/TF_OG.member_fraction_all_candidates.tsv'))
put(tf,file.path(root,'results/TF_OG.all_TF_comparison_conservation.before_TF80.tsv.gz'))
tf<-merge(tf,fractions,by='orthogroup',all.x=TRUE,sort=FALSE)[passes_TF80==TRUE]
setorder(tf,comparison_order,orthogroup)
put(tf,file.path(root,'results/TF_OG.comparison_conservation.tsv.gz'))
put(tf[is_up_selected|is_down_selected],file.path(root,'results/TF_OG.selected_associations.tsv'))
put(tf[,.(n_tf_ogs=.N,n_up_selected=sum(is_up_selected),n_down_selected=sum(is_down_selected),
 n_up_and_down_selected=sum(is_up_selected&is_down_selected),
 n_OGs_with_within_species_direction_conflict=sum(n_species_with_both_directions>0)),
 by=.(comparison_id,n_species,selected_min_species,n_original_units,n_display_original_units,inference_scope)],
 file.path(root,'results/comparison_candidate_counts.tsv'))
put(all[,.(n_all_OG=.N,n_up_selected=sum(is_up_selected),n_down_selected=sum(is_down_selected)),
 by=.(comparison_id,n_species,selected_min_species)],file.path(root,'results/all_OG_candidate_counts.tsv'))
cols<-fread(file.path(root,'results/plot_columns.tsv'))
make_data<-function(links,cc,pdir,direction){
 links<-copy(links)
 links[,`:=`(support=if(direction=='up')n_up_species else n_down_species,
  original_unit_support=if(direction=='up')n_up_units else n_down_units,
  fraction=(if(direction=='up')n_up_species else n_down_species)/n_species,
  directional_effect=if(direction=='up')median_log2FC else -median_log2FC)]
 setorder(links,orthogroup,-support,-fraction,-directional_effect,comparison_order)
 rr<-links[!duplicated(orthogroup)]
 rr[,`:=`(primary_comparison=comparison_id,primary_cell_type=comparison_id,plot_direction=direction)]
 setorder(rr,comparison_order,-fraction,-directional_effect,orthogroup);rr[,row_order:=seq_len(.N)]
 put(rr,file.path(pdir,'row_annotations.tsv'));put(links,file.path(pdir,'all_comparison_associations.tsv'));put(cc,file.path(pdir,'column_annotations.tsv'))
 for(field in c('plot_log2FC','marker_state','plot_value_source')){
  dd<-raw[orthogroup%in%rr$orthogroup & column_id%in%cc$column_id]
  if(nrow(rr)){
   mm<-dcast(dd,orthogroup~column_id,value.var=field)[match(rr$orthogroup,orthogroup)]
   setcolorder(mm,c('orthogroup',cc$column_id))
  }else mm<-as.data.table(setNames(replicate(nrow(cc)+1,character(),simplify=FALSE),c('orthogroup',cc$column_id)))
  put(mm,file.path(pdir,paste0(field,'_matrix.tsv')))
 }
}
for(direction in c('up','down')){
 flag<-paste0('is_',direction,'_selected')
 make_data(tf[get(flag)==TRUE & eligible_for_plot==TRUE],cols,
  file.path(root,'results/plot_data/majority_all_celltypes',direction),direction)
}
manual<-fread(file.path(broad,'results/legacy_manual_gene_og_mapping.tsv'))
versions<-list()
for(direction in c('up','down')){
 rr<-fread(file.path(root,'results/plot_data/majority_all_celltypes',direction,'row_annotations.tsv'))
 manual[,(paste0('speciesvotes_',direction)):=orthogroup%in%rr$orthogroup]
 versions[[direction]]<-data.table(version='majority_all_celltypes',direction=direction,n_OGs=nrow(rr),n_columns=nrow(cols))
}
put(manual,file.path(root,'results/legacy_manual_gene_og_audit.tsv'))
put(rbindlist(versions),file.path(root,'results/version_summary.tsv'))
cat('SELECTION COMPLETE\n');print(rbindlist(versions))
