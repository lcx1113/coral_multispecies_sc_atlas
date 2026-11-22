#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(data.table);library(Matrix);library(digest)})
setDTthreads(2L)
arg<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
root<-normalizePath(file.path(dirname(arg),'..'))
for(d in c('common_og_scores','qc','logs')) dir.create(file.path(root,d),recursive=TRUE,showWarnings=FALSE)
species<-c('Hcoe','Xesp','Pdam','Spin','Fspp','Pcru','Plut','Gfas','Adig','Amur','Ahya')
types<-c('calicoblast','epidermis','gastrodermis');mods<-c('M3','M4');sets<-'narrative44'
write_tsv<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',na='NA')}
defs<-lapply(sets,function(s)fread(file.path(root,'inputs',paste0(s,'_definitions.tsv'))))
names(defs)<-sets
audit<-fread(file.path(root,'inputs/species_manifest.tsv'))
target_ogs<-sort(unique(defs$narrative44$og))
stopifnot(length(target_ogs)==44L)
args<-commandArgs(trailingOnly=TRUE)
stopifnot(all(args %in% c(paste0('--species=',species),'--summarize-only')))
chosen<-sub('^--species=','',args[grepl('^--species=',args)])
stopifnot(length(chosen)<=1L)

score_species<-function(sp){
  message('Loading ',sp,' ',Sys.time())
  ar<-audit[object_species==sp];stopifnot(nrow(ar)==1L)
  sha<-digest(file=ar$object_path,algo='sha256')
  obj<-readRDS(ar$object_path)
  data<-GetAssayData(obj,assay='SCT',layer='data');meta<-obj[[]]
  stopifnot(identical(colnames(data),rownames(meta)),ncol(data)==ar$n_total_cells,nrow(data)==ar$n_pool_features)
  use<-as.character(meta$cell_type) %in% types
  ids<-rownames(meta)[use];ct<-as.character(meta$cell_type[use]);pool<-rownames(data)
  stopifnot(length(ids)==ar$n_selected_cells,all(table(factor(ct,levels=types))==as.integer(ar[,.(n_calicoblast,n_epidermis,n_gastrodermis)])))
  X<-data[,ids,drop=FALSE]
  w<-fread(file.path(root,'inputs/shared_controls',paste0('OG_feature_weights.',sp,'.tsv.gz')))[og %in% target_ogs]
  features<-sort(unique(w$SCT_feature_id));ogs<-sort(unique(w$og))
  draws<-fread(file.path(root,'inputs/shared_controls',paste0('control_draws.',sp,'.tsv.gz')))[SCT_feature_id %in% features]
  stopifnot(all(features %in% pool),all(draws$control_SCT_feature_id %in% pool),
    setequal(features,draws$SCT_feature_id),all(draws[, .N,by=SCT_feature_id]$N==100L),
    !anyDuplicated(draws[,.(SCT_feature_id,control_SCT_feature_id)]),!anyDuplicated(w[,.(og,SCT_feature_id)]),
    all(abs(w[,sum(within_OG_weight),by=og]$V1-1)<1e-12))
  cm<-sparseMatrix(i=match(draws$SCT_feature_id,features),j=match(draws$control_SCT_feature_id,pool),x=rep(.01,nrow(draws)),dims=c(length(features),length(pool)))
  contrasts<-X[features,,drop=FALSE]-cm%*%X
  om<-sparseMatrix(i=match(w$og,ogs),j=match(w$SCT_feature_id,features),x=w$within_OG_weight,dims=c(length(ogs),length(features)))
  os<-om%*%contrasts;rownames(os)<-ogs;colnames(os)<-ids
  stopifnot(all(is.finite(os)))
  saveRDS(list(scores=os,cell_type=setNames(ct,ids),object_species=sp),file.path(root,'common_og_scores',paste0('OG_scores.',sp,'.rds')))
  status44<-fread(file.path(root,'inputs/narrative44',paste0('module_og_status.',sp,'.tsv')))[module %in% mods]
  for(set in sets){
    dest<-file.path(root,set,'scores/OG_equal/per_species');scores<-list();coverage<-list();ogsum<-list()
    for(mod in mods){
      declared<-defs[[set]][module==mod,og];available<-intersect(declared,ogs);K<-length(available);stopifnot(K>0L)
      label<-defs[[set]][module==mod,unique(module_label)];stopifnot(length(label)==1L)
      value<-as.numeric(Matrix::colMeans(os[available,,drop=FALSE]))
      scores[[mod]]<-data.table(object_species=sp,module_set=set,analysis_set='main',scheme='prose',module=mod,module_label=label,
        cell_id=ids,cell_type=ct,score=value,score_minus_gastrodermis_median=value-median(value[ct=='gastrodermis']))
      st<-status44[og %in% declared];mw<-w[og %in% available]
      stopifnot(setequal(st$og,declared),sum(st$n_SCT_features>0)==K)
      coverage[[mod]]<-data.table(object_species=sp,module_set=set,analysis_set='main',scheme='prose',module=mod,module_label=label,
        n_defined_OG=length(declared),n_OG_with_Broccoli_member=sum(st$n_native_member_genes>0),n_OG_with_SCT_feature=K,
        n_unique_SCT_features=uniqueN(mw$SCT_feature_id),coverage_fraction=K/length(declared),missing_OG=paste(setdiff(declared,available),collapse=';'),
        nbin_used=10L,ctrl=100L,seed=1L,assay='SCT',layer='data',control_source='frozen_marker308_332OG_union',
        n_calicoblast=sum(ct=='calicoblast'),n_epidermis=sum(ct=='epidermis'),n_gastrodermis=sum(ct=='gastrodermis'))
      ogsum[[mod]]<-rbindlist(lapply(types,function(t){z<-os[available,ct==t,drop=FALSE]
        data.table(object_species=sp,module_set=set,module=mod,module_label=label,og=available,cell_type=t,n_cells=ncol(z),
                   mean_score=Matrix::rowMeans(z),median_score=apply(z,1,median))}))
    }
    results<-rbindlist(scores)
    write_tsv(results,file.path(dest,paste0('cell_scores.',sp,'.tsv.gz')))
    write_tsv(rbindlist(coverage),file.path(dest,paste0('module_coverage.',sp,'.tsv')))
    write_tsv(rbindlist(ogsum),file.path(dest,paste0('og_score_by_celltype.',sp,'.tsv.gz')))
    write_tsv(status44[og %in% defs[[set]]$og],file.path(dest,paste0('module_og_status.',sp,'.tsv')))
    write_tsv(w[og %in% defs[[set]]$og],file.path(dest,paste0('OG_feature_weights.',sp,'.tsv.gz')))
  }
  write_tsv(data.table(object_species=sp,object_path=ar$object_path,object_sha256=sha,n_total_cells=ncol(data),n_pool_features=nrow(data),
    n_selected_cells=length(ids),n_calicoblast=sum(ct=='calicoblast'),n_epidermis=sum(ct=='epidermis'),n_gastrodermis=sum(ct=='gastrodermis')),
    file.path(root,'qc',paste0('object_audit.',sp,'.tsv')))
  message('Completed ',sp,' ',Sys.time())
}
if(!'--summarize-only' %in% args)for(sp in if(length(chosen))chosen else species){score_species(sp);gc(verbose=FALSE)}
if(length(chosen))quit(save='no',status=0)
aud<-rbindlist(lapply(species,function(sp)fread(file.path(root,'qc',paste0('object_audit.',sp,'.tsv')))))
stopifnot(sum(aud$n_selected_cells)==63295L);write_tsv(aud,file.path(root,'qc/object_audit.tsv'))
for(set in sets){
  dest<-file.path(root,set,'scores/OG_equal')
  combine<-function(stem,suffix='.tsv')rbindlist(lapply(species,function(sp)fread(file.path(dest,'per_species',paste0(stem,'.',sp,suffix)))))
  dat<-combine('cell_scores','.tsv.gz');cov<-combine('module_coverage');ogsum<-combine('og_score_by_celltype','.tsv.gz')
  stopifnot(nrow(dat)==126590L,!anyDuplicated(dat[,.(object_species,module,cell_id)]),all(is.finite(dat$score)))
  sm<-dat[,.(n_cells=.N,mean=mean(score),median=median(score),q25=as.numeric(quantile(score,.25)),q75=as.numeric(quantile(score,.75)),
    sd=sd(score),min_score=min(score),max_score=max(score)),by=.(object_species,module_set,analysis_set,scheme,module,module_label,cell_type)]
  stopifnot(nrow(sm)==66L,nrow(cov)==22L)
  write_tsv(dat,file.path(dest,'cell_scores.tsv.gz'));write_tsv(sm,file.path(dest,'cell_type_score_summary.tsv'))
  write_tsv(cov,file.path(dest,'module_coverage.tsv'));write_tsv(ogsum,file.path(dest,'og_score_by_celltype.tsv.gz'))
  wide<-dcast(sm,object_species+module_set+module+module_label~cell_type,value.var=c('n_cells','mean','median','q25','q75'))
  wide[,`:=`(median_difference_calicoblast_minus_epidermis=median_calicoblast-median_epidermis,
    median_difference_calicoblast_minus_gastrodermis=median_calicoblast-median_gastrodermis,
    median_difference_epidermis_minus_gastrodermis=median_epidermis-median_gastrodermis,
    mean_difference_calicoblast_minus_epidermis=mean_calicoblast-mean_epidermis,
    mean_difference_calicoblast_minus_gastrodermis=mean_calicoblast-mean_gastrodermis)]
  write_tsv(wide,file.path(dest,'within_species_contrasts.tsv'))
  contrib<-dcast(ogsum,object_species+module+og~cell_type,value.var='mean_score')
  contrib<-merge(contrib,cov[,.(object_species,module,n_OG_with_SCT_feature)],by=c('object_species','module'))
  contrib[,`:=`(contribution_calicoblast_minus_epidermis=(calicoblast-epidermis)/n_OG_with_SCT_feature,
                contribution_calicoblast_minus_gastrodermis=(calicoblast-gastrodermis)/n_OG_with_SCT_feature)]
  cc<-merge(contrib[,.(ce=sum(contribution_calicoblast_minus_epidermis),cg=sum(contribution_calicoblast_minus_gastrodermis)),by=.(object_species,module)],wide,by=c('object_species','module'))
  stopifnot(max(abs(cc$ce-cc$mean_difference_calicoblast_minus_epidermis),abs(cc$cg-cc$mean_difference_calicoblast_minus_gastrodermis))<1e-10)
  write_tsv(contrib,file.path(dest,'og_score_contributions.tsv.gz'))
}
capture.output(sessionInfo(),file=file.path(root,'logs/sessionInfo_scoring.txt'))
writeLines(c('PASS: SCT dimensions and exact three target-cell labels match the species manifest',
  'PASS: 44-OG scores use frozen per-feature controls and within-OG weights; no random sampling',
  'PASS: 126590 scores and 66 group summaries; OG mean contributions sum to module mean contrasts'),file.path(root,'qc/SCORE_VALIDATION.txt'))
message('ALL SCORING COMPLETE')
