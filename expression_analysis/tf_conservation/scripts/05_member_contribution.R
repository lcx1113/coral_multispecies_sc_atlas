#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite);library(digest)})
setDTthreads(4)
args<-commandArgs(TRUE);root<-normalizePath(args[1]);spi<-args[2];cfg<-fromJSON(file.path(root,'config/analysis_config.json'))
put<-function(x,p){dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);fwrite(x,p,sep='\t',na='NA',quote=FALSE,compress=if(grepl('gz$',p))'gzip' else 'none')}
ogsum<-fread(file.path(root,'results/TF_OG.comparison_conservation.tsv.gz'))
chosen<-sort(unique(ogsum[is_up_selected | is_down_selected,orthogroup]),method='radix')
manifest<-fread(file.path(cfg$source_root,'config/species_manifest.tsv'))
rnafile<-manifest[object_species==spi,seurat_file][1]
ogfile<-file.path(cfg$source_root,'intermediate/02_og_sum_objects',paste0(spi,'.OG_SUM.seurat.rds'))
cat('READING',spi,'native RNA and OG_SUM\n');flush.console()
seu<-readRDS(rnafile);obj<-readRDS(ogfile)
rna<-LayerData(seu,assay='RNA',layer='counts');ogs<-LayerData(obj,assay='OG_SUM',layer='counts')
stopifnot(identical(colnames(rna),colnames(ogs)),identical(as.character(seu$cell_type),as.character(obj$cell_type)))
mapfile<-cfg$gene_og_mapping
map<-fread(mapfile,select=c('gene','object_species','gene_seurat','orthogroup','is_orphan'))[object_species==spi & orthogroup%in%chosen]
map<-map[gene_seurat%in%rownames(rna)]
evidence<-fread(file.path(root,'results/tf_og_member_evidence.tsv'))
flags<-evidence[,.(has_tf=any(is_tf%in%TRUE),has_false=any(is_tf%in%FALSE),has_unknown=any(is.na(is_tf))),by=.(orthogroup,gene)]
map<-merge(map,flags,by=c('orthogroup','gene'),all.x=TRUE,sort=FALSE)
member<-map[,.(canonical_genes=paste(sort(unique(gene),method='radix'),collapse=';'),
  tf_marked=any(has_tf%in%TRUE),annotation_conflict=any(has_tf%in%TRUE)&any(has_false%in%TRUE),
  annotation_unknown=all(is.na(has_tf)) | any(has_unknown%in%TRUE)),by=.(orthogroup,gene_seurat)]
features<-sort(unique(member$gene_seurat),method='radix');ogids<-sort(unique(member$orthogroup),method='radix')
membership<-sparseMatrix(i=match(member$orthogroup,ogids),j=match(member$gene_seurat,features),x=rep(1,nrow(member)),
                        dims=c(length(ogids),length(features)),dimnames=list(ogids,features))
gene_counts<-rna[features,,drop=FALSE]
rebuilt<-membership%*%gene_counts
available<-intersect(ogids,rownames(ogs))
raw_error<-if(length(available))max(abs(rebuilt[available,,drop=FALSE]-ogs[available,,drop=FALSE])) else 0
stopifnot(raw_error==0)
scaling<-10000/Matrix::colSums(ogs);scaling[!is.finite(scaling)]<-0
normalized_genes<-gene_counts%*%Diagonal(x=scaling)
linear<-LayerData(obj,assay='OG_SUM',layer='data');linear@x<-expm1(linear@x)
labels<-as.character(obj$cell_type);targets<-cfg$celltype_order[cfg$celltype_order%in%unique(labels)]
results<-list();gresults<-list();max_error<-0
for(ct in targets){
  ii<-which(labels==ct);gm<-Matrix::rowMeans(normalized_genes[,ii,drop=FALSE])
  d<-copy(member);d[,mean_normalized_contribution:=as.numeric(gm[match(gene_seurat,features)])]
  d[,`:=`(object_species=spi,cell_type=ct)]
  z<-d[,.(mean_all_members=sum(mean_normalized_contribution),
    mean_tf_marked_members=sum(mean_normalized_contribution[tf_marked]),
    n_mapped_features=.N,n_tf_marked_features=sum(tf_marked),n_conflicting_features=sum(annotation_conflict)),by=orthogroup]
  setorder(d,orthogroup,-mean_normalized_contribution,gene_seurat)
  top<-d[!duplicated(orthogroup),.(orthogroup,dominant_gene_seurat=gene_seurat,dominant_canonical_genes=canonical_genes,
                                  dominant_gene_tf_marked=tf_marked,dominant_gene_contribution=mean_normalized_contribution)]
  z<-merge(z,top,by='orthogroup',all.x=TRUE)
  expected<-data.table(orthogroup=available,mean_OG_SUM=as.numeric(Matrix::rowMeans(linear[available,ii,drop=FALSE])))
  z<-merge(z,expected,by='orthogroup',all.x=TRUE)
  err<-max(abs(z$mean_all_members-z$mean_OG_SUM),na.rm=TRUE);max_error<-max(max_error,err)
  z[,`:=`(tf_marked_contribution_fraction=fifelse(mean_all_members>0,mean_tf_marked_members/mean_all_members,NA_real_),
            dominant_gene_fraction=fifelse(mean_all_members>0,dominant_gene_contribution/mean_all_members,NA_real_),
            object_species=spi,cell_type=ct)]
  z[mean_all_members==0,`:=`(dominant_gene_seurat=NA_character_,dominant_canonical_genes=NA_character_,dominant_gene_tf_marked=NA)]
  d<-merge(d,z[,.(orthogroup,mean_all_members)],by='orthogroup',all.x=TRUE)
  d[,contribution_fraction:=fifelse(mean_all_members>0,mean_normalized_contribution/mean_all_members,NA_real_)]
  results[[ct]]<-z;gresults[[ct]]<-d
}
stopifnot(max_error<1e-8)
out<-file.path(root,'results/member_contribution',spi)
put(rbindlist(results),file.path(out,'OG_TF_member_contribution.tsv.gz'))
put(rbindlist(gresults),file.path(out,'gene_contributions.tsv.gz'))
put(data.table(object_species=spi,n_candidate_ogs=length(chosen),n_reconstructed_ogs=length(available),
               raw_UMI_reconstruction_max_error=raw_error,normalized_mean_reconstruction_max_error=max_error,
               source_RNA=rnafile,source_RNA_sha256=digest(file=rnafile,algo='sha256',serialize=FALSE)),file.path(out,'validation.tsv'))
cat('MEMBER CONTRIBUTIONS COMPLETE',spi,'max error',max_error,'\n')
