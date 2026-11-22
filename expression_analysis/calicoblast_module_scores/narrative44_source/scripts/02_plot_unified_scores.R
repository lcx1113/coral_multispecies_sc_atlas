#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
file_arg<-sub("^--file=","",commandArgs()[grepl("^--file=",commandArgs())][1])
root<-normalizePath(file.path(dirname(file_arg),".."))
audit<-fread(file.path(root,"inputs/species_manifest.tsv"))
for(module_set in "narrative44"){
  out<-file.path(root,module_set,"figures")
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  input<-file.path(root,module_set,"scores/OG_equal")
figure_width_mm<-190
figure_height_mm<-110
panel_width_mm<-190
panel_height_mm<-68
violin_width<-0.27
cell_offset<-0.28
violin_outline_linewidth<-0.45
iqr_bar_linewidth<-0.70
median_point_size<-1.20
axis_linewidth<-0.45
species<-c("Hcoe","Xesp","Pdam","Spin","Fspp","Pcru","Plut","Gfas","Adig","Amur","Ahya")
types<-c("calicoblast","epidermis","gastrodermis")
palette<-c(gastrodermis="wheat1",epidermis="cadetblue2",calicoblast="mediumpurple1")
labels<-expression(italic("H. coerulea"),italic(Xenia)~"sp.",italic("P. damicornis"),
  italic("S. pistillata"),italic(Fungia)~"sp.",italic("P. crustacea"),italic("P. lutea"),
  italic("G. fascicularis"),italic("A. digitifera"),italic("A. muricata"),italic("A. hyacinthus"))
dat<-rbindlist(lapply(species,function(sp)fread(file.path(input,"per_species",paste0("cell_scores.",sp,".tsv.gz")))[analysis_set=="main"]))
manifest<-fread(file.path(input,"module_manifest.tsv"))
summary<-fread(file.path(input,"cell_type_score_summary.tsv"))[analysis_set=="main"]
coverage<-fread(file.path(input,"module_coverage.tsv"))[analysis_set=="main"]
stopifnot(nrow(dat)==2L*sum(audit$n_selected_cells),all(is.finite(dat$score)),uniqueN(dat[,.(object_species,module,cell_id)])==nrow(dat),nrow(summary)==66L)
stats<-dat[,.(n_cells=.N,mean=mean(score),median=median(score),q25=as.numeric(quantile(score,.25)),
 q75=as.numeric(quantile(score,.75)),min_score=min(score),max_score=max(score)),by=.(object_species,module,cell_type)]
check<-merge(stats,summary,by=c("object_species","module","cell_type"),suffixes=c("_plot","_score"))
error<-max(abs(check$mean_plot-check$mean_score),abs(check$median_plot-check$median_score))
stopifnot(error<1e-12,all(check$n_cells_plot==check$n_cells_score))
fwrite(stats,file.path(out,"figures_plot_group_summary.tsv"),sep="\t")
bandwidths<-list()
for(mod in c("M3","M4")){
  bw<-dat[module==mod,.(bw=bw.nrd0(score)),by=.(object_species,cell_type)]
  bandwidths[[mod]]<-data.table(module=mod,bandwidth=median(bw$bw[is.finite(bw$bw)&bw$bw>0]),
     rule="median of 33 group bw.nrd0; shared within module")
}
bandwidths<-rbindlist(bandwidths)
for(layout in c("no_species_gap","with_species_gap")){
  dest<-file.path(out,layout);dir.create(file.path(dest,"panels"),recursive=TRUE,showWarnings=FALSE)
  positions<-setNames(if(layout=="no_species_gap")seq_along(species) else c(1,2,seq(3.5,11.5,by=1)),species)
  offsets<-setNames(c(-cell_offset,0,cell_offset),types)
  d<-copy(dat);s<-copy(stats)
  d[,`:=`(cell_type=factor(cell_type,levels=types),plot_x=unname(positions[object_species])+unname(offsets[cell_type]),
     violin_group=paste(object_species,cell_type,sep="__"))]
  s[,plot_x:=unname(positions[object_species])+unname(offsets[cell_type])]
  panels<-list()
  for(i in 1:2){
    mod<-c("M3","M4")[i];dd<-d[module==mod];ss<-s[module==mod]
    title<-paste0(LETTERS[i],"   ",sub("Ca2+","Ca\u00b2\u207a",manifest[module==mod,module_label],fixed=TRUE))
    p<-ggplot(dd,aes(plot_x,score,group=violin_group,fill=cell_type))+
      geom_hline(yintercept=0,color="#777777",linewidth=.25,linetype="dashed")+
      geom_violin(scale="width",trim=TRUE,bw=bandwidths[module==mod,bandwidth],width=violin_width,alpha=1,color="#555555",linewidth=violin_outline_linewidth)+
      geom_linerange(data=ss,aes(x=plot_x,ymin=q25,ymax=q75),inherit.aes=FALSE,color="#303030",linewidth=iqr_bar_linewidth)+
      geom_point(data=ss,aes(x=plot_x,y=median),inherit.aes=FALSE,color="#303030",size=median_point_size)+
      scale_fill_manual(values=palette,breaks=types,labels=c("Calicoblast","Epidermis","Gastrodermis"),drop=FALSE)+
      scale_x_continuous(breaks=unname(positions),labels=labels,
        limits=if(layout=="no_species_gap")c(.55,11.5) else c(.55,12),expand=expansion(mult=0))+
      scale_y_continuous(breaks=scales::breaks_pretty(n=4),expand=expansion(mult=c(.035,.055)))+
      labs(title=title,x=NULL,y="OG-equal module score",fill=NULL)+
      theme_classic(base_size=8,base_family="Liberation Sans")+
      theme(plot.title=element_text(size=9,face="bold",hjust=0,margin=margin(b=3)),plot.title.position="panel",
        axis.text.x=element_text(angle=45,hjust=1,vjust=1),axis.title.y=element_text(size=8,margin=margin(r=5)),
        axis.text=element_text(color="#333333",size=7.5),axis.line=element_line(linewidth=axis_linewidth),
        axis.ticks=element_line(linewidth=axis_linewidth),axis.ticks.length=grid::unit(1.2,"mm"),legend.position="bottom",
        legend.text=element_text(size=8),legend.key.width=grid::unit(4,"mm"),legend.key.height=grid::unit(3,"mm"),
        plot.margin=margin(3,5,3,3))+
      guides(fill=guide_legend(override.aes=list(color="#555555",alpha=1)))
    if(layout=="with_species_gap")p<-p+geom_vline(xintercept=2.75,color="#C7C7C7",linewidth=.25,linetype="dotted")
    ggsave(file.path(dest,"panels",paste0(mod,"_raw_score_violin.pdf")),p,width=panel_width_mm,height=panel_height_mm,units="mm",device=cairo_pdf)
    if(i<2)p<-p+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    panels[[i]]<-p
  }
  figure<-wrap_plots(panels,ncol=1,guides="collect")&theme(legend.position="bottom")
  stem<-file.path(dest,paste0("calicoblast_",module_set,"_OG_equal_",layout,"_wider_violins_bold"))
  ggsave(paste0(stem,".pdf"),figure,width=figure_width_mm,height=figure_height_mm,units="mm",device=cairo_pdf)
  ggsave(paste0(stem,".png"),figure,width=figure_width_mm,height=figure_height_mm,units="mm",dpi=600,bg="white")
  ggsave(paste0(stem,".tiff"),figure,width=figure_width_mm,height=figure_height_mm,units="mm",dpi=600,compression="lzw",bg="white")
  fwrite(bandwidths,file.path(dest,"density_bandwidths.tsv"),sep="\t")
  fwrite(stats,file.path(dest,"plot_group_summary.tsv"),sep="\t")
  fwrite(coverage,file.path(dest,"module_OG_feature_coverage.tsv"),sep="\t")
  fwrite(data.table(object_species=species,x_position=unname(positions)),file.path(dest,"species_plot_order.tsv"),sep="\t")
  fwrite(data.table(cell_type=types,R_color_name=unname(palette[types])),file.path(dest,"cell_type_palette.tsv"),sep="\t")
}


layout_parameters<-data.table(figure_width_mm=figure_width_mm,figure_height_mm=figure_height_mm,
  panel_width_mm=panel_width_mm,panel_height_mm=panel_height_mm,
  violin_width=violin_width,cell_offset=cell_offset,
  violin_outline_linewidth=violin_outline_linewidth,iqr_bar_linewidth=iqr_bar_linewidth,
  median_point_size=median_point_size,axis_linewidth=axis_linewidth,
  font_sizes_changed=FALSE,density_bandwidth_rule_changed=FALSE,
  y_axis_rule="Full observed score range; separate scale for each module")
fwrite(layout_parameters,file.path(out,"layout_parameters.tsv"),sep="\t")
writeLines(c(paste0("PASS: ",nrow(dat)," raw OG-equal cell x module scores; ",nrow(stats)," violin groups"),
 paste0("PASS: summary reproduction max error ",format(error,digits=6)),
 "PASS: M3/M4 two-row 190 x 110 mm canvas; original 190 x 68 mm standalone panels; exact final bold style",
 "PASS: shared bandwidth recalculated from updated groups using original median-of-33-group rule",
 "PASS: both species layouts; vector PDF, 600 dpi PNG/TIFF and individual module PDFs"),
 file.path(out,"VALIDATION.txt"))
message("PLOTTING COMPLETE: ",module_set,"; summary error ",error)
}

dir.create(file.path(root,"logs"),recursive=TRUE,showWarnings=FALSE)
capture.output(sessionInfo(),file=file.path(root,"logs/sessionInfo_plotting.txt"))
