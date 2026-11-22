labels_for<-function(r){
 stopifnot(all(!is.na(r$display_annotation)),all(nzchar(r$display_annotation)))
 lab<-r$display_annotation
 paste(r$orthogroup,ifelse(nchar(lab)>68,paste0(substr(lab,1,65),'...'),lab),sep=' | ')
}
render_transposed<-function(mat,states,r,c,prefix,title,overview=FALSE,width_inches,height_inches){
 nr<-nrow(mat);nc<-ncol(mat)
 stopifnot(nr==nrow(c),nc==nrow(r),identical(rownames(mat),c$column_id),identical(colnames(mat),r$orthogroup))
 OG_group<-r$primary_cell_type
 combo_group<-unname(celltype_to_group[c$cell_type])
 stopifnot(!anyNA(combo_group),all(OG_group%in%names(group_colors)),all(combo_group%in%names(group_colors)))
 source_species<-species_labels[match(c$object_species,object_species)]
 species_expr<-as.expression(lapply(source_species$figure_label,function(s){
  if(s=='Xenia sp.')return(quote(italic(Xenia)~sp.))
  if(s=='Fungia sp.')return(quote(italic(Fungia)~sp.))
  bquote(italic(.(s)))
 }))
 strip_size<-unit(layout$annotation_strip_mm,'mm')
 make_row_strip<-function(){
  rowAnnotation(Cell_type=anno_simple(combo_group,col=group_colors,border=FALSE,gp=gpar(col=layout$group_border_color,lwd=layout$group_border_lwd_pt,lty=1)),
   show_annotation_name=FALSE,show_legend=FALSE,simple_anno_size=strip_size)
 }
 top<-HeatmapAnnotation(
  conservation=anno_text(paste0(r$support,'/',r$n_species),rot=layout$horizontal_axis_text_rotation_degrees,
   gp=gpar(fontsize=if(overview)6 else 8),just='right',location=unit(0,'npc'),
   height=unit(layout$conservation_annotation_height_mm,'mm'),show_name=TRUE),
  Cell_type=anno_simple(OG_group,col=group_colors,border=FALSE,gp=gpar(col=layout$group_border_color,lwd=layout$group_border_lwd_pt,lty=1)),
  show_annotation_name=c(conservation=TRUE,Cell_type=FALSE),show_legend=FALSE,
  annotation_name_side='left',annotation_name_rot=layout$horizontal_axis_text_rotation_degrees,annotation_name_gp=gpar(fontsize=7),
  simple_anno_size=strip_size,gap=unit(.8,'mm'))
 bottom<-HeatmapAnnotation(Cell_type=anno_simple(OG_group,col=group_colors,border=FALSE,gp=gpar(col=layout$group_border_color,lwd=layout$group_border_lwd_pt,lty=1)),
  show_annotation_name=FALSE,show_legend=FALSE,simple_anno_size=strip_size)
 rowsplit<-if('plot_column_group'%in%names(c))factor(c$plot_column_group,levels=unique(c$plot_column_group)) else
  if(overview)factor(c$cell_type,levels=unique(c$cell_type)) else NULL
 row_titles<-if(is.null(rowsplit))unique(c$cell_type) else levels(rowsplit)
 colsplit<-factor(OG_group,levels=unique(OG_group))
 parent_mark<-states%in%c('up','down');dim(parent_mark)<-dim(states)
 mark<-states%in%'up';dim(mark)<-dim(states)
 mark<-mark & !is.na(mat) & mat>0
 dimnames(mark)<-dimnames(mat)
 h<-Heatmap(mat,name='OG_SUM log2FC',col=color_fun,na_col=layout$na_color,
  cluster_rows=FALSE,cluster_columns=FALSE,cluster_row_slices=FALSE,cluster_column_slices=FALSE,
  show_row_names=TRUE,row_names_side='right',row_labels=species_expr,
  row_names_gp=gpar(fontsize=if(overview)6 else 9),row_names_rot=0,row_names_max_width=unit(36,'mm'),
  show_column_names=TRUE,column_names_side='bottom',column_labels=labels_for(r),
  column_names_rot=layout$horizontal_axis_text_rotation_degrees,column_names_gp=gpar(fontsize=if(overview)6 else 8),
  column_names_max_height=unit(if(overview)115 else 125,'mm'),
  row_split=rowsplit,column_split=colsplit,row_title=row_titles,row_title_side='left',
  row_title_gp=gpar(fontsize=7),row_title_rot=0,
  column_title=unique(OG_group),column_title_side='top',column_title_rot=layout$horizontal_axis_text_rotation_degrees,column_title_gp=gpar(fontsize=7),
  top_annotation=top,bottom_annotation=bottom,left_annotation=make_row_strip(),right_annotation=make_row_strip(),
  rect_gp=gpar(col=NA),border=TRUE,border_gp=gpar(col=layout$group_border_color,lwd=layout$group_border_lwd_pt),
  row_gap=unit(layout$group_gap_mm,'mm'),column_gap=unit(layout$group_gap_mm,'mm'),use_raster=FALSE,
  heatmap_legend_param=list(at=c(-3,-2,-1,0,1,2,3),title_gp=gpar(fontsize=8),labels_gp=gpar(fontsize=7)),
  layer_fun=function(j,i,x,y,w,h,fill){
   k<-mark[cbind(i,j)]
   if(any(k))grid.points(x[k],y[k],pch=16,size=unit(if(overview)layout$marker_size_overview_mm else layout$marker_size_comparison_mm,'mm'),gp=gpar(col='#252525'))
  })
 strip_gps<-lapply(list(h@top_annotation,h@bottom_annotation,h@left_annotation,h@right_annotation),function(a)a@anno_list$Cell_type@fun@var_env$gp)
 stopifnot(all(vapply(strip_gps,function(g)identical(g$col,layout$group_border_color) && g$lwd==h@matrix_param$border_gp$lwd && g$lty==1,logical(1))))
 # ComplexHeatmap normalizes -90 to 270 for titles but retains a left anchor.
 # Anchor the end of top titles at their lower edge so all text stays above the matrix.
 h@column_title_param$just<-c(1,0.5)
 extra<-Legend(labels=c('Upregulated marker','Missing expression value'),type='points',pch=c(16,15),
  legend_gp=gpar(col=c('#252525',layout$na_color)),labels_gp=gpar(fontsize=7))
 drawit<-function(){
  draw(h,heatmap_legend_side='right',annotation_legend_side='right',heatmap_legend_list=list(extra),
   column_title=title,column_title_gp=gpar(fontsize=11),padding=unit(c(5,5,5,5),'mm'))
 }
 marker_audit_dir<-file.path(root,'results/marker_display',figure_key)
 put(data.table(column_id=rownames(mat),as.data.table(mark)),file.path(marker_audit_dir,'marker_display_mask.tsv'))
 indices<-which(mark,arr.ind=TRUE)
 put(data.table(row_index=indices[,1],column_index=indices[,2],column_id=rownames(mat)[indices[,1]],
  orthogroup=colnames(mat)[indices[,2]],plot_log2FC=mat[indices],marker_state=states[indices]),
  file.path(marker_audit_dir,'rendered_points.tsv'))
 dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
 pdf(paste0(prefix,'.pdf'),width=width_inches,height=height_inches,useDingbats=FALSE)
 drawn<-drawit()
 actual_row_order<-as.integer(unlist(row_order(drawn),use.names=FALSE))
 actual_column_order<-as.integer(unlist(column_order(drawn),use.names=FALSE))
 stopifnot(identical(actual_row_order,seq_len(nr)),identical(actual_column_order,seq_len(nc)))
 required_OG_label_height_mm<-as.numeric(convertWidth(ComplexHeatmap::max_text_width(
  h@column_names_param$labels,gp=gpar(fontsize=if(overview)6 else 8)),'mm',valueOnly=TRUE))
 allocated_OG_label_height_mm<-as.numeric(convertHeight(h@column_names_param$max_height,'mm',valueOnly=TRUE))
 stopifnot(allocated_OG_label_height_mm>=required_OG_label_height_mm)
 required_conservation_name_height_mm<-as.numeric(convertWidth(grobWidth(textGrob('conservation',gp=gpar(fontsize=7))),'mm',valueOnly=TRUE))
 allocated_conservation_height_mm<-as.numeric(convertHeight(h@top_annotation@anno_list$conservation@height,'mm',valueOnly=TRUE))
 stopifnot(allocated_conservation_height_mm>=required_conservation_name_height_mm)
 dev.off()
 if(overview){
  preview_dpi<-min(160,3000/width_inches,6500/height_inches)
  png(paste0(prefix,'.preview.png'),width=round(width_inches*preview_dpi),
   height=round(height_inches*preview_dpi),res=preview_dpi,type='cairo')
  drawit();dev.off()
 }
 ra<-copy(c)
 ra[,`:=`(plot_row_order=seq_len(.N),figure_label=source_species$figure_label,celltype_group=combo_group,
  species_render_expression=vapply(as.list(h@row_names_param$labels),function(z)paste(deparse(z),collapse=' '),character(1)),
  left_strip_color=unname(group_colors[combo_group]),right_strip_color=unname(group_colors[combo_group]))]
 put(ra,file.path(root,'results/layout_annotations',figure_key,'row_display_annotations.tsv'))
 ca<-copy(r)
 ca[,`:=`(plot_column_order=seq_len(.N),conservation=paste0(support,'/',n_species),
  top_strip_color=unname(group_colors[OG_group]),bottom_strip_color=unname(group_colors[OG_group]),
  TF_display_label=h@column_names_param$labels)]
 put(ca,file.path(root,'results/layout_annotations',figure_key,'column_display_annotations.tsv'))
 geometry<-data.table(actual_row_order_is_parent_column_order=identical(actual_row_order,seq_len(nr)),
  actual_column_order_is_parent_row_order=identical(actual_column_order,seq_len(nc)),
  new_n_rows=nr,new_n_columns=nc,marker_size_mm=if(overview)layout$marker_size_overview_mm else layout$marker_size_comparison_mm,
  annotation_cell_borders=TRUE,annotation_border_color=layout$group_border_color,annotation_border_lwd_pt=layout$group_border_lwd_pt,annotation_border_lty='solid',rendered_dot_count=sum(mark),parent_dot_count=sum(parent_mark),
  rendered_downregulated_dot_count=sum(mark & (states%in%'down')),
  all_rendered_dots_upregulated=all(states[mark]=='up'),all_rendered_dots_positive=all(mat[mark]>0),
  required_OG_label_height_mm=required_OG_label_height_mm,
  allocated_OG_label_height_mm=allocated_OG_label_height_mm,legacy_truncation_enabled=TRUE,annotation_character_limit=68L,
  row_gap_mm=as.numeric(h@matrix_param$row_gap),column_gap_mm=as.numeric(h@matrix_param$column_gap),
  border=as.logical(h@matrix_param$border),border_lwd_pt=h@matrix_param$border_gp$lwd,
  internal_grid_removed=is.na(h@matrix_param$gp$col),
  top_annotation_names=paste(names(h@top_annotation@anno_list),collapse=';'),
  conservation_name_displayed=isTRUE(h@top_annotation@anno_list$conservation@name_param$show),
  conservation_text_rot=h@top_annotation@anno_list$conservation@fun@var_env$rot,
  conservation_name_rot=h@top_annotation@anno_list$conservation@name_param$rot,
  required_conservation_name_height_mm=required_conservation_name_height_mm,
  allocated_conservation_height_mm=allocated_conservation_height_mm,
  bottom_annotation_names=paste(names(h@bottom_annotation@anno_list),collapse=';'),
  left_annotation_names=paste(names(h@left_annotation@anno_list),collapse=';'),
  right_annotation_names=paste(names(h@right_annotation@anno_list),collapse=';'),
  row_names_side=h@row_names_param$side,column_names_side=h@column_names_param$side,
  column_names_rot=h@column_names_param$rot,column_title_rot=h@column_title_param$rot,row_title_rot=h@row_title_param$rot,
  bottom_OG_text_just=h@column_names_param$anno@var_env$just,
  top_group_title_horizontal_just=h@column_title_param$just[1],
  row_cluster_enabled=FALSE,column_cluster_enabled=FALSE,axes_reversed=FALSE)
 put(geometry,file.path(root,'results/layout_annotations',figure_key,'render_geometry.tsv'))
 figures[[length(figures)+1]]<<-data.table(figure=paste0(prefix,'.pdf'),
  n_rows=nr,n_columns=nc,n_OGs=ncol(mat),n_species_celltype_combinations=nrow(mat),missing_cells=sum(is.na(mat)),marker_cells=sum(mark),parent_marker_cells=sum(parent_mark),removed_down_marker_cells=sum(parent_mark)-sum(mark),
  width_inches=width_inches,height_inches=height_inches,scale_limit=3,row_zscore=FALSE,axes_transposed=TRUE,
  horizontal_group_gap_mm=layout$group_gap_mm,vertical_group_gap_mm=layout$group_gap_mm,
  border_lwd_pt=layout$group_border_lwd_pt,annotation_strip_mm=layout$annotation_strip_mm,
  top_celltype_strip=TRUE,bottom_celltype_strip=TRUE,left_celltype_strip=TRUE,right_celltype_strip=TRUE,
  species_color_annotation=FALSE,conservation_location='top')
}
