################################################################################
# formatting
################################################################################
# shorten colnmaes
shorten_names<-function(list_in){
  shortened_names=sub("", '',list_in)
  return(shortened_names)
}

################################################################################
# DE
################################################################################
deg_comparison<-function(cntrl_in,treat_in){
  contras=c(treat_in,cntrl_in)
  
  # subset sample info for contrasts
  full_contrast=paste0(cntrl_in,"-",treat_in)
  sub_df=subset(groups_df,group %in% contras)
  sampleinfo= data.frame(sub_df$group,row.names = sub_df$sampleid)
  colnames(sampleinfo)=c("condition")
  sampleinfo$condition=as.factor(sampleinfo$condition)
  
  # read in Raw count file
  # example: DEG_KO-CRISPR_53_without_IFNb_0.5_0.5/RawCountFile_RSEM_genes_filtered.txt
  fpath=paste0(input_dir,"DEG_",treat_in,"-",
               cntrl_in,"_0.5_0.5/RawCountFile_RSEM_genes_filtered.txt")
  x=read.csv(fpath,sep="\t")
  rownames(x)=x$symbol
  x=x[,c(2:ncol(x))]
  
  # run DESEQ
  ddsHTSeq<-DESeqDataSetFromMatrix(countData=x,colData=sampleinfo, design=~condition)
  dds<-DESeq(ddsHTSeq)
  
  # prep df
  mfc=c()
  mpval=c()
  i=1
  
  # pull results for conditions
  res<-results(dds,contrast=c("condition",as.character(contras[1]),as.character(contras[2])))
  res1=as.data.frame(res)
  write.table(res1,
              file=paste(output_dir,"DESeq2_",contras[1],"-",
                         contras[2],"_DEG_allgenes_res1.txt",
                         sep=""),
              sep="\t",col.names=NA) 
  restmp=res1
  
  # calc fc, create inital df
  restmp$FoldChange <- ifelse(restmp$log2FoldChange<0, -1/(2^restmp$log2FoldChange), 2^restmp$log2FoldChange)
  mfc=cbind(mfc,restmp$FoldChange)
  mpval=cbind(mpval,restmp$pvalue)
  
  x=rownames(restmp)
  ensID=apply(array(as.character(x)),1,function(z) unlist(strsplit(z, "\\|"))[1])
  gene=apply(array(as.character(x)),1,function(z) unlist(strsplit(z, "\\|"))[2])
  restmp=cbind(ensID,gene,restmp)
  
  #remove rownames for df merging downstream
  deseq2out=restmp
  deseq2out$X=rownames(deseq2out)
  
  #subselect df, recalc new values
  deseq2out=deseq2out[,which(names(deseq2out) %in% c("X", "gene","log2FoldChange","pvalue"))]
  deseq2out$fc=2^deseq2out$log2FoldChange
  down_reg=deseq2out$log2FoldChange<0
  deseq2out$fc[down_reg]=-1/deseq2out$fc[down_reg]
  
  # pull final values
  deseq2out=deseq2out[,c("X","gene","fc","log2FoldChange","pvalue")]
  colnames(deseq2out)=c("ensid_gene","gene","fc","log2fc","pvalue")
  deseq2out$fdr=p.adjust(deseq2out$pvalue,method='fdr',n=length(deseq2out$pvalue))
  deseq2out$gsea_ranking_score=-log10(deseq2out$pvalue)*sign(deseq2out$log2fc)
  write.table(deseq2out,file=paste(output_dir,"DESeq2_",contras[1],"-",
                                   contras[2],"_DEG_allgenes.txt",sep=""),
              row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  return(deseq2out)
}

deg_group_add<-function(cntrl_in,treat_in){
  # run DESE2
  res_comp<-deg_comparison(cntrl_in,treat_in)
  
  # rename cols with sampleid
  tmp_deg=as.data.frame(res_comp)
  colnames(tmp_deg)=paste0(treat_in,"_",colnames(tmp_deg))
}

################################################################################
# heatmaps
################################################################################
# creates gene df of top n_in significant genes per sample
create_sig_gene_df<-function(cntrl_in,treat_in,n_in){
  #set contrasat
  contras=c(treat_in,cntrl_in)
  
  #read in path
  source_path=paste0(input_dir,"DEG_",contras[1],"-",contras[2],"_0.5_0.5/")
  res1=read.csv(paste0(output_dir,"DESeq2_",contras[1],"-", contras[2],"_DEG_allgenes_res1.txt"),sep="\t")
  
  # subset for sig genes
  sub_res1=subset(res1,(log2FoldChange>log2fc_cutoff) | (log2FoldChange<log2fc_cutoff))
  sub_res1=subset(sub_res1,padj<padj_cutoff)
  
  # subset for top X genes
  sub_top_fc=sub_res1[order(sub_res1$padj),][c(1:n_in),]
  
  #select cols, rename
  sub_top_fc=sub_top_fc[,c("X","log2FoldChange","padj")]
  colnames(sub_top_fc)=c("genes",
                         paste0(contras[1],"_vs_",contras[2],"--log2FoldChange"),
                         paste0(contras[1],"_vs_",contras[2],"--padj"))
  
  return(sub_top_fc)
}

# fills in df for other samples
fillin_sig_gene_df<-function(df_in){
  # add rownames, subset
  rownames(df_in)=df_in$genes
  df_in <- subset(df_in, select = -c(genes))
  
  # for each column, search through and fill in any NA's
  for (colid in colnames(df_in)){
    # pull all NA values for this column
    df_in[df_in=='NA'] <- NA
    missing_rows=rownames(df_in[is.na(df_in[,colid]),])
    
    #read in path
    contras=strsplit(strsplit(colid,"--")[[1]],"_vs_")
    source_path=paste0(input_dir,"DEG_",contras[[1]][1],"-",contras[[1]][2],"_0.5_0.5/")
    res1=read.csv(paste0(output_dir,"DESeq2_",
                         contras[[1]][1],"-", contras[[1]][2],"_DEG_allgenes_res1.txt"),sep="\t")
    rownames(res1)=res1$X
    
    # fill in missing information
    for (rowid in missing_rows){
      df_in[rowid,colid]=res1[rowid,"log2FoldChange"]
    }
  }
  
  return(df_in)
}

# creates heatmap formatted df with only log2fc values
create_heatmap_df<-function(df_in,extra_filter=""){
  df_out=dplyr::select(df_in,contains("log"))
  
  # cleanup col names
  output_list=sub('--log2FoldChange', '',colnames(df_out))
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    output_list=sub(ct[2], '',output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # use extra filter
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_out)=output_list
  return(df_out)
}

# creates df of log2fc and pvalues
create_output_df<-function(df_in,n_in,extra_filter=""){
  output_list=colnames(df_in)
  
  # split treat from control
  cntrl=strsplit(output_list,"_vs_")
  for (ct in cntrl){
    metric=strsplit(ct[2],"--")[[1]][2]
    output_list=sub(ct[2], paste0("_",metric),output_list)
    output_list=sub('_vs_', '',output_list)
  }
  
  # cleanup col names
  output_list=sub('log2FoldChange', 'log2FC',output_list)
  
  # use extra filter
  extra_filter="_without_IFNb"
  if (extra_filter != ""){
    output_list=sub(extra_filter, '',output_list)
  }
  
  colnames(df_in)=output_list
  
  # round df
  for (colid in colnames(df_in)){
    df_in[,colid]=signif(df_in[,colid], digits=3)
  }
  
  # if number of pathways is less than n_in, adjust
  if (nrow(df_in)<n_in){
    n_in=nrow(df_in)
  }
  
  # print out table
  caption_title=paste0("Expression values for ",n_in," genes")
  p = DT::datatable(df_in, extensions = 'Responsive', 
                    caption=htmltools::tags$caption(paste0(caption_title) ,
                                                    style="color:gray; font-size: 18px" ))
  return(p)
}

# Overwrites the pheatmap defaults
draw_colnames_45 <- function (coln, gaps, ...) {
  "Overwrites body of pheatmap:::draw_colnames, customizing it my liking"
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)
}

# creates heatmap
generate_heat_map<-function(df_in,show_names="ON"){
  
  ####################
  # formatting
  #####################
  # Overwrite pheatmaps default draw_colnames with new version
  assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap")) 
  
  # Heatmap Color Gradients 
  paletteLength <- 1000
  mycolors <- colorRampPalette(c("blue","white","red"), interpolate = "linear")(paletteLength)
  
  ####################
  # metadata
  ####################
  # Creating Dataframe to map samplenames to groups
  meta = groups_df
  groups <- data.frame(as.factor(meta$group))
  colnames(groups) <- "Groups"
  rownames(groups) <- meta$sampleid
  
  # Creating Group Column Annotation Colors
  columnColors <- c("lightpink","lightblue","orange","purple","red","green")
  names(columnColors) <- unique(groups$Groups)
  anno_colors <- list(Groups = columnColors)
  
  # set title
  title_in=paste0("Significant Genes (N=",nrow(df_in),")")
  
  ####################
  # function
  ####################
  if (show_names=="OFF"){
    pheatmap(df_in, 
             scale = "none", main=title_in,
             cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
             border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = FALSE)
  } else{
    pheatmap(df_in, 
             scale = "none", main=title_in,
             cellwidth = 30, fontsize = 12, fontsize_row = 7, fontsize_col = 8, color = mycolors, 
             border_color = "NA",cluster_cols=F,annotation_colors = anno_colors, show_rownames = TRUE)
  }
  
}

################################################################################
# volcano plots
################################################################################
generate_volcano_plots<-function(cntrl_in,treat_in,type_in){
  # read in res from DEG merge
  contras=c(treat_in,cntrl_in)
  fpath=paste0(output_dir,"DESeq2_",contras[1],"-", contras[2],"_DEG_allgenes_res1.txt")
  res1=read.csv(fpath,sep="\t")
 
  
  # Volcano Plots
  if (type_in=="pvalue"){
    log_pval=-log10(res1$pvalue)
    y_title="-Log10 pvalue"
  } else{
    ## logfc and FDR
    log_pval=-log10(res1$padj)
    y_title="-Log10 FDR"
  }
  log_FC=res1$log2FoldChange
  Significant=rep("1_NotSignificant",length(log_FC))
  Significant[which(res1$pvalue<padj_cutoff & abs(res1$log2FoldChange)>=log2fc_cutoff)]=paste0("3_LogFC_and_",type_in)
  Significant[which(res1$pvalue<padj_cutoff & abs(res1$log2FoldChange)<log2fc_cutoff)]=paste0("2b_",type_in,"_Only")
  Significant[which(res1$pvalue>=padj_cutoff & abs(res1$log2FoldChange)>=log2fc_cutoff)]="2a_LogFC_Only"
  gene=res1$X
  volcano_data=as.data.frame(cbind(gene,log_FC,log_pval,Significant))
  p <- plot_ly(data = volcano_data, x = log_FC, y = log_pval, text = gene,
               mode = "markers", 
               color = Significant) %>% layout(title =paste0(contras[1]," vs. ", contras[2]),
                                               xaxis=list(title="Fold Change",
                                                          range =c(-5,5),
                                                          tickvals=c(-5,-4,-3,-2,-1,0,1,2,3,4,5),
                                                          ticktext=c('-32','-16','-8','-4','-2',
                                                                     '1','2','4','8','16','32')),
                                               yaxis=list(title=y_title,range =c(0,15)))
  return(p)
}

################################################################################
# functions ORA/GSEA
################################################################################
capture_entrezids<-function(input_df){
  sep_df=input_df %>%
    separate(ensid_gene,sep="[|]",c("ENSEMBL","SYMBOL"))%>%
    separate(ENSEMBL,sep="[.]",c("ENSEMBL","ID"))
  
  # search for ENTREZ by ENSEMBL
  gene_df_e <- bitr(sep_df[,c("ENSEMBL")], fromType = "ENSEMBL",
                    toType = c("ENTREZID"),
                    OrgDb = org.Hs.eg.db)
  
  # if any genes did not map, try by SYMBOL
  missing_df=subset(sep_df,ENSEMBL %ni% gene_df_e$ENSEMBL)
  if (nrow(missing_df)>0){
    gene_df_s <- bitr(sep_df[,c("SYMBOL")], fromType = "SYMBOL",
                      toType = c("ENTREZID"),
                      OrgDb = org.Hs.eg.db)
    
    # merge dfs together
    gene_output=merge(gene_df_e,gene_df_s,all=TRUE,by="ENTREZID")
    
  }
  
  # remove duplicated ENSEMBL ID's that have a gene symbol
  dup_list=gene_output[duplicated(gene_output$ENSEMBL),]
  filter_list=dup_list[is.na(dup_list$SYMBOL),]$ENTREZID
  gene_output2=subset(gene_output, ENTREZID %ni% filter_list)
  
  # fill in missing gene symbols
  filter_list=gene_output2[is.na(gene_output2$SYMBOL),]
  for (i in rownames(filter_list)){
    eid=filter_list[i,"ENSEMBL"]
    #pull symbol from deg
    gene_output2[i,"SYMBOL"]=unique(subset(sep_df,ENSEMBL==eid)$SYMBOL)
  }
  
  # fill in missing ENSEMBL
  filter_list=gene_output2[is.na(gene_output2$ENSEMBL),]
  for (i in rownames(filter_list)){
    sym=filter_list[i,"SYMBOL"]
    #pull symbol from deg
    gene_output2[i,"ENSEMBL"]=subset(sep_df,SYMBOL==sym)$ENSEMBL
  }
  
  #remove dups
  output_df=gene_output2[!duplicated(gene_output2[,c("ENSEMBL")]),]
  
  # rename cols to match deg
  colnames(output_df)=c("ENTREZID","ENSEMBL","gene")
  return(output_df)
}

#set genelist for GSEA
deg2geneList<-function(deg,t2g){
  # create refs of entrez:genes
  gene_ref_db=capture_entrezids(deg)
  
  #add entrezs to deg df
  deg_anno_df=merge(gene_ref_db,deg,by="gene")
  
  # create genelist
  gsea_genelist=deg_anno_df$log2fc
  
  if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    names(gsea_genelist)=as.character(deg_anno_df$gene)
  } else{
    names(gsea_genelist)=as.character(deg_anno_df$ENTREZID)
  }
  gsea_genelist=sort(gsea_genelist,decreasing=TRUE)
  
  return(gsea_genelist)
}

# set the annotation dbs
db_lookup<-function(t2g){
  # generate gene lists for C1
  # generate gene lists for C2 with subtypes biocarta, kegg, reactome, wiki
  # generate gene lists for C5 with subtypes MF, BP, CC
  # generate gene lists for Hallmark
  if (t2g=="C1"){
    db_out=msigdbr(species = species_in, category = "C1") %>% 
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C2:BIOCARTA"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "BIOCARTA") %>% 
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C2:KEGG"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "KEGG") %>% 
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C2:REACTOME"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "REACTOME") %>%
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C2:WIKIPATHWAYS"){
    db_out=msigdbr(species = species_in, category = "C2", subcategory = "WIKIPATHWAYS") %>%
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C5:MF"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:MF") %>%
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C5:BP"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:BP") %>%
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="C5:CC"){
    db_out=msigdbr(species = species_in,  category = "C5", subcategory = "GO:CC") %>%
      dplyr::select(gs_name,gene_symbol)
  } else if (t2g=="H"){
    db_out=msigdbr(species = species_in, category = "H") %>% 
      dplyr::select(gs_name,gene_symbol)  
  } else{
    print("DB does not exist. Please review")
  }
  
  return(db_out)
}

# create a small legend
addSmallLegend <- function(myPlot, pointSize = 2, textSize = 5, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# plot ORA
ora_plus_plot <- function(gl,t2g,contrast_in,n_show=3){
  # pull the DB
  pulled_db=db_lookup(t2g)
  
  # run ORA
  result=enricher(gene=gl, TERM2GENE=pulled_db, pvalueCutoff = padj_cutoff)
  resultdf=as.data.frame(result)
  
  # write out pathways file
  ttl_abbrev=sub(" ","_",sub(":","_",t2g))
  fpath=paste0(output_dir,"ORA_",contrast_in[1],"-",contrast_in[2],"_table_",ttl_abbrev,".txt")
  write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
  
  # create dotplot if pathways are sig
  if(nrow(resultdf)==0){
    p1 = ggparagraph( paste0("\n\n\n No Sig Results for ", "ORA:",t2g,"\n-",
                             contrast_in[1],"-",contrast_in[2]), 
                      color = NULL, size = 20, face = "bold", 
                      family = NULL, lineheight = NULL)
  } else{
    p1 = dotplot(result,
                 title=paste0(contrast_in[1],"-",contrast_in[2],"\nORA:",t2g),
                 font.size = 6, showCategory=n_show)
  }
  # add small legend
  pf = addSmallLegend(p1)
  return(pf)
}

# plot GSEA
gsea_plus_plot <- function(gl,t2g,contrast_in,select_flag="OFF"){
  # shorthand species
  if (species_in=="Homo sapiens"){
    species_short="hsa"
    org_db="org.Hs.eg.db"
  } else if (species_in=="Mus musculus"){
    species_short="mmu"
    org_db="org.Mm.eg.db"
  } else{
    print (paste0(species_in,": species not found"))
  }
  
  # run GSEA
  if (t2g=="C2:KEGG"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    result=gseKEGG(geneList=gl,
                   pvalueCutoff=padj_cutoff,
                   eps=0,
                   pAdjustMethod="BH", 
                   organism=species_short,
                   verbose=FALSE)
    
  } else if (t2g=="C2:REACTOME"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    result=gsePathway(gene=gl, 
                      pvalueCutoff = padj_cutoff,
                      eps=0,
                      pAdjustMethod = "BH", 
                      verbose = FALSE)
    
  } else if (t2g=="C2:WIKIPATHWAYS"){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html
    result=gseWP(gene=gl, 
                 pvalueCutoff = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 organism=species_in)
    
  } else if ((t2g=="C5:MF") | (t2g=="C5:BP") | (t2g=="C5:CC")){
    #https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
    ont_id=strsplit(t2g,":")[[1]][2]
    result=gseGO(geneList=gl,
                 pvalueCutof = padj_cutoff,
                 eps=0,
                 pAdjustMethod = "BH",
                 OrgDb= get(org_db),
                 ont=ont_id,
                 verbose= FALSE)
    
  } else if ((t2g=="C1") | (t2g=="C2:BIOCARTA") | (t2g=="H")){
    pulled_db=db_lookup(t2g)
    result=GSEA(geneList=gl,
                pvalueCutoff = padj_cutoff,
                eps=0,
                pAdjustMethod = "BH",
                TERM2GENE = pulled_db)
  } else {
    print (paste0(t2g, ": DB selected is not valid"))
  }
  
  # breakpoint for the select_pathway function versus standard analysis
  if (select_flag=="ON"){ 
    return(result)
  } else{
    # save datatable
    resultdf=as.data.frame(result)
    ttl_abbrev=sub(" ","_",sub(":","_",t2g))
    fpath=paste0(output_dir,"GSEA_",contrast_in[1],"-",contrast_in[2],"_table_",ttl_abbrev,".txt")
    write.table(resultdf,file=fpath,quote=FALSE,sep="\t",row.names = FALSE,col.names = TRUE)
    
    # create dot plots for all DB, ridgeplots for specific DB's
    if(nrow(result)==0){
      pf = ggparagraph( paste0("\n\n\n No Sig Results for GSEA:",t2g,
                               "\n-",contrast_in[1]), 
                        size = 20, face = "bold")
    } else{
      p1 = dotplot(result,
                   title=paste0(contrast_in,"\nGSEA:",t2g),
                   font.size = 6, showCategory=2, split=".sign",orderBy="p.adjust") +
        facet_grid(.~.sign)
      p2 = ridgeplot(result, label_format = 30, showCategory = 4, orderBy="p.adjust") +
        labs(x = "Enrichment distribution for top 5 pathways") + 
        theme(text = element_text(size=6),
              axis.text.x = element_text(size=6),
              axis.text.y = element_text(size=5.5))
      pcol <- cowplot::plot_grid(
        p1 + theme(legend.position="none"),
        p2 + theme(legend.position="none"),
        nrow = 2
      )
      legend<-get_legend(p1)
      pf=cowplot::plot_grid(pcol,legend,rel_widths = c(3, .4))
    }
    return(pf)
  }
}

# save plots
print_save_plots<-function(plot_list,contrast_id,type_in){
  # set figure caption letter
  p=1
  
  # ORA plots need to be merged two plots per figure
  if (type_in == "ORA"){
    for (i in seq(from=1,to=length(plot_list),by=2)){
      #if it's the final image, no merging needed
      if ((i+1)>length(plot_list)){
        pf=cowplot::plot_grid(plot_list[[i]],
                              ncol=1, labels=LETTERS[p])
      } else{
        p1=plot_list[[i]]
        p2=plot_list[[i+1]]
        pf=cowplot::plot_grid(p1,p2,ncol=1,
                              labels=LETTERS[p])
      }
      print(pf)
      ggsave(filename = paste0(output_dir, type_in,"_",
                               contrast_id[1],"-",contrast_id[2], 
                               "_dotplot_",LETTERS[p],".png"),
             height = 8.90, width = 12.80, device = "png", plot = pf)
      # increase counter
      p=p+1
    }
  } else{
    # GSEA plots are already merged, one plot into two, just print
    for (i in seq(from=1,to=length(plot_list))){
      pf = cowplot::plot_grid(plot_list[[i]],
                              ncol=1, labels=LETTERS[p])
      
      #print and save
      print(pf)
      ggsave(filename = paste0(output_dir, type_in,"_",
                               contrast_id[1],"-",contrast_id[2], 
                               "_dotplot_",LETTERS[p],".png"),
             height = 8.90, width = 12.80, device = "png", plot = pf)
      # increase counter
      p=p+1
    }
  }
}

# create output dt that summarized pathways
create_dts<-function(type_in,t2g,contras,n_in,db_list){
  
  # read in datatable created during GSEA/ORA plotting
  merged_df=data.frame()
  for (t2g in db_list){
    ttl_abbrev=sub(" ","_",sub(":","_",t2g))
    fpath=paste0(output_dir,type_in,"_",contras[1],"-",contras[2],"_table_",ttl_abbrev,".txt")
    tmp_df=read.csv(fpath,sep="\t")
    
    # if there are no signifcant pathways, skip
    if (nrow(tmp_df)==0){ next }
    tmp_df$anno=t2g
    
    # merge
    merged_df=rbind(merged_df,tmp_df)
  }
  
  # split leading_edge col into three
  if (nrow(merged_df)>0){
    if (type_in=="GSEA"){
      merged_df=separate(
        merged_df,
        leading_edge,
        c("percent_included","percent_list","percent_signal"),
        sep=", "
      )
      # remove value= in cols
      merged_df$percent_included=sub("tags=","",merged_df$percent_included)
      merged_df$percent_list=sub("list=","",merged_df$percent_list)
      merged_df$percent_signal=sub("signal=","",merged_df$percent_signal)
      
      # round cols
      col_list=c("p.adjust","enrichmentScore","NES")
      for (colid in col_list){
        merged_df[,colid]=signif(merged_df[,colid], digits=3)
      }
      
      # select cols
      output_df=merged_df[,c("anno","ID","Description","setSize","percent_included","p.adjust",
                             "enrichmentScore","NES","core_enrichment")]
      colnames(output_df)=c("anno_db","ID","Desc","total_genes","percent_included","p.adj",
                            "enrichmentScore","NES","genes_included")
    } else{
      # pull values for genes included, genes in set
      merged_df=separate(
        merged_df,
        BgRatio,
        c("set_total","n_bgratio"),
        sep="/"
      )
      # calculate percent genes included
      merged_df$set_total=as.numeric(merged_df$set_total)
      merged_df$percent_included=(merged_df$Count/merged_df$set_total)*100
      
      #round cols
      col_list=c("pvalue","p.adjust","qvalue","percent_included")
      for (colid in col_list){
        merged_df[,colid]=signif(merged_df[,colid], digits=3)
      }
      merged_df$percent_included=paste0(merged_df$percent_included,"%")
      
      # select cols
      output_df=merged_df[,c("anno","ID","set_total","percent_included","p.adjust","geneID")]
      colnames(output_df)=c("anno_db","ID","total_genes","percent_included","p.adj","genes_included")
    }
    # sort by p.adjust
    output_df=output_df[order(output_df$p.adj),]
    
    # if number of pathways is less than n_in, adjust
    if (nrow(output_df)<n_in){
      n_in=nrow(output_df)
    }
    
    # create DT
    caption_title=paste0("Top ", n_in, " Pathways for all annotation databases (", type_in, ")\n")
    p <- DT::datatable(output_df[1:n_in,], extensions = 'Responsive', 
                       caption=htmltools::tags$caption(paste0(caption_title, contras[1],"_vs_", contras[2]) ,
                                                       style="color:gray; font-size: 18px" ),rownames=F)
  } else{
    p="No Significant Genes"
  }
  return(p)
}

# main function
main_gsea_ora_function<-function(cntrl_in,treat_in,db_list,top_path_value,ORA_flag="",GSEA_flag=""){
  
  # set contrast
  contras=c(treat_in,cntrl_in)
  
  # read in deg
  deg_file = paste0(input_dir, "DEG_", contras[1],"-",contras[2],"_0.5_0.5/",analysis_type,
                    "_DEG_",contras[1],"-",contras[2],"_all_genes.txt")
  deg=read.csv(deg_file,header=TRUE,sep="\t")
  
  # run ORA
  o=list()
  if (ORA_flag=="ON"){
    # create ORA genelist
    siggenes=subset(deg,fdr <= padj_cutoff)
    siggenes=subset(deg,(fc < -log2fc_cutoff) | (fc > log2fc_cutoff))
    sigGeneList=siggenes$gene
    
    # for each annotation db, run ORA, save plots
    i=1
    for (db_id in db_list){
      o[[i]]=ora_plus_plot(gl=sigGeneList,t2g=db_id,contrast_in=contras)
      i=i+1
    }
    
    # print,save plots
    print_save_plots(o,contras,"ORA")
    
    # print DT
    create_dts("ORA",t2g,contras,top_path_value,db_list)
  }
  
  # run GSEA
  g=list()
  if (GSEA_flag=="ON"){
    
    # for each annotation db, create GSEA gene list - must be done within loop
    # to handle differences in naming list 
    i=1
    for (db_id in db_list){
      # create GSEA genelist
      gsea_genelist=deg2geneList(deg,t2g=db_id)
      
      # run GSEA, save plots
      g[[i]]=gsea_plus_plot(gl=gsea_genelist,t2g=db_id,contrast_in=contras)
      i=i+1
    }
    
    # print,save plots
    print_save_plots(g,contras,"GSEA")
    
    # print DT
    create_dts("GSEA",t2g,contras,top_path_value,db_list)
  }
}

######################################################################
# functions secondary pathway analysis
######################################################################
#create heatmap and DT
generate_heat_map_select<-function(select_deg,contras){
  # prep df for heatmap generation
  rownames(select_deg)=select_deg$ensid_gene
  log_name=paste0(contras[1],"_vs_",contras[2],"--log2FoldChange")
  p_name=paste0(contras[1],"_vs_",contras[2],"--padj")
  select_deg=select_deg[,c("log2fc","fdr")]
  colnames(select_deg)=c(log_name,p_name)
  
  # create heatmap df of sig genes in pathway list
  sig_df=subset(select_deg, get(log_name) > log2fc_cutoff | get(log_name) < -log2fc_cutoff)
  sig_df=subset(sig_df, get(p_name) < padj_cutoff )
  heat_df=create_heatmap_df(sig_df,"")
  
  # run functions
  generate_heat_map(heat_df)
  p = create_output_df(select_deg,nrow(select_deg),"")
  return(p)
}

# create gseaplot
gsea_plus_plots_select<-function(deg,t2g,path_id,contras){
  
  # create GSEA genelist
  gsea_genelist=deg2geneList(deg,t2g=t2g)
  
  # run GSEA, save plots
  result=gsea_plus_plot(gl=gsea_genelist,t2g=t2g,contrast_in=contras,select_flag="ON")
  
  #get rowname of pathway
  result_df=as.data.frame(result)
  rownames(result_df) <- NULL
  row_id=as.numeric(rownames(subset(result_df,ID==path_id))[[1]])
  path_desc=subset(result_df,ID==path_id)$Description[[1]]
  
  # create plot
  p1=gseaplot(result, by = "all", title = paste0(path_id,": ",path_desc),
              geneSetID =row_id)
  pf=p1&theme(text = element_text(size=8),
              axis.text.x = element_text(size=8),
              axis.text.y = element_text(size=8),
              axis.text =element_text(size=8),
              axis.title=element_text(size=8,face="bold")) 
  print(pf)
}

main_selectpath_function<-function(cntrl_in,treat_in,type_in,t2g,path_id){
  t2g="C1"
  type_in="GSEA"
  path_id="chr18q21"
  cntrl_in=cntrl
  treat_in=treatment
  
  # set contrast
  contras=c(treat_in,cntrl_in)
  
  # read in datatable created during GSEA/ORA plotting, get list of gene names in pathway selected
  ttl_abbrev=sub(" ","_",sub(":","_",t2g))
  fpath=paste0(output_dir,type_in,"_",contras[1],"-",contras[2],"_table_",ttl_abbrev,".txt")
  path_df=read.csv(fpath,sep="\t")
  
  #check pathway exists
  if (nrow(subset(path_df,ID==path_id))==0){
    stop(paste0("The selected pathway (",path_id,") does not exist in the annotation database (",t2g,
                "). Please select a valid combination"))
  }
  
  # create gene list from pathway
  genes_in_pathway=strsplit(subset(path_df,ID==path_id)$core_enrichment,"/")[[1]]
  
  # read in created deg
  fpath=paste0(output_dir,"DESeq2_",contras[1],"-",contras[2],"_DEG_allgenes.txt")
  deg=read.csv(fpath,header=TRUE,sep="\t")
  
  # subset deg for genes,
  #convert ENTREZID if necessary
  special_dbs=c("C1","C2:BIOCARTA","H")
  if (t2g %ni% special_dbs){
    gene_ref_db = capture_entrezids(deg)
    genes_in_pathway=subset(gene_ref_db,ENTREZID %in% genes_in_pathway)$gene
  } 
  select_deg=subset(deg, gene %in% genes_in_pathway)
  
  # create heatmaps
  p = generate_heat_map_select(select_deg,contras)
  
  # create gseaPlot
  gsea_plus_plots_select(deg,t2g,path_id,contras)
  
  # return DT
  return(p)
}