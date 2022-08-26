

library(ggplot2)
library(dplyr)
library(DESeq2)


# make deseq dataset dds for nuclear fractionation


make_deseq <- function(cell="A549",
                       compare = c("IS","IP"),
                       row.names="gene_id",
                       data=gene_data,
                       countcut =1,
                       sampleNcut=2,
                       ddsName="dds/A549_I.dds")
{
  selected_sample <- lib_info %>% dplyr::filter(CellType == cell & condition %in% compare ) # should take care for column information
  
  countData <- lapply(selected_sample$Sample, function(x){
    data.x <- subset.data.frame(x = data,select  = colnames(data) ==x)
  }) %>% bind_cols() %>% as.data.frame 
  row.names(countData) <- data[[row.names]]
  
  countData <-na.omit(countData)
  keep <- rowSums(edgeR::cpm(countData)>countcut) >= sampleNcut
  
  countData.kept <- data.frame(countData[keep,],row.names = row.names(countData)[keep])
  
  colData <- data.frame(BC=selected_sample$Sample,condition=selected_sample$condition)
  dds <- DESeqDataSetFromMatrix(countData = countData.kept,
                                colData   = colData,
                                design    = ~condition) %>% DESeq
  
  rownames(dds) <- row.names(countData.kept) 
  saveRDS(dds,file = ddsName)
  
}

# make deseq results from deseq dds
results_deseq <- function(dds="dds/A549_I.dds",contrast=c("condition","IP","IS")){
  library(DESeq2)
  deseq_d <- results(object = readRDS(dds),contrast = contrast)
  res <- deseq_d@listData %>%  as.data.frame(row.names = row.names(deseq_d))
}



get_DESeq_results <- function(dds = "dds/A549_I.dds",contrast =c("condition","IP","IS"),padj_cut = 0.05,log2FoldChange_cut=1,outdir="results/deseq/",name="A549_I.gz"){
  
  results <- results_deseq(dds = dds,contrast =contrast ) %>% tibble::rownames_to_column(var = "gene_symbol") %>% mutate(direction = ifelse(padj< padj_cut, ifelse(log2FoldChange> log2FoldChange_cut,"Chr", ifelse(log2FoldChange < 0-log2FoldChange_cut,"Soluble","NC")),"NC" ))

  data.table::fwrite(results,file = paste0(outdir,"DESeqResult_",name))
  return(results)
  
  
}


# useful functions
hash <- function(key,value){
  hash <- value
  names(hash) <- key
  return( hash)
}
parse_ratio <- function(x) {
  y= as.numeric(sub("/\\d+", "", x))/as.numeric(sub("\\d+/", "", x))
} 

Venn.diagram2<- function (x, filename, height = 12, width = 12, resolution = 500, 
                          imagetype = "pdf", units = "px", compression = "lzw", na = "stop", 
                          main = NULL, sub = NULL, main.pos = c(0.5, 1.05), main.fontface = "plain", 
                          main.fontfamily = "serif", main.col = "black", main.cex = 1, 
                          main.just = c(0.5, 1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                          sub.fontfamily = "serif", sub.col = "black", sub.cex = 1, 
                          sub.just = c(0.5, 1), category.names = names(x), force.unique = TRUE, 
                          print.mode = "raw", sigdigs = 3, direct.area = FALSE, area.vector = 0, 
                          hyper.test = FALSE, total.population = NULL, lower.tail = TRUE, 
                          ...) 
{
  library(VennDiagram)
  time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
  
  out.list = as.list(sys.call())
  out.list[[1]] <- NULL
  out.string = capture.output(out.list)
  flog.info(out.string, name = "VennDiagramLogger")
  if (direct.area) {
    if (1 == length(area.vector)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = area.vector[1], 
                                                 category = list.names, ind = FALSE, ...)
    }
    if (3 == length(area.vector)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = area.vector[1], 
                                                   area2 = area.vector[2], cross.area = area.vector[3], 
                                                   category = category.names, ind = FALSE, print.mode = print.mode, 
                                                   sigdigs = sigdigs, ...)
    }
    if (7 == length(area.vector)) {
      grob.list <- VennDiagram::draw.triple.venn(area1 = 0, 
                                                 area2 = 0, area3 = 0, n12 = 0, n23 = 0, n13 = 0, 
                                                 n123 = 0, category = category.names, ind = FALSE, 
                                                 list.order = 1:3, print.mode = print.mode, sigdigs = sigdigs, 
                                                 area.vector = area.vector, direct.area = TRUE, 
                                                 ...)
    }
    if (15 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quad.venn(area1 = 0, 
                                               area2 = 0, area3 = 0, area4 = 0, n12 = 0, n13 = 0, 
                                               n14 = 0, n23 = 0, n24 = 0, n34 = 0, n123 = 0, 
                                               n124 = 0, n134 = 0, n234 = 0, n1234 = 0, category = category.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               area.vector = area.vector, direct.area = TRUE, 
                                               ...)
    }
    if (31 == length(area.vector)) {
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = 0, 
                                                    area2 = 0, area3 = 0, area4 = 0, area5 = 0, n12 = 0, 
                                                    n13 = 0, n14 = 0, n15 = 0, n23 = 0, n24 = 0, 
                                                    n25 = 0, n34 = 0, n35 = 0, n45 = 0, n123 = 0, 
                                                    n124 = 0, n125 = 0, n134 = 0, n135 = 0, n145 = 0, 
                                                    n234 = 0, n235 = 0, n245 = 0, n345 = 0, n1234 = 0, 
                                                    n1235 = 0, n1245 = 0, n1345 = 0, n2345 = 0, n12345 = 0, 
                                                    category = category.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, area.vector = area.vector, 
                                                    direct.area = TRUE, ...)
    }
  }
  else {
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    if ("none" == na) {
      x <- x
    }
    else if ("stop" == na) {
      for (i in 1:length(x)) {
        if (any(is.na(x[[i]]))) {
          flog.error("NAs in dataset", call. = FALSE, 
                     name = "VennDiagramLogger")
          stop("NAs in dataset", call. = FALSE)
        }
      }
    }
    else if ("remove" == na) {
      for (i in 1:length(x)) {
        x[[i]] <- x[[i]][!is.na(x[[i]])]
      }
    }
    else {
      flog.error("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"", 
                 name = "VennDiagramLogger")
      stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
    }
    if (0 == length(x) | length(x) > 5) {
      flog.error("Incorrect number of elements.", call. = FALSE, 
                 name = "VennDiagramLogger")
      stop("Incorrect number of elements.", call. = FALSE)
    }
    if (1 == length(x)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                 category = list.names, ind = FALSE, ...)
    }
    else if (2 == length(x)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = length(x[[1]]), 
                                                   area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                                         x[[2]])), category = category.names, ind = FALSE, 
                                                   print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (3 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      list.names <- category.names
      nab <- intersect(A, B)
      nbc <- intersect(B, C)
      nac <- intersect(A, C)
      nabc <- intersect(nab, C)
      grob.list <- VennDiagram::draw.triple.venn(area1 = length(A), 
                                                 area2 = length(B), area3 = length(C), n12 = length(nab), 
                                                 n23 = length(nbc), n13 = length(nac), n123 = length(nabc), 
                                                 category = list.names, ind = FALSE, list.order = 1:3, 
                                                 print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else if (4 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n34 <- intersect(C, D)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n134 <- intersect(n13, D)
      n234 <- intersect(n23, D)
      n1234 <- intersect(n123, D)
      grob.list <- VennDiagram::draw.quad.venn(area1 = length(A), 
                                               area2 = length(B), area3 = length(C), area4 = length(D), 
                                               n12 = length(n12), n13 = length(n13), n14 = length(n14), 
                                               n23 = length(n23), n24 = length(n24), n34 = length(n34), 
                                               n123 = length(n123), n124 = length(n124), n134 = length(n134), 
                                               n234 = length(n234), n1234 = length(n1234), category = list.names, 
                                               ind = FALSE, print.mode = print.mode, sigdigs = sigdigs, 
                                               ...)
    }
    else if (5 == length(x)) {
      A <- x[[1]]
      B <- x[[2]]
      C <- x[[3]]
      D <- x[[4]]
      E <- x[[5]]
      list.names <- category.names
      n12 <- intersect(A, B)
      n13 <- intersect(A, C)
      n14 <- intersect(A, D)
      n15 <- intersect(A, E)
      n23 <- intersect(B, C)
      n24 <- intersect(B, D)
      n25 <- intersect(B, E)
      n34 <- intersect(C, D)
      n35 <- intersect(C, E)
      n45 <- intersect(D, E)
      n123 <- intersect(n12, C)
      n124 <- intersect(n12, D)
      n125 <- intersect(n12, E)
      n134 <- intersect(n13, D)
      n135 <- intersect(n13, E)
      n145 <- intersect(n14, E)
      n234 <- intersect(n23, D)
      n235 <- intersect(n23, E)
      n245 <- intersect(n24, E)
      n345 <- intersect(n34, E)
      n1234 <- intersect(n123, D)
      n1235 <- intersect(n123, E)
      n1245 <- intersect(n124, E)
      n1345 <- intersect(n134, E)
      n2345 <- intersect(n234, E)
      n12345 <- intersect(n1234, E)
      grob.list <- VennDiagram::draw.quintuple.venn(area1 = length(A), 
                                                    area2 = length(B), area3 = length(C), area4 = length(D), 
                                                    area5 = length(E), n12 = length(n12), n13 = length(n13), 
                                                    n14 = length(n14), n15 = length(n15), n23 = length(n23), 
                                                    n24 = length(n24), n25 = length(n25), n34 = length(n34), 
                                                    n35 = length(n35), n45 = length(n45), n123 = length(n123), 
                                                    n124 = length(n124), n125 = length(n125), n134 = length(n134), 
                                                    n135 = length(n135), n145 = length(n145), n234 = length(n234), 
                                                    n235 = length(n235), n245 = length(n245), n345 = length(n345), 
                                                    n1234 = length(n1234), n1235 = length(n1235), 
                                                    n1245 = length(n1245), n1345 = length(n1345), 
                                                    n2345 = length(n2345), n12345 = length(n12345), 
                                                    category = list.names, ind = FALSE, print.mode = print.mode, 
                                                    sigdigs = sigdigs, ...)
    }
    else {
      flog.error("Invalid size of input object", name = "VennDiagramLogger")
      stop("Invalid size of input object")
    }
  }
  if (length(x) == 2 & !is.null(total.population) & hyper.test) {
    fisher.mat <- matrix(c(
      total.population-length(union(x[[1]],x[[2]])),
      length(setdiff(x[[1]],x[[2]])),
      length(setdiff(x[[2]],x[[1]])),
      length(intersect(x[[1]],x[[2]]))),nrow=2
    )
    t <- fisher.test(fisher.mat)
    
    if (is.null(sub)) {
      sub = paste0("Hyper.test p = ", signif(t[[1]], digits = 4), "\n OR= ",signif(t[[3]],digits = 4))
    }
    else {
      sub = paste0(sub, ",Hyper.test p = ", signif(t[[1]], digits = 4), "\n OR= ",signif(t[[3]],digits = 4))
    }
  }
  if (!is.null(sub)) {
    grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
                           fontface = sub.fontface, fontfamily = sub.fontfamily, 
                           col = sub.col, cex = sub.cex)
  }
  if (!is.null(main)) {
    grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
                           fontface = main.fontface, fontfamily = main.fontfamily, 
                           col = main.col, cex = main.cex)
  }
  if (!is.null(filename)) {
    current.type <- getOption("bitmapType")
    # if (length(grep("Darwin", Sys.info()["sysname"]))) {
    #   options(bitmapType = "quartz")
    # }
    # else {
    #   options(bitmapType = "cairo")
    # }
    if ("tiff" == imagetype) {
      tiff(filename = filename, height = height, width = width, 
           units = units, res = resolution, compression = compression)
    }
    else if ("pdf" == imagetype) {
      pdf(file =  filename, height = height, width = width)
    }
    else if ("png" == imagetype) {
      png(filename = filename, height = height, width = width, 
          units = units, res = resolution)
    }
    else if ("svg" == imagetype) {
      svg(filename = filename, height = height, width = width)
    }
    else {
      flog.error("You have misspelled your 'imagetype', please try again", 
                 name = "VennDiagramLogger")
      stop("You have misspelled your 'imagetype', please try again")
    }
    grid.draw(grob.list)
    dev.off()
    options(bitmapType = current.type)
    return(1)
  }
  return(grob.list)
}

plot_GO_fortified <-function(formula_res,showCategory = 10){
  gdata <- fortify(formula_res,showCategory = showCategory) %>% filter(pvalue < 0.05) %>% 
    mutate(log.FoldEnrichment =log(parse_ratio(GeneRatio)/ parse_ratio(BgRatio)))
  
  g <- ggplot(data=gdata,aes(x=log.FoldEnrichment,Description))+
    geom_point(aes(size=GeneRatio,color=pvalue))+
    #  geom_text(aes(label=geneID,))+
    #  annotate(geom = "text",label=formula_res_hubgene_LRS_GO$geneID)+
    facet_grid(.~factor(Cluster))+
    scale_color_gradient(low = "#ff0000",high = "#0000ff")+
    labs(x="log(FoldEnrichment)",y="GO Term",size="GeneRatio",color="p.adjust")+
    scale_y_discrete(labels = function(x) stringr::str_wrap(x,width = 40))+
    scale_x_continuous(limits = c(0,max(gdata$log.FoldEnrichment)))+
    theme_bw()+theme(axis.text.x = element_text(angle = 60,hjust = 1),panel.grid.major.x   = element_blank())
  
}



# for hyper tests 
odds.test <- function(data=Count,group="direction",term="locus_group"){
  term.Levels <- levels(factor(data[[term]]))
  group.Levels <- levels(factor(data[[group]]))
  
  Test <-  lapply(group.Levels,function(Group){
    T<- lapply(term.Levels,function(Term){
      S1 = data$n[data[[group]]==Group & data[[term]]==Term] %>% sum()
      S2 = data$n[data[[group]]!=Group & data[[term]]==Term] %>% sum()
      S3 = data$n[data[[group]]==Group & data[[term]]!=Term] %>% sum()
      S4 = data$n[data[[group]]!=Group & data[[term]]!=Term] %>% sum()
      
      t <- fisher.test(x = matrix(data = c(S1,S2,S3,S4),nrow = 2))
      
      res <- data.frame(Group=Group,Term=Term,p=t$p.value,OR=t$estimate,method=t$method,alternative=t$alternative,null.value=t$null.value)
      
    }) %>% bind_rows()
  }) %>% bind_rows()
  n <- nrow(Test)
  
  Test <- Test %>% mutate(bonferroni=p.adjust(p = p,method = "bonferroni",n =n ))
  
}


# annotation 
anno_gene <- function(data=cheRNA_A549,gene_anno=gene_tmap){
  data <- data %>% dplyr::mutate(transcript_id = hash(key = gene_anno$qry_gene_id,value = gene_anno$major_iso_id)[gene_id],
                                 num_exons = hash(key = gene_anno$qry_gene_id,value = gene_anno$num_exons)[gene_id],
                                 transcript_len = hash(key = gene_anno$qry_gene_id,value = gene_anno$len)[gene_id]
  )
}


anno_tx_id <-function(data=cheRNA_A549,gene_anno=gene_tmap){
  data <- data %>% dplyr::mutate(transcript_id = hash(key = gene_anno$ref_id,value = gene_anno$major_iso_id)[gene_id],
                                 num_exons = hash(key = gene_anno$ref_id,value = gene_anno$num_exons)[gene_id],
                                 transcript_len = hash(key = gene_anno$ref_id,value = gene_anno$len)[gene_id]
  )
}