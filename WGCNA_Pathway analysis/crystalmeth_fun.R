# Function for applying missMethyl GO and KEGG analysis
# Creates output tables, adapted by Dr. Roy Lardenoije and Dr. Ehsan Pishva

crystalmeth <- function(sig.cpg, 
                        all.cpg, 
                        collection = "GO", 
                        array.type = c("450K","EPIC")){
  
  ## Check and load required packages
  pkgs <- c("missMethyl", "goseq", "annotate", "org.Hs.eg.db", 
            "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "clusterProfiler",
            "KEGGREST", "reshape2","IlluminaHumanMethylation450kanno.ilmn12.hg19")
  for(i in 1:length(pkgs)){
    suppressPackageStartupMessages(OK <- requireNamespace(pkgs[i], quietly = TRUE))
    if(!OK) stop(pkgs[i]," package required but not not installed (or can't be loaded)")
  }
  
  library(missMethyl, quietly = TRUE)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, quietly = TRUE)
  library(clusterProfiler, quietly = TRUE)
  library(reshape2, quietly = TRUE)
  library(KEGGREST, quietly = TRUE)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  if(collection == "GO"){
    ## Run gometh
    results <- missMethyl::gometh(sig.cpg = sig.cpg, 
                                  all.cpg = all.cpg, 
                                  collection = "GO", 
                                  array.type = array.type,
                                  plot.bias = FALSE,
                                  prior.prob = TRUE)
    
    
    ## Select terms with 10 - 2000 genes and adjust FDR
    results <- results[results$N >= 10 & results$N <= 2000, ]
    results$FDR <- p.adjust(results$P.DE, method = "fdr")
    
    
    ## Add gene names in term and test list
    # Find gene names mapped to CpG sites
    flat.u <-missMethyl:::.getFlatAnnotation(array.type)
    out <- getMappedEntrezIDs(sig.cpg = sig.cpg, 
                              all.cpg = all.cpg,
                              array.type = array.type)
    de <- out$sig.eg
    universe <- out$universe
    
    # Construct dataframe with gene IDs mapped to GO terms
    egGO2ALLEGS <- getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
    m <- match(AnnotationDbi::Lkeys(egGO2ALLEGS), universe, 0L)
    universe <- universe[m]
    AnnotationDbi::Lkeys(egGO2ALLEGS) <- universe
    EG.GO <- AnnotationDbi::toTable(egGO2ALLEGS)
    d <- duplicated(EG.GO[ , c("gene_id", "go_id", "Ontology")])
    EG.GO <- EG.GO[!d, ]
    
    for(i in 1:length(de)){
      names(de)[i] <- flat.u[flat.u$entrezid == de[i], "symbol"][1]
    }
    
    # Add column with genes in term and test list
    results$Genes <- NA
    for(i in 1:nrow(results)){
      results$Genes[i] <- paste(intersect(unique(names(de)),
                                          unique(flat.u[flat.u$entrezid %in% EG.GO[EG.GO$go_id == rownames(results)[i], "gene_id"], "symbol"])), collapse = "; ")
    }
  }
  
  
  if(collection == "KEGG"){
    ## Load KEGG data
    hsaKegg <- download_KEGG('hsa')
    kegg <- hsaKegg$KEGGPATHID2EXTID
    k <- dcast(kegg, kegg$to ~ kegg$from)
    k <- data.frame(k[2:ncol(k)])
    
    ## Run gsameth
    results <- data.frame(missMethyl::gsameth(sig.cpg = sig.cpg, 
                                              all.cpg = all.cpg, 
                                              collection = k, 
                                              prior.prob = TRUE, 
                                              plot.bias = FALSE))
    
    ## Add pathway description
    desc <- ""
    for(i in 1:nrow(results)){
      hsa <- rownames(results)[i] 
      xxx <- keggGet(hsa)
      desc[i] <- xxx[[1]]$NAME
    }
    results <- data.frame(desc, results)
    
    ## Clean up description column
    results$desc <- gsub(" \\- Homo sapiens \\(human\\)", "", results$desc)
    results$desc <- gsub(" \\- Homo sapiens \\(human\\)", "", results$desc)
    
    ## Select pathways with > 10 genes and adjust FDR
    results <- results[results$N >= 10, ]
    results$FDR <- p.adjust(results$P.DE, method = "fdr")
    
    ## Add gene names in term and test list
    # Find gene names mapped to CpG sites
    flat.u <-missMethyl:::.getFlatAnnotation(array.type)
    out <- getMappedEntrezIDs(sig.cpg = sig.cpg, 
                              all.cpg = all.cpg,
                              array.type = array.type)
    de <- out$sig.eg
    
    # Construct dataframe with gene IDs mapped to KEGG pathways
    d <- duplicated(kegg[ , c("from", "to")])
    kegg <- kegg[!d, ]
    
    
    for(i in 1:length(de)){
      names(de)[i] <- flat.u[flat.u$entrezid == de[i], "symbol"][1]
    }
    
    # Add column with genes in term and test list
    results$Genes <- NA
    for(i in 1:nrow(results)){
      results$Genes[i] <- paste(intersect(unique(names(de)),
                                          unique(flat.u[flat.u$entrezid %in% kegg[kegg$from == rownames(results)[i], "to"], "symbol"])), collapse = "; ")
    }
  }
  
  
  ## Order by p-value
  results <- results[order(results$P.DE), ]
  results
}
