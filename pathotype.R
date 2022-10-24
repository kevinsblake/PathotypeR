
amrfinder_process <- function(indir){

  filenames <- list.files(indir, full.names=TRUE)
  filenames2 <- list.files(indir)
  samplelist <- gsub("_out.tsv","",filenames2) 
  
  mydfs <- lapply(filenames, function(df){
                        tryCatch(read.table(df,
                        sep="\t",
                        comment.char="",
                        check.names = FALSE,
                        quote="",
                        header=TRUE,
                        stringsAsFactors=FALSE),
             
             # ignore if there is an error like an empty file
             error=function(e) NULL)
  })
  
  names(mydfs) <- samplelist
  df <- plyr::ldply(mydfs,.id="sample")
  
  return(df)
  
}

patho_pred <- function(indir, all=FALSE){
  
  amrfinder_output <- amrfinder_process(indir)
  
  filenames2 <- list.files(indir)
  sample <- gsub("_out.tsv","",filenames2)
  num_sampl <- length(sample)
  
  main <- data.frame(sample)
  
  # Rename `Gene symbol`s
  amrfinder_fil <- amrfinder_output %>%
    dplyr::mutate(Gene_name = `Gene symbol`) %>%
    dplyr::mutate(Gene_name = recode(Gene_name, 
                                     "stxA1a" = "stx1", "stxA1c" = "stx1", "stxB1a" = "stx1", "stxB1c" = "stx1",
                                     "stxA2b" = "stx2", "stxA2c" = "stx2", "stxA2d" = "stx2", "stxA2f" = "stx2", "stxB2a" = "stx2", "stxB2b" = "stx2", "stxB2c" = "stx2", "stxB2f" = "stx2", 
                                     "ipaH1" = "ipaH", "ipaH2" = "ipaH", "ipaH3" = "ipaH", "ipaH4" = "ipaH", "ipaH5" = "ipaH"))
  
  # Filter just VFs and remove multiple occurence of a VF per sample
  output_fil <- amrfinder_fil %>% 
    dplyr::filter(`Element type` == "VIRULENCE") %>%
    dplyr::group_by(sample, Gene_name) %>%
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(Gene_name = stringr::str_replace_all(Gene_name, "[/-]", "_")) #necessary?
  
  # Calculate VF prevalence        
  vf_prev <- output_fil %>% 
    dplyr::group_by(Gene_name) %>% 
    dplyr::summarise(Cnt = n(), 
                     prev = round(Cnt/num_sampl, 4)) %>%
    dplyr::ungroup()
  
  # VF presence/absence (1 = present, 0 = absent)
  sample_vfs <- output_fil %>%
    filter(Gene_name %in% vf_prev$Gene_name) %>% 
    inner_join(main, by="sample") %>%
    group_by(sample, Gene_name) %>%
    summarise(Cnt = n()) %>% 
    reshape2::dcast(data = ., sample ~ Gene_name, value.var = "Cnt", fill=0, fun.aggregate = sum) %>% 
    mutate_if(is.numeric, ~1 * (. > 0))
  
    vfcols <- colnames(sample_vfs[-1])
  
    # Count total number VFs per sample
    main <-  main %>% 
      inner_join(sample_vfs, by ="sample") %>%  
      mutate(total_vfs = rowSums(select(., vfcols), na.rm = TRUE)) %>%
      dplyr::select(sample, total_vfs, everything())
    
    # Predict pathotype based on VFs present
    vf.pred <- main %>% dplyr::select(sample, vfcols) %>% 
      mutate(
        STEC = case_when(stx1 == 1 & eae == 0 ~ 1,
                         stx2 == 1 & eae == 0 ~ 1,
                         TRUE ~ 0),
        EPEC = case_when(eae == 1 & stx1 == 0 & stx2 == 0 ~ 1,
                         bfpA == 1 & stx1 == 0 & stx2 == 0 ~ 1,
                         TRUE ~ 0),
        EHEC = case_when(stx1 == 1 & eae == 1 ~ 1,
                         stx2 == 1 & eae == 1 ~ 1,
                         TRUE ~ 0),
        EAEC = case_when(aggR == 1 ~ 1,
                         aaiC == 1 ~ 1,
                         aatA == 1 ~ 1,
                         TRUE ~ 0),
        ETEC = case_when(ltcA == 1 ~ 1,
                         sta1 == 1 ~ 1, #astA
                         TRUE ~ 0),
        DAEC = case_when(afaC == 1 ~ 1,
                         TRUE ~ 0),
        EIEC = case_when(ipaH == 1 ~ 1, #ipaD
                         TRUE ~ 0),
        none = case_when(stx1 == 0 & stx2 == 0 & eae == 0 & bfpA == 0 & aggR == 0 & aaiC == 0 & aatA == 0 &
                         ltcA == 0 & sta1 == 0 & afaC == 0 & ipaH == 0 ~ 1, 
                         TRUE ~ 0)) %>% 
      dplyr::select(sample, 
                    EAEC, ETEC, DAEC, EPEC, STEC, EHEC, EIEC, none)
    
    pathcols <- c("EAEC", "EPEC", "STEC", "ETEC", "EHEC", "EIEC", "DAEC", "none")
    
    tmp.df <- vf.pred %>% 
      dplyr::select(sample, all_of(pathcols)) %>%
      tidyr::gather(key="Pathotype", value="flag", pathcols) %>%
      filter(flag == 1) %>% 
      group_by(sample) %>% 
      summarise(pred = paste(Pathotype, collapse = "-"))
    
    patho_output <- main %>% 
      inner_join(tmp.df, by ="sample") %>%  
      dplyr::select(sample, pred)
    
    if (all) {
      patho_output <- main %>% 
        inner_join(tmp.df, by ="sample") %>%  
        dplyr::select(sample, pred, total_vfs, everything())
    }

    return(patho_output)
    
}

