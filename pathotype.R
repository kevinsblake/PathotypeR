
amrfinder_process <- function(indir, suffix=".tsv"){

  filenames <- list.files(indir, full.names=TRUE)
  filenames2 <- list.files(indir)
  samplelist <- gsub(suffix,"",filenames2) 
  
  mydfs <- lapply(filenames, function(df){
                        tryCatch(read.table(df,
                        sep="\t",
                        comment.char="",
                        check.names = FALSE,
                        quote="",
                        header=TRUE,
                        stringsAsFactors=FALSE),
             
             error=function(e) NULL)
  })
  
  names(mydfs) <- samplelist
  df <- plyr::ldply(mydfs,.id="sample")
  
  return(df)
  
}


pathotypeR <- function(indir, suffix=".tsv", output=c("patho_pred", "patho_prev", "vf_pres", "vf_prev")){
  
  amrfinder_output <- amrfinder_process(indir)
  
  filenames2 <- list.files(indir)
  sample <- gsub(suffix,"",filenames2)
  num_sampl <- length(sample)
  
  main <- data.frame(sample)

  
  # Rename `Gene symbol`s
  amrfinder_fil <- amrfinder_output %>%
    mutate(Gene_name = `Gene symbol`) %>%
    mutate(Gene_name = recode(Gene_name, 
                                     "stxA1a" = "stx1", "stxA1c" = "stx1", "stxB1a" = "stx1", "stxB1c" = "stx1",
                                     "stxA2b" = "stx2", "stxA2c" = "stx2", "stxA2d" = "stx2", "stxA2f" = "stx2", "stxB2a" = "stx2", "stxB2b" = "stx2", "stxB2c" = "stx2", "stxB2f" = "stx2", 
                                     "ipaH1" = "ipaH", "ipaH2" = "ipaH", "ipaH3" = "ipaH", "ipaH4" = "ipaH", "ipaH5" = "ipaH"))
  
  # Filter just VFs and remove multiple occurence of a VF per sample
  output_fil <- amrfinder_fil %>% 
    filter(`Element type` == "VIRULENCE") %>%
    group_by(sample, Gene_name) %>%
    filter(row_number() == 1) %>% 
    ungroup() %>% 
    mutate(Gene_name = stringr::str_replace_all(Gene_name, "[/-]", "_"))
  
  # Calculate VF prevalence        
  vf_prev <- output_fil %>% 
    group_by(Gene_name) %>% 
    summarise(Cnt = n(), 
                     prev = round(Cnt/num_sampl, 4)) %>%
    ungroup()
  
  # VF presence/absence (1 = present, 0 = absent)
  sample_vfs <- output_fil %>%
    filter(Gene_name %in% vf_prev$Gene_name) %>% 
    inner_join(main, by="sample") %>%
    group_by(sample, Gene_name) %>%
    summarise(Cnt = n()) %>% 
    reshape2::dcast(data = ., sample ~ Gene_name, value.var = "Cnt", fill=0, fun.aggregate = sum) %>% 
    mutate_if(is.numeric, ~1 * (. > 0)) %>%
    add_cols(c("stx1", "stx2", "eae", "bfpA", "aggR", "aaiC", "aatA", "ltcA", "sta1", "afaC", "afaE", "ipaH"))
  
    vfcols <- colnames(sample_vfs[-1])
  
    # Count total number VFs per sample
    main <-  main %>% 
      inner_join(sample_vfs, by ="sample") %>%  
      rowwise %>% mutate(total_vfs = sum(c_across(where(is.numeric)))) %>%
      select(sample, total_vfs, everything()) 
    
    # Predict pathotype based on VFs present
    vf.pred <- main %>% select(sample, any_of(vfcols)) %>% 
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
                         sta1 == 1 ~ 1,
                         TRUE ~ 0),
        DAEC = case_when(afaC == 1 ~ 1,
                         afaE == 1 ~ 1,
                         TRUE ~ 0),
        EIEC = case_when(ipaH == 1 ~ 1,
                         TRUE ~ 0),
        none = case_when(stx1 == 0 & stx2 == 0 & eae == 0 & bfpA == 0 & aggR == 0 & aaiC == 0 & aatA == 0 &
                         ltcA == 0 & sta1 == 0 & afaC == 0 & afaE == 0 & ipaH == 0 ~ 1, 
                         TRUE ~ 0)) %>% 
      select(sample, 
                    EAEC, ETEC, DAEC, EPEC, STEC, EHEC, EIEC, none)
    
    pathcols <- c("EAEC", "EPEC", "STEC", "ETEC", "EHEC", "EIEC", "DAEC", "none")
    
    tmp.df <- vf.pred %>% 
      tidyr::gather(key="Pathotype", value="flag", any_of(pathcols)) %>%
      filter(flag == 1) %>% 
      group_by(sample) %>% 
      summarise(pred = paste(Pathotype, collapse = "-"))
    
    # Predict pathotype
    patho_pred <- main %>% 
      inner_join(tmp.df, by ="sample") %>%  
      select(sample, total_vfs, pred)
    
    # Calculate pathotype prevalence        
    patho_prev <- tmp.df %>% 
      group_by(pred) %>% 
      summarise(Cnt = n(), 
                       prev = round(Cnt/num_sampl, 4)) %>%
      ungroup()

    if (missing(output)){
      output <- "patho_pred"
    }
    
    if (output == "patho_pred") {
      patho_output <- patho_pred
    }
    
    if (output == "vf_pres") {
      patho_output <- sample_vfs
    }
    
    if (output == "vf_prev") {
      patho_output <- vf_prev
    }
    
    if (output == "patho_prev") {
      patho_output <- patho_prev
    }

    return(patho_output)
    
}

add_cols <- function(df, cols) {
  add <- cols[!cols %in% names(df)]
  if(length(add) != 0) df[add] <- 0
  return(df)
}

