read_maf <- function(file) {
  df <- read.delim(file, comment.char="#")
    df <- df[,c(1,9,16)]
      df2 <- df[df$Variant_Classification %in% c("Missense_Mutation",
                                               "3'UTR", 
                                               "Nonsense_Mutation",
                                               "5'UTR", 
                                               "Frame_Shift_Del",
                                               "Frame_Shift_Ins"),]
    df2<-df2 %>% 
      group_by(Hugo_Symbol) %>% 
      filter(n() == 1)
   colnames(df2)[3] <- file
    df2$Variant_Classification <- NULL
    df2[2] <-1 
    return(df2) 
    }  
  

final_df <- purrr::map_dfr(files, read_maf)
final_df[is.na(final_df)] <- 0

df2 <- final_df %>% group_by(Hugo_Symbol) %>% summarise_all(funs(sum))

df2$sum <- rowSums(df2[,2:68])
strict.df <- df2[which(df2$sum >= 2),]
