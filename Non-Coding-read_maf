read_maf2 <- function(file) {
  df <- read.delim(file, comment.char="#")
  df <- df[,c(1,9,16)]
  df2 <- df[df$Variant_Classification %in% c("Intron",
                                             "IGR"),]
  df2<-df2 %>% 
    group_by(Hugo_Symbol) %>% 
    filter(n() == 1)
  colnames(df2)[3] <- file
  df2$Variant_Classification <- NULL
  df2[2] <-1 
  return(df2) 
}  

final_df2 <- purrr::map_dfr(files, read_maf2)

final_df2[is.na(final_df2)] <- 0

df3 <- final_df2 %>% group_by(Hugo_Symbol) %>% summarise_all(funs(sum))

df3$sum <- rowSums(df3[,2:68])
strict.df2 <- df3[which(df3$sum >= 2),]
