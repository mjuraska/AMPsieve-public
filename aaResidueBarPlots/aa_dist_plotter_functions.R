library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(seqinr)

#' Function that builds the foundation panels for the final output figure
#'
#' @param sequenceData A dataframe with a column for each observed sequence and a column for the associated pubids
#' @param sites A vector of numerical values indicating the sites that should be analyzed
#' @param map A dataframe with a column for the hxb2 position and another column with the associated alignmnet position
#' @param seqMaster A dataframe with a column for each unique pubid in \code{sequenceData}, and a column with the associated treatment group
#' 
#' @return Returns a list of ggplot objects, each with a collection of vertical bars associated with each given position. 
#' The height of each color is dependant on the proportion of that amino acid in the total sequences for the associated pubid.
#' 
dist.plot.builder <- function(sequenceData, sites, map, seqMaster){
  
  # vector with all possible colors
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # colors <- c('white', col_vector[1:18])
  # names(colors) <- c('-', 'D', 'S', 'N', 'T', 'E', 'G', 'Q', 'P', 'I', 'A', 'V', 'C', 'M', 'K', 'L', 'F', 'R', 'H')
  
  # order seqMaster and the sequenceData so that they line up with one another
  seqMaster <- seqMaster[order(seqMaster$pub_id), ]
  sequenceData <- sequenceData[order(sequenceData$pubid), ]
  
  plots <- lapply(sites, function(hxb2){
    
    # use hxb2 mapping to find the correct position
    pos <-  subset(map, hxb2Pos==hxb2)$posNum
    
    # find all residues at the current position
    # siteAA <- tapply(sequenceData$sequence, sequenceData$pubid, function(x){substr(x, pos, pos)})
    siteAA <- sequenceData %>%
      group_by(pubid) %>%
      summarise(res=substr(sequence, pos, pos)) %>%
      ungroup()
    
    # order unique AAs by their prevalence
    # AAs <- names(sort(table(unlist(siteAA))))
    AAs <- names(sort(table(siteAA$res)))
    
    # create 'df' for the current position
    df <- data.frame(pubid = rep(unique(sequenceData$pubid), length(AAs)),
                     residue = rep(rev(AAs), each=length(unique(sequenceData$pubid))),
                     tx = rep(seqMaster$tx, length(AAs)))
    
    # populate the 'prop' column with a proportion of each acid for each subject
    # df$prop <- sapply(1:length(df$pubid), function(i){
    #   mean(unlist(siteAA[df$pubid[i]])==df$residue[i])
    # })
    df$prop <- sapply(1:length(df$pubid), function(i){
      siteAA.pubid <- filter(siteAA, pubid==df$pubid[i])
      mean(siteAA.pubid$res==df$residue[i])
    })
    
    # spread 'df' out so each AA has a column
    # 'spread' is equivalent to pivot_wider(names_from=residue, values_from=prop)
    df <- df %>% spread(residue, prop)
    
    # round each proportion for a lower resolution sort
    df1 <- df %>% mutate_if(is.numeric, round, digits=1)
    
    # order 'df1' by the most prevalent AA, suborder by next most prevalent, and so on
    # always put the '-' in front
    if('-' %in% AAs) { AAs <- c(AAs[!grepl("-", AAs)], '-') }
    df1 <- df1[do.call('order', df1[rev(AAs)]), ]
    
    # relevel 'pubid' by the order of 'pubid' in 'df1'
    df$pubid <- factor(df$pubid, levels=rev(unique(df1$pubid)))
    
    # pivot 'df' to its original structure for plotting
    # 'gather' equivalent to pivot_longer
    df <- df %>% gather(residue, prop, -pubid, -tx)
    
    # relevel the amino acids by their prevalence, so they're stacked in this order
    if('-' %in% AAs) { AAs <- c(AAs[!grepl("-", AAs)], '-') }
    df$residue <- factor(df$residue, levels=AAs)
    
    # color each residue so that gap is first and white, N is tomato3, and P is skyblue3
    colors <- rep(NA, length(AAs))
    names(colors) <- c(AAs[AAs != '-'], AAs[AAs == '-'])
    brewerColors <- col_vector[1 : length(AAs[!(AAs %in% c("-", "N", "P"))])]
    colors[!(AAs %in% c('-', 'N', 'P'))] <- rev(brewerColors)
    colors[AAs == '-'] <- 'white'
    colors[AAs == 'N'] <- 'tomato3'
    colors[AAs == 'P'] <- 'skyblue3'
    
    # create the AA residue bar plot with the color palette chosen earlier
    p <- ggplot(df, aes(x=pubid, y=prop, fill=forcats::fct_rev(residue))) +
      geom_col(width=0.7) +
      ylab(hxb2) +
      scale_fill_manual(values=colors, guide=guide_legend(title=NULL, reverse=TRUE, ncol=3)) +
      facet_grid(. ~ tx, scales="free", space="fixed") +
      theme_void() +
      theme(axis.title.y=element_text(angle=90), strip.text=element_blank(), legend.key.size=unit(0.25, 'cm'), legend.text = element_text(size=8.5),
            plot.margin=unit(c(0,0.75,0,0.1),'cm'), plot.background = element_rect(fill="grey90", color='grey90'),
            legend.justification = 'top')
    
    return(p)
  }) # end position
  
  return(plots)
}

#' A function which gathers sequences spread across several files, such as the sequence data associated with the AA alignment for HVTN 704 and 703
#'
#' @param path The path to the .fasta files containing sequence data
#' 
#' @return A data frame, with columns indicating sampleid, pubid, and genetic sequence. One row is given for each observed sequence, not each unique sequence.
#' 
gatherSequences <- function(path){
  
  files <- list.files(path, pattern="uncollapsed")
  
  if(any(!grepl('.fasta', files))) {
    print("There seem to be non-fasta files in the given directory. They will be ignored.")
    files <- files[grepl('.fasta', files)]
  }
  
  # Populate the data frame with sequence data
  df <- data.frame()
  
  for(i in 1:length(files)){
    
    tempData <- read.fasta(paste0(path, "/", files[i]))
    tempDf <- data.frame(id=names(tempData), sequence.name=names(tempData), sequence=sapply(tempData, function(x){paste(x, collapse="")}))
    rownames(tempDf) <- NULL
    df <- rbind(df, tempDf)
    
  }
  
  # delete rows with HXB2
  df <- subset(df, id!="HXB2_translation")
  
  # replace last two characters "_2" with "-2" so that the number of separators "_" remains constant
  # for each value of df$id
  tmp <- substring(df$id, nchar(df$id) - 1)
  substring(df$id, nchar(df$id) - 1) <- ifelse(tmp=="_2", "-2", tmp)

  # Separate id into individual columns
  df <- separate(data=df, col=id, into=c("study", "pubid", "visit", "region", "seqID"), sep="_")
  # Unite several columns to form the sampleid according to the Mullins Lab
  df <- unite(data=df, col=sampleid, c("study", "pubid", "visit", "region"), sep="_")
  df$sampleid <- substr(df$sampleid, 2, nchar(df$sampleid))
  
  df$pubid <- substr(df$sampleid, 1, 8)
  df$sequence <- toupper(df$sequence)
  
  # Remove incomplete entries (HXB2 translations, for example)
  return(na.omit(df))
}

# a version of 'gatherSequences' for a single fasta file, from which a single variant (mf, ms, ls)
# needs to be filtered
gatherSequences2 <- function(path, filename, variant){
  if (!grepl(".fasta", filename)){ stop("The supplied 'filename' must have the file extension '.fasta'.") }
  
  # populate the data frame with sequence data
  df <- read.fasta(file.path(path, filename))
  df <- data.frame(id=names(df), sequence.name=names(df), sequence=sapply(df, function(x){ paste(x, collapse="") }))
  rownames(df) <- NULL
  
  # delete rows with HXB2 and keep only a single variant
  df <- filter(df, !grepl("HXB2", id), grepl(variant, id))
  
  # replace last two characters "_2" with "-2" so that the number of separators "_" remains constant
  # for each value of df$id
  # tmp <- substring(df$id, nchar(df$id) - 1)
  # substring(df$id, nchar(df$id) - 1) <- ifelse(tmp=="_2", "-2", tmp)
  
  # Separate id into individual columns
  df <- separate(data=df, col=id, into=c("study", "pubid", "visit", "region", "seqID", "seqType"), sep="_")
  # Unite several columns to form the sampleid according to the Mullins Lab
  df <- unite(data=df, col=sampleid, c("study", "pubid", "visit", "region"), sep="_")
  df$sampleid <- substr(df$sampleid, 2, nchar(df$sampleid))
  
  df$pubid <- substr(df$sampleid, 1, 8)
  df$sequence <- toupper(df$sequence)
  
  # Remove incomplete entries (HXB2 translations, for example)
  return(na.omit(df))
}