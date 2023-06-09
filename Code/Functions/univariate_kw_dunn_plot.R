univariate_kw_dunn_plot <- function(y
                                    , Letter_height = 5400
                                    , y_lab = NULL
                                    , p_type = c("adjusted", "unadjusted")){
  
  # Load required packages
  if(!require(multcompView)){library(multcompView)}
  if(!require(data.table)){library(data.table)}
  if(!require(ggplot2)){library(ggplot2)}
  
  
  #print("y and Treatment must be part of a data table called dt")
  #y <- get(y)
  y_var <- dt[, ..y]
  dt[, y_col := y_var]
  #return(dt)
  # univariate_kw_dunn_plot("Species_richness")

  # Kruskal wallis test
  kw <- kruskal.test(y_col ~ Treatment
                        , data = dt)
  #   return(kw)
  # }
  
  # Dunn's test
  dunn <- dunnTest(y_col ~ Treatment
                        , data = dt
                        , method = "bh")

  # Dunn comparison letters. Might need checking
  if(p_type == "adjusted"){a <- dunn$res$P.adj} else{a <- dunn$res$P.unadj}
  
  names(a) <- dunn$res$Comparison
  names(a) <- gsub("-Year", ";Year", names(a))
  a  <-  na.omit(a)
  a <- multcompLetters(a)
  names(a$Letters) <- gsub(";Year", "-Year", names(a$Letters))
  names(a$Letters) <- gsub("\\ ", "", names(a$Letters))
  
  if(p_type == "unadjusted"){
    dunn_letters <- sapply(unique(names(a$Letters)), FUN = function(x){
      wrk <- a$Letters[names(a$Letters) == x]
      out <- unique(wrk[which(nchar(wrk) == min(nchar(wrk)))])
      return(out)
    })
  } else{
    dunn_letters <- unique(a$Letters)
  }
  
  
  
  DunnLetters <- data.table(Letters = dunn_letters
                            , Treatment = unique(names(a$Letters))
                            , Y = rbindlist(tapply(dt$y_col
                                                                , dt$Treatment
                                                                , mean_ci))$ymax
                            )
  #  return(list(kw, dunn,  DunnLetters))
  # }
  # univariate_kw_dunn_plot("Species_richness")
  
  p <- ggplot(dt, aes(x = Treatment, y = y_col
                      , colour = Treatment)) + 
    geom_point(colour = 'grey50', position = position_jitter(0.15)) +
    stat_summary(fun.data = mean_ci, size = 1.1) +
    theme_classic() + 
    scale_colour_manual(breaks = levels(dt$Treatment)
                        , values = BlackoutPalette) +
    scale_x_discrete(breaks=levels(dt$Treatment)) +
    guides(colour = 'none') +
    ylab(y_lab) +
    geom_text(data = DunnLetters, aes(label = Letters
                                      , x = Treatment
                                      , y = Letter_height)
              , colour = 'black') +
    theme(axis.text = element_text(size = rel(1)))
  
  ef <- dt[, .(Median = as.double(median(y_col)), IQR = as.double(IQR(y_col)))
           , by = Treatment]
  
    return(list(Kruskal_wallis = kw
                , Dunn = dunn
                , Dunn_letters = DunnLetters
                , Plot = p
                , Effect = ef))
}

# Test this
# univariate_kw_dunn_plot("Species_richness", y_lab = "Species richness")
  
