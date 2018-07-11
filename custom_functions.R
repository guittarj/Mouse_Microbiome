# Custom functions for microbiome analyses. Created Feb 28, 2017.

bac_extract <- function(x, trait) {
  
  # for each species query result, filter blanks and pick best value
  for(i in 1:length(tax_bacdat)) {
    names(x[[i]]) <- tax_bacdat[[i]]
    x[[i]] <- x[[i]][sapply(x[[i]], length) > 0]
    x[[i]] <- sapply(x[[i]], pick_best)
  }
  
  # unlist, clean, return df
  x <- unlist(x)
  x[x == 'NULL'] <- NA
  x_names <- strsplit(names(x), '.', fixed = TRUE)
  x <- data.frame(
    Genus = sapply(x_names, function(y) y[[1]]),
    Species = sapply(x_names, function(y) y[[2]]),
    trait = trait,
    val = x
  )
  x <- x[!is.na(x$val), ]
  x$val <- as.character(x$val)
  return(x)
}

bac_search <- function (x) {
  # needs list of genera to search for
  
  #just use the first page of results
  page <- '1'
  taxon <- 'https://bacdive.dsmz.de/api/bacdive/taxon/'
  
  for (i in 1:length(x)) {
    
    url_genus <- URLencode(paste0(taxon, x[i], '/?page=', page, '&format=json'))
    response <- getURL(url_genus, userpwd="guittarj@msu.edu:twkYkcJbxCQzEvwkUi", httpauth = 1L)
    jsondata <- fromJSON(response)
    urls <- unlist(jsondata$results)
    
    if (length(urls) > 0) {
      genus <- x[i] 
      bacdat[[genus]] <- list()
      for (j in 1:length(urls)) {
        urlx <- URLencode(paste0(urls[j],'?&format=json'))
        response <- getURL(urlx, userpwd="guittarj@msu.edu:twkYkcJbxCQzEvwkUi", httpauth = 1L)
        jsondata <- fromJSON(response)
        bacdat[[genus]][[j]] <- jsondata
      }
    }
    print(paste("Downloaded", length(urls), "entries for", x[i]))
  }
}


bac_search_ID <- function (x, bacdat = list()) {
  
  #just use the first page of results
  page <- '1'
  bacURL <- 'https://bacdive.dsmz.de/api/bacdive/'
  
  for (i in x) {
    
    url_bacID <- URLencode(paste0(bacURL, 'bacdive_id/', i, '/?page=', page, '&format=json'))
    response <- getURL(url_bacID, userpwd="guittarj@msu.edu:twkYkcJbxCQzEvwkUi", httpauth = 1L)
    bacdat[[i]] <- fromJSON(response)
    
    prog <- match(i, x)
    if(match(i, x) %% 500 == 0) {
      
      #Print status
      print(paste("Downloaded", match(i, x), "of", length(x), "entries"))
      
      #save progress
      saveRDS(bacdat, file = paste0('bacdat', match(i, x), '.RDS'))
    }
  
  }
  
  return(bacdat)

}


loadpax <- function(pkg){
  # (1) checks package installation, (2) installs them if not, then (3) loads them
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}


grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  # A function that lets you put multiple plots together and have a single shared legend (from the first plot)
  # from hadley on the internet...
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#suppress warnings as numeric
swan <- function (x) suppressWarnings(as.numeric(x))

# print list of loaded functions
print(data.frame(Custom_Functions = 
  c('bac_search: search bacdat',
  	'bac_extract: extract data from bacdat object',
  	'loadpax: install+load multiple packages',
    'grid_arrange_shared_legend: Multiple plots, one legend',
    'multiplot: multiple plots',
    'swan: a function to suppress warnings, convert to numeric')))
