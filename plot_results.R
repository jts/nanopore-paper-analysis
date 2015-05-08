# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, cols) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)

    numPlots = length(plots)

    # Make the panel
    plotCols = cols                          # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols

    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
            curCol = (i-1) %% plotCols + 1
            print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
}

plot_count_correlation <- function(filename) {
    
    require(ggplot2)
    require(reshape2)

    data <- read.table(filename, header=T)

    # rename columns to refer to reference and assembly
    names(data)[names(data) == "NC_000913.fna.gi.556503834.ref.NC_000913.3..fwd"] <- "reference"
    names(data)[names(data) == "circular_draft.fasta.1.rc"] <- "draft"
    names(data)[names(data) == "circular_polished.fasta.1.rc"] <- "polished"
    print(head(data))
    
    # reorder by reference count descending
    data$kmer = reorder(data$kmer, -data$reference)
    
    # Extract cases where the reference and draft assembly count differ by more than a factor
    F = 1.5
    w = 0.5 # width of bars

    data$description = ifelse(data$reference / data$draft > F, "outlier", "normal") 
    data$kmer_label = ifelse(data$description == "outlier", data$kmer, "") 
 
    print(cor(data$reference, data$draft))
    print(cor(data$reference, data$polished))
   
    #
    # Draft assembly plots
    #
    p1 <- ggplot(data, aes(reference, draft)) + 
          geom_point(aes(color=description)) +
          #geom_text(aes(label=kmer_label)) +
          scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
          xlab("5-mer count in reference") + 
          ylab("5-mer count in draft assembly") +
          theme_bw(base_size=13)
    
    # Over represented in reference
    outliers = subset(data, description == "outlier")

    md <- melt(outliers[,c('kmer', 'draft', 'reference')], id.vars=1)
    p2 <- ggplot(md, aes(kmer, value)) + geom_bar(aes(fill=variable), width = w, position = position_dodge(width = w), stat="identity") + ylab("count") + theme_bw(base_size=13)
   
    p3 <- ggplot(data, aes(reference, polished)) + 
          geom_point(aes(color=description)) +
          #geom_text(aes(label=kmer_label)) +
          scale_colour_manual(values=c("black", "red"), guide=FALSE) + 
          xlab("5-mer count in reference") + 
          ylab("5-mer count in polished assembly") + theme_bw(base_size=13)

    mp <- melt(outliers[,c('kmer', 'polished', 'reference')], id.vars=1)
    p4 <- ggplot(mp, aes(kmer, value)) + geom_bar(aes(fill=variable), width = w, position = position_dodge(width = w), stat="identity") + ylab("count") + theme_bw(base_size=13)

    pdf("figure_kmer_counts.pdf", 16, 16)
    multiplot(p1, p2, p3, p4, cols=2, rows=2)
    dev.off() 
    png("figure_kmer_counts.png", 1024, 1024)
    multiplot(p1, p2, p3, p4, cols=2, rows=2)
    dev.off()  
}

args = commandArgs(TRUE)
plot_count_correlation(args[1])

