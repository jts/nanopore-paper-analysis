require(ggplot2)
means = c(60.3, 40.6, 52.2, 54.1, 61.5, 72.7, 49.4, 45.2)
duration = c(521, 112, 356, 51, 291, 15, 141, 301)
sd = c(0.7, 1.0, 2.0, 1.2, 1.0, 1.5, 1.5, 1.2)

data = data.frame(c( rnorm(duration[1], means[1], sd[1]), 
                     rnorm(duration[2], means[2], sd[2]), 
                     rnorm(duration[3], means[3], sd[3]), 
                     rnorm(duration[4], means[4], sd[4]), 
                     rnorm(duration[5], means[5], sd[5]),
                     rnorm(duration[6], means[6], sd[6]),
                     rnorm(duration[7], means[7], sd[7]),
                     rnorm(duration[8], means[8], sd[8])))

colnames(data) <- "samples"
data$index <- 1:nrow(data)

p <- ggplot(data, aes(index / 1000, samples)) + 
    geom_point() + 
    ylim(0, 100) + 
    ylab("Current (pA)") + 
    xlab("time (s)")

# draw lines
last_end = 0
for (i in seq(1, length(duration))) {
    start = last_end
    end = last_end + duration[i] / 1000
    m = means[i]

    p <- p + geom_segment( x = start, xend = end, y = m, yend = m, col="red" )
    last_end <- end

    s <- sprintf("%d & %.1f & %.1f & %.3f", i, m, sd[i], duration[i] / 1000)
    print(s)
}
 #   p + geom_segment( aes( x = duration[1] / 1000, xend = duration[2] / 1000, y = means[2], yend = means[2], col="red" ) )
ggsave("figures/simulation.pdf")
