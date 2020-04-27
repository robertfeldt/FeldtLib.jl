using RCall

function ggplot2_distribution(values::AbstractVector, filename;
    title = "Distribution",
    xaxistitle = "Cases",
    types = [:histogram, :density],
    yaxistitle = "Count",
    xlog10 = false)
    values = convert(Vector{Float64}, values)
    R"library(ggplot2)"
    R"library(ggthemes)"
    @rput values
    R"df <- data.frame(dummytoensurerightsize = values)"
    @rput xaxistitle
    R"df[, xaxistitle] <- values"
    @rput yaxistitle
    @rput title
    extras = xlog10 ? "+ scale_x_log10()" : ""
    R"p <- ggplot(df, aes_string(x = xaxistitle)) + theme_tufte() + xlab(xaxistitle) + ggtitle(title) $extras"
    basefilename, ext = splitext(filename)
    if in(:histogram, types)
        R"phist <- p + geom_histogram(binwidth=.5, color='black', fill='blue', alpha=0.5) + ylab(yaxistitle)"
        fname = basefilename * "_histogram" * ext
        @rput fname
        R"ggsave(fname, phist)"
    end
    if in(:density, types)
        R"pdens <- p + geom_density(fill='blue', alpha=0.5)"
        fname = basefilename * "_density" * ext
        @rput fname
        R"ggsave(fname, pdens)"
    end
end
