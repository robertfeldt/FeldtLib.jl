module FeldtLib

export setgooglemapsapikey!, swedish_postnr2position
export ggplot2_distribution

#include("ranks.jl")
#include("baumgartner_weis_schindler.jl")
#include("fitall.jl")
include("cacheddict.jl")
include("geocoding.jl")
include("geonames.jl")
include("ggplot2_graphs.jl")

end