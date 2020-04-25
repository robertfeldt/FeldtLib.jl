using Serialization

# To map Swedish postal numbers to Long and Lat we can use the export dumps
# of Geonames.org available here:
const GeonamesPostalDumpBaseUrl = "https://download.geonames.org/export/zip/"
# However, the resolution of this data is poor and most postal numbers are just
# mapped to the city that is closest to them.

function get_latest_geonames_postal_country_dump(country = "SE", redownload = false)
    filename = country * ".zip"
    url = GeonamesPostalDumpBaseUrl * filename
    destdir = joinpath(dirname(@__FILE__()), "..", "data")
    destfile = joinpath(destdir, filename)
    resultfile = joinpath(destdir, country * ".txt")
    if redownload || !isfile(resultfile)
        run(`curl -o $destfile $url`)
        if isfile(destfile)
            cd(destdir) do
                run(`unzip -u $destfile`)
            end
            if isfile(resultfile)
                cd(destdir) do
                    rm(destfile)
                    rm("readme.txt")
                end
                return resultfile
            else
                error("Couldn't unpack file $destfile")
            end
        else
            error("Couldn't download from $url")
        end
    end
    return resultfile
end

function parse_geonames_country_dump(rawdatafile)
    d = Dict{String,Tuple{Float64, Float64, String, String, String}}()
    for line in readlines(rawdatafile)
        parts = split(line, r"\s+")
        postalnumber = replace(strip(parts[2] * parts[3]), r"\s+" => "")
        # Index from the back for lat and long since the number of names
        # in between sometimes varies
        lat = parse(Float64, parts[end-2])
        long = parse(Float64, parts[end-1])
        d[postalnumber] = (lat, long, parts[4], parts[5], parts[6])
    end
    d
end

SwedishPostnrDictCache = nothing # This is parsed on first access below

function swedish_postnr_dict(redownload = false)
    rawdatafile = get_latest_geonames_postal_country_dump("SE", redownload)
    jdfile = rawdatafile * ".juliadata"
    global SwedishPostalNumbersDictCache
    if !isfile(jdfile)
        SwedishPostalNumbersDictCache = parse_geonames_country_dump(rawdatafile)
        serialize(jdfile, SwedishPostalNumbersDictCache)
        return SwedishPostalNumbersDictCache
    else
        try
            return deserialize(jdfile)
        catch _err
            dict = parse_geonames_country_dump(rawdatafile)
            serialize(dict, jdfile)
            return dict
        end
    end
end

function swedish_postnr2position(postnr::String, redownload = false)
    postnr = replace(strip(postnr), r"\s+" => "")
    return swedish_postnr_dict(redownload)[postnr]
end

swedish_postnr2position(postnr::Int, redownload = false) =
    swedish_postnr2position(string(postnr), redownload)