using Serialization
using JSON
using HTTP
using ZipFile

const DataDir = joinpath(dirname(@__FILE__()), "..", "data")
if !isdir(DataDir)
    mkdir(DataDir)
end

function http_get_binary_file(url, destfile)
    r = HTTP.request("GET", url)
    open(destfile, "w") do fh
        print(fh, String(r.body))
    end
end

function unzip(filename, writefiles = [])
    r = ZipFile.Reader(filename)
    for f in r.files
        if length(writefiles) == 0 || in(f.name, writefiles)
            open(f.name, "w") do fh
                write(fh, read(f))
            end
        end
    end
end

# To map Swedish postal numbers to Long and Lat we can use the export dumps
# of Geonames.org available here:
const GeonamesPostalDumpBaseUrl = "https://download.geonames.org/export/zip/"
# However, the resolution of this data is poor and most postal numbers are just
# mapped to the city that is closest to them.

function get_latest_geonames_postal_country_dump(country = "SE", redownload = false)
    filename = country * ".zip"
    url = GeonamesPostalDumpBaseUrl * filename
    destdir = DataDir
    if !isdir(destdir)
        mkpath(destdir)
    end
    destfile = joinpath(destdir, filename)
    resultfile = joinpath(destdir, country * ".txt")
    if redownload || !isfile(resultfile)
        http_get_binary_file(url, destfile)
        if isfile(destfile)
            cd(destdir) do
                unzip(destfile)
            end
            if isfile(resultfile)
                cd(destdir) do
                    isfile(destfile) && rm(destfile)
                    isfile("readme.txt") && rm("readme.txt")
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

GeonamesSwedishPostnrDictCache = nothing # This is parsed on first access below

function swedish_postnr_dict(redownload = false)
    rawdatafile = get_latest_geonames_postal_country_dump("SE", redownload)
    jdfile = rawdatafile * ".juliadata"
    global GeonamesSwedishPostnrDictCache
    if !isfile(jdfile)
        GeonamesSwedishPostnrDictCache = parse_geonames_country_dump(rawdatafile)
        serialize(jdfile, GeonamesSwedishPostnrDictCache)
        return GeonamesSwedishPostnrDictCache
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

function geonames_swedish_postnr2position(postnr::String, redownload = false)
    postnr = replace(strip(postnr), r"\s+" => "")
    return swedish_postnr_dict(redownload)[postnr]
end

# Google maps api string:
# https://maps.googleapis.com/maps/api/geocode/json?address=417+56+G%C3%B6teborg&key=YOUR_API_KEY
mutable struct GoogleMapsApi
    key::String
    mindelay::Float64
    lasttime::Float64
    GoogleMapsApi(key::String, mindelay::Float64 = 0.2) = 
        new(key, mindelay, 0.0)
end

function addapikey(g::GoogleMapsApi, prefix, postfix = "")
    prefix * "&key=$(g.key)" * postfix
end

function get(g::GoogleMapsApi, url::String)
    t = time()
    if t - g.lasttime < g.mindelay
        sleep(g.mindelay - (t - g.lasttime))
    end
    g.lasttime = time()
    HTTP.get(url)
end

function getjson(g::GoogleMapsApi, url::String)
    searchurl = addapikey(g, url)
    JSON.parse(String(get(g, searchurl).body))
end

DefaultGoogleMapsApi = nothing

function setgooglemapsapikey!(key::String, defaultmindelay = 0.5)
    HTTP.setuseragent!("Mozilla/5.0 (Windows NT 10.0; rv:68.0) Gecko/20100101 Firefox/68.0")
    global DefaultGoogleMapsApi
    DefaultGoogleMapsApi = GoogleMapsApi(key, defaultmindelay)
end

removews(s) = replace(s, r"\s+" => "")
replacews(s, c="+") = replace(strip(s), r"\s+" => c)

function format_postnummer_for_search(nr::String)
    nr = removews(nr)
    occursin(r"\d{5}", nr) ? (nr[1:3] * "+" * nr[4:5]) : nr
end

function geocode(g::GoogleMapsApi, searchterms...)
    sts = map(replacews, searchterms)
    s = join(sts, "+")
    url = "https://maps.googleapis.com/maps/api/geocode/json?address=$(s)"
    js = getjson(g, url)
    if !isnothing(js) && haskey(js, "results") && (length(js["results"]) >= 1)
        return first(js["results"])
    else
        return nothing
    end
end

function geocode(searchterms...)
    global DefaultGoogleMapsApi
    if isnothing(DefaultGoogleMapsApi)
        error("You must first set the Google Maps API key with setgooglemapsapikey!")
    end
    geocode(DefaultGoogleMapsApi, searchterms...)
end

const SwePostnrDictCache = joinpath(DataDir, "SwePostnrDict.juliadata")
const SwePostnrDict = CachedDict(SwePostnrDictCache) do
    Dict{Union{String,Int},Tuple{Float64, Float64, String}}()
end

function swedish_postnr2position(postnr::String, addr::String = "";
    redownload = false, usegeonames = false)
    shortform = replace(strip(postnr), r"\s+" => "")
    global SwePostnrDict
    if !redownload && haskey(SwePostnrDict, shortform)
        return SwePostnrDict[shortform]
    end
    searchform = format_postnummer_for_search(shortform)
    res = if usegeonames
        lat, long, c1, c2, c3 = geonames_swedish_postnr2position(shortform)
        geocode(searchform, c1, c2, c3)
    else
        geocode(searchform, addr)
    end
    isnothing(res) && return nothing
    loc = res["geometry"]["location"]
    lat = loc["lat"]
    lng = loc["lng"]
    entry = (lat, lng, res["formatted_address"])
    update!(SwePostnrDict, entry, postnr)
    update!(SwePostnrDict, entry, shortform)
    if length(shortform) == 5
        longform = shortform[1:3] * " " * shortform[4:5]
        update!(SwePostnrDict, entry, longform)
    end
    asint = parse(Int, shortform)
    update!(SwePostnrDict, entry, asint, true)
    return entry
end

const AddressDictCache = joinpath(DataDir, "AddressDict.juliadata")
const AddressDict = CachedDict(AddressDictCache) do
    Dict{Union{String,Int},Tuple{Float64, Float64, String}}()
end

address2position(a::String) = address2position([a])

normalizeaddress(elems) = 
    join(map(s -> replace(lowercase(strip(s)), r"\s{2}" => " "), elems), ", ")

function address2position(addressParts::Vector{String}; redownload = false)
    key = normalizeaddress(addressParts)
    global AddressDict
    if !redownload && haskey(AddressDict, key)
        return AddressDict[key]
    end
    res = geocode(addressParts...)
    isnothing(res) && return nothing
    loc = res["geometry"]["location"]
    lat = loc["lat"]
    lng = loc["lng"]
    entry = (lat, lng, res["formatted_address"])
    update!(AddressDict, entry, key)
    return entry
end
