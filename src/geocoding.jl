using HTTP
using JSON

export get_lat_lon, vgregion_lat_lon_bounds

#nominatim requires as user agent and will block you otherwise
HTTP.setuseragent!("Mozilla/5.0 (Windows NT 10.0; rv:68.0) Gecko/20100101 Firefox/68.0")

const NominatimSearchPrefix = "https://nominatim.openstreetmap.org/search?q="

# Nominatim will block you if you access more than once per second so keep track 
# of time.
struct NominatimSearcher
    lasttime::Float64
    NominatimSearcher() = new(0.0)
end

const DefaultSearcher = NominatimSearcher()

"""
    get_lat_lon(address::String)

Returns (lat,lon) coordinate for address (typically in middle of street)
as well as additional data returned from Nominatim.
"""
function get_lat_lon(address::String, n::NominatimSearcher = DefaultSearcher)
    addressFormatted = join(split(address),"+")
    searchurl = NominatimSearchPrefix * addressFormatted * "&format=json&limit=1"

    if time() - n.lasttime < 1.0
        sleep(1.05 - (time() - n.lasttime))
        n.lasttime = time()
    end

    rawdata = HTTP.get(searchurl)

    metadata = JSON.parse(String(rawdata.body))[1]
    coordinates = parse(Float64, metadata["lat"]), parse(Float64, metadata["lon"])

    return coordinates, metadata
end

# Yttergränser för Västra Götaland från Google Maps 2020-04-22:
const VGRegion = Dict(
    :south => (57.140190, 13.127416), # Gislaved V
    :east  => (58.687169, 14.714468), # Karlsborg N
    :north => (59.260856, 12.208255), # Årjäng SO
    :west  => (58.887149, 10.963946)  # väster om Nordkoster
)

function vgregion_lat_lon_bounds(delta = 0.10)
    lats = map(first, values(VGRegion))
    lons = map(t -> t[2], values(VGRegion))

    # Boundingbox for ggmap given as lowerleftlon, lowerleftlat, upperightlon, upperrightlat
    return (minimum(lons)-delta, minimum(lats)-delta, maximum(lons)+delta, maximum(lats)+delta)
end
