using FeldtLib: get_lat_lon, vgregion_lat_lon_bounds

@testset "geocoding" begin

@testset "test one known address" begin
    bgc, bgd = get_lat_lon("Bratteråsgatan, Sweden")
    # Manually taken from Nominatim 2020-04-23:
    @test bgc[1] == 57.7021997
    @test bgc[2] == 11.9131308
end

@testset "check that address of region Västra Götaland is within its bounding box" begin
    vgrcoord, vgrd = get_lat_lon("Västra Götaland, Sweden")
    lowleftlon, lowleftlat, uprightlon, uprightlat = vgregion_lat_lon_bounds()
    vgrlat, vgrlon = vgrcoord
    @test lowleftlat < vgrlat < uprightlat
    @test lowleftlon < vgrlon < uprightlon
end

end