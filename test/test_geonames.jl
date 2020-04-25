using FeldtLib: swedish_postnr2position

@testset "postal numbers" begin

# Can't run this since you need to manually set your Google Maps API key...
@testset "test one known postal number" begin
    #lat, lng = swedish_postnr2position("417 56")
    #@test lat ~ 57.7121456
    #@test lat ~ 11.9484268
end

end