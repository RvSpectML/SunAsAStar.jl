using SunAsAStar
using Test
using Dates

@testset "SunAsAStar.jl" begin
    # Write your tests here.
    using SunAsAStar

    @testset "Common" begin
        @test_nowarn get_obs_loc(:WIYN)
        @test_nowarn get_obs_loc(:NEID)
        @test_nowarn get_obs_loc(:DCT)
        @test_nowarn get_obs_loc(:EXPRES)
        @test_nowarn get_obs_loc(:LAPALMA)
        @test_nowarn get_obs_loc(:HARPSN)
        @test_nowarn get_lst(Dates.datetime2julian(DateTime(2021,1,1)), obs=:WIYN )
        @test_nowarn get_lst([Dates.datetime2julian(DateTime(2020,1,1)), Dates.datetime2julian(DateTime(2021,1,1))], obs=:WIYN )
    end

    t = Dates.datetime2julian(DateTime(2021,1,1))
    @testset "Solar rotation" begin
        @test_nowarn SolarRotation.calc_Δfwhm_solar_rotation(t, obs=:WIYN)

    end
    @testset "Differential extinction" begin
        @test_nowarn DifferentialExtinction.calc_Δv_diff_extinction(t, obs=:WIYN)
        tlist = collect(t .+ range(0.0,stop=1.0,step=0.05))
        output = DifferentialExtinction.compare_Horizons_Astropy(tlist, :WIYN, max_airmass=4.0)
        # TODO: Improve agreement.  Likely due to astropy breaking refraction=true option
        @test maximum(abs.(output.delta_vr_old .- output.delta_vr_new)) < 0.01
    end

    include("test_solar_rotation.jl")
end
