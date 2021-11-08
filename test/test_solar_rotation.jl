#=
if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("SunAsAStar")   end
if occursin(r"SunAsAStar", pwd())   cd("test")   end
 using Pkg
 Pkg.activate(".")
=#

using SunAsAStar
using Test
#using CSV
using DataFrames

make_plots = false
if make_plots
   using Plots
end

@testset "Compare to Collier-Cameron table" begin

   HARPSN_df = SunAsAStar.CSV.read("data/Sun_harpsn_qualflag.csv",DataFrame)
   jd = HARPSN_df.JD .+ 2400000
   HARPSN_df.Delta = HARPSN_df.FWHMobs.^2 .- HARPSN_df.FWHMsid.^2

   # Test scalar
   @test SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs(jd[1], obs=:HARPSN) ≈ -0.5967647214702758

   # Test on array of times
   Delta = SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs.(jd[1:50:end],obs=:HARPSN)

#Delta = SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs(jd,obs=:HARPSN)
#jdsmall = copy(jd[1:50:end])
#Delta = SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs(jdsmall,obs=:HARPSN)

   if make_plots
      using Plots
      plot(jd,HARPSN_df.Delta, label="ACC")
      scatter!(jd[1:50:end],Delta, ms=2,label="Julia")
      xlabel!("Time (d)")
      ylabel!("FWHM_obs^2 - FWHM_sid^2")

   end

   if make_plots
      plot(jd[1:50:end],Delta.-HARPSN_df.Delta[1:50:end], label="Julia-ACC")
      xlabel!("Time (d)")
      ylabel!("Residuals in Δ(FWHM^2)")
   end


   if make_plots && false
      plot(jd[1:50:end],Delta./HARPSN_df.Delta[1:50:end], label="Julia/ACC")
      xlabel!("Time (d)")
      ylabel!("Ratio of Δ(FWHM^2)'s")
   end
   @test maximum(abs.(HARPSN_df.Delta[1:50:end].-Delta)) < 0.01
end
