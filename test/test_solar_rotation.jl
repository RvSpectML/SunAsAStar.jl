if occursin(r"RvSpectMLEcoSystem$", pwd())   cd("SunAsAStar")   end
 if occursin(r"SunAsAStar$", pwd())   cd("test")   end
 using Pkg
 Pkg.activate(".")

using SunAsAStar
using CSV, DataFrames

HARPSN_df = CSV.read("../data/Sun_harpsn_qualflag.csv",DataFrame)
jd = HARPSN_df.JD .+ 2400000
HARPSN_df.Delta = HARPSN_df.FWHMobs.^2 .- HARPSN_df.FWHMsid.^2

# Test scalar
Delta = SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs(jd[1], obs=:HARPSN)

# Test on array of times
Delta = SolarRotation.CalculateFWHMDifference_SolarRotation_from_obs.(jd[1:50:end],obs=:HARPSN)

using Plots
plot(jd,HARPSN_df.Delta, label="ACC")
plot!(jd[1:50:end],Delta, label="Julia")
xlabel!("Time (d)")
ylabel!("FWHM_obs^2 - FWHM_sid^2")

plot(jd[1:50:end],Delta.-HARPSN_df.Delta[1:50:end], label="Julia-ACC")
xlabel!("Time (d)")
ylabel!("Residuals in Δ(FWHM^2)")

#=
plot(jd[1:50:end],Delta./HARPSN_df.Delta[1:50:end], label="Julia/ACC")
xlabel!("Time (d)")
ylabel!("Ratio of Δ(FWHM^2)'s")
=#
