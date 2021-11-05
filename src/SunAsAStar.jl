module SunAsAStar

using Dates, LinearAlgebra

# For computing effects of apparent solar rotation rate
include("solar_rotation.jl")
export SolarRotation
export CalculateFWHMDifference_SolarRotation_from_long_lat_alt, CalculateFWHMDifference_SolarRotation_from_obs
export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial

# For computing differential extinction
include("diff_solar_extinction.jl")
export DifferentialExtinction

end
