__precompile__() # this module is safe to precompile
module SunAsAStar

include("common.jl")
export get_obs_loc, get_lst, get_solar_ra_dec, get_solar_hour_angle, get_earth_sun_dist, get_solar_info

# For interface to JPL Horizons & caching to disk
include("horizons.jl")
export Horizons

# For computing differential extinction
include("diff_solar_extinction.jl")
export DifferentialExtinction
export calc_Δv_diff_extinction

# For computing effects of apparent solar rotation rate
include("solar_rotation.jl")
export SolarRotation
export calc_Δfwhm_solar_rotation
#export CalculateFWHMDifference_SolarRotation_from_obs, CalculateFWHMDifference_SolarRotation_from_long_lat_alt
#export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial

# Temporary workaround since need to pin CSV
using CSV
export CSV

end
