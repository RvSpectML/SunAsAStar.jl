"""
   DifferentialExtinction
   Original Python code (solar_diffex.py) by Andrea Lin
   Ported to Julia by Eric Ford
   Based on Collier-Cameron et al. (2019) MNRAS 487, 1082.  doi:10.1093/mnras/stz1215
   See also Astronomical Algorithms (Meeus, 1991) p151 & p178
"""
__precompile__() # this module is safe to precompile
module Horizons

using PyCall
using ..SunAsAStar

using Dates
using DataFrames, CSV

# Import Python package for querying JPL Horizons database
const JplHorizons = PyNULL()

function __init__()
	copy!(JplHorizons , pyimport("astroquery.jplhorizons") )
end


function query_horizons(;epochs::Union{Vector{Float64},Dict}, obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	if typeof(epochs) <: Dict
		@assert haskey(epochs,"start")
		@assert haskey(epochs,"stop")
		@assert haskey(epochs,"step")
	elseif typeof(epochs) <: Vector{Float64}
		@assert 1 <= length(epochs) <= 10000
	end
	@assert 1 < max_airmass <= 10.0
	loc = get_obs_loc(obs)
	# Sometimes this fails.  Why?
	obj = JplHorizons.Horizons(id="10", location=loc, epochs=epochs, id_type="majorbody")
	if verbose println("# obj = ", obj) end
	#eph = obj.ephemerides(quantities="1,4,8,17,20,42", airmass_lessthan=max_airmass, refraction=true)  # WARNING: refraction=true stopped working! Why?
	eph = obj.ephemerides(quantities="1,4,8,17,20,42", airmass_lessthan=max_airmass )  # WARNING: refraction=true stopped working! Why?
	if verbose  println("# eph = ", eph)  end
	return eph
end


function cache_horizons_output_year(year::Int64, obs::Symbol)
	fn = "horizons_" * string(obs) * "_" * string(year) * ".csv"
	filename = joinpath(dirname(pathof(SunAsAStar)),"..","data","horizons",fn)
	if isfile(filename) && fielsize(filename) > 0
		return
	end
	df = DataFrame(:jd=>Float64[], :sol_npole_ang=>Float64[], :sol_dist_au=>Float64[])
	epoch_start = string(year) * "-01-01 00:00:00"
	epoch_stop  = string(year+1) * "-01-01 00:00:00"
	step = "6h"
	epochs = Dict("start"=>epoch_start, "stop"=>epoch_stop, "step"=>step)
	eph = query_horizons(epochs=epochs, obs=obs, max_airmass=10.0, verbose=true)
	df = DataFrame(:jd=>eph.columns["datetime_jd"], :sol_npole_ang=>eph.columns["NPole_ang"], :sol_dist_au=>eph.columns["delta"])
	CSV.write(filename,df)
	return df
end

function read_horizons_output(jd::Union{Float64,Vector{Float64}}, obs::Symbol)
	#@assert obs ∈ valid_obs
	years = unique(year.(Dates.julian2datetime.(jd)))
	@assert length(years) >= 1
	cache = Dict{Int64,Any}()
	for year in years
		cache[year] = read_horizons_output_year(year, obs=obs)
		#=
		fn = "horizons_" * string(obs) * "_" * string(year) * ".csv"
		filename = joinpath(dirname(pathof(SunAsAStar)),"..","data","horizons",fn)
		if isfile(filename)
			cache[year] = CSV.read(filename, DataFrame)
		else
			cache[year] = cache_horizons_output_year(year, obs)
		end
		cache[year].sol_npole_ang[cache[year].sol_npole_ang .> 180] .-= 360
		=#
	end
	if length(years) == 1
		return cache[first(keys(cache))]
	else
		sorted_keys = sort(collect(keys(cache)))
		output = cache[first(sorted_keys)]
		for k in sorted_keys[2:end]
			append!(output,cache[k])
		end
		return output
	end
end

function read_horizons_output_year(year::Integer ; obs::Symbol)
	fn = "horizons_" * string(obs) * "_" * string(year) * ".csv"
	filename = joinpath(dirname(pathof(SunAsAStar)),"..","data","horizons",fn)
	if isfile(filename) && filesize(filename) > 0
		cache = CSV.read(filename, DataFrame)
	else
		cache = cache_horizons_output_year(year, obs)
	end
	cache.sol_npole_ang[cache.sol_npole_ang .> 180] .-= 360
	return cache
	#first(keys(cache))
end

function get_horizons_output(jd::Float64; obs::Symbol)
  df = read_horizons_output(jd, obs)
  idxlo = searchsortedlast(df.jd, jd)
  @assert 1 <= idxlo < size(df,1)
  idxhi = idxlo + 1
  result = Dict{Symbol,Float64}()
  num = jd-df.jd[idxlo]
  denom = df.jd[idxhi]-df.jd[idxlo]
  whi = num/denom
  wlo = 1-whi
  for k in Symbol.(names(df))
	if k == :jd continue end
  	result[k] = df[idxlo,k] * wlo +  df[idxhi,k] * whi
  end
  result[:jd] = jd
  return result
end

#function get_horizons_output(jd::Vector{Float64}; obs::Symbol)
#  get_horizons_output.(jd,obs)
  #=
  df = read_horizons_output(jd, obs)
  idxlo = map(t->searchsortedlast(df.jd, t),jd)
  mask = idxlo .>= size(df,1)
  if length(mask) >= 1
  	println("# Maximum Δt = ", maximum(jd[mask] .- df.jd[end]) )
  	idxlo[mask] .= size(df,1)-1
end
  @assert all(1 .<= idxlo .<= size(df,1))
  result = DataFrame(:jd=>jd)
  whi = (jd.-df.jd[idxlo]) ./ (df.jd[idxlo.+1].-df.jd[idxlo])
  #wlo = 1.-whi
  for k in  Symbol.(names(df))
	if k == :jd continue end
  	result[!,k] = df[idxlo,k] .* (1.0.-whi) .+  df[idxlo.+1,k] .* whi
  end
  return result
  =#
#end

export query_horizons, get_horizons_output

end # module Horizons
