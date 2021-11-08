# Intiailize Python packages
using PyCall
const AstropyCoordinates = PyNULL()
const AstropyTime = PyNULL()
const Barycorrpy = PyNULL()

function __init__()
	copy!(AstropyCoordinates , pyimport("astropy.coordinates") )
	copy!(AstropyTime , pyimport("astropy.time") )
	copy!(Barycorrpy , pyimport("barycorrpy") )
	@assert AstropyCoordinates != PyNULL()
	@assert AstropyTime != PyNULL()
	@assert Barycorrpy != PyNULL()
end

# Constants
AU = 149597870.7  # km
R_sun = 695700.0  # km # should we use better, but non-standard value, e.g, 695660?
valid_obs = [:WIYN, :NEID, :LAPALMA, :HARPSN, :DCT, :EXPRES]

#using Dates, LinearAlgebra
#using DataFrames

""" get_obs_loc(obs::Symbol)
 Returns a Dict with long & lat (degrees) and elevation (km)
 Warning: Currently only has info for :WIYN and :HARPSN.
"""
function get_obs_loc(obs::Symbol)
	valid_obs = [:WIYN, :NEID, :LAPALMA, :HARPSN, :DCT, :EXPRES]
	@assert obs ∈ valid_obs
	if obs == :LAPALMA || obs == :HARPSN
		loc = Dict("lon"=> -17.88905, "lat"=> 28.754, "elevation"=> 2.3872)
	elseif obs == :WIYN || obs == :NEID
		loc = Dict("lon"=> -111.600562, "lat"=> 31.958092, "elevation"=> 2.091)
	elseif obs == :DCT || obs == :EXPRES
		loc = Dict("lon"=> -111.421944, "lat"=> 34.744444, "elevation"=> 2.219)
	else
		@error("Don't have coordinates for obs = " * string(obs))
	end
	return loc
end

function get_lst(jd::Union{Real,Vector{<:Real}}; obs::Symbol)
   loc = get_obs_loc(obs)
   pyloc = AstropyCoordinates.EarthLocation.from_geodetic(loc["lon"], loc["lat"], height=loc["elevation"])
   TimeObs = AstropyTime.Time(jd, format="jd", scale="utc", location=pyloc)
   LST = TimeObs.sidereal_time("mean")
   if typeof(jd) <: Real
	   return first(LST)
   else
	   return LST
   end
end

function get_solar_ra_dec(jd::Union{Real,Vector{<:Real}}; obs::Symbol)
	jlloc = get_obs_loc(obs)
    pyloc = AstropyCoordinates.EarthLocation.from_geodetic(jlloc["lon"], jlloc["lat"], height=jlloc["elevation"]) #*1000)
    TimeObs = AstropyTime.Time(jd, format="jd", scale="utc", location=pyloc)
	SolarObj = AstropyCoordinates.get_sun(TimeObs)
	RA, Dec = SolarObj.ra, SolarObj.dec
	return (ra=vec(RA), dec=vec(Dec))
end

function get_solar_hour_angle(jd::Union{Real,Vector{<:Real}}; obs::Symbol)
	lst = get_lst(jd, obs=obs)
	ra = get_solar_ra_dec(jd, obs=obs).ra
	ha = lst .- ra./15
	ha[ha.<-12] .+= 24
	ha[ha.>12] .-= 24
	return ha
end

function get_earth_sun_dist(jd::Union{Real,Vector{<:Real}}; obs::Symbol)
	jlloc = get_obs_loc(obs)
    loc = AstropyCoordinates.EarthLocation.from_geodetic(jlloc["lon"], jlloc["lat"], height=jlloc["elevation"])
	ephemeris = "de430"
	JDUTC = AstropyTime.Time(jd , format="jd", scale="utc")
	# Convert times to obtain TDB and TT
	JDTDB = JDUTC.tdb
	JDTT = JDUTC.tt

	get_body_barycentric_posvel = AstropyCoordinates.get_body_barycentric_posvel
	earth_geo = get_body_barycentric_posvel("earth", JDTDB, ephemeris=ephemeris)[1].xyz * 1000.0 # [km->m]
	solar_ephem = get_body_barycentric_posvel("sun", JDTDB, ephemeris=ephemeris)[1].xyz * 1000.0 # [km->m]
	r_pint, v_pint = Barycorrpy.PINT_erfautils.gcrs_posvel_from_itrf(loc, JDUTC, JDTT) # [m]
	diff_vector = (earth_geo + r_pint') -solar_ephem #[m]
	dist = vec(sqrt.(sum(diff_vector.^2,dims=1)))
	AU = 149597870.7  # km
	dist = dist/(AU*1000)
	if typeof(jd) <:Real
		return first(dist)
	else
		return dist
	end
end

function get_solar_info(jd::Union{Real,Vector{<:Real}}; obs::Symbol)
	jlloc = get_obs_loc(obs)
    pyloc = AstropyCoordinates.EarthLocation.from_geodetic(jlloc["lon"], jlloc["lat"], height=jlloc["elevation"]) #*1000)
    TimeObs = AstropyTime.Time(jd, format="jd", scale="utc", location=pyloc)
	SolarObj = AstropyCoordinates.get_sun(TimeObs)
	frame_obsnight = AstropyCoordinates.AltAz(obstime=TimeObs, location=pyloc)
	SolarAltAz = SolarObj.transform_to(frame_obsnight)
	RA, Dec = SolarObj.ra, SolarObj.dec
	Alt, Az = SolarAltAz.alt, SolarAltAz.az
	Airmass = SolarAltAz.secz
	LST = TimeObs.sidereal_time("mean")

	ephemeris = "de430"
	JDUTC = AstropyTime.Time(jd , format="jd", scale="utc")
	# Convert times to obtain TDB and TT
	JDTDB = JDUTC.tdb
	JDTT = JDUTC.tt
	r_pint, v_pint = Barycorrpy.PINT_erfautils.gcrs_posvel_from_itrf(pyloc, JDUTC, JDTT) # [m]
	get_body_barycentric_posvel = AstropyCoordinates.get_body_barycentric_posvel
	earth_geo = get_body_barycentric_posvel("earth", JDTDB, ephemeris=ephemeris)[1].xyz * 1000.0 # [km->m]
	solar_ephem = get_body_barycentric_posvel("sun", JDTDB, ephemeris=ephemeris)[1].xyz * 1000.0 # [km->m]
	diff_vector = (earth_geo + r_pint') -solar_ephem #[m]
	dist = sqrt.(sum(diff_vector.^2,dims=1))
	AU = 149597870.7  # km
	dist = dist./(AU*1000)

	if typeof(jd) <: Real
		result = Dict(:ra=>first(RA), :dec=>first(Dec), :alt=>first(Alt), :az=>first(Az), :airmass=>first(Airmass), :lst=>first(LST), :sol_dist_au=>first(dist) )
	else
		#result = DataFrame(:ra=>RA, :dec=>Dec, :alt=>Alt, :az=>Az, :airmass=>Airmass, :lst=>LST, :sol_dist_au=>dist )
		result = (;ra=RA, dec=Dec, alt=Alt, az=Az, airmass=Airmass, lst=LST, sol_dist_au=dist )
	end
	return result
end

#export get_obs_loc, get_lst, get_solar_hour_angle, get_solar_ra_dec, get_earth_sun_dist, get_solar_info
