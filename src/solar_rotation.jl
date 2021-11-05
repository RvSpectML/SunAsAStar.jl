__precompile__() # this module is safe to precompile
module SolarRotation
using PyCall


const AstropyCoordinates = PyNULL()
const AstropyTime = PyNULL()
const Barycorrpy = PyNULL()
const SciPySpatialTransform = PyNULL()

function __init__()
	copy!(AstropyCoordinates , pyimport("astropy.coordinates") )
	copy!(AstropyTime , pyimport("astropy.time") )
	copy!(Barycorrpy , pyimport("barycorrpy") )
	# Uncomment below and remove pyimport from CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic if start using it
	#copy!(SciPySpatialTransform , pyimport("scipy.spatial.transform") )
end


using Dates, LinearAlgebra
using DataFrames

"""
  `CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic(jd; loc )`
Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
Inputs:
	jd: Julian Date (scalar)
	loc: Astropy Earth Location object. https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html

Output:
	Delta:  F_obs**2 - F_sid**2 [(km/s)^2]

Based on Colier Cameron et al. (2019)
Uses ses JPL ephemeris rather than assuming a Keplerian orbit
Performs calculations in ecliptic reference frame.

Python version by Shubham Kanodia March 2021 (with input from Jason Wright and Eric Ford)
Translated to Julia by Eric Ford March 9, 2021.
"""
function CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic(jd::Real; loc::PyObject )
#function CalculateFWHMDifference_SolarRotationGen(loc::PyObject, jd::AVR ) where { T<:Real, AVR<:AbstractVector{T} } # JDUTC::PyObject)
	# Functions and constants from AstroPy or SciPy that we'll use multiple times
	get_body_barycentric_posvel = AstropyCoordinates.get_body_barycentric_posvel
	get_body_barycentric = AstropyCoordinates.get_body_barycentric
	speed_of_light = Barycorrpy.PhysicalConstants.c
	R_sun = 695700  # km # should we use better, but non-standard value, e.g, 695660?

	SciPySpatialTransform = pyimport("scipy.spatial.transform")
	SciPyRotation = SciPySpatialTransform.Rotation

	JDUTC = AstropyTime.Time(jd , format="jd", scale="utc")
	# Convert times to obtain TDB and TT
	JDTDB = JDUTC.tdb
	JDTT = JDUTC.tt

	################################
	######EARTH EPHEMERIS ############
	################################
	##### NUTATION, PRECESSION, ETC. #####

	# Observatory position wrt Geocenter
	r_pint, v_pint = Barycorrpy.PINT_erfautils.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)
	r_eci = r_pint'  # [m]
	v_eci = v_pint'  # [m/s]

	##### EPHEMERIDES #####
	ephemeris = "de430"
	earth_geo = get_body_barycentric_posvel("earth", JDTDB, ephemeris=ephemeris) # [km]
	r_geo = earth_geo[1].xyz*1000.0 # [m]
	v_geo = earth_geo[2].xyz*1000.0/86400.0  # [m/s]

	PosVector_EarthSSB = r_eci + r_geo # [m]MethodError

	# Relativistic Addition of Velocities
	VelVector_EarthSSB = (v_eci + v_geo) / (1.0 + dot(v_eci,v_geo)/speed_of_light^2) # [m/s]

	################################
	######SOLAR EPHEMERIS ############
	################################
	solar_ephem = get_body_barycentric_posvel("sun", JDTDB, ephemeris=ephemeris)
	PosVector_SolSSB = solar_ephem[1].xyz*1000. #[m]
	VelVector_SolSSB = solar_ephem[2].xyz*1000.0/86400.0  # [m/s]

	################################
	####EQUATORIAL COORD VECTORS ####
	################################
	PosVector_EarthSol = PosVector_EarthSSB - PosVector_SolSSB
	PosMag_EarthSol = norm(PosVector_EarthSol)

	VelVector_EarthSol = (VelVector_EarthSSB - VelVector_SolSSB) / (1.0 + dot(VelVector_SolSSB,VelVector_EarthSSB)/speed_of_light^2)

	#OmegaVector = map(i->cross(PosVector_EarthSol[:,i], VelVector_EarthSol[:,i]) ./ dot(PosVector_EarthSol[:,i],PosVector_EarthSol[:,i]), 1:size(PosVector_EarthSol,2))
	OmegaVector = cross(vec(PosVector_EarthSol), vec(VelVector_EarthSol)) ./ dot(PosVector_EarthSol,PosVector_EarthSol)
	################################
	################################

	#DeltaCentury = map(i->(Dates.year(JDUTC.datetime) .- 2000)./100, 1:length(jd) )
	DeltaCentury = (Dates.year(JDUTC.datetime) - 2000)/100

	# https://www2.mps.mpg.de/homes/fraenz/systems/systems3art.pdf

	#Eqn 14
	SolarInclination = 7.25
	SolarLongitude = 75.76 + 1.397*DeltaCentury
	EclipticEpsilon = 23.44

	# Need to perform Extrinsic Euler Rotations to rotate from equatorial to ecliptic
	#REcliptic = R.from_euler("X",  -EclipticEpsilon, degrees=True).as_matrix()
	REcliptic = SciPyRotation.from_euler("X",  -EclipticEpsilon, degrees=true).as_dcm()
	# may need to change to .as_matrix() with SciPy 1.4?

	# Intrinsic rotation to go from solar axis to ecliptic
	#RObliquity = R.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_matrix()
	RObliquity = SciPyRotation.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_dcm()

	################################
	####ROTATED COORD VECTORS ####
	################################

	#RotatedPositionVector = map(i->REcliptic * PosVector_EarthSol[:,i], 1:size(PosVector_EarthSol,2) )
	#RotatedPositionHat = RotatedPositionVector ./ norm.(RotatedPositionVector)
	RotatedPositionVector = REcliptic * PosVector_EarthSol
	RotatedPositionHat = RotatedPositionVector / norm(RotatedPositionVector)

	#OmegaEarthVectorRotated = map(i->REcliptic * OmegaVector[i],1:size(OmegaVector,2))
	OmegaEarthVectorRotated = REcliptic * OmegaVector
	################################
	################################

	OmegaSolVector = [0.0, 0.0,  2.972 *1e-6]
	#OmegaSolHat = [0.0 , 0.0 , 1.0]

	# Rotated Solar rotation vector to ecliptic plane
	#OmegaSolVector = map(i->RObliquity * OmegaSolVector[:,i], 1:size(OmegaSolVector,2) )
	#OmegaSolHat = OmegaSolVector./norm.(OmegaSolVector)
	OmegaSolVector = RObliquity * OmegaSolVector
	OmegaSolHat = OmegaSolVector/norm(OmegaSolVector)

	sini = sqrt(1 - dot(OmegaSolHat,RotatedPositionHat)^2)
	Gamma = 1.04
	DeltaOmega = OmegaSolVector .- OmegaEarthVectorRotated

	Delta = ((Gamma* R_sun)^2) * (dot(DeltaOmega,DeltaOmega)*sini^2 - dot(OmegaSolVector,OmegaSolVector))
end

"""
   `CalculateFWHMDifference_SolarRotation_from_loc_Equatorial(jd; loc )`
Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
Inputs:
  jd: Julian Date (scalar)
  loc: Astropy Earth Location object. https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html

Output:
  Delta:  F_obs**2 - F_sid**2 [(km/s)^2]

Based on Colier Cameron et al. (2019)
Uses ses JPL ephemeris rather than assuming a Keplerian orbit.
Performs calculations in equatorial reference frame.

Python version by Shubham Kanodia (with input from Jason Wright and Eric Ford) March 2021
Translated to Julia by Eric Ford March 9, 2021.
"""
function CalculateFWHMDifference_SolarRotation_from_loc_Equatorial(jd::Real; loc::PyObject )
#function CalculateFWHMDifference_SolarRotationGen(loc::PyObject, jd::AVR ) where { T<:Real, AVR<:AbstractVector{T} } # JDUTC::PyObject)
	# Functions and constants from AstroPy or SciPy that we'll use multiple times
	get_body_barycentric_posvel = AstropyCoordinates.get_body_barycentric_posvel
	get_body_barycentric = AstropyCoordinates.get_body_barycentric
	speed_of_light = Barycorrpy.PhysicalConstants.c
	R_sun = 695700  # km # should we use better, but non-standard value, e.g, 695660?

	JDUTC = AstropyTime.Time(jd , format="jd", scale="utc")
	#return JDUTC
	# Convert times to obtain TDB and TT
	JDTDB = JDUTC.tdb
	JDTT = JDUTC.tt

	################################
	######EARTH EPHEMERIS ############
	################################
	##### NUTATION, PRECESSION, ETC. #####

	# Observatory position wrt Geocenter
	r_pint, v_pint = Barycorrpy.PINT_erfautils.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)
	r_eci = r_pint'  # [m]
	v_eci = v_pint'  # [m/s]

	##### EPHEMERIDES #####
	ephemeris = "de430"
	earth_geo = get_body_barycentric_posvel("earth", JDTDB, ephemeris=ephemeris) # [km]
	r_geo = earth_geo[1].xyz*1000.0 # [m]
	v_geo = earth_geo[2].xyz*1000.0/86400.0  # [m/s]

	PosVector_EarthSSB = r_eci + r_geo # [m]

	# Relativistic Addition of Velocities
	VelVector_EarthSSB = (v_eci + v_geo) / (1.0 + dot(v_eci,v_geo)/speed_of_light^2) # [m/s]

	################################
	######SOLAR EPHEMERIS ############
	################################
	solar_ephem = get_body_barycentric_posvel("sun", JDTDB, ephemeris=ephemeris)
	PosVector_SolSSB = solar_ephem[1].xyz*1000. #[m]
	VelVector_SolSSB = solar_ephem[2].xyz*1000.0/86400.0  # [m/s]

	################################
	####EQUATORIAL COORD VECTORS ####
	################################
	PosVector_EarthSol = PosVector_EarthSSB - PosVector_SolSSB
	PosMag_EarthSol = norm(PosVector_EarthSol)
	PosHat_EarthSol = PosVector_EarthSol ./ PosMag_EarthSol
	VelVector_EarthSol = (VelVector_EarthSSB - VelVector_SolSSB) / (1.0 + dot(VelVector_SolSSB,VelVector_EarthSSB)/speed_of_light^2)

	#OmegaVector = map(i->cross(PosVector_EarthSol[:,i], VelVector_EarthSol[:,i]) ./ dot(PosVector_EarthSol[:,i],PosVector_EarthSol[:,i]), 1:size(PosVector_EarthSol,2))
	OmegaVector = cross(vec(PosVector_EarthSol), vec(VelVector_EarthSol)) ./ dot(PosVector_EarthSol,PosVector_EarthSol)

	################################
	################################
	# Rotation Axis of the Sun with respect to Equatorial System
	# https://www2.mps.mpg.de/homes/fraenz/systems/systems3art.pdf

	#Eqn 13
	alpha = 286.13 # Rotate along the z axis
	dec = 63.87 # Rotate along the x axis
	Theta = alpha * π/180
	Phi = dec * π/180
	EclipticEpsilon = 23.44
	SolarInclination = 7.25
	SolarLongitude = 75.76

	# Need to perform Extrinsic Euler Rotations
	# Transform to solar rotation axis
	REquatorial = [cos(Theta)*cos(Phi), sin(Theta)*cos(Phi), sin(Phi)]

	OmegaSolVector = [0.0, 0.0,  2.972 *1e-6]
	#OmegaSolHat = [0.0 , 0.0 , 1.0]

	OmegaSolHat = REquatorial
	OmegaSolVector = REquatorial * norm(OmegaSolVector)
	#OmegaSolHat = OmegaSolVector / norm(OmegaSolVector)

	sini = sqrt(1 - dot(OmegaSolHat,PosHat_EarthSol)^2)
	Gamma = 1.04
	DeltaOmega = OmegaSolVector .- OmegaVector

	Delta = ((Gamma* R_sun)^2) * (dot(DeltaOmega,DeltaOmega)*sini^2 - dot(OmegaSolVector,OmegaSolVector))

end

"""
	`CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd; long, lat, alt)`
Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
	Inputs:
		jd: Julian Date (scalar)
		long: longitude (degrees, west negative)
		lat: latitude (degrees)
		alt: altitude (m)
"""
function CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd::Real; long::Real, lat::Real, alt::Real )
	loc = AstropyCoordinates.EarthLocation.from_geodetic(long, lat, height=alt)
	delta = CalculateFWHMDifference_SolarRotation_from_loc_Equatorial(jd, loc=loc)
	#delta = CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic(jd, loc=loc)
	return delta
end


"""
	`CalculateFWHMDifference_SolarRotation_from_obs(jd; obs)`
Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
Inputs:
	jd: Julian Date (scalar)
	obs: observatory name (symbol)
Currently has HARPS-N and NEID/WIYN coordinates.
"""
function CalculateFWHMDifference_SolarRotation_from_obs(jd::Real; obs::Symbol )
	loc = get_obs_loc(obs)
	#=
	if obs == :HARPSN
		longi = -17.88905
		lat = 28.754
		alt = 2387.2
	elseif  obs == :WIYN
		longi = -111.600562
		lat =  31.958092
		alt =  2091.0
	else
		@error("Don't have coordinates for obs = " * string(obs))
	end
	=#
	CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd, long=loc["lon"], lat=loc["lat"], alt=loc["elevation"])
end

#=
function get_wiyn_loc()
	obsname = "WIYN"
	longi = -111.600562
	lat =  31.958092
	alt =  2091.0
	AstropyCoordinates.EarthLocation.from_geodetic(longi, lat, height=alt)
end

function get_harpsn_loc()
	obsname = "HARPS-N"
	longi = -17.88905
	lat = 28.754
	alt = 2387.2
	loc = AstropyCoordinates.EarthLocation.from_geodetic(longi, lat, height=alt)
end
=#

export CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic, CalculateFWHMDifference_SolarRotation_from_loc_Equatorial
export CalculateFWHMDifference_SolarRotation_from_long_lat_alt, CalculateFWHMDifference_SolarRotation_from_obs

#export get_wiyn_loc, get_harpsn_loc

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
	if type(jd) <:Real
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
		result = DataFrame(:ra=>RA, :dec=>Dec, :alt=>Alt, :az=>Az, :airmass=>Airmass, :lst=>LST, :sol_dist_au=>dist )
	end
	return result
end

export get_obs_loc, get_lst, get_solar_hour_angle, get_solar_ra_dec, get_earth_sun_dist, get_solar_info

end # module SolarRotation
