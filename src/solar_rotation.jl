__precompile__() # this module is safe to precompile
module SolarRotation
using PyCall

using ..SunAsAStar
using ..SunAsAStar: AstropyCoordinates, AstropyTime, Barycorrpy

using Dates, LinearAlgebra
using DataFrames

const SciPySpatialTransform = PyNULL()

function __init__()
	# Uncomment below and remove pyimport from CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic if start using it
	copy!(SciPySpatialTransform , pyimport("scipy.spatial.transform") )
	@assert SciPySpatialTransform != PyNULL()
end


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
#function CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic(jd::Union{Real,Vector{<:Real}}; loc::PyObject )
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
	#if typeof(jd) <: Real
		OmegaVector = cross(vec(PosVector_EarthSol), vec(VelVector_EarthSol)) ./ dot(PosVector_EarthSol,PosVector_EarthSol)
	#else
	#	OmegaVector = map(i->cross(PosVector_EarthSol[:,i], VelVector_EarthSol[:,i]) ./ dot(PosVector_EarthSol[:,i],PosVector_EarthSol[:,i]), 1:size(PosVector_EarthSol,2))
	#end
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
	#REcliptic = R.from_euler("X",  -EclipticEpsilon, degrees=true).as_matrix()
	REcliptic = SciPyRotation.from_euler("X",  -EclipticEpsilon, degrees=true).as_matrix()
	#REcliptic = SciPyRotation.from_euler("X",  -EclipticEpsilon, degrees=true).as_dcm()
	# may need to change to .as_matrix() with SciPy 1.4?

	# Intrinsic rotation to go from solar axis to ecliptic
	#RObliquity = R.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_matrix()
	RObliquity = SciPyRotation.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_matrix()
	#RObliquity = SciPyRotation.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_dcm()

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


function CalculateFWHMDifference_SolarRotation_from_loc_Ecliptic(jd::Vector{<:Real}; loc::PyObject )
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

	OmegaVector = map(i->cross(PosVector_EarthSol[:,i], VelVector_EarthSol[:,i]) ./ dot(PosVector_EarthSol[:,i],PosVector_EarthSol[:,i]), 1:size(PosVector_EarthSol,2))
	################################
	################################

	DeltaCentury = map(i->((JDUTC.datetime[i].year) .- 2000)./100, 1:length(jd) )
	#DeltaCentury = map(i->(Dates.year(JDUTC.datetime) .- 2000)./100, 1:length(jd) )
	#DeltaCentury = (Dates.year(JDUTC.datetime) - 2000)/100

	# https://www2.mps.mpg.de/homes/fraenz/systems/systems3art.pdf

	#Eqn 14
	SolarInclination = 7.25
	SolarLongitude = 75.76 .+ 1.397*DeltaCentury
	EclipticEpsilon = 23.44

	# Need to perform Extrinsic Euler Rotations to rotate from equatorial to ecliptic
	#REcliptic = R.from_euler("X",  -EclipticEpsilon, degrees=true).as_matrix()
	#REcliptic = SciPyRotation.from_euler("X",  -EclipticEpsilon, degrees=true).as_dcm()
	REcliptic = SciPyRotation.from_euler("X",  -EclipticEpsilon, degrees=true).as_dcm()
	# may need to change to .as_matrix() with SciPy 1.4?

	# Intrinsic rotation to go from solar axis to ecliptic
	#RObliquity = R.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_matrix()
	#RObliquity = SciPyRotation.from_euler("xz",  [SolarInclination, SolarLongitude], degrees=true).as_dcm()
	RObliquity = map(i->SciPyRotation.from_euler("xz",  [SolarInclination, SolarLongitude[i]], degrees=true).as_dcm(),1:length(SolarLongitude))

	################################
	####ROTATED COORD VECTORS ####
	################################

	RotatedPositionVector = map(i->REcliptic * PosVector_EarthSol[:,i], 1:size(PosVector_EarthSol,2) )
	RotatedPositionHat = RotatedPositionVector ./ norm.(RotatedPositionVector)
	#RotatedPositionVector = REcliptic * PosVector_EarthSol
	#RotatedPositionHat = RotatedPositionVector / norm(RotatedPositionVector)

	OmegaEarthVectorRotated = map(i->REcliptic * OmegaVector[i],1:size(OmegaVector,2))
	#OmegaEarthVectorRotated = REcliptic * OmegaVector
	################################
	################################

	OmegaSolVector = [0.0, 0.0,  2.972 *1e-6]
	#OmegaSolHat = [0.0 , 0.0 , 1.0]

	# Rotated Solar rotation vector to ecliptic plane
	OmegaSolVector = map(i->RObliquity[i] * OmegaSolVector[:,i], 1:size(OmegaSolVector,2) )
	OmegaSolHat = OmegaSolVector./norm.(OmegaSolVector)
	#OmegaSolVector = RObliquity * OmegaSolVector
	#OmegaSolHat = OmegaSolVector/norm(OmegaSolVector)

	sini = sqrt.(1 .- dot.(OmegaSolHat,RotatedPositionHat).^2)
	Gamma = 1.04
	DeltaOmega = OmegaSolVector .- OmegaEarthVectorRotated

	Delta = ((Gamma* R_sun)^2) * (dot.(DeltaOmega,DeltaOmega).*sini.^2 .- dot.(OmegaSolVector,OmegaSolVector))
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

#=
# Still needs work due to pain of Python objects
function CalculateFWHMDifference_SolarRotation_from_loc_Equatorial(jd::Vector{<:Real}; loc::PyObject )
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

	OmegaVector = map(i->cross(PosVector_EarthSol[:,i], VelVector_EarthSol[:,i]) ./ dot(PosVector_EarthSol[:,i],PosVector_EarthSol[:,i]), 1:size(PosVector_EarthSol,2))
	#OmegaVector = cross(vec(PosVector_EarthSol), vec(VelVector_EarthSol)) ./ dot(PosVector_EarthSol,PosVector_EarthSol)

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

	println("size(PosHat_EarthSol) = ", size(PosHat_EarthSol))
	println("size(first(PosHat_EarthSol)) = ", size(first(PosHat_EarthSol)))
	sini = map(i->sqrt(1 - dot(OmegaSolHat,PosHat_EarthSol[i])^2), 1:length(PosHat_EarthSol))
	Gamma = 1.04
	DeltaOmega = map(i->OmegaSolVector .- OmegaVector[i],1:length(OmegaVector))
	println("# size(DeltaOmega) = ", size(DeltaOmega))
	println("# size(OmegaSolVector) = ", size(OmegaSolVector))
	println("# size(sini) = ", size(sini))

	Delta = ((Gamma* R_sun)^2) .* (dot.(DeltaOmega,DeltaOmega).*sini^.2 .- dot.(OmegaSolVector,OmegaSolVector))

end
=#

"""
	`CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd; long, lat, alt)`
Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
	Inputs:
		jd: Julian Date (scalar)
		long: longitude (degrees, west negative)
		lat: latitude (degrees)
		alt: altitude (m)
"""
function CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd::Union{Real,Vector{<:Real}}; long::Real, lat::Real, alt::Real )
#function CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd::Real; long::Real, lat::Real, alt::Real )
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
function CalculateFWHMDifference_SolarRotation_from_obs(jd::Union{Real,Vector{<:Real}}; obs::Symbol )
#function CalculateFWHMDifference_SolarRotation_from_obs(jd::Real; obs::Symbol )
	loc = get_obs_loc(obs)
	CalculateFWHMDifference_SolarRotation_from_long_lat_alt(jd, long=loc["lon"], lat=loc["lat"], alt=loc["elevation"])
end

calc_Δfwhm_solar_rotation(jd::Union{Real,Vector{<:Real}}; obs::Symbol ) = CalculateFWHMDifference_SolarRotation_from_obs(jd, obs=obs)
export calc_Δfwhm_solar_rotation

end # module SolarRotation
