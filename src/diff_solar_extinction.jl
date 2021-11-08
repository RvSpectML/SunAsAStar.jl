__precompile__() # this module is safe to precompile
"""
   DifferentialExtinction
   Original Python code (solar_diffex.py) by Andrea Lin
   Ported to Julia by Eric Ford
   Based on Collier-Cameron et al. (2019) MNRAS 487, 1082.  doi:10.1093/mnras/stz1215
   See also Astronomical Algorithms (Meeus, 1991) p151 & p178
"""
module DifferentialExtinction
using ..SunAsAStar

using Dates
using DataFrames, CSV


import ..SunAsAStar: AU, R_sun, valid_obs

f_z(alt) = π/2-alt

function f_airmass(z)
	secz = sec(z)
	return secz * (1.0 - 0.0012*(secz^2 - 1.0))
end

default_k_extinction = 0.157482  # Fit for WIYN, revisit at some point?
"""fractional brightness diff (from sun center)"""
function f_epsilon(z_center, dist_sun; k = default_k_extinction ) #
        z_upper = z_center - asin(R_sun/dist_sun)
        z_lower = z_center + asin(R_sun/dist_sun)
        delta_airmass = f_airmass(z_upper) - f_airmass(z_lower)
        #=
		 k_ex = 0.157482  # WARNING: Hardwired for WIYN.
		 X = [ones(length(eph.columns["airmass"])) eph.columns["airmass"]]
		 y =  eph.columns["magextinct"]
		 coeff = X \y
		 k_ex = coeff[2]
		 =#
		return k * delta_airmass/2.0 * log(10.)/2.5
end

"parallactic angle"
function f_q(H_sun, z_center, lat)
	#lat = deg2rad(wiyn["lat"]) #* u.deg # observer's latitude # no longer hard-wired
	return asin(sin(H_sun) * cos(lat)/sin(z_center))
end

f_phi(pole_angle, q) = (pole_angle - q)

function f_pole_angle_ecl(alpha_sun)
	e_tilt = deg2rad(23.4) # ecliptic tilt
    asin(sin(e_tilt) * cos(alpha_sun))
end

f_psi(pole_angle_ecl, q) = (pole_angle_ecl - q)

 """
 f_B0(JD)
  see Astronomical Algorithms (Meeus, 1991) p151 & p178
 """
 function f_B0(JD)

    T = (JD - 2451545.) / 36525.
    L0 = deg2rad(280.46645 + 36000.76983*T + 0.0003032*T^2)
    MA = deg2rad(357.52910 + 35999.05030*T - 0.0001559*T^2 - 0.00000048*T^3)
    sun_center = (1.914600 - 0.004817*T - 0.000014*T^2) * sin(MA) +
    			 (0.019993 - 0.000101*T) * sin(2*MA) +
    			 0.000290 * sin(3*MA)
	sun_center = deg2rad(sun_center)
    sun_true_lon = L0 + sun_center
    lam = sun_true_lon - deg2rad(0.00569) - deg2rad(0.00478 * sin(deg2rad(125.04 - 1934.136*T)))
    kappa = deg2rad(73.6667 + 1.3958333 * (JD - 2396758.) / 36525.)
    return asin(sin(lam - kappa) * sin(deg2rad(7.25)))
end

function diff_extinct_corr(JD, alt, HA, RA, pole_angle, dist_sun, lat)
    z = f_z(deg2rad(alt))
    epsilon = f_epsilon(z, dist_sun)
    q = f_q(deg2rad(HA*15), z, deg2rad(lat))

    phi = (f_phi(deg2rad(pole_angle), q))
    psi = f_psi(f_pole_angle_ecl(deg2rad(RA)), q)

    inc_sol_axis = f_B0(JD)

	ld_coeff = 0.6 # limb darkening coeff (linear)

	# sidereal diff rotation (Snodgrass & Ulrich 1990)
	sr_A = 2.972e-6 #* u.rad / u.s
	sr_B = -0.484e-6 #* u.rad / u.s
	sr_C = -0.361e-6 #* u.rad / u.s

	omega_orb = (2*π) / (365.25 * 24. * 3600.) #* u.rad / u.s # omega of Earth around Sun

    diff_rot_fac = sr_A*(7*ld_coeff - 15.) / (20*(3-ld_coeff)) +
      sr_B*(19 *ld_coeff - 35 + cos(inc_sol_axis)^2 * (3*ld_coeff - 35)) / (280*(3-ld_coeff)) +
      sr_C*(187*ld_coeff - 315 - cos(inc_sol_axis)^2 * (630 - 118*ld_coeff) -
      105 *cos(inc_sol_axis)^4 * (ld_coeff-1.)) / (6720*(3-ld_coeff))

    delta_vr = (epsilon * sin(phi) * R_sun * sin(inc_sol_axis) * diff_rot_fac -
       (7*ld_coeff - 15) * epsilon * sin(psi) / (20*(3-ld_coeff)) * R_sun * omega_orb)
	delta_vr *= 1e3 # km/s -> m/s.  WARNING: Andrea's python code returns cm/s

    return delta_vr
end

"""
   calc_Δv_diff_extinction(bjd; obs, ...)
   calc_Δv_diff_extinction(start, stop; obs, dt, dt_unit, ...)
Estimate apparent Δv (m/s) in solar observations due to differential extinction.
"""
function calc_Δv_diff_extinction end

function calc_Δv_diff_extinction(epoch::Float64; obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	calc_Δv_diff_extinction_astropy(epoch; obs=obs, max_airmass=max_airmass, verbose=verbose)
end

function calc_Δv_diff_extinction_horizons(;epochs::Union{Vector{Float64},Dict}, obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	loc = get_obs_loc(obs)
	eph = Horizons.query_horizons(epochs=epochs,obs=obs,max_airmass=max_airmass, verbose=verbose)
	delta_vr = diff_extinct_corr.(eph.columns["datetime_jd"], eph.columns["EL"],
	 	eph.columns["hour_angle"], eph.columns["RA"],
		eph.columns["NPole_ang"], eph.columns["delta"].*AU, loc["lat"])
	return (;jd=eph.columns["datetime_jd"], Δv=delta_vr)
end

function calc_Δv_diff_extinction_astropy(epoch::Float64; obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	@assert 1 < max_airmass <= 10.0
	sol = get_solar_info(epoch, obs=obs)
	ha = sol[:lst] .- sol[:ra]./15
	idx_visible = 0 < sol[:airmass] <= max_airmass
	#if !idx_visible return (jd=epoch, Δv = 0.0) end
	if !idx_visible return Δv = 0.0 end
	horizons = Horizons.get_horizons_output(epoch,obs=obs)
	#delta = SolarRotation.get_earth_sun_dist(epoch, obs=obs)
	loc = get_obs_loc(obs)
	delta_vr = diff_extinct_corr(epoch, sol[:alt], ha, sol[:ra], horizons[:sol_npole_ang], sol[:sol_dist_au].*AU, loc["lat"]) #eph["NPole_ang"], eph["delta"])
	return delta_vr
end

function calc_Δv_diff_extinction_astropy(epochs::Vector{Float64}, obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	@assert 1 <= length(epochs) <= 10000
	@assert 1 < max_airmass <= 10.0
	sol = get_solar_info(epochs, obs=obs)
	ha = sol.lst .- sol.ra./15
	idx_visible = 0 .< sol.airmass .<= max_airmass
	horizons = Horizons.get_horizons_output.(epochs[idx_visible],obs=obs)
	delta = get_earth_sun_dist(epochs[idx_visible], obs=obs)  # CHECK:  OK TO REMOVE DOT?
	#delta = sun.sol_dist_au
	loc = get_obs_loc(obs)
	#println(" jd = ", epochs, " delta_h = ", horizons.sol_dist_au, " delta_ap = ", delta)
	#println(" diff = ", horizons.sol_dist_au.-delta)
	delta_vr = diff_extinct_corr.(epochs[idx_visible], sol.alt[idx_visible], ha[idx_visible], sol.ra[idx_visible], horizons.sol_npole_ang, horizons.sol_dist_au.*AU, loc["lat"]) #eph["NPole_ang"], eph["delta"])
	return delta_vr
end

function compare_Horizons_Astropy(epochs::Vector{Float64}, obs::Symbol; max_airmass::Float64 = 10.0, verbose::Bool = false)
	if verbose println("# compare_Horizons_Astropy: ") end
	loc = get_obs_loc(obs)
	eph = Horizons.query_horizons(epochs=epochs,obs=obs,max_airmass=max_airmass, verbose=verbose)
	delta_vr = diff_extinct_corr.(eph.columns["datetime_jd"], eph.columns["EL"],
	 	eph.columns["hour_angle"], eph.columns["RA"],
		eph.columns["NPole_ang"], eph.columns["delta"].*AU, loc["lat"])

	sol = get_solar_info(epochs, obs=obs)
	ha = sol.lst .- sol.ra./15
	idx_visible = 0 .< sol.airmass .<= max_airmass
	delta = get_earth_sun_dist.(epochs[idx_visible], obs=obs)

	if verbose
		println("Horizons airmass =", eph.columns["airmass"])
		println("Astropy  airmass =", sol.airmass[idx_visible])
		#println("Δjd = ", eph.columns["datetime_jd"], " vs ", epochs[idx_visible])
		println("Δjd = ", eph.columns["datetime_jd"].-epochs[idx_visible])
		println("Δalt = ", eph.columns["EL"].- sol.alt[idx_visible])
		println("Δha = ", eph.columns["hour_angle"].- ha[idx_visible])
		println("ΔRA = ", eph.columns["RA"].- sol.ra[idx_visible])
	end
	horizons = Horizons.get_horizons_output.(epochs[idx_visible],obs=obs)
	delta_vr_new = diff_extinct_corr.(epochs[idx_visible], sol.alt[idx_visible], ha[idx_visible], sol.ra[idx_visible], zeros(length(epochs[idx_visible])), fill(AU,length(epochs[idx_visible])), loc["lat"]) #eph["NPole_ang"], eph["delta"])

	if verbose
		println("delta_vr_old = ", delta_vr, "\ndelta_vr_new = ", delta_vr_new)
		println("delta vr = ", delta_vr.- delta_vr_new)
	end
	return (eph=eph, sol=sol, delta_vr_old=delta_vr, delta_vr_new=delta_vr_new)
end

#=
function calc_Δv_diff_extinction(bjd::Real; obs::Symbol, max_airmass::Float64 = 10.0, verbose::Bool = false)
	epochs = [bjd]
	(bjd_out, Δv) = calc_Δv_diff_extinction(epochs=epochs, obs=obs, max_airmass=max_airmass, verbose=verbose)
	return first(Δv)
end
=#

#=
function calc_Δv_diff_extinction(start_date::String, stop_date::String; obs::Symbol, dt::Integer = 1, dt_unit::String = "m", max_airmass::Float64 = 10.0, verbose::Bool = false)
	@assert 1 <= dt <= 24*60
	@assert dt_unit ∈ ["m", "h", "d"]
	@assert 1 < max_airmass <= 10.0
	epochs = Dict("start"=>start_date, "stop"=>stop_date, "step"=>string(dt) * dt_unit)
	#(bjd_out, Δv) = calc_Δv_diff_extinction_horizons(epochs=epochs, obs=obs, max_airmass=max_airmass, verbose=verbose)
	# Doesn't work yet
	(bjd_out, Δv) = calc_Δv_diff_extinction_astropy(epochs, obs=obs, max_airmass=max_airmass, verbose=verbose)

	return (bjd=bjd_out, Δv=Δv)
end
=#

#=
function calc_Δv_diff_extinction(start_bjd::Real, stop_bjd::Real; obs::Symbol, dt::Integer = 1, dt_unit::String = "m", max_airmass::Float64 = 10.0, verbose::Bool=false)
	start_date = string(Date(Dates.julian2datetime(start_bjd)))
	stop_date = string(Date(Dates.julian2datetime(stop_bjd)))
	calc_Δv_diff_extinction(start_date, stop_date, obs=obs, dt=dt, dt_unit=dt_unit, max_airmass=max_airmass, verbose=verbose)
end

function calc_Δv_diff_extinction_demo(; verbose::Bool = false)
	start_date = "2021-01-01"
	stop_date = "2021-02-11"
	calc_Δv_diff_extinction(start_date, stop_date, obs=:WIYN, verbose=verbose)
end
=#

export calc_Δv_diff_extinction

end # module DifferentialExtinction
