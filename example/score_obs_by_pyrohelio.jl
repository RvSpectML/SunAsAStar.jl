### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 723b4419-4413-431a-9a2b-5a501cbd2607
begin
	import Pkg;
    Pkg.activate(joinpath(homedir(),"Code","RvSpectMLEcoSystem","NeidSolarScripts"))
	begin 
	using CSV, JLD2, FileIO
	using DataFrames, Query, DataValues
	using Dates
	using NeidSolarScripts
	using Polynomials
	using Statistics, NaNMath
	using PlutoUI
	using Plots, ColorSchemes
end
end

# ╔═╡ 67b3208b-428b-4f48-94e5-3cb11beee67d
begin
	path = "/mnt/data_simons/NEID/pyrohelio"
	files = readdir(path)
end

# ╔═╡ b2c1942e-7d8a-4b5e-afba-d42aff7f6471
begin
	df = DataFrame()
	for f in files
		if !occursin("neid_ljpyrohelio_chv0",f) continue end
		df_tmp = CSV.read(joinpath(path,f), header=[:time,:val1,:val2], DataFrame) |> @map({time=DateTime(_.time[1:23]),flux=_.val2}) |> DataFrame
		append!(df,df_tmp)
	end
	size(df)
end

# ╔═╡ b677c0ee-0c58-4a48-8d19-d1c54f7771bb
days_to_keep = 360

# ╔═╡ 90172c7e-3339-49da-bd3c-27a7fb0dd3f0
df_by_date = df |> @take(days_to_keep*24*60*60) |> @mutate( date=Date(_.time #=-Dates.Hour(13)=#) ) |> @groupby(_.date) |> @map({ date=_.date, time=_.time, flux=_.flux }) |> DataFrame;

# ╔═╡ 003ff73f-aad3-40bf-a673-7c28091e51e4
begin
	datestr = "20210312"
	manifest = CSV.read(joinpath(path,datestr,"manifest.csv"), DataFrame)
end

# ╔═╡ 8b662f93-b270-4ce4-ace0-ac3f3bb1dc20
first(df_by_date[1,:date])

# ╔═╡ 1276658b-97fc-4859-84e2-8df7cd43d1f6
m = match(r"neidL1_(\d{8}T\d{6})\.fits",manifest.Filename[3])

# ╔═╡ 1cffcbe4-fa66-4cbb-835b-4ee8e29b87b5
idx_date = findfirst(isequal(Date(m.captures[1][1:8],"YYYYmmdd")),map(i->first(df_by_date[i,:date]),1:size(df_by_date,1)))

# ╔═╡ e4dc87e4-f7d4-462a-a24f-d22ee036ba9b
files

# ╔═╡ 24ff107b-f337-4667-a09d-4e0342b085de
pwd()

# ╔═╡ 34f82d2f-4740-4422-9a17-e35702b51556
manifest[!,:solar_hour_angle] = SolarRotation.get_solar_hour_angle(manifest[!,:bjd], obs=:WIYN)
#manifest[!,:solar_hour_angle] = map(i->SolarRotation.get_solar_hour_angle(manifest[i,:bjd], obs=:WIYN),1:size(manifest,1))

# ╔═╡ 451d154c-29ad-4862-894a-d96dc97e8f12
names(manifest)

# ╔═╡ 7fde57bb-8bb9-499b-89f5-3f31fc5050ce
df_out = select(manifest, :Filename => (x->basename.(x))=>:filename, :bjd, :pyrheliometer_flux_over_model, :pyrheliometer_rms_flux, :alt_sun, :solar_hour_angle, :Δv_diff_ext, :Δfwhm²)

# ╔═╡ a013b123-05c8-4054-bc69-eedab7da28cc
CSV.write(datestr * "_solar_info.csv", df_out)

# ╔═╡ cb47c02c-cd8f-4a05-a3f4-d7baaf6afb8d
size(manifest)

# ╔═╡ b5afbe2c-b3c6-4825-9804-9df9d22c9a56
df_by_date[idx_date,:time]|>collect

# ╔═╡ a83c43ba-fc15-4756-bebb-53d29635db5c
 SolarRotation.get_solar_info(manifest.bjd[1], obs=:WIYN)

# ╔═╡ f5bfccb9-109b-48ae-ab4c-39b72ff070c8
SolarRotation.get_solar_hour_angle(manifest.bjd[1], obs=:WIYN)

# ╔═╡ fe8a34e8-afc6-11eb-1ff5-dfac1585fd04
#=
begin
	local plt3 = plot(legend=:none)
	local pal = palette(:berlin, length(df_plt))
	for day_idx in idx_clear_days[rms_list.<=0.015]
		if size(df_plt[day_idx],1) > 1
		local t = Time.(df_plt[day_idx][:,:times]).-Hour(7)
		local model = flux_model.(df_plt[day_idx][!,:airmass],df_plt[day_idx][!,:sol_dist])
		local mask = (0.95 .<= df_plt[day_idx][!,:flux]./model .<= 1.05)
		if sum(mask) >= 4
			local coeff_daily = fit(Dates.value.(t[mask].-t[1]),df_plt[day_idx][mask,:flux]./model[mask],2)
	
			model .*= coeff_daily.(Dates.value.(t.-t[1]))
			local norm = quantile(df_plt[day_idx][mask,:flux]./model[mask],0.75)
			model .*= norm
			local mask2 = (0.995 .<= df_plt[day_idx][:,:flux]./model .<= 1.005)
			#plot!(plt3, df_plt[day_idx][!,:airmass], model, color=pal[day_idx])
			scatter!(plt3,df_plt[day_idx][!,:airmass],df_plt[day_idx][mask,:flux]./(model[mask]) ,ms=1.0, label=:none) 
			scatter!(plt3,df_plt[day_idx][!,:airmass],df_plt[day_idx][mask2,:flux]./(model[mask2]) ,ms=2.5, label=:none) 
		end
		end
		ylims!(plt3,0.98,1.015)
	end
	plt3
end	
=#

# ╔═╡ d56c42dc-66b9-4730-99a6-fdb18debf798
function get_idx_in_pyro_data(d_start::DateTime, tvec::AVD; exp_time::Real = 55) where {AVD <: AbstractVector{DateTime} }
	idx_lo = searchsortedfirst(tvec,d_start)
	d_stop = julian2datetime(datetime2julian(d_start)+exp_time/(24*60*60))
	idx_hi = searchsortedlast(tvec,d_stop)
	if idx_lo >1 idx_lo -= 1 end
	if idx_hi <length(tvec) idx_hi += 1 end
	return (idx = idx_lo:idx_hi, time_start = tvec[idx_lo], time_stop=tvec[idx_hi])
end

# ╔═╡ fff4ec82-2dd3-4d56-9faf-4b309bfee108
get_idx_in_pyro_data(DateTime(m.captures[1],"YYYYmmddTHHMMSS"), df_by_date[idx_date,:time] )

# ╔═╡ d9f7d500-7f43-4217-9946-7c20647d3940
idx_sec = get_idx_in_pyro_data(julian2datetime(manifest[3,:bjd]), (df_by_date[idx_date,:time]) )

# ╔═╡ 3c9725b1-ea31-40ab-b61c-90e1f07c23ea
function get_alt_airmass(d::DateTime)
	t = datetime2julian(d)
	return get_alt_airmass(t)
end

# ╔═╡ 72c08181-2a6f-4475-abfa-e89ca5a4965e
function get_alt_airmass(t::Real)
	result = SolarRotation.get_solar_info(t, obs=:WIYN)
	return (alt = result[:alt], airmass = result[:airmass], sol_dist = result[:sol_dist_au])
end

# ╔═╡ 51c3761e-0a85-4b37-90b0-499538120f4c
df_obs = DataFrame(get_alt_airmass.(df_by_date[idx_date,:time][idx_sec]))

# ╔═╡ 86c3d69a-57d2-4828-bb75-f9b88f6036a5
function flux_model(airmass::Real, dist::Real = 1)
	poly = Polynomial([7.061980163705281, -0.07628424166629734])
	exp(poly(airmass))/ dist^2
end

# ╔═╡ cf5b8699-6957-4067-a4a1-f32780ef2644
function flux_model_obs(dt::DateTime; exp_time::Real = 55)
	d = Date(dt)
	idx_date = findfirst(isequal(d),map(i->first(df_by_date[i,:date]),1:size(df_by_date,1)))
	(idx_sec, tlo, thi) = get_idx_in_pyro_data(dt, dropna(df_by_date[idx_date,:time]), exp_time=exp_time )
	#return (1.0, length(idx_sec))
	#(alts, airmasses, soldists) = get_alt_airmass.(df_by_date[idx_date,:time][idx_sec])
	#fm = flux_model.(df_obs.airmass,df_obs.sol_dist)
	df_obs_start = get_alt_airmass(df_by_date[idx_date,:time][first(idx_sec)])
	df_obs_stop  = get_alt_airmass(df_by_date[idx_date,:time][last(idx_sec)])
	flux_start = flux_model(df_obs_start[:airmass],df_obs_start[:sol_dist])
	flux_stop  = flux_model(df_obs_stop[:airmass],df_obs_stop[:sol_dist])
	flux_mean = 0.5*(flux_start+flux_stop)
	#flux_slope = (flux_stop-flux_start)/Dates.value(Dates.Second(thi-tlo))
	#return df_by_date[idx_date,:flux][idx_sec] |> collect
	#return dropna(df_by_date[idx_date,:flux])[idx_sec], std(df_by_date[idx_date,:flux][idx_sec] |> collect)
	flux_ratio = dropna(df_by_date[idx_date,:flux]) 
	flux_ratio = flux_ratio[idx_sec]./flux_mean
	mean_flux_ratio = mean(flux_ratio)
	#flux_ratio = dropna(df_by_date[idx_date,:flux])[idx_sec] ./flux_mean
	#flux_ratio = mean(dropna(df_by_date[idx_date,:flux][idx_sec])) ./flux_mean
	#return flux_ratio #, std(flux_ratio) )
	#NaNMath.std(flux_ratio)
	return (ave = mean_flux_ratio, rms = std(flux_ratio) ) 
	#slope = flux_slope, idx_date=idx_date, idx_sec=idx_sec, flux_ratio = flux_ratio)
end

# ╔═╡ 5af436fe-bb4a-4069-bb82-42259d7d57c0
mean(flux_model.(df_obs.airmass,df_obs.sol_dist))

# ╔═╡ 23bd0a0b-2ef3-4acd-8ab6-d61f8cdd9c2b
df_by_date[188

# ╔═╡ 63de83e7-91c8-41aa-91de-51f280534c67
function flux_model_obs(fn::String)
	m = match(r"neidL1_(\d+T\d+)\.fits",basename(fn))
	dt = DateTime(m.captures[1],"yyyymmddTHHMMSS")
	flux_model_obs(dt)	
end

# ╔═╡ 8f484f4a-c194-474e-825c-4301aae07e92
function normalized_flux_obs(df::DataFrame)
	flux_ratios = zeros(size(df,1))
	std_flux_ratios = zeros(size(df,1))
	#idx_rngs = fill(0:1,size(df,1))
	for (i,row) in enumerate(eachrow(df))
		#flux_ratio = flux_model_obs.(df[!,:Filename])
		(mean_flux_ratio, std_flux_ratio) = flux_model_obs(row.Filename)
		flux_ratios[i] = mean_flux_ratio
		std_flux_ratios[i] = std_flux_ratio
	end
	
		local mask = (0.95 .<= flux_ratios .<= 1.05)
		if sum(mask) >= 4
			coeff_daily = fit((df[mask,:bjd].-manifest[1,:bjd]),dropna(flux_ratios[mask]),2)
			flux_ratios ./= coeff_daily.(df[:,:bjd].-df[1,:bjd])
			flux_ratios ./= quantile(dropna(flux_ratios[mask]),0.75)
		end
	return (flux_ratios, std_flux_ratios)
end

# ╔═╡ b671c156-4ea7-4d4b-93fa-69f632e3da18
for (i,f) in enumerate(files)
	m = match(r"(\d{8})$", f) 
	if m == nothing continue end
	datestr = m.captures[1]
	manifest_fn = joinpath(path,datestr,"manifest.csv")
	println("# Processing ",manifest_fn)
	if !isfile(manifest_fn)
		println("# Can't find ", manifest_fn)
		continue
	end
	output_fn = joinpath(path,datestr * "_solar_info.csv")
	if filesize(output_fn) > 0 continue end
	manifest = CSV.read(manifest_fn, DataFrame)
	println("# Read ", size(manifest,1), " rows.")
	if size(manifest,1) == 0 continue end
	(fm, fs) = normalized_flux_obs(manifest)
	manifest[!,:pyrheliometer_flux_over_model] = fm
	manifest[!,:pyrheliometer_rms_flux] = fs
	manifest[!,:solar_hour_angle] = SolarRotation.get_solar_hour_angle(manifest[!,:bjd], obs=:WIYN)
	df_out = select(manifest, :Filename => (x->basename.(x))=>:filename, :bjd, :pyrheliometer_flux_over_model, :pyrheliometer_rms_flux, :alt_sun, :solar_hour_angle, :Δv_diff_ext, :Δfwhm²)

	CSV.write(output_fn, df_out)
	#if i>3 break end
end


# ╔═╡ 89beb337-444f-4953-a6ed-a51c75d22cd9
(fm, fs) = normalized_flux_obs(manifest)

# ╔═╡ 0a6793ad-99e4-4d2a-9ec6-a82f2e1278da
scatter(manifest[!,:bjd],fm)

# ╔═╡ e34fbb33-05c0-48bf-b4b8-3b5d74d77d57
scatter(fm,fs,xlims=(0.9,1.02),ylims=(0.,0.1))

# ╔═╡ 4471d7dd-f26d-4a0d-9f37-de86d736579d
manifest[!,:pyrheliometer_flux_over_model] = fm

# ╔═╡ b96e3c2e-87e4-4311-9887-f5692d5e2259
manifest[!,:pyrheliometer_rms_flux] = fs

# ╔═╡ b4ba194d-8123-4ddf-b80b-f79470926d6a
size(fs)

# ╔═╡ 657c2fbc-5d68-4549-93ac-aabee96a4304
#DateTime(Date(now()),Time(Dates.hour(now()),Dates.minute(now()),Dates.second(now())+55))
#(Date(now()) + 
	Dates.Second(Time(0,30,46)-Time(0,23,46))
#Dates.minute(now())

# ╔═╡ c9e6b8fe-df4d-4c56-85f2-1e8d9b38a70b
df_by_date[idx_date,:flux][idx_sec]./1000

# ╔═╡ Cell order:
# ╠═723b4419-4413-431a-9a2b-5a501cbd2607
# ╠═67b3208b-428b-4f48-94e5-3cb11beee67d
# ╠═b2c1942e-7d8a-4b5e-afba-d42aff7f6471
# ╠═b677c0ee-0c58-4a48-8d19-d1c54f7771bb
# ╠═90172c7e-3339-49da-bd3c-27a7fb0dd3f0
# ╠═003ff73f-aad3-40bf-a673-7c28091e51e4
# ╠═8b662f93-b270-4ce4-ace0-ac3f3bb1dc20
# ╠═1cffcbe4-fa66-4cbb-835b-4ee8e29b87b5
# ╠═1276658b-97fc-4859-84e2-8df7cd43d1f6
# ╠═fff4ec82-2dd3-4d56-9faf-4b309bfee108
# ╠═d9f7d500-7f43-4217-9946-7c20647d3940
# ╠═0a6793ad-99e4-4d2a-9ec6-a82f2e1278da
# ╠═e34fbb33-05c0-48bf-b4b8-3b5d74d77d57
# ╠═b671c156-4ea7-4d4b-93fa-69f632e3da18
# ╠═e4dc87e4-f7d4-462a-a24f-d22ee036ba9b
# ╠═24ff107b-f337-4667-a09d-4e0342b085de
# ╠═89beb337-444f-4953-a6ed-a51c75d22cd9
# ╠═4471d7dd-f26d-4a0d-9f37-de86d736579d
# ╠═b96e3c2e-87e4-4311-9887-f5692d5e2259
# ╠═34f82d2f-4740-4422-9a17-e35702b51556
# ╠═451d154c-29ad-4862-894a-d96dc97e8f12
# ╠═7fde57bb-8bb9-499b-89f5-3f31fc5050ce
# ╠═a013b123-05c8-4054-bc69-eedab7da28cc
# ╠═b4ba194d-8123-4ddf-b80b-f79470926d6a
# ╠═cb47c02c-cd8f-4a05-a3f4-d7baaf6afb8d
# ╠═8f484f4a-c194-474e-825c-4301aae07e92
# ╠═cf5b8699-6957-4067-a4a1-f32780ef2644
# ╠═51c3761e-0a85-4b37-90b0-499538120f4c
# ╠═5af436fe-bb4a-4069-bb82-42259d7d57c0
# ╠═b5afbe2c-b3c6-4825-9804-9df9d22c9a56
# ╠═a83c43ba-fc15-4756-bebb-53d29635db5c
# ╠═f5bfccb9-109b-48ae-ab4c-39b72ff070c8
# ╠═fe8a34e8-afc6-11eb-1ff5-dfac1585fd04
# ╠═d56c42dc-66b9-4730-99a6-fdb18debf798
# ╠═3c9725b1-ea31-40ab-b61c-90e1f07c23ea
# ╠═72c08181-2a6f-4475-abfa-e89ca5a4965e
# ╠═86c3d69a-57d2-4828-bb75-f9b88f6036a5
# ╠═23bd0a0b-2ef3-4acd-8ab6-d61f8cdd9c2b
# ╠═63de83e7-91c8-41aa-91de-51f280534c67
# ╠═657c2fbc-5d68-4549-93ac-aabee96a4304
# ╠═c9e6b8fe-df4d-4c56-85f2-1e8d9b38a70b
