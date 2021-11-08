ENV["PYTHON"] = ""
using PyCall
using Conda

pipcmd = joinpath(Conda.PYTHONDIR,"pip")

#Conda.add("astroquery.jplhorizons")
Conda.add("astroquery", channel="conda-forge")
#run(`$pipcmd install --pre astroquery`)
#Conda.add("astroquery.jplhorizons", channel="conda-forge")

run(`$pipcmd install barycorrpy`)

Conda.add("astropy")
#Conda.add("astropy.coordinates")
#Conda.add("astropy.time")

Conda.add("scipy")
#Conda.add("scipy.spatial.transform")

#run(`$pipcmd install pyneid`)
