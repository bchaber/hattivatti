# run with $ julia --track-allocation profile.jl poiseuille.jl
using Profile
include(first(ARGS))

step()
Profile.clear_malloc_data()
step()
