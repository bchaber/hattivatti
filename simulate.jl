using ProgressMeter
include(first(ARGS))

init()
@showprogress for t=1:NT
    step()
end
evaluate()