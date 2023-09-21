#
# Copyright [2023] [Simon Bernier]
#
using ITensorInfiniteMPS
using ITensors

include(
  joinpath(
    pkgdir(ITensorInfiniteMPS), "examples", "vumps", "src", "vumps_subspace_expansion.jl"
  ),
)

##############################################################################
# VUMPS parameters
#

maxdim = 50 # Maximum bond dimension
cutoff = 1e-6 # Singular value cutoff when increasing the bond dimension
max_vumps_iters = 200 # Maximum number of iterations of the VUMPS algorithm at each bond dimension
vumps_tol = 1e-5
outer_iters = 5 # Number of times to increase the bond dimension
localham_type = MPO # or ITensor
conserve_qns = true
eager = true

##############################################################################
# CODE BELOW HERE DOES NOT NEED TO BE MODIFIED
#

N = 3 # Number of sites in the unit cell

initstate(n) = isodd(n) ? "↑" : "↓"
s = infsiteinds("Electron", N; initstate, conserve_qns)
ψ = InfMPS(s, initstate)

function ITensorInfiniteMPS.unit_cell_terms(::Model"extendedHopping"; t, g=0.0)
    opsum = OpSum()
    # nearest neighbour hopping for a unit cell
    opsum += -t, "Cdagup", 1, "Cup", 2
    opsum += -t, "Cdagup", 2, "Cup", 1
    opsum += -t, "Cdagdn", 1, "Cdn", 2
    opsum += -t, "Cdagdn", 2, "Cdn", 1
    opsum += -t, "Cdagup", 2, "Cup", 3
    opsum += -t, "Cdagup", 3, "Cup", 2
    opsum += -t, "Cdagdn", 2, "Cdn", 3
    opsum += -t, "Cdagdn", 3, "Cdn", 2
    # next-nearest neighbour hopping
    opsum += -t, "Cdagup", 1, "Cup", 3
    opsum += -t, "Cdagup", 3, "Cup", 1
    opsum += -t, "Cdagdn", 1, "Cdn", 3
    opsum += -t, "Cdagdn", 3, "Cdn", 1
    # three site g-interaction
    opsum += t*g, "Cdagup", 1, "Ntot", 2, "Cup", 3
    opsum += t*g, "Cdagup", 3, "Ntot", 2, "Cup", 1
    opsum += t*g, "Cdagdn", 1, "Ntot", 2, "Cdn", 3
    opsum += t*g, "Cdagdn", 3, "Ntot", 2, "Cdn", 1
    return opsum
end
model = Model("extendedHopping")

# Form the Hamiltonian
H = InfiniteSum{MPO}(model, s; t=1.0, g=0.0)

# Check translational invariance
println("\nCheck translation invariance of the initial VUMPS state")
@show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

outputlevel = 1
vumps_kwargs = (tol=vumps_tol, maxiter=max_vumps_iters, outputlevel, eager)
subspace_expansion_kwargs = (cutoff=cutoff, maxdim=maxdim)

ψ = vumps_subspace_expansion(H, ψ; outer_iters, subspace_expansion_kwargs, vumps_kwargs)

# Check translational invariance
println()
println("==============================================================")
println()
println("\nCheck translation invariance of the final VUMPS state")
@show norm(contract(ψ.AL[1:N]..., ψ.C[N]) - contract(ψ.C[0], ψ.AR[1:N]...))

Sz = [expect(ψ, "Sz", n) for n in 1:N]

energy_infinite = expect(ψ, H)
@show energy_infinite

## using JLD2
## jldsave("infmps.jld2"; ψ)