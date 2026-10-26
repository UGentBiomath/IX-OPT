### A Pluto.jl notebook ###
# v0.20.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ eee0fd64-96b6-4753-b262-9c9c4bfd81a1
using PlutoUI; TableOfContents(depth=4, include_definitions=false)

# ╔═╡ 051f3f5c-a493-4e22-b77c-45222096cdd1
using PlutoUI: combine

# ╔═╡ 9b1b7dbe-5fbc-4426-9d91-664a67938d34
using ModelingToolkit, MethodOfLines, DomainSets, OrdinaryDiffEq

# ╔═╡ 794e1837-5b56-42a1-af76-c9b566523ec4
using Plots, Plots.PlotMeasures, LaTeXStrings, StatsPlots

# ╔═╡ c21a3c39-d908-4444-9ebe-63f4569a0b31
using NumericalIntegration, Statistics, LinearAlgebra

# ╔═╡ 2862b00a-a238-4448-a7e5-e5f121e4b0d8
using Optim

# ╔═╡ 14ddcaa9-2d26-40ea-a1cd-a8f1c144d272
using FiniteDiff

# ╔═╡ 63850d62-693a-4651-851e-3e4665a950d1
using Interpolations

# ╔═╡ 3844c32a-1212-421e-a25c-82a859aebb10
using Distributions

# ╔═╡ 00dd236d-444d-4227-9068-8d1cdb26f0dd
using Printf

# ╔═╡ d833b79d-ca76-43a9-b3fc-a3e5e6a989dd
using PlutoTeachingTools

# ╔═╡ 5519a8d7-4abe-449d-83f4-d50ecce71494
using Random

# ╔═╡ 9da02a78-6118-47f6-aa32-0157e2cad312
using GlobalSensitivity, QuasiMonteCarlo

# ╔═╡ 55e3366b-e1e8-44b2-ad02-44146ce0d839
using Turing

# ╔═╡ b8d2511d-0865-4d5a-999e-bf40b2314dd8
md"""
# Libraries
"""

# ╔═╡ 34266a3a-004e-495c-8eac-3df81d7f2232
cp = theme_palette(:default)

# ╔═╡ 581c7094-6e13-4118-a3a7-6d4f9a2abff4
md"""
# IX multicomponent model
"""

# ╔═╡ 8b071f14-cfcf-4ba6-a65d-85962b6dedc2
md"""
## Parameters
"""

# ╔═╡ 57dec26e-3e14-45ec-83b3-34567fc16414
md"""
##### Components
"""

# ╔═╡ 3354bf0d-ff7e-4f47-9b07-c9dd7995ac52
c = comp = (:H, :Na, :K, :Mg, :Ca)#, :Ba)

# ╔═╡ 9a763df8-dfd8-42b3-9519-e2401f08bd7e
N::Int = length(comp)

# ╔═╡ 00503082-7ea5-486a-b121-5e130c2fde4e
lab = ["H", "Na", "K", "Mg", "Ca"]

# ╔═╡ 791bdb1e-cea6-4f94-8bce-eb50fc38ddfb
col = vcat(:grey, cp[1:N-1])

# ╔═╡ 78e93282-dcd5-4517-98f7-86059c31a551
md"""
##### Independent variables
"""

# ╔═╡ aa4be22a-5dc8-4592-b4fa-b8fbe53447d7
tmax = 2000.

# ╔═╡ cfed0866-324d-4f74-9aae-83966de26969
Lmax = 1.

# ╔═╡ 0668d48e-bad2-4d1d-88a5-d0c80e6bfc06
@parameters t x

# ╔═╡ 565221c8-15a9-4482-bcf4-885009acbe41
domains = [t ∈ Interval(0.0, tmax), x ∈ Interval(0.0, Lmax)]

# ╔═╡ 9ce4958a-914a-412a-a283-2f93d1912b9f
md"""
##### Dependent variables
"""

# ╔═╡ 79f86646-ecf8-47cd-972c-08b1dccd4268
@variables C(..)[1:N] [description = "Liquid concentration C(t,x)"]

# ╔═╡ a03fc089-0c67-4e52-a283-b8553c074f53
@variables q(..)[1:N] [description = "Solid concentration q(t,x)"]

# ╔═╡ 5faf1dfa-4b4c-44db-8159-b8629a1346d6
dvs = [C(t,x), q(t,x)]

# ╔═╡ ad35a9f3-f94b-4bb4-aae9-56f67a50f51a
C |> typeof

# ╔═╡ 19d1c246-925b-4373-90ab-2f415416a5b6
C.metadata

# ╔═╡ ed5ff484-b0c2-42dc-b91e-ae5cc19cca38
md"""
##### Parameters
"""

# ╔═╡ 75b26c32-0c4f-4583-8635-acdcfb8ddcfc
# Column dimensions
@parameters L D A V

# ╔═╡ 5ef51740-4097-47be-b0ed-559bbeceac42
column = Dict([
		L => 0.6,   		# 600 	mm = 0.6 m – height of bed in the column
		D => 0.0264,  		# 26.4 	mm = 0.0264 m
		A => π/4 * D^2,
		V => A*L, 			# 383 	mL
	])

# ╔═╡ aa0d5203-09c4-40fc-9792-20a56ebe1d4a
# Operational
@parameters Q ε u τ

# ╔═╡ 81d0e150-ae98-4998-bc1f-c58e2c58e40b
operational = Dict([
		Q => 1.2e-4/60, 	# [m³/s] 120 mL/min = 2 mL/s = 2e-6 m3/s
		ε => 0.37, 			# [-] bed porosity (L effluent/L bed)
		u => Q/(A*ε), 		# [m/s] fluid effective velocity
		τ => L/u, 			# [s] column residence time (empty bed contact time)
	])

# ╔═╡ 148119ed-49c8-4106-9c41-60788d3ff7d6
# Fluid
@parameters ρf μ

# ╔═╡ 215cf20b-c2b5-4e80-9584-9b9510d1df7a
fluid = Dict([
		ρf => 1e3,
		μ  => 1e-3,
	])

# ╔═╡ 03178e82-d523-486c-ac5e-49ed20871dca
# Resin
@parameters ρ rₚ qmax

# ╔═╡ 1966ab11-75b6-4855-906e-16af2cd77f8d
resin = Dict([
		ρ  	 => 800., 		# [g/L] resin density
		rₚ 	 => 6.3e-4/2, 	# [m] resin particle radius
		qmax => 3.0, 		# [eq/L resin] maximum resin capacity, ~1.8 eq/L bed
	])

# ╔═╡ 21954982-4d00-410a-a1e4-ba247a90d744
# Initial conditions
@parameters (C0::Any)(..)[1:N] [description = "C0 – initial conditions"]

# ╔═╡ 6890ad0b-a50a-48de-95eb-d294bfffd6de
@parameters (q0::Any)(..)[1:N]

# ╔═╡ 1afa3007-1c4d-420d-a564-97954b703f66
# Boundary conditions
@parameters Cx0[1:N] [description = "Cx0 – Dirichlet boundary conditions"]

# ╔═╡ eb29fafd-b527-4ec1-a76e-292df208e5aa
@parameters qx0[1:N] [description = "qx0 – Dirichlet boundary conditions"]

# ╔═╡ bbe1895a-cbce-4a0a-a7e6-147b2361b499
# Dimensionless
@parameters Cref qref 		# [eq/L] reference concentrations

# ╔═╡ fa5e1464-735d-4999-9893-0e50f4360f6f
@parameters Re Dg 			# [-] dimensionless numbers

# ╔═╡ 9accd87d-68dc-4b73-8eda-193255022ada
@parameters kp ki

# ╔═╡ 7ee93e74-65c5-48c9-b456-6186bc81c71c
dimensionless = Dict([
		Cref => Cx0[2], 		# [eq/L] max. concentration of ion in the Sliquid
		# Cref => sum([Cx0[i] for i in 1:N]),
		qref => qmax,			# [eq/L] max. capacity of counterion in the resin
		Re => 2rₚ*u*ρf/μ, 		# [-] Reynolds Number
		Dg => (1-ε)/ε*qref/Cref,# [-] partition ratio
		# Dg => (1-ε)/ε*(1/Cref)
	])

# ╔═╡ 3e07fe4c-131e-4d34-bbdd-836972049a16
# Electric charge
@parameters Z[1:N] [description = "Zi – Electric charge"]

# ╔═╡ 8dc31ff5-525e-40ac-8e50-2906c89eefff
# Molar mass
@parameters M[1:N] [description = "Mi – Molar mass"]

# ╔═╡ 481221d6-960d-4f37-920c-7c8031e0c6b8
# Diffusivities in water
@parameters Dz[1:N] [description = "Dz - Diffusivity coefficient"]

# ╔═╡ f3fb5bb6-4500-4174-9c4a-636817625ce3
# Peclet number
@parameters Pe[1:N] [description = "Pe - Peclet number"]

# ╔═╡ 9da4312a-00eb-400e-ba3d-39f03a25a9e2
# Schmidt number
@parameters Sc[1:N] [description = "Sc - Schmidt number"]

# ╔═╡ 45679967-8e1b-4f3c-af73-41939bf99956
# Sherwood number
@parameters Sh[1:N] [description = "Sh - Sherwood number"]

# ╔═╡ 0ee7fee1-5a45-45b1-ba8d-944a31414ad5
# Diffusivities in resin
@parameters Dp[1:N] [description = "Dp - Diffusivity in resin"]

# ╔═╡ d3c12be8-bfd3-4312-adbc-50dad911b61b
# Separation factors
@parameters R[1:N] [description = "Ri - Separation factor"]

# ╔═╡ cffe894d-d84e-4673-8c94-6d739410841a
# Selectivity coefficients
@parameters K[1:N] [description = "Ki – Selectivity coefficient"]

# ╔═╡ ddd71bc0-a30a-4114-8d5b-5400c796a343
# Diffusivities in water [m²/s]
Di = Dict([
	:H  => 9.3e-9,
	:Na => 1.33e-9,
	:K  => 1.96e-9,
	:Mg => 0.705e-9,
	:Ca => 0.793e-9,
	:Ba => 0.85e-9,
])

# ╔═╡ d9878823-74d3-4a41-a043-89c9f772d4d6
# Selectivity coefficients
Ki = Dict([
	:H  => 1.0, 
    :Na => 1.56, 
    :K  => 2.28, 
    :Mg => 2.94, 
    :Ca => 4.06,
	:Ba => 9.06,
])

# ╔═╡ cd46bc3f-334f-4dc5-859a-b581de1ebdcc
# Electric charge
Zi = Dict([
	:H  => 1, 
    :Na => 1, 
    :K  => 1, 
    :Mg => 2, 
    :Ca => 2,
	:Ba => 2,
])

# ╔═╡ 66f76310-6861-45c7-8fb0-3bc5879d9f45
# Resistance coefficients
Ri = Dict([
	:H  => 10., 
    :Na => 10., 
    :K  => 10., 
    :Mg => 100., 
    :Ca => 100.,
	:Ba => 100.,
])

# ╔═╡ 9e8d5262-e0e4-40ef-a8be-5c38cc650747
# Molar mass [mg/mmol]
Mi = Dict([
	:H  =>   1.0079, 
    :Na => 	22.9897, 
    :K  => 	39.0983, 
    :Mg => 	24.305, 
    :Ca => 	40.078,
	:Ba => 137.33,
])

# ╔═╡ 71ea145c-23fa-409f-a235-354110004312
md"""
Minimum equilibrium concentration considered – 1e-5 ~
"""

# ╔═╡ 6a455ed0-8147-49a2-a729-93a1749e0f47
zero_eq = 1e-5

# ╔═╡ 215e8b8d-b113-4349-8522-8ebf395e54cd
# Initial conditions in the liquid phase [eq/L]
C₀ = Dict([
	:H  => 1e-7, 
    :Na => zero_eq, 
    :K  => zero_eq, 
    :Mg => zero_eq, 
    :Ca => zero_eq, 
	:Ba => 1e-9, 
])

# ╔═╡ e9d33fe7-3bb0-4ac7-909c-a4aef8e8e72f
# Initial conditions in the solid phase [eq/L]
q₀ = Dict([c[i] => 0. for i in 1:N]); q₀[:H] = 1.0; q₀

# ╔═╡ 65caa2ef-2e28-49f8-8e02-55a144abca8b
# Boundary conditions in the liquid phase – left side - mass [mg/L]
Cₓ₀mass = Dict{Symbol, Float64}([
	:H  => 1e-7, 
    :Na => 100.0, 
    :K  => 100.0, 
    :Mg => 100.0,
    :Ca => 100.0,
	:Ba => 1e-3,
])

# ╔═╡ 2e6b1b04-294a-4395-8e69-0348f4442a5d
Cₓ₀=Dict([i=>Cₓ₀mass[i]/(Mi[i]*1e3)*Zi[i] for i in c[2:5]]); Cₓ₀[:H] = 1e-7 ; Cₓ₀[:Ba] = 1e-7; Cₓ₀

# ╔═╡ 0c894ba4-0753-4327-854a-cf4093a9020e
# Boundary conditions in the solid phase – left side [eq/L]
qₓ₀ = Dict([
	:H  => 0., 
    :Na => 0., 
    :K  => 0., 
    :Mg => 0., 
    :Ca => 0., 
	:Ba => 0., 
])

# ╔═╡ b7e24c21-0d85-4c24-8000-4ba700cbd87e
# Boundary conditions in the liquid phase – right side [eq/L]
Cₓ₁ = Dict([
	:H  => 0., 
    :Na => 0., 
    :K  => 0., 
    :Mg => 0., 
    :Ca => 0., 
	:Ba => 0., 
])

# ╔═╡ b4b71052-9609-4fb1-a541-e9e21919f275
# Boundary conditions in the solid phase – right side [eq/L]
qₓ₁ = Dict([
	:H  => 0., 
    :Na => 0., 
    :K  => 0., 
    :Mg => 0., 
    :Ca => 0., 
	:Ba => 0., 
])

# ╔═╡ 886f1479-c180-4f0f-bd22-1c1632382c27
species = Dict([
		Z  => [Zi[i] for i in c],
		M  => [Mi[i] for i in c],
		Dz => [Di[i] for i in c],
		Pe => [u*L/Di[i] for i in c],
		Sc => [μ/ρf/Di[i] for i in c],
		Sh => [(1.09/ε)*(Re*Sc[i])^0.33 for i in 1:N],
		Dp => [Di[c[i]]/R[i] for i in 1:N],
		R => [Ri[c[i]] for i in 1:N],
		Cx0 => [Cₓ₀[i] for i in c],
		qx0 => [qₓ₀[i] for i in c], 
		K  => [Ki[i] for i in c],
	])

# ╔═╡ b33073d7-7676-4184-a426-0cb19bd75c7d
IC = Dict([C0 => [C₀[i] for i in c], q0 => [q₀[i] for i in c]])

# ╔═╡ 3b2bd4d5-1733-472e-ab2e-bf756488209a
# Merge all parameters
pars = merge(column, operational, fluid, resin, dimensionless, species, IC)

# ╔═╡ 398f3af2-fe1f-4f87-976b-5e9b2f041d58
pars[C0] |> typeof

# ╔═╡ c94b8d18-d346-412a-ab1d-5cc7768b9d88
# parameter set
ps = [L, D, A, V, Q, ε, u, τ, ρf, μ, ρ, rₚ, qmax, Cx0, qx0, Cref, qref, Re, Sc, Sh, Dg, K, Z, M, Dz, Pe, Dp, R, C0, q0]

# ╔═╡ 8f918e3e-edf7-4e77-80db-0248a18b4983
md"""
## Boundary conditions
"""

# ╔═╡ 344e091a-fff8-47c0-8bfc-85c33c8bfcb5
begin
	# Differentials
		Dt = Differential(t)
		Dx = Differential(x)
		Dxx = Differential(x)^2
end;

# ╔═╡ d2de844f-4124-4cf4-80c7-52300ded90d3
begin
	bcs = Array{Symbolics.Equation}(undef, 2N + 2N) # 2 ICs + 2x2 BCs
		for (i, j) in zip(1:N, c)
			# Initial conditions
			bcs[i] 		= C(0.0, x)[i] ~ C₀[j] 
			bcs[i + N] 	= q(0.0, x)[i] ~ q₀[j] 
			# Boundary conditions - Dirichlet for Cₓ₀, Neumann for the rest
			bcs[i + 2N] = C(t, 0.0)[i] ~ Cx0[i]/Cref
			bcs[i + 3N] = Dx(q(t, 0.0)[i]) ~ qx0[i]
		end
	bcs
end

# ╔═╡ 1abc008b-4a60-44fc-b624-52c4840aebbb
md"""
## Model equations
"""

# ╔═╡ 3fe71419-b2e4-4705-b5fd-655c8e8c69df
md"""
##### Equilibrium isotherm
"""

# ╔═╡ 8975d13c-bf27-4a78-812f-21122ff45aab
function qᵢ′(i; N=N)
	
	# denominator
	den = sum([K[j]*C(t,x)[j] for j in 1:N])
	
	# numerator
	num = K[i]*C(t, x)[i]
	
	return num/den 
end

# ╔═╡ 01039316-bdd0-4b08-a75a-ad0d63b781b2
qᵢ′(2)

# ╔═╡ cd0d9c44-9299-4c9d-a6d2-2f197d21c201
qᵢ′(2; N=2)

# ╔═╡ 77e04104-139c-4be3-862c-5d2655c7bdca
md"""
### Averaged mobility, $k_f = f(C_i, D_i)$
"""

# ╔═╡ 303dab3b-b49e-4cc2-9c32-22d34c23a9a0
function Dpi()
	Dpi = sum([Dp[i]*q(t,x)[i] for i in 1:N])/sum([q(t,x)[j] for j in 1:N])
end

# ╔═╡ 6d9c39d7-b725-4722-bec7-684173c8e3d5
Dpi()

# ╔═╡ 72b1485f-2d7f-4a06-adce-4e05c37aced7
begin
	eqs_same = Array{Symbolics.Equation}(undef, 2N)
	for i in 1:N
		eqs_same[i] = Dt(C(t,x)[i]) ~  (1/Pe[i])*Dxx(C(t,x)[i]) +
									  - Dx(C(t,x)[i]) - Dg*Dt(q(t,x)[i]) 
									
		eqs_same[N+i] = Dt(q(t,x)[i]) ~ 15/rₚ^2*Dpi()*τ*(qᵢ′(i) - q(t,x)[i])
	end
	eqs_same
end

# ╔═╡ e6ec0b6c-ae86-406f-8a18-b629ab39dda6
bcs_same = let
	bcs = Array{Symbolics.Equation}(undef, 2N + 2N) # 2 ICs + 2x2 BCs
		for (i, j) in zip(1:N, c)
			# Initial conditions
			bcs[i] 		= C(0.0, x)[i] ~ C₀[j]
			bcs[i + N] 	= q(0.0, x)[i] ~ q₀[j]
			# Boundary conditions - Dirichlet for Cₓ₀, Neumann for the rest
			bcs[i + 2N] = C(t, 0.0)[i] ~ Cx0[i]/Cref
			# bcs[i + 2N] = Dx(C(t, 0.0)[i]) ~ Pe[i]*(C(t, 0.0)[i] - Cx0[i]/Cref)
			bcs[i + 3N] = Dx(q(t, 0.0)[i]) ~ qx0[i]
		end
	bcs
end

# ╔═╡ ca9830b5-0b56-4ec7-a2a0-4841fd7f30c7
@named pdesys_same = PDESystem(eqs=eqs_same, bcs=bcs_same, domain=domains, 
							   ivs=[t, x], dvs=dvs, ps=ps, defaults=pars)

# ╔═╡ cb90e8da-6d8a-49b7-9677-2dbbcf0745e2
md"""
### Full conservation, $k_f = f(C_i, D_i, z_i)$

Helfferich (1976)
"""

# ╔═╡ b3fadf2c-4a39-4dcc-9f6a-77a38584c000
md"""
##### Mobility
"""

# ╔═╡ 0eeaec2e-28ab-4826-973d-8a650d31373a
function Dij(i, j; n=N)
	den = sum([Dp[j] * q(t,x)[j] for j in 1:N]) #*Z[j]
	if i == j
		num =  Dp[i] * (Dp[i] - Dp[n]) * q(t,x)[i]
		return Dp[i] - num/den
	else
		num = -Dp[i] * (Dp[j] - Dp[n]) * q(t,x)[i]
		return num/den
	end
end

# ╔═╡ 96891cfa-bd75-43e3-8a42-41b019d2fbc8
Dij(1, 1) |> simplify

# ╔═╡ e320c41b-a06c-4953-b258-de6ffa983a6a
Dij(1, 2) |> simplify

# ╔═╡ 34c81e22-e726-4183-b10a-f9ea21e0acc1
Dij(1, 5)

# ╔═╡ a0660c1b-5497-43bb-8da4-d32b6c2cb759
Dij(5, 5)

# ╔═╡ 90ea266e-c8c8-4848-9b9e-58e867c71b23
md"""
##### Rate equation
"""

# ╔═╡ fd909090-e1be-42ee-ba88-45c37793dfbd
begin
	eqs_full = Array{Symbolics.Equation}(undef, 2N)
	for i in 1:N
		eqs_full[i]   = Dt(C(t,x)[i]) ~ (1/Pe[i]) * Dxx(C(t,x)[i]) +
										- Dx(C(t,x)[i]) - Dg * Dt(q(t,x)[i]) 
									
		eqs_full[N+i] = Dt(q(t,x)[i]) ~ 
					15/rₚ^2 * τ * sum([Dij(i,j) * (qᵢ′(j) - q(t,x)[j]) for j in 1:N])
	end
	eqs_full
end

# ╔═╡ e5238e11-583b-492c-9cc9-79d6c48583d0
begin
	eqs_full_nodiff = Array{Symbolics.Equation}(undef, 2N)
	for i in 1:N
		eqs_full_nodiff[i] 	 = Dt(C(t,x)[i]) ~ - Dx(C(t,x)[i]) - Dg*Dt(q(t,x)[i]) 
									
		eqs_full_nodiff[N+i] = Dt(q(t,x)[i]) ~ 
					15/rₚ^2 * τ * sum([Dij(i,j) * (qᵢ′(j) - q(t,x)[j]) for j in 1:N])
	end
	eqs_full_nodiff
end

# ╔═╡ b2786e24-3119-42b9-9b09-eddb2c14c80d
md"""
##### Problem
"""

# ╔═╡ fc389dab-2d23-4923-845a-d895399e57da
@named pdesys_full = PDESystem(eqs=eqs_full, bcs=bcs, domain=domains, 
							   ivs=[t, x], dvs=dvs, ps=ps, defaults=pars);

# ╔═╡ bb99d2a3-dc1e-4418-add5-47bebca9124a
pdesys_full |> typeof

# ╔═╡ 05ff089a-6d32-482b-85d2-c985feea65fd
md"""
## Discretization
"""

# ╔═╡ 60499827-399c-4053-aa56-731f20012ad5
function discretize_problem(pdesys::PDESystem; dx=0.025, grid_align=edge_align)
	"""
	Function to discretize a set of PDEs with the MOL. Default grid size is Np = 22 points, with two additional ghost points outside of the domain at distance dx/2.
	"""
	
	discretization = MOLFiniteDifference([x => dx], t, grid_align=grid_align, advection_scheme=UpwindScheme(), approx_order=2)
	
	problem = discretize(pdesys, discretization)
end

# ╔═╡ aa96d243-3dd3-4a0d-bb57-c9494cb35b8e
md"""
##### Average mobility
"""

# ╔═╡ e7bcaaac-64c6-456b-b52b-f6006bf05f37
# @time prob22S = discretize_problem(pdesys_same, dx=0.05, grid_align=edge_align);

# ╔═╡ e7d5f938-b446-4e71-93ae-ec77a55f18ed
@time prob42S = discretize_problem(pdesys_same, dx=0.025);

# ╔═╡ 659597cf-5258-4997-80de-1a866988a3b4
md"""
##### Discretization analysis
"""

# ╔═╡ 4bc7ce04-9f4e-4779-88fc-823d3ffb6ecb
bcs_disc = let
	bcs = Array{Symbolics.Equation}(undef, 2N + 4N) # 2 ICs + 2x2 BCs
		for (i, j) in zip(1:N, c)
			# Initial conditions
			bcs[i] 		= C(0.0, x)[i] ~ C₀[j]
			bcs[i + N] 	= q(0.0, x)[i] ~ q₀[j]
			# Boundary conditions - Dirichlet for Cₓ₀, Neumann for the rest
			bcs[i + 2N] = C(t, 0.0)[i] ~ Cx0[i]/Cref
			# bcs[i + 2N] = Dx(C(t, 0.0)[i]) ~ Pe[i]*(C(t, 0.0)[i] - Cx0[i]/Cref)
			bcs[i + 3N] = Dx(q(t, 0.0)[i]) ~ 0.0
			bcs[i + 4N] = Dx(C(t, 1.0)[i]) ~ 0.0
			bcs[i + 5N] = Dx(q(t, 1.0)[i]) ~ 0.0
		end
	bcs
end

# ╔═╡ f98fc0ba-72c7-44a5-9b15-0b9b30e6b13e
@named pdesys_disc = PDESystem(eqs=eqs_same, bcs=bcs_disc, domain=domains, 
							   ivs=[t, x], dvs=dvs, ps=ps, defaults=pars);

# ╔═╡ ef6206f7-305e-4199-b023-97fa55107dc4
@time prob10C = discretize_problem(pdesys_disc, dx=0.1, grid_align=center_align);

# ╔═╡ 37698f81-a6a2-4843-91b2-09eb8ef902e9
@time prob10E = discretize_problem(pdesys_disc, dx=0.1, grid_align=edge_align);

# ╔═╡ 754f07f9-9fed-430e-804b-d40a292dded4
@time solC = solve(prob10C, Rodas5P(), saveat=1, progress=true, progress_steps=1);

# ╔═╡ 37f89836-caf0-466a-abcf-d14041cc559c
@time solE = solve(prob10E, Rodas5P(), saveat=1, progress=true, progress_steps=1);

# ╔═╡ 88a4f5ba-9857-45fc-87bd-9da63b1da4af
solC[x]

# ╔═╡ 48360ff2-9923-41a4-bac1-448c7f5a6131
prob10C.problem_type.discretespace.grid[x]

# ╔═╡ 0d4433db-9a67-4a5a-bd40-4dbfa7bf0f43
prob10C.problem_type.discretespace.axies[x]

# ╔═╡ 40ce2635-e918-455c-889f-4ab43a2dfa15
solE[x]

# ╔═╡ a8de4251-3319-429f-aace-c18adccfb9ff
prob10E.problem_type.discretespace.grid[x]

# ╔═╡ 0bf5f346-4046-433d-bbcc-7fb04ac0ef0d
prob10E.problem_type.discretespace.axies[x]

# ╔═╡ 6542cd88-7816-4b78-a6ab-8e76967d98ff
let
	Nt = length(solC[t])
	Nt = 401
	m = 34
	n = length(1:m:Nt)
	col = permutedims(cp[1:n])
	Cmax = 3.3

	# Profiles
	p1 = plot(ylims=(0, Cmax), 
			  legendposition=:topright, legendcolumns=1)
	plot!([], [[], []]; lab=["center" "edge"], c=cp[1],
		  ls=[:solid :dot], lw=[1 2], m=[:circle :utriangle])

	# edge-align
	plot!(solE[x], [solE[C(t,x)[2]][i, :] for i in 1:m:Nt]; 
		  lab="", ls=:dot, c=col, lw=2)
	scatter!(solE[x], [solE(i, :, dv=C(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:utriangle, c=col)
	scatter!(solC[x], [solE(i, solC[x], dv=C(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:utriangle, c=col)

	# center-align
	plot!(solC[x], [solC[C(t,x)[2]][i, :] for i in 1:m:Nt]; lab="", c=col)
	scatter!(solC[x], [solC(i, :, dv=C(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:circle, c=col)
	x_ext = solE[x][2:end-1]
	scatter!(x_ext, [solC(i, x_ext, dv=C(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:circle, c=col) 
	

	# Breakthrough
	p2 = plot(xlims=(0, Nt-1), ylims=(0, Cmax), 
			  legendposition=:topright, legendcolumns=2)
	# plot!([], [[]]; lab=["center" "edge"], 
	# 	  ls=[:solid :dot], lw=[1 2], m=[:circle :utriangle])

	
	# inlet
	# center-align
	# plot!(solC[t][1:Nt], solC[C(t,x)[2]][1:Nt, 1]; 
	# 	  lab="inlet (actual)", lw=1, c=cp[1])
	plot!(solC[t][1:Nt], solC(0:Nt-1, 0, dv=C(t,x)[2]); 
		  lab="center (0.00)", lw=2, c=cp[1], ls=:dot)
	
	# edge-align
	plot!(solE[t][1:Nt], solE(0:Nt-1, 0, dv=C(t,x)[2]); 
		  lab="edge (0.00)", lw=2, c=cp[4], ls=:dot)
	plot!(solE[t][1:Nt], solE[C(t,x)[2]][1:Nt, 1]; 
		  lab="edge (-0.05)", lw=2, c=cp[2])
	plot!(solE[t][1:Nt], solE[C(t,x)[2]][1:Nt, 2]; 
		  lab="edge (0.05)", lw=2, c=cp[3])
	

	# breakthrough
	# center-align
	# plot!(solC[t][1:Nt], solC[C(t,x)[2]][1:Nt, end]; 
	# 	  lab="center (actual)", lw=3, c=cp[1])
	plot!(solC[t][1:Nt], solC(0:Nt-1, 1, dv=C(t,x)[2]); 
		  lab="center (1.00)", lw=3, c=cp[1], ls=:dot)
	scatter!(solC[t][1:m:Nt], [solC[C(t,x)[2]][1:m:Nt, end]]; lab="", c=cp[1])
	scatter!(solC[t][1:m:Nt], [solC(0:m:Nt-1, 1, dv=C(t,x)[2])]; 
			 lab="", c=cp[1], m=:diamond)

	# edge-align
	plot!(solE[t][1:Nt], solE(0:Nt-1, 1, dv=C(t,x)[2]); 
		  lab="edge (1.00)", c=cp[4], lw=3, ls=:dashdot)
	plot!(solE[t][1:Nt], solE[C(t,x)[2]][1:Nt, end]; 
		  lab="edge (1.05)", c=cp[2], lw=3, ls=:dot)
	plot!(solE[t][1:Nt], solE[C(t,x)[2]][1:Nt, end-1]; 
		  lab="edge (0.95)", c=cp[3], lw=3, ls=:dash)
	
	scatter!(solE[t][1:m:Nt], solE[C(t,x)[2]][1:m:Nt, end]; 
			 lab="", c=cp[2], m=:utriangle)
	scatter!(solE[t][1:m:Nt], solE[C(t,x)[2]][1:m:Nt, end-1]; 
			 lab="", c=cp[3], m=:dtriangle)
	scatter!(solE[t][1:m:Nt], solE(0:m:Nt-1, 1, dv=C(t,x)[2]); 
		     lab="", c=cp[4], m=:diamond)

	
	# Accumulated in solid [eq.]
	Nt = 500
	Vres = prob10C.ps[V]*(1-prob10C.ps[ε])*1e3
	BV = solC[t][1:Nt]*prob10C.ps[τ]*prob10C.ps[Q]*1e3
	Vmax = BV[Nt]

	qe =  prob10C.ps[qmax]*prob10C.ps[K[2]]*Cₓ₀[c[2]] / 
			sum([Cₓ₀[c[j]]*prob10C.ps[K[j]] for j in 1:N])
	
	p3 = plot(xlabel="Volume [L]", title="Solid", 
			  xlims=(0, Vmax), ylims=(-0.001, Vmax*Cₓ₀[:Na]))
	
	# inlet
	plot!(BV, BV*Cₓ₀[:Na]; lab="", ls=:dashdot, c=:black)

	# equilibrium
	hline!([qe*Vres]; lab="", c=:blue, ls=:dot)

	# center-align
	plot!(BV, prob10C.ps[qmax]*Vres*
		  integrate(solC[x], solC(solC[t][1:Nt], :, dv=q(t,x)[2]); dims=1);
		  lab="", c=cp[1], ls=:solid, lw=2)
	# plot!(BV, prob10C.ps[qmax]*Vres*
	# 	  integrate(solC[x], solC[q(t,x)[2]][1:Nt, :]; dims=1);
	# 	  lab="", c=cp[1], ls=:solid, lw=2)

	# edge-align
	plot!(BV, prob10E.ps[qmax]*Vres*
		  integrate(solE[x], solE(solE[t][1:Nt], :, dv=q(t,x)[2]); dims=1);
		  lab="", c=cp[2], ls=:solid, lw=2)
	plot!(BV, prob10E.ps[qmax]*Vres*
		  integrate(solE[x][2:end], 
					solE(solE[t][1:Nt], solE[x][2:end], dv=q(t,x)[2]); dims=1);
		  lab="", c=cp[2], ls=:dot, lw=2)
	plot!(BV, prob10E.ps[qmax]*Vres*
		  integrate(solE[x], solE[q(t,x)[2]][1:Nt, :]; dims=1);
		  lab="", c=cp[3], ls=:dash, lw=2)
	plot!(BV, prob10E.ps[qmax]*Vres*
		  integrate(solE[x][2:end], solE[q(t,x)[2]][1:Nt, 2:end]; dims=1);
		  lab="", c=cp[4], ls=:dashdot, lw=2)

	
	# Accumulated in liquid [eq.]
	p4 = plot(xlabel="Volume [L]", title="Liquid", 
			  xlims=(0, Vmax), ylims=(-0.001, Vmax*Cₓ₀[:Na]))
	
	# inlet
	plot!(BV, BV*Cₓ₀[:Na]; lab="", ls=:dashdot, c=:black)

	# equilibrium
	plot!(BV, BV*Cₓ₀[:Na] .- qe*Vres; lab="", ls=:dot, c=:blue, lw=1.5)

	# verification
	plot!(BV, BV*Cₓ₀[:Na] .- prob10E.ps[qmax]*Vres*
		  integrate(solE[x][2:end], solE[q(t,x)[2]][1:Nt, 2:end]; dims=1); 
		  lab="", ls=:solid, c=:black)

	# center-align
	plot!(BV[2:Nt], Cₓ₀[:Na]*
		  [integrate(BV[1:i], solC(solC[t][1:i], 1, dv=C(t,x)[2])) for i in 2:Nt];
		  lab="", c=cp[1], ls=:solid, lw=2)

	# edge-align
	plot!(BV[2:Nt], Cₓ₀[:Na]*
		  [integrate(BV[1:i], solE(solE[t][1:i], 1, dv=C(t,x)[2])) for i in 2:Nt];
		  lab="", c=cp[2], ls=:dot, lw=2)
	plot!(BV[2:Nt], Cₓ₀[:Na]*
		  [integrate(BV[1:i], solE[C(t,x)[2]][1:i, end]) for i in 2:Nt];
		  lab="", c=cp[3], ls=:dot, lw=2)
	plot!(BV[2:Nt], Cₓ₀[:Na]*
		  [integrate(BV[1:i], solE[C(t,x)[2]][1:i, end-1]) for i in 2:Nt];
		  lab="", c=cp[4], ls=:dot, lw=2)	

	l = @layout [a b; c{0.3h}; d{0.3h}]
	figure = plot!(p1, p2, p3, p4, layout=l, dpi=300, size=(600, 800), 
				   suptitle="Center vs. edge discretization") 
	figure |> x -> aside(x; v_offset=1400)
end

# ╔═╡ 090782ac-3e17-453b-953d-911615b6a45d
let
	Nt = length(solC[t])
	Nt = 401
	m = 50
	n = length(1:m:Nt)
	col = permutedims(cp[1:n])
	Cmax = 1.0

	# Profiles
	p1 = plot(ylims=(0, Cmax), 
			  legendposition=:topright, legendcolumns=1)
	plot!([], [[], []]; lab=["center" "edge"], c=cp[1],
		  ls=[:solid :dot], lw=[1 2], m=[:circle :utriangle])

	# edge-align
	plot!(solE[x], [solE[q(t,x)[2]][i, :] for i in 1:m:Nt]; 
		  lab="", ls=:dot, c=col, lw=2)
	scatter!(solE[x], [solE(i, :, dv=q(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:utriangle, c=col)
	scatter!(solC[x], [solE(i, solC[x], dv=q(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:utriangle, c=col)

	# center-align
	plot!(solC[x], [solC[q(t,x)[2]][i, :] for i in 1:m:Nt]; lab="", c=col)
	scatter!(solC[x], [solC(i, :, dv=q(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:circle, c=col)
	x_ext = solE[x][2:end-1]
	scatter!(x_ext, [solC(i, x_ext, dv=q(t,x)[2]) for i in 0:m:Nt-1];
		  lab="", lw=2, ls=:dot, m=:circle, c=col) 

	# sum of profiles
	scatter!(solC[x], 
			 [sum([solC[q(t,x)[i]][j, :] for i in 1:N]) for j in 1:m:Nt]; 
			 lab="", lw=2, c=cp[1], m=:circle)
	scatter!(solC[x], 
			 [sum([solC(j, :, dv=q(t,x)[i]) for i in 1:N]) for j in 1:m:Nt]; 
			 lab="", lw=2, c=cp[1], m=:diamond)
	
	scatter!(solE[x], 
			 [sum([solE[q(t,x)[i]][j, :] for i in 1:N]) for j in 1:m:Nt]; 
			 lab="", lw=2, c=cp[2], m=:utriangle)
	scatter!(solE[x], 
			 [sum([solE(j, :, dv=q(t,x)[i]) for i in 1:N]) for j in 1:m:Nt]; 
			 lab="", lw=2, c=cp[2], m=:dtriangle)
	

	# Breakthrough
	p2 = plot(xlims=(0, Nt-1), ylims=(0, Cmax), 
			  legendposition=:topright, legendcolumns=2)
	# plot!([], [[]]; lab=["center" "edge"], 
	# 	  ls=[:solid :dot], lw=[1 2], m=[:circle :utriangle])

	
	# inlet
	# center-align
	# plot!(solC[t][1:Nt], solC[C(t,x)[2]][1:Nt, 1]; 
	# 	  lab="inlet (actual)", lw=1, c=cp[1])
	plot!(solC[t][1:Nt], solC(0:Nt-1, 0, dv=q(t,x)[2]); 
		  lab="center (0.00)", lw=2, c=cp[1], ls=:dot)
	
	# edge-align
	plot!(solE[t][1:Nt], solE(0:Nt-1, 0, dv=q(t,x)[2]); 
		  lab="edge (0.00)", lw=2, c=cp[4], ls=:dot)
	plot!(solE[t][1:Nt], solE[q(t,x)[2]][1:Nt, 1]; 
		  lab="edge (-0.05)", lw=2, c=cp[2])
	plot!(solE[t][1:Nt], solE[q(t,x)[2]][1:Nt, 2]; 
		  lab="edge (0.05)", lw=2, c=cp[3])
	

	# breakthrough
	# center-align
	# plot!(solC[t][1:Nt], solC[C(t,x)[2]][1:Nt, end]; 
	# 	  lab="center (actual)", lw=3, c=cp[1])
	plot!(solC[t][1:Nt], solC(0:Nt-1, 1, dv=q(t,x)[2]); 
		  lab="center (1.00)", lw=3, c=cp[1], ls=:dot)
	scatter!(solC[t][1:m:Nt], [solC[q(t,x)[2]][1:m:Nt, end]]; lab="", c=cp[1])
	scatter!(solC[t][1:m:Nt], [solC(0:m:Nt-1, 1, dv=q(t,x)[2])]; 
			 lab="", c=cp[1], m=:diamond)

	# edge-align
	plot!(solE[t][1:Nt], solE(0:Nt-1, 1, dv=q(t,x)[2]); 
		  lab="edge (1.00)", c=cp[4], lw=3, ls=:dashdot)
	plot!(solE[t][1:Nt], solE[q(t,x)[2]][1:Nt, end]; 
		  lab="edge (1.05)", c=cp[2], lw=3, ls=:dot)
	plot!(solE[t][1:Nt], solE[q(t,x)[2]][1:Nt, end-1]; 
		  lab="edge (0.95)", c=cp[3], lw=3, ls=:dash)
	
	scatter!(solE[t][1:m:Nt], solE[q(t,x)[2]][1:m:Nt, end]; 
			 lab="", c=cp[2], m=:utriangle)
	scatter!(solE[t][1:m:Nt], solE[q(t,x)[2]][1:m:Nt, end-1]; 
			 lab="", c=cp[3], m=:dtriangle)
	scatter!(solE[t][1:m:Nt], solE(0:m:Nt-1, 1, dv=q(t,x)[2]); 
		     lab="", c=cp[4], m=:diamond)

	figure = plot!(p1, p2, dpi=300, size=(600, 400), 
				   suptitle="Center vs. edge discretization") 
	figure |> x -> aside(x; v_offset=200)
end

# ╔═╡ a8bb2a58-08f3-4495-a22e-f4fed33a936e
md"""
##### Full conservation
"""

# ╔═╡ c8b83cee-caa9-4087-8487-aa93a39c0e9b
# @time prob22H = discretize_problem(pdesys_full, dx=0.05, grid_align=edge_align);

# ╔═╡ b50e0b99-434c-415d-9639-afdf5e8726d3
# grid_22 = prob22H.problem_type.discretespace.grid[x]

# ╔═╡ b9aadf8f-5a77-46e5-986a-bb4c2a314335
# @time prob42H = discretize_problem(pdesys_full, dx=0.025, grid_align=edge_align);

# ╔═╡ 410d04f4-8158-41a1-aadc-59277cf260dd
md"""
## Functions
"""

# ╔═╡ 460c29de-39bf-4dcc-afb0-be20d117860a
md"""
##### Parameters
"""

# ╔═╡ f174d441-2934-478e-beeb-3b2e50c796ca
# show default numerical value for symbol – substitute dependencies
function val!(p; params=params, times=2)
	temp = copy(params)
	for i in 1:times
		for k in keys(temp)
			if k isa Symbolics.Num
				temp[k] = substitute(temp[k], temp) |> Symbolics.eval |> Symbolics.value
				# @assert temp[k] isa Symbolics.Num "Value must be a Float"
			elseif k isa Symbolics.Arr{Num, 1}
	        	for j in 1:length(temp[k])
					temp[k][j] = substitute(temp[k][j], temp) |> Symbolics.eval
					# @assert temp[k][j] isa Symbolics.Num "Value must be a Float"
				end
			end
	    end
	end
	out = Dict{Any, Any}()
	for k in keys(temp)
		out[k] = temp[k] |> Symbolics.value
    end
	out[p]
	
	# out = Symbolics.value(temp[p])
	# out = float(temp[p])
	# temp[p] |> Symbolics.value #|> float
end

# ╔═╡ ad0e7db9-e293-4aa6-8e33-683a9ed05440
# substitute all dependencies and provide all numerical values in a dictionary
function value!(params::Dict; times=3, verbose=false)
    temp = deepcopy(params)
	for time in 1:times # repeat substitution
	    for k in keys(temp)
			if k isa Symbolics.Num
				temp[k] = substitute(temp[k], temp) |> Symbolics.eval
				# @assert temp[k] isa Real
			elseif k isa Symbolics.Arr{Num, 1}
	        	for j in 1:length(temp[k])
					temp[k][j] = substitute(temp[k][j], temp) |> Symbolics.eval
					# @assert temp[k][j] isa Symbolics.Num
				end
			end
	    end
	end
	verbose && println("Dependencies updated.") 
	out = Dict{Any, Any}()
	for time in 1:times 
		for k in keys(temp)
			out[k] = temp[k] |> Symbolics.value
    	end
	end
    out
end

# ╔═╡ e6862433-a14d-41b1-8f0f-970c47f6f14f
pars |> deepcopy |> value!

# ╔═╡ df52ae4c-97d4-48f7-9a5d-ef281b26577d
# update provided parameters for the given problem before running simulation
function update!(params, prob, pairs::Vararg{Pair{<:Union{Symbolics.Num, Symbolics.Arr{Symbolics.Num, 1}},<:Union{Any,Vector{Any}}}}; verbose=false)
	
	temp = deepcopy(params)

	# update parameters
	verbose && println("Updates:", pairs) 
	verbose && println("Updated:") 
	for (key, value) in pairs
		# println(value |> typeof)
		temp[key] = value
		verbose && println(key, "=>", temp[key]) 
	end

	# update all dependencies
	out = value!(temp; verbose=verbose)

	# update all dependencies in problem
	for k in keys(out)
		if prob.ps[k] isa Float64
			prob.ps[k] = out[k]
		elseif prob.ps[k] isa SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}
			[prob.ps[k][i] = Symbolics.value(out[k][i]) for i in 1:length(out[k])]
		end
	end
	out
end

# ╔═╡ 1fcd177b-dee2-4e56-8817-fa57bc4eea0c
md"""
### Breakthrough detection
"""

# ╔═╡ 569fa407-e152-4c5b-aa8c-dd696ac2326d
function BTdetection(solved; threshold=0.05, flip=false, tmin=10)
	if flip == false
		index = findfirst(x -> x > threshold, solved[C(t,x)[2]][tmin:end,end]) + tmin - 1
	elseif flip == true
		index = findfirst(x -> -x > -threshold, solved[C(t,x)[2]][tmin:end,end]) + tmin - 1
	end
	tload = solved[t][index]
	return tload, index
end

# ╔═╡ 44f56d11-a5a1-444e-826b-a86515ce78d9
# let
# 	threshold = 0.2
# 	tload, index = BTdetection(solreg; threshold=threshold, flip=false)
# 	plot(solreg[C(t,x)[2]][:,end]; lab="elution")
# 	solreg[C(t,x)[2]][index,end] |> println
# 	vline!([index]; lab=:none, ls=:dash, c=:black, xlims=(0, 500))
# 	hline!([threshold]; lab=L"$C/C_0$ = %$threshold", ls=:dashdot, c=:red)
# end

# ╔═╡ 3a359f7e-919e-41aa-94e6-c1004fbe0a62
md"""
##### pH slope detection
"""

# ╔═╡ 67b555d4-edac-4731-b167-06e2265f1dcc
function pHdetection(sol; slope=0.5, tmin=10, flip=false, Cref=Cₓ₀[:Na])
	# function to detect changes in pH, returns time point and time vector index
	# tmin |> println
	pH  = -log10.(Cref * sol[C(t,x)[1]][:, end])
	# pH |> length |> println
	dpH = diff(pH)./diff(sol[t]) # diff removes one element from pH
	# dpH |> length |> println
	# dpH[tmin:end] |> length |> println
	# TODO rename tmin to skip? (means x elements from vector are skipped)
	
	if flip == false
		index = findfirst(x->x >= slope, dpH[tmin:end]) + tmin
	elseif flip == true
		index = findfirst(x->-x <= slope, dpH[tmin:end]) + tmin
	end
	# index |> println 
	# dpH[index] |> println
	tload = sol[t][index]
	# tload |> println 
	return tload, index
end

# ╔═╡ 509a586c-0be9-4a11-aec3-ac3e71cf29d7
md"""
##### Plotting
"""

# ╔═╡ f596e2e6-81a1-4fd6-9a2e-d9ce5b03ec3a
function plot_data(data_t, data_y, 
				   pairs::Vararg{Pair{<:Union{Symbolics.Num, Symbolics.Arr{Symbolics.Num, 1}},<:Union{Any,Vector{Any}}}}; 
				   Vmax = 60., pHmax = 6., Cmax = 4.55, gray_Mg=5, gray_Ca=4, 
				   prob=prob22H, title="")

	# Plot options
	xlims = (0, Vmax)
	ylims = (-0.15, Cmax)
	ylims_pH = (-0.5, pHmax) 
	yticks_pH = 0:2:pHmax
	label=["Na", "K", "Mg", "Ca"]
	marker=[:square, :circle, :diamond, :utriangle]

	# Update model parameters
	update!(pars, prob, pairs...)

	# Convert data units
	data_BV = data_t*60*prob.ps[Q]*1e3
	
	# pH
	p1 = plot(xlims=xlims, ylims=ylims_pH, ylabel="pH", yticks=yticks_pH, 
			  legendposition=:topleft)
	plot!(data_BV, data_y[:,1]; lab="pH", c=:grey, marker=:circle, lw=1, grid=:off)
	plot!([], []; lab="ΔpH/Δt", c=:grey, ls=:dash, lw=1, marker=:xcross, 
		  markerstrokewidth=2)
	plot!(twinx(), data_BV[2:end], diff(data_y[:,1]) ./ diff(data_BV); 
		lab=:none, c=:grey, ls=:dash, lw=1, marker=:xcross, markerstrokewidth=2,
		xlims=(0, Vmax), ylims=ylims_pH./2, yticks=yticks_pH./2, ylabel="ΔpH/Δt", grid=:off, minorgrid=:off, minorticks=1)
	gridlines = 0:2:pHmax
	for g in gridlines
	    hline!([g], color=:gray, alpha=0.1, linestyle=:solid, lab=:none)
	end
	
	# ylims = (-0.15, 4.55)
	# Breakthrough
	p2 = plot(xlims=xlims, ylims=ylims, grid=:off,
			  xlabel="Treated effluent volume [L]", 
			  ylabel="Normalized concentration", legendposition=:topleft)
	plot!(twinx(); xlims=xlims, ylims=ylims.*Cₓ₀mass[:Na], grid=:off, 
		  ylabel="Concentration [mg/L]")
	
	for i in 2:5
		plot!(data_BV, data_y[:,i]./Cₓ₀mass[c[i]]; 
			  c=cp[i-1], lab=label[i-1], marker=marker[i-1])
	end
	# gray out last data points for Mg
	r = gray_Mg-1
	for i in 4
		plot!(data_BV[end-r:end], data_y[end-r:end,i]/Cₓ₀mass[c[i]]; 
			  lab="", marker=marker[i-1], c=:lightgrey, markerstrokecolor=cp[i-1], 
			  linecolor=cp[i-1])
	end
	# gray out last data points for Ca
	r = gray_Ca-1
	for i in 5
		plot!(data_BV[end-r:end], data_y[end-r:end,i]/Cₓ₀mass[c[i]]; 
			  lab="", marker=marker[i-1], c=:lightgrey, markerstrokecolor=cp[i-1], 
			  linecolor=cp[i-1])
	end
	
	hline!([1.]; lab=:none, lw=1, ls=:dash, c=:black)
	gridlines = 0:1:Cmax
	for g in gridlines
	    hline!([g], color=:gray, alpha=0.1, linestyle=:solid, lab=:none)
	end

	# Plot figure
	plot(p1, p2, layout=grid(2, 1, heights=[0.3, 0.7]), dpi=300, size=(680, 400), 
		 suptitle=title)
end

# ╔═╡ 7e1c7fa6-e329-417d-95e2-ba2896a05cd2
# plot_data(t4, y4, Q=>105/60/1e6; Vmax=70., Cmax=5, title="Calibration dataset, pH 4")

# ╔═╡ ae7704b0-ad6a-41b0-9f05-67bababe40d6
# plot_model(t4, y4, 
# 		   Cx0  => [Cₓ₀[i] for i in c], 
# 		   Q 	=> 105/60*1e-6, 
# 		   qmax => 2.6, 
# 		   K 	=> [Ki[i] for i in c], 
# 		   R 	=> vcat(10., [Ri[i] for i in c[2:N]]),
# 		   ; prob=prob42H, Vmax=70., units="L")

# ╔═╡ 47bbd2ab-0829-4136-84de-f98f7b3408a7
md"""
##### Regeneration
"""

# ╔═╡ 9f1329e5-04b0-4395-983a-8dad2cf4bef1
function sliders(parameters::Vector, widths::Vector, defaults::Vector)
	
	return confirm(combine() do Child
		
		inputs = [
			md""" $(name): $(
				Child(name, Slider(width, show_value=true, default=default))
			)"""
			
			for (name, width, default) in zip(parameters, widths, defaults)
		]
		
		md"""
		$(inputs)
		"""
	end)
end

# ╔═╡ d7f657c3-e540-4e2d-933a-5a3534957015
function updateIC!(sol, prob, tstep; reverse=false)
	# Np = length(sol[x]) - 2
	# Nx = sol.ivdomain[2]
	Nx = prob.problem_type.discretespace.axies[x]; #Nx |> println
	Np = length(Nx)
	# i  = findfirst(x -> x >= tstep, sol[t])  # index for time in seconds
	# reverse ? begin idcs = Np+1:-1:2 end : idcs = 2:Np+1
	reverse ? begin idcs = Np:-1:1 end : idcs = 1:Np

	temp = deepcopy(prob)

	# For each node, q and C alternated for each component, from last to first
	for j in 1:N
		temp.u0[@. 2Np*(N-j)+2*(1:Np)-1] = sol(tstep, Nx[idcs], dv=q(t,x)[j])
		# temp.u0[@. 2Np*(N-j)+2*(1:Np)-1] = sol[q(t,x)[j]][i, idcs]
	end
	
	return temp
end

# ╔═╡ 33e1d1ad-85f4-48ad-933b-a85f2203c300
md"""
## Experimental data
"""

# ╔═╡ 67c6d1c1-ba68-4966-a3cb-f23d38a9049a
md"""
##### Calibration, pH 4
"""

# ╔═╡ df88c1dc-ad21-4ef4-b41e-416ee3543a4a
t4, y4 = data_t4, data_y4 = begin

	t4 = [10, 50, 100, 140, 197, 217, 237, 257, 262, 267, 272, 292, 312, 332, 352, 372, 397, 447, 507, 537, 567] # [min] 
	pH4 = [2., 2., 1.97, 1.94, 1.87, 1.91, 2.07, 3.14, 3.35, 3.55, 3.66, 3.82, 3.88, 3.91, 3.95, 4.04, 4.02, 4.04, 4.04, 4.06, 4.07] 
	Na4 = [0.22, 0.21, 0.22, 0.22, 7.56, 38.24, 293.04, 455.20, 446.10, 464.72, 428.96, 286.46, 208.76, 173.30, 151.37, 111.46, 109.83, 108.26, 108.61, 111.19, 110.60] #./ mean([108.61, 111.19, 110.60]) * Cₓ₀mass[:Na]
	K4 = [2.76, 2.70, 2.68, 2.70, 2.67, 0., 0., 29.82, 41.92, 63.94, 92.90, 306.86, 479.74, 426.37, 329.25, 135.54, 111.53, 111.89, 109.22, 110.38, 110.25] #./ mean([109.22, 110.38, 110.25]) * Cₓ₀mass[:K]
	Mg4 = [0.12, 0.11, 0.09, 0.10, 0.33, 0., 0., 0., 0., 0., 0., 0., 7.92, 38.12, 81.68, 153.26, 168.06, 168.40, 165.44, 165.42, 162.64] #./ mean([165.44, 165.42, 162.64]) * Cₓ₀mass[:Mg]
	Ca4 = [0.11, 0.09, 0.08, 0.08, 0.36, 0., 0., 0., 0., 0., 0., 0., 0., 9.27, 19.71, 37.44, 40.87, 42.54, 41.45, 44.86, 46.39] #./ mean([ 41.45, 44.86, 46.39]) * Cₓ₀mass[:Ca]
	y4 = hcat(pH4, Na4, K4, Mg4, Ca4)

	t4, y4
end

# ╔═╡ 9ce03269-5734-40cf-a87f-7b9896d93bdd
md"""
##### Loading, pH 6
"""

# ╔═╡ a9853f0c-6e44-42b2-96bb-24fc22164df0
t6, y6 = data_t6, data_y6 = begin

	t6 = [30, 90, 130, 150, 170, 190, 210, 230, 235, 255, 275, 295, 305, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400, 405, 410, 415, 418] # [min] 
	pH6 = [2.05, 2.03, 2.02, 2.01, 2.03, 2.03, 2.03, 2.03, 2.02, 2.02, 2.02, 2.02, 2.03, 2.08, 2.13, 2.21, 2.37, 2.53, 2.68, 2.88, 3.21, 3.48, 3.76, 4.06, 4.36, 4.62, 4.78, 4.87, 4.87, 4.91, 4.92, 4.92, 4.89, 4.89, 4.92]
	Na6 = [0., 0., 0., 0., 0., 0.02, 0., 0.25, 0.38, 0., 2.07, 34.21, 69.40, 121.80, 162.28, 219.50, 303.12, 365.50, 373.02, 393.28, 433.10, 428.43, 432.42, 406.54, 387.56, 361.75, 328.25, 318.22, 298.80, 278.21, 269.30, 258.66, 248.63, 235.50, 230.73]
	K6 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.50, 0.90, 1.66, 3.08, 5.02, 9.86, 19.88, 24.34, 31.73, 50.87, 70.38, 98.31, 125.98, 205.69, 203.80, 222.26, 278.66, 310.05, 336.90, 383.16, 400.25, 414.00, 419.71, 421.32]
	Ca6 = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.69, 1.18, 0., 2.48, 0., 0., 0., 0., 0., 4.32, 0., 0., 0., 0., 0., 0., 6.55, 6.02, 5.35, 6.37, 0., 0., 5.27, 0.]	
	Mg6 = [0.66, 1.02, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4.43, 7.4, 7.01]
	y6 = hcat(pH6, Na6, K6, Mg6, Ca6)

	t6, y6
end

# ╔═╡ 8eb84797-2805-4f94-911d-fedf2ca10fc6
md"""
##### Elution, 0.01 HCl
"""

# ╔═╡ 4e58cf1e-623a-4296-b3fa-a7dd4375c481
t6e, y6e = data_t6e, data_y6e = begin

	t6e = [0, 10, 40, 70, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330] # [min] 
	pH6e = [4.88, 4.88, 4.87, 4.84, 4.83, 4.8, 4.79, 4.78, 4.68, 3.87, 3.21, 3.04, 2.91, 2.94, 2.79, 2.75, 2.7, 2.66, 2.63, 2.62, 2.6, 2.6, 2.57, 2.57, 2.55, 2.54, 2.54, 2.54]
	Na6e = [0, 68.78, 68.78, 113.48, 110.04, 95.14, 88.26, 88.26, 111.19, 83.68, 58.46, 61.90, 63.05, 60.75, 63.05, 57.31, 47.00, 51.58, 43.56, 42.41, 34.39, 34.39, 32.10, 29.80, 29.80, 33.24, 28.66, 25.22]
	K6e = [0, 119.21, 202.89, 185.70, 163.92, 135.26, 114.63, 114.63, 131.82, 96.29, 77.95, 65.34, 72.22, 68.78, 66.48, 61.90, 48.14, 56.17, 49.29, 47.00, 45.85, 43.56, 37.83, 36.68, 42.41, 49.29, 33.24, 32.10]
	Mg6e = zeros(length(t6e))
	Ca6e = zeros(length(t6e))
	y6e = hcat(pH6e, Na6e, K6e, Mg6e, Ca6e)

	t6e, y6e
end

# ╔═╡ 4b8afbf2-6bc7-43b4-b5f1-0bd3e7854698
md"""
##### Elution, 5% HCl
"""

# ╔═╡ f3b30039-70c0-4f8e-85f2-d94bb9a160f3
t6e2, y6e2 = data_t6e2, data_y6e2 = begin

	t6e2 = [0, 5., 10., 15., 35., 40.] # [min] 
	pH6e2 = [7, 0.125, 1.2, 0.86, 0.78, 0.72]
	Na6e2 = [0, 0., 0., 1080., 810., 450.]
	K6e2 = [0, 0., 0., 1470., 2300., 3030.]
	Ca6e2 = [0, 0., 0., 7050., 5850., 4230.]
	Mg6e2 = [0, 0., 0., 8950., 7320., 5510.]
	y6e2 = hcat(pH6e2, Na6e2, K6e2, Mg6e2, Ca6e2)

	t6e2, y6e2
end

# ╔═╡ 54b90a93-992c-42bf-8672-61905b874c47
md"""
# Scenario analysis
"""

# ╔═╡ 4242c3d6-20e4-4a27-9998-6498dafaffc0
md"""

Phases of operation:
1. Loading
2. Regeneration

Backwashes and rinses are often automated and will not be considered.



**Stopping criteria**:
1. Fixed time: in volume or time, often calculated from equivalents and resin capacity
2. pH level:
3. Na level (%, fixed value)
4. pH slope

Model-predictive control → stop at % capacity, model errors

**Performance metrics**:
- Resin capacity (remaining/restored)
- Cumulative ion mass (breakthrough/eluted/adsorbed)
- Regenerant spent (volume-concentration/equivalents)

and these for combinations of loading-regeneration cycles with different criteria.

"""

# ╔═╡ dd1a1903-872f-42cb-b138-ebc044c40c14
md"""
## Influent loading
"""

# ╔═╡ a0f7e2f3-2693-4241-a4c7-9af5617e04dd
begin
	Cₓ₀[:H] = 10.0^(-6); 	Cₓ₀mass[:Na] = 100.0; 	Cₓ₀mass[:K] = 100.0;
							Cₓ₀mass[:Ca] = 100.0; 	Cₓ₀mass[:Mg] = 100.0;
	[Cₓ₀[i] = Cₓ₀mass[i]/(Mi[i]*1e3)*Zi[i] for i in c[2:5]] # [eq/L]
	Cₓ₀
end

# ╔═╡ bff5109f-88ad-4269-8ce5-4a9810f4053e
p0 = [4.3, 1.4, 2.2, 10, 10.5]

# ╔═╡ f1dbaf46-580f-4511-b012-6efb7f582230
sol = let
	prob = prob42S
	update!(pars, prob,	Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
			qmax => p0[1], K => vcat(1.0, p0[2:N]), 
	)
	@time sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)
end;

# ╔═╡ ce14146a-f9dc-425b-bc5d-dc628f33dc56
let
	# sol = sol
	threshold = 0.05
	tload, index = BTdetection(sol; threshold=threshold)
	plot(sol[C(t,x)[2]][:,end]; lab="breakthrough")
	vline!([index]; lab=:none, ls=:dash, c=:black, xlims=(0, 500))
	hline!([threshold]; lab=L"$C/C_0$ = %$threshold", ls=:dashdot, c=:red)
end

# ╔═╡ e2973d09-defd-42cc-93f4-e55e38a3e05a
let
	prob = prob42S
	# sol = sol
	slope = 0.02 # equivalent to ~5% breakthrough Na
	tmin = 4
	tload, index = pHdetection(sol; slope=slope, tmin=tmin, Cref=prob.ps[Cref])
	BV = sol[t][tmin:end]
	pH = -log10.(sol[C(t,x)[1]][tmin:end,end] * prob.ps[Cref])
	dpH = diff(pH)./diff(BV)
	Na = sol[C(t,x)[2]][:,end]
	println(dpH[index-tmin+1])
	
	p1 = plot(BV, pH; lab="pH", lw=2, c=:grey, ls=:dash, ylims=(0, 6))
	hline!([pH[index]]; lab=:none, ls=:dash, c=:black, lw=1)
	scatter!([tload],[pH[index]]; marker=:cross, c=:red, 
			 lab="pH = $(round(pH[index];digits=1))")
	plot!(Na[:,end]; lab="Na", c=cp[1], lw=1, legendposition=:topleft)
	vline!([tload]; lab=:none, ls=:dash, c=:black)
	scatter!([tload],[Na[index]]; marker=:cross, c=:black, 
			 lab="$(round(Na[index]*100;digits=1))%")
	hline!([Na[index]]; c=:blue, ls=:dashdot, lab=:none)
	
	
	p2 = plot(BV[2:end], dpH; lab=:none, ylims=(0, 0.10), yticks=0:0.02:0.1)
	hline!([slope]; lab="slope = $slope", ls=:dashdot, c=:red)
	vline!([tload]; lab=:none, ls=:dash, c=:black)
	scatter!([tload],[slope]; marker=:cross, c=:black, lab=:none) 			# slope
	scatter!([tload],[dpH[index-tmin+1]]; marker=:cross, c=:red, lab=:none) # actual
	plot(p1, p2, layout=(2,1), dpi=300, size=(680, 400), xlims=(0, 500))
end

# ╔═╡ b0b839a1-1bd4-43cc-9f43-325189278f5f
function plot_model(data_t, data_y, 
				   	pairs::Vararg{Pair{<:Union{Symbolics.Num, Symbolics.Arr{Symbolics.Num, 1}},<:Union{Any,Vector{Any}}}}; 
				   	Vmax = 60., pHmax = 6., Cmax = 4.5, gray_Mg=6, gray_Ca=4, 
				   	prob=prob42H, solver=Rodas5P(), progress=true,
					filename="", title="", units="L", BT=0.1)

	# Simulate model
	update!(pars, prob, pairs...)
	@time sol = solve(prob, solver, saveat=1, progress=progress, progress_steps=1)

	# Convert data units
	data_BV = data_t*60*prob.ps[Q]*1e3
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3
	if units == "L"
		xlabel = "Treated effluent volume [L]"
		step = 10
	elseif units == "BV"
		xlabel = "Bed volumes"
		step = 20
		Vmax = 200.
		data_BV /= prob.ps[V]*1e3
		BV /= prob.ps[V]*1e3
	end

	# Plot options
	xlims = (0, Vmax)
	xticks = 0:step:Vmax
	ylims = (-0.15, Cmax)
	ylims_pH = (-0.5, pHmax) 
	yticks_pH = 0:2:pHmax
	label=["Na", "K", "Mg", "Ca"]
	marker=[:square, :circle, :diamond, :utriangle]

	# Store solution
	Ci = zeros(length(sol[t]), N)
	Ci[:,1] = -log10.(prob.ps[Cref] * sol[C(t,x)[1]][:,end])
	[Ci[:,i] = prob.ps[Cref] * sol[C(t,x)[i]][:,end]/Cₓ₀[c[i]] for i in 2:5]

	# pH
	p1 = plot(xlims=xlims, ylims=ylims_pH, ylabel="pH", yticks=yticks_pH, 
			  legendposition=:topleft, xticks=xticks, minorticks=2)
	gridlines = 0:2:pHmax
	for g in gridlines
	    hline!([g], color=:gray, alpha=0.1, linestyle=:solid, lab=:none)
	end
	# model
	plot!(BV[10:end], Ci[10:end,1]; ls=:solid, c=:grey, lw=2.5, lab=:none)
	plot!(twinx(), BV[12:end], diff(Ci[11:end,1]) ./ diff(BV[11:end]); 
		  ls=:dot, c=:grey, lw=2.5, lab=:none, grid=:on, ylabel="ΔpH/Δt",
		  xlims=xlims, ylims=ylims_pH./2, yticks=yticks_pH./2, 
		  legendposition=:bottomright)
	# data
	plot!(data_BV, data_y[:,1]; lab="pH", c=:grey, marker=:circle, lw=1, grid=:off)
	plot!([], []; lab="ΔpH/Δt", c=:black, ls=:dash, lw=1, marker=:xcross, 
		  markerstrokewidth=2)
	plot!(twinx(), data_BV[2:end], diff(data_y[:,1]) ./ diff(data_BV); 
		lab=:none, c=:black, ls=:dash, lw=1, marker=:xcross, markerstrokewidth=2,
		xlims=(0, Vmax), ylims=ylims_pH./2, yticks=yticks_pH./2, ylabel="ΔpH/Δt", grid=:off, xticks=xticks, minorticks=1)

	# Breakthrough
	p2 = plot(xlims=xlims, ylims=ylims, yticks=0:1:Cmax, grid=:off,
			  xlabel=xlabel, xticks=xticks, minorticks=2,
			  ylabel="Normalized concentration", legendposition=:topleft)
	plot!(twinx(); xlims=xlims, ylims=(-15, Cmax*Cₓ₀mass[:Na]), yticks=0:100:100Cmax,
		  ylabel="Concentration [mg/L]", grid=:off, minorticks=2)
	gridlines = 0:1:Cmax
	for g in gridlines
	    hline!([g], color=:gray, alpha=0.1, linestyle=:solid, lab=:none)
	end
	# equilibrium
	hline!(p2, [1.]; lab=:none, c=:black, ls=:dot, lw=0.5)
	# model
	for i in 2:5
		plot!(BV, Ci[:, i]; lab=:none, lw=2.5, ls=:solid, c=cp[i-1])
	end
	# data
	for i in 2:5
		plot!(data_BV, data_y[:,i]./Cₓ₀mass[c[i]]; 
			  c=cp[i-1], lab=label[i-1], marker=marker[i-1])
	end
	# gray out last data points for Mg
	r = gray_Mg-1
	for i in 4
		plot!(data_BV[end-r:end], data_y[end-r:end,i]/Cₓ₀mass[c[i]]; 
			  lab="", marker=marker[i-1], c=:lightgrey, markerstrokecolor=cp[i-1], 
			  linecolor=cp[i-1])
	end
	# gray out last data points for Ca
	r = gray_Ca-1
	for i in 5
		plot!(data_BV[end-r:end], data_y[end-r:end,i]/Cₓ₀mass[c[i]]; 
			  lab="", marker=marker[i-1], c=:lightgrey, markerstrokecolor=cp[i-1], 
			  linecolor=cp[i-1])
	end

	p = plot(p1, p2, layout=grid(2, 1, heights=[0.3, 0.7]), 
			 dpi=300, size=(680, 400), suptitle=title)

	# breakthrough detection
	idx = findfirst(x -> x >= BT, Ci[10:end,2])+10-1
	vline!(p, [BV[idx]]; lab="", ls=:dash, c=:black, lw=0.75, subplot = 4)
	BTlabel = L"$t_{%$(round(Int,100BT))\%} = %$(round(Int,BV[idx]))\,%$units$"
	annotate!(p, 0.95*BV[idx], 0.95*Cmax, (BTlabel, 9, :right); subplot = 4)

	savefig(filename*".pdf")
	plot!()
end

# ╔═╡ c07242d9-5015-4cce-a2d7-24d34200b6d1
md"""
## Stopping criteria
"""

# ╔═╡ 86754b3c-a9bf-4138-a85a-aae4a3e7c78d
md"""
### 1. Fixed-time operation
"""

# ╔═╡ 220c7e50-022e-45fc-a79f-d44596e89b52
function fixed_load(; prob=prob42S)

	update!(pars, prob,	Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
			qmax => p0[1], K => vcat(1.0, p0[2:N]), 
	)
	@time sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

	return sol
end

# ╔═╡ 1464fa94-e1c6-4626-89fe-e28d10155d16
md"""
Treated volume: $(@bindname BVstop confirm(Slider(30:1:40, default=30, show_value=true))) Current: $BVstop L
"""

# ╔═╡ e10178ba-d810-408c-a65c-0d206280e541
md"""
### 2. pH control
"""

# ╔═╡ 8f74a1f1-7fa3-4672-ae77-51c904fe2cc2
md"""
For a noticeable change in pH of 1 decimal digit, the Na breakthrough is already significant: 1.4 g Na at 42.2 L, and a breakthrough concentration of 95 mg/L, close to the initial concentration.

A resolution of 0.05 pH units gets these figures down to 54 mg/L and 777 mg at a volume of 41.1 L.
"""

# ╔═╡ cc953d89-ebf3-4812-a9b0-f9e73e40f611
md"""
pH change: $(@bind pHchange Select([0.01, 0.05, 0.1], default=0.05)) units
"""

# ╔═╡ 8b6f1c83-cc1d-466f-94e0-9f04ba0e7cbb
# let
# 	prob = prob42S; 

# 	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
# 	BVmax = 70
# 	Nt = length(BV)

# 	# pH
# 	pHmax = 8
# 	idxmin = 3
# 	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
# 	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
	
# 	p1 = plot(xlims=(0, BVmax), xlabel="", ylims=(0, pHmax), ylabel="pH")
# 	plot!(BV[idxmin:end], pH[idxmin:end]; 
# 		  lab="", c=col[1], lw=3)

# 	# initial pH
# 	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
# 	prob.ps[qmax]*(1-prob.ps[ε])*prob.ps[V]*1e3
# 	hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
# 	annotate!(2, 1.5*pH0, ("initial pH: $(round(pH0;digits=2))", 9, :left))

# 	# first detection of pH change
# 	# pHchange = (round(10*(pH0 - floor(pH0));digits=1) + 1.0)/10 + 1.0
# 	hline!([pH0 + pHchange]; lab="", c=cp[N+2], ls=:dashdot, lw=3)
# 	idxpH = findfirst(x->x>=pH0 + pHchange, pH[idxmin:end]) + idxmin - 1
# 	vline!([BV[idxpH]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	annotate!(0.95*BV[idxpH], 0.95*pHmax, 
# 			  ("pH $(round(pH[idxpH];digits=2)) ($pHchange pH units) detected at $(round(BV[idxpH];digits=1)) L", 9, :right))

# 	# Breakthrough
# 	Cmax = 500
# 	p2 = plot(xlims=(0, BVmax), xlabel="Volume [L]", ylims=(0, Cmax),
# 			  ylabel="Breakthrough concentration [mg/L]")
# 	plot!(BV, Na; lab=lab[2], col=col[2], lw=3)

# 	hline!([Cₓ₀[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=1, ls=:dot)

# 	# breakthrough detection
# 	NaBT = Na[idxpH]
# 	NaBTint = integrate(BV[1:idxpH], Na[1:idxpH])

# 	vline!([BV[idxpH]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	hline!([NaBT]; lab="", c=cp[N+2], lw=3, ls=:dashdot)
# 	annotate!(0.95*BV[idxpH], 2.5*NaBT, 
# 			  ("Na breakthrough at $(round(BV[idxpH];digits=1)) L: $(round(Int, NaBTint)) mg ($(round(NaBT;digits=1)) mg/L)", 9, :right))
	
# 	figure = plot(p1, p2, layout=grid(2, 1, heights=[0.4, 0.6]), 
# 				  dpi = 300, size=(680, 600), 
# 				  suptitle="pH change breakthrough stop control scenario")
# 	figure |> aside
# end

# ╔═╡ e8c5ecfd-78a8-451a-94da-d26a12dbbe3a
# let
# 	solx = prob42S.problem_type.discretespace.axies[x]
# 	plot(sol[x], sol[C(t,x)[2]][501, :]; lw=2)
# 	scatter!(solx, sol(500, solx, dv=C(t,x)[2]))
# 	scatter!([1], [sol(500, 1, dv=C(t,x)[2])])
# 	scatter!([sol[x][end]], [sol[C(t,x)[2]][500, end]])
# end

# ╔═╡ b14e753a-7c41-48d5-88c3-5f2dac364e47
md"""
### 3. Na level
"""

# ╔═╡ 60125c5d-4997-43fd-a48e-f06373ec54b0
md"""
For very low thresholds, the difference in pH measured at breakthrough is undetectable.

Depending on the resolution and update time of the Na analyzer, the detection can be swift.

The lowest range of detection is 10$^{-6}$ M or ~0.02 ppm Na (23 ppb). The major issue of ion selective electrodes is the interference with other ions such as H+, K+ and Li+. If Na is the first ion to break through the column, and in absence of Li ions, sodium analyzers are an alternative to pH sensors.
"""

# ╔═╡ f08b716a-cbbe-4c1d-b9b0-03e415a9d1ed
md""" Threshold for Na breakthrough: 
$(@bindname Na_threshold Select([0.02 => "20 ppb", 0.1 => "100 ppb", 1 => "1 ppm", 10 => "10 ppm", 100 => "100 ppm"]))
"""

# ╔═╡ 2bf1ffa1-d687-4cca-8072-98908b0e4a4c
let
	prob = prob42S; 

	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
	BVmax = 70
	Nt = length(BV)

	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref] # [eq/L]

	# breakthrough detection in ppb!
	maxNa = Na_threshold

	idxmin = 3
	idxNa = findfirst(x -> x >= maxNa, Na[idxmin:Nt]*Mi[:Na]*1e3) + idxmin - 1
	NaBT = Na[idxNa]*Mi[:Na]*1e3
	NaBTint = integrate(BV[1:idxNa], Na[1:idxNa]*Mi[:Na]*1e3)

	# pH
	pHmax = 8
	
	p1 = plot(xlims=(0, BVmax), xlabel="", ylims=(0, pHmax), ylabel="pH")
	plot!(BV[idxmin:end], pH[idxmin:end]; 
		  lab="", c=col[1], lw=3)

	# initial pH
	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
	prob.ps[qmax]*(1-prob.ps[ε])*prob.ps[V]*1e3
	hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
	annotate!(2, 1.5*pH0, ("initial pH: $(round(pH0;digits=1))", 9, :left))

	vline!([BV[idxNa]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([pH[idxNa]]; lab="", c=cp[N+2], lw=3, ls=:dashdot)
	scatter!([BV[idxNa]], [pH[idxNa]]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	annotate!([0.95BV[idxNa]], [0.95*pHmax], ("pH detected at $(round(BV[idxNa];digits=1)) L: $(round(pH[idxNa];sigdigits=4))", 9, :right))
	

	# Breakthrough
	Cmax = 500
	p2 = plot(xlims=(0, BVmax), xlabel="Volume [L]", ylims=(0, Cmax),
			  ylabel="Breakthrough concentration [mg/L]")
	plot!(BV, [prob.ps[Cref]*Mi[c[i]]/Zi[c[i]]*1e3*
		  sol(:, 1, dv=(C(t,x)[i])) for i in 2:N];
		  lab=permutedims(lab[2:N]), col=permutedims(col[2:N]), lw=3)

	hline!([Cₓ₀[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=1, ls=:dot)
	
	vline!([BV[idxNa]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([NaBT]; lab="", c=cp[N+2], lw=3, ls=:dashdot)
	scatter!([BV[idxNa]], [Na[idxNa]*Mi[:Na]*1e3]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

	annotate!(0.95*BV[idxNa], 25 + NaBT, 
			  ("BT at $(round(BV[idxNa];digits=1)) L: $(round(NaBTint;digits=2)) mg ($(round(NaBT;digits=2)) mg/L)", 9, :right))

	# Resin capacity
	p3 = plot()

	solx = prob.problem_type.discretespace.axies[x]
	plot!(BV, 
		  [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]* 100 ;
		  lab="", c=:grey, lw=3, grid=:off, title="Resin capacity", 
		  xlabel="Volume [L]", xlims=(0, BVmax),
		  ylabel="Resin capacity, %", ylims=(-0.5, 101), )
	capacity = 1 .- [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]
	plot!(BV, 100capacity; lab="", c=:darkred, lw=3, ls=:dot)
	vline!([BV[idxNa]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	annotate!(0.98*BV[idxNa], 95, 
			  ("saturation at breakthrough: $(round(100capacity[idxNa];digits=1))%", 9, :right))
	
	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				  dpi = 300, size=(680, 800), 
				  suptitle="pH change breakthrough stop control scenario")
	figure |> x -> aside(x; v_offset=0)
end

# ╔═╡ d37bdba9-ac5b-4023-8aa2-8ca70300549d
md"""
### 4. pH slope
"""

# ╔═╡ 584c3a98-9bd7-4d46-9b6c-a098694e835d
md"""
Assuming a resolution of 0.01 pH units and a sampling frequency of 1 Hz (1 sample per second), electrochemical probes used in industry can detect pH changes at a rate of 0.01 pH units s$^{-1}$.

For a measurement accuracy of ±0.05 pH units, and for an average response time of 3 s (~0.3 Hz), a more realistic lower bound for pH change could be ~0.015 pH units s$^{-1}$.

Based on these assumptions and a slope of 0.01, the breakthrough detection reaches 23 mg/L Na, or equivalently 23 ppm or 0.001 M. This is much higher than the best sodium analyzers.

Another important issue is the daily calibration of the pH probes to avoid sensor drift.
"""

# ╔═╡ 3fb3da00-ee78-4339-a458-b6f947f6900f
let
	# pH with Gaussian noise
	
	α = 0.05
	Random.seed!(1234);
	μ = 1.7
	σ = 0.1/2
	y = μ .+ σ .* randn(10000) 			# pH: 1.7 ± 0.1
	
	p975 = quantile(y, 1-α/2)
	p025 = quantile(y, α/2)

	count(<(p025), y) |> println

	figure = plot(dpi=300, ylims=(0, 0.25), legend=false) 
	histogram!(y; bins=range(1.5, 1.9, 20+1), normalize=:probability)
	vline!([mean(y)]; lw=3, ls=:dash)
	vline!([p975]; lw=2, ls=:dot, c=:black)
	vline!([p025]; lw=2, ls=:dot, c=:black)
	annotate!(1.005mean(y), 0.2, ("mean = $(round(mean(y);digits=4))", :left, 10))
	annotate!(1.005p975, 0.1, ("97.5% = $(round(p975;digits=4))", :left, 10))
	annotate!(0.995p025, 0.1, ("2.5% = $(round(p025;digits=4))", :right, 10))
	figure |> aside
end

# ╔═╡ 82c98042-7e70-4243-b8b4-30143322a077
function causal_moving_average(y, w)
    n = length(y)
    y_filt = similar(y, Float64)

    halfw = div(w, 2)

    for i in eachindex(y)
        i1 = max(1, i - halfw)
        i2 = min(n, i + halfw)

        y_filt[i] = mean(y[i1:i2])
    end

    return y_filt
end

# ╔═╡ b2ca0f67-05b7-4a94-a679-5a649d002fc7
function filter(y, w1)
    n = length(y)
    y_filt = similar(y, Float64)

    for i in eachindex(y)
        i1 = max(1, i - w1 + 1)
        i2 = i

        y_filt[i] = mean(y[i1:i2])
    end

    return y_filt
end

# ╔═╡ d4b5e8f6-f288-47de-b05f-377181ce3326
let
	prob = prob42S; 

	BV = sol[t]*prob.ps[τ]/60#*prob.ps[Q]*1e3 # [L]  [min]
	Nt = length(BV)
	tmax = 500
	ΔpH = pHchange
	σ = pHchange/2
	w1 = 3
	Random.seed!(1234);
	
	pHmax = 8; 	dpHmin = -0.0008; 	dpHmax = 0.2001
	idxmin = 3
	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	dpH = diff(pH)./diff(BV)
		
	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3 # [mg/L]
	K = sol(:, 1, dv=C(t,x)[3])*prob.ps[Cref]*Mi[:K]*1e3 # [mg/L]
	Mg = sol(:, 1, dv=C(t,x)[4])*prob.ps[Cref]*Mi[:Mg]/Zi[:Mg]*1e3 # [mg/L]
	Naint = [integrate(BV[1:j], Na[1:j]) for j in 2:Nt] # [mg]

	# pH
	p1 = plot(xlims=(0, tmax), xlabel="Time [min]", ylims=(0, pHmax), ylabel="pH")
	plot!(BV[idxmin:end], pH[idxmin:end]; 
		  lab="pH", c=col[1], lw=3, legendposition=:topright)	
	# plot!([], [], lab="ΔpH/Δt", c=col[1], ls=:dashdot, lw=3, grid=:off)

	# initial pH
	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N])) # pH = 1.7
	# hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
	annotate!(10, 0.5*pH0, ("initial pH: $(round(pH0;sigdigits=3))", 9, :left))

	# noisy measurement
	noisy_pH = pH[idxmin:end] .+ σ.*randn(length(pH[idxmin:end]))
	
	# first detection of pH change
	idxpH = findfirst(x -> x >= pH0+ΔpH, filter(noisy_pH, w1)) + idxmin - 1
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2.5)
	hline!([pH0+ΔpH]; lab="", c=cp[N+1], ls=:dot, lw=2.5)
	scatter!([BV[idxpH]], [pH0+ΔpH]; 
			 lab="", m=:xcross, markersize=5, c=:black, markerstrokewidth=2)
	annotate!(0.955*BV[idxpH], 0.5*pH0, ("pH: $(round(pH0+ΔpH;sigdigits=3)) at $(round(Int, BV[idxpH])) min ($(round(BV[idxpH]*60*prob.ps[Q]*1e3;digits=1)) L)", 9, :right))	

	# inset
	plot!(BV, pH; inset = (1, bbox(0.075, 0.05, 0.5+ΔpH, 0.5)), subplot=2,
			 xlims=(280, 380), ylims=(1.6, 1.9), yticks=1:0.1:2, ls=:solid,
			 lab="", c=col[1], lw=3, box=:box, grid=:off,
		 	 legendposition=:bottomright, legendfontsize=8)
	# noisy measurement
	plot!(BV[idxmin:end], noisy_pH; 
		  subplot=2, c=col[1], ls=:dot, lw=2, lab="pH probe")
	# moving average filter
	plot!(BV[idxmin:end], filter(noisy_pH, w1); 
		  subplot=2, lw=2, c=:darkred, lab="average")
	hline!([pH0+ΔpH]; lab="", c=cp[N+1], ls=:dot, lw=2.5, subplot=2)
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2.5, subplot=2)
	# scatter!([BV[idxpH]], [pH0+ΔpH]; subplot=2,
	# 		 lab="", m=:xcross, markersize=9, c=:black, markerstrokewidth=2)

	# pH slope
	# plot!(twinx(), BV[idxmin:end-1], dpH[idxmin:end];
	# 	  lab="", c=col[1], ls=:dashdot, lw=3, grid=:off,
	# 	  xlims=(0, tmax), ylims=(dpHmin, dpHmax), ylabel="ΔpH/Δt")

	# Breakthrough
	Cmax = 500
	p2 = plot(xlims=(0, tmax), xlabel="Time [min]", ylims=(-0.5, Cmax),
			  ylabel="Concentration [mg/L]", grid=:off, legendposition=:topleft)
	plot!(BV[1:idxpH], Na[1:idxpH]; lab=lab[2], c=cp[1], lw=3)
	plot!(BV[idxpH:end], Na[idxpH:end]; lab="", c=cp[1], lw=3, ls=:dot, alpha=0.9)
	plot!(BV[1:idxpH], Na[1:idxpH]; fillrange=0, c=cp[1], alpha=0.3, lab="")
	plot!(BV[1:end], K[1:end]; lab=lab[3], c=cp[2], lw=3, ls=:dot, alpha=0.9)

	# breakthrough detection
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2.5)
	hline!([Na[idxpH]]; lab="", c=cp[N+1], lw=2.5, ls=:dot)
	scatter!([BV[idxpH]], [Na[idxpH]]; 
			 lab="", m=:xcross, markersize=5, c=:black, markerstrokewidth=2)
	annotate!(0.95*BV[idxpH], 2*Na[idxpH],
			  ("Breakthrough at $(round(Int, BV[idxpH])) min: $(round(Na[idxpH];digits=1)) mg/L (total: $(round(Int, Naint[idxpH-1])) mg)", 9, :right))
	
	figure = plot(p1, p2, layout=grid(2, 1, heights=[0.5, 0.5]), 
				  dpi = 300, size=(680, 600), 
				  suptitle="Breakthrough stop with pH change detection")
	figure |> x-> aside(x; v_offset=0)
end

# ╔═╡ d8ed6fd2-0ae6-4b8a-889b-87ce71580edb
let
	# derivative of pH: ΔpH/Δt
	
	α = 0.05
	Random.seed!(1234);
	μ = 1.7
	σ = 0.1
	x = filter(μ .+ σ .* randn(10000), 20) 			# pH: 1.7 ± 0.1
	y = diff(x)./diff(1:1:length(x))
	
	p975 = quantile(y, 1-α/2)
	p025 = quantile(y, α/2)

	count(<(p025), y) |> println

	figure = plot(dpi=300, ylims=(0, 0.25), legend=false) 
	histogram!(y; bins=range(-0.03, 0.03, 20+1), normalize=:probability)
	vline!([mean(y)]; lw=3, ls=:dash)
	vline!([p975]; lw=2, ls=:dot, c=:black)
	vline!([p025]; lw=2, ls=:dot, c=:black)
	annotate!(1.005mean(y), 0.2, ("mean = $(round(mean(y);digits=4))", :left, 10))
	annotate!(1.005p975, 0.1, ("97.5% = $(round(p975;digits=4))", :left, 10))
	annotate!(0.995p025, 0.1, ("2.5% = $(round(p025;digits=4))", :right, 10))
	figure |> aside
end

# ╔═╡ 97b01366-9ff8-40c3-b6dd-fee9510eb58f
let
	Random.seed!(1234);
	y = randn(10)
	w1 = 2
	
	plot()
	plot!(y; m=:circle, lab="measured", lw=2, ls=:dot)
	# plot!(causal_moving_average(y, w1), lab="mov. average", lw=2)
	plot!(filter(y, w1), lab="filter, w=$w1", lw=2)
	plot!(filter(y, 10), lab="filter, w=10", lw=2)
end

# ╔═╡ 3d573f46-8def-4a6d-83b7-aa577996b079
@bind dpHslope confirm(Slider(0.001:0.0005:0.02, default=0.01, show_value=true))

# ╔═╡ aaeeca5e-08f6-479d-a300-02c4132752df
let
	prob = prob42S; 

	BV = sol[t]*prob.ps[τ]/60#*prob.ps[Q]*1e3 # [L]  [min]
	Nt = length(BV)
	tmax = 500
	ΔpH = 0.1
	σ = 0.1/2
	w1 = 5
	w2 = 15
	Random.seed!(1233);
	SF = 1.2
	
	pHmax = 8
	imin = 3
	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	dpH = diff(pH)./diff(BV)
		
	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3 # [mg/L]
	K = sol(:, 1, dv=C(t,x)[3])*prob.ps[Cref]*Mi[:K]*1e3 # [mg/L]
	Mg = sol(:, 1, dv=C(t,x)[4])*prob.ps[Cref]*Mi[:Mg]/Zi[:Mg]*1e3 # [mg/L]
	Naint = [integrate(BV[1:j], Na[1:j]) for j in 2:Nt] # [mg]

	# (a) pH
	p1 = plot(xlims=(0, tmax), xlabel="Time [min]", legendposition=:topright,
			  ylims=(0, pHmax), ylabel="pH", grid=:off, title="(a) pH control: ΔpH = $ΔpH, moving average ($w1 min)")
	plot!(BV[imin:end], pH[imin:end]; lab="pH", c=col[1], lw=3)	

	# initial pH
	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N])) # pH = 1.7
	annotate!(10, 0.5*pH0, ("initial pH: $(round(pH0;digits=2))", 9, :left))

	# noisy pH measurement (upsampled)
	# pHint = LinearInterpolation(BV, pH)
	# BVup = 0:(3/60):BV[end]
	# pHup = pHint(BVup)
	# noisy_pH_up = pHup .+ σ.*randn(length(pHup))

	# noisy pH measurement
	noisy_pH = pH[imin:end] .+ σ.*randn(length(pH[imin:end]))
	noisy_dpH = filter(diff(noisy_pH)./diff(BV[imin:end]), w2)
	plot!(BV[imin:end], noisy_pH; c=:darkblue, ls=:dot, lw=1, lab="")
	
	# first detection of pH change
	idxpH = findfirst(x -> x >= pH0+ΔpH, filter(noisy_pH, w1)) + imin - 1
	
	hline!([filter(noisy_pH, w1)[idxpH-imin+1]]; 
		   lab="", c=cp[N+1], ls=:dot, lw=3)
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=3)
	scatter!([BV[idxpH]], [pH0+ΔpH]; 
			 lab="", m=:xcross, markersize=6, c=:black, markerstrokewidth=2)
	annotate!(0.955*BV[idxpH], 0.5*pH0, ("pH = $(round(filter(noisy_pH, w1)[idxpH-imin+1];digits=2)) at $(round(Int,BV[idxpH])) min ($(round(BV[idxpH]*60*prob.ps[Q]*1e3;digits=1)) L)", 9, :right))	

	# inset
	plot!(BV[imin:end], pH[imin:end]; 
		  inset = (1, bbox(0.07, 0.02, 0.5+ΔpH, 0.5)), subplot=2, 
		  xlims=(300, 380), ylims=(1.6, 1.9), yticks=1.6:0.1:1.9, 
		  lab="", c=col[1], lw=3, box=:box, grid=:off, ylabel="",
		  legendposition=:bottomright, legendfontsize=8)
	
	# noisy measurement
	plot!(BV[imin:end], noisy_pH; 
		  subplot=2, c=:blue, ls=:dot, lw=1.5, lab="pH probe")
	# moving average filter
	plot!(BV[imin:end], filter(noisy_pH, w1); 
		  subplot=2, lw=2.5, c=:darkred, lab="average")

	hline!([pH0+ΔpH]; subplot=2, lab="", ls=:dot, c=:black, lw=0.75)
	hline!([filter(noisy_pH, w1)[idxpH-imin+1]]; 
		   lab="", c=cp[N+1], ls=:dot, lw=2.5, subplot=2)
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2.5, subplot=2)

	# up-sample
	# plot!(BVup, noisy_pH_up; subplot=2, lab="", c=:orange, lw=2)
	# plot!(BVup, filter(noisy_pH_up, 3w); 
	# 	  subplot=2, lw=2, c=:darkred, lab="0.3 Hz")

	
	# (b) pH slope
	dpHmin = -0.02; 	dpHmax = 0.06
	p2 = plot(yticks=dpHmin:(dpHmax/3):dpHmax, xlabel="Time [min]",
			  xlims=(0, tmax), ylims=(dpHmin, dpHmax), ylabel="ΔpH/Δt", 
			  grid=:off, minorgrid=:off, minorticks=2, legend=:bottomright,
			  title="(b) pH slope control: ΔpH/Δt = $dpHslope min⁻¹, $w2 min,  $(round(Int, 100*(SF-1)))% SF")
	
	plot!(BV[imin+1:end], dpH[imin:end]; lab="ΔpH/Δt", c=:green, lw=3)
	# hline!([0.01]; lab="", ls=:dot, lw=0.75, c=:black)
	plot!(BV[imin+1:end], noisy_dpH; c=:darkblue, ls=:dot, lw=1, lab="")

	# pH slope detection
	idxpH2 = findfirst(x -> x >= SF*dpHslope, noisy_dpH[3imin:end]) + 3imin + 2

	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2)
	vline!([BV[idxpH2]]; c=cp[N+2], lw=3, lab="", ls=:dot)
	hline!([dpH[idxpH2]]; c=cp[N+2], lw=3, lab="", ls=:dot)
	scatter!([BV[idxpH2]], [dpH[idxpH2]]; 
			 lab="", m=:star, markersize=9, c=:black, markerstrokewidth=2)
	annotate!(0.95*BV[idxpH2], -0.015, ("ΔpH/Δt = $(round(noisy_dpH[idxpH2-imin];sigdigits=2)) at $(round(Int,BV[idxpH2])) min ($(round(BV[idxpH2]*60*prob.ps[Q]*1e3;digits=1)) L)", 9, :right))	
	
	# inset
	plot!(BV[2:end], dpH; inset = (1, bbox(0.09, 0.02, 0.55, 0.5)), 
		  subplot=2, xlims=(300, 380), ylims=(-0.01, 0.02),yticks=-0.01:0.01:0.02,
		  lab="", c=:green, lw=3, box=:box, grid=:off, ylabel="",
		  legendposition=:bottomright, legendfontsize=8)
	
	hline!([0.01]; subplot=2, lab="", ls=:dot, lw=0.75, c=:black)
	plot!(BV[imin+1:end], noisy_dpH; 
		  subplot=2, c=:darkblue, ls=:solid, lw=2, lab="average")
	vline!([BV[idxpH2]]; c=cp[N+2], lw=2.5, lab="", ls=:dot, subplot=2)
	hline!([noisy_dpH[idxpH2-imin]]; 
		   c=cp[N+2], lw=2.5, lab="", ls=:dot, subplot=2)
	
	# plot!(BVup[1:end-1], noisy_dpH_up; 
		  # subplot=2, c=:orange, ls=:solid, lw=2, lab="pH probe")
	# noisy_dpH_up = diff(moving_average(noisy_pH_up, 2w))./diff(BVup)
	# moving average filter
	# plot!(BV[idxmin:end], moving_average(noisy_pH, w); 
	# 	  subplot=2, lw=2, c=:darkred, lab="mov. average"

	
	# (c) Breakthrough
	Cmax = 500
	p3 = plot(xlims=(0, tmax), xlabel="Time [min]", ylims=(-0.5, Cmax),
			  ylabel="Concentration [mg/L]", grid=:off, legendposition=:topleft)
	
	plot!(BV[1:idxpH], Na[1:idxpH]; lab=lab[2], c=col[2], lw=3)
	plot!(BV[idxpH:end], K[idxpH:end]; lab=lab[3], c=col[3], lw=3, ls=:dash, alpha=0.9)
	plot!(BV[idxpH:end], Na[idxpH:end]; lab="", c=col[2], lw=3, alpha=0.7)
	
	hline!([Cₓ₀[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=0.75, ls=:dot)
	
	plot!(BV[1:idxpH], Na[1:idxpH]; fillrange=0, c=col[2], alpha=0.2, lab="")
	plot!(BV[1:idxpH2], Na[1:idxpH2]; fillrange=0, c=col[2], alpha=0.5, lab="")
	
	# plot!(BV[1:end], Mg[1:end]; lab=lab[4], c=col[4], lw=3, ls=:dot, alpha=0.9)
	

	# breakthrough detection
	vline!([BV[idxpH]]; lab="", c=cp[N+1], ls=:dot, lw=2)
	# hline!([Na[idxpH]]; lab="", c=cp[N+1], lw=2, ls=:dot)
	annotate!(0.95*BV[idxpH2], 1.5*prob.ps[Cref]*Mi[:Na]*1e3,
			  ("Na breakthrough at $(round(Int,BV[idxpH2])) min: $(round(Na[idxpH2];digits=1)) mg/L ($(round(Int, Naint[idxpH2-1])) mg)", 9, :right))

	vline!([BV[idxpH2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	# hline!([Na[idxpH2]]; lab="", c=cp[N+2], lw=3, ls=:dot)
	# scatter!([BV[idxpH2]], [Na[idxpH2]];
	# 		 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	# annotate!(0.95*BV[idxpH2], 0.6Cmax, 
	# 		  ("BT at $(round(BV[idxpH2];digits=1)) L: $(round(Int, Na[idxpH2])) mg ($(round(Naint[idxpH2-1];digits=1)) mg/L)", 9, :right))
	
	
	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.4, 0.4, 0.2]), 
				  dpi = 300, size=(740, 900), #xlabel="Time [min]",
				  suptitle="pH change breakthrough stop control scenario")
	figure |> x-> aside(x; v_offset=0)
end

# ╔═╡ c5dabff0-fe48-42c9-b292-d3880fd2df1f
md"""
## Consecutive cycles
"""

# ╔═╡ 66466e20-66d0-47dc-81fc-f2123d9bce09
md"""
Assuming quick rinse or backwash between cycles, selecting different stopping criteria for combined cycles of loading and regeneration will produce various performance results.


"""

# ╔═╡ c8df5818-a945-4f40-a5db-db866b7397f5
Cₓ₀reg = Dict([i => 0. for i in c]); Cₓ₀reg[:H] = 1.4; Cₓ₀reg

# ╔═╡ 39f3f4a9-5513-4973-80b1-f07c152f7d28
C₀reg = Dict([c[j] => zero_eq for j in 1:N]); C₀reg[:H]=1e-7; C₀reg

# ╔═╡ eaa12a3d-23e6-41d6-9e9c-913b9b2c0b82
bcsreg = let
	bcs = Array{Symbolics.Equation}(undef, 2N + 2N) # 2 ICs + 2 BCs
		for (i, j) in zip(1:N, c)
			# Initial conditions
			bcs[i] 		= C(0.0, x)[i] ~ C₀reg[j]
			bcs[i + N] 	= q(0.0, x)[i] ~ q₀[j]
			# Boundary conditions - Danckwerts for fluid, Neumann for the solid
			bcs[i + 2N] = C(t, 0.0)[i] ~ Cx0[i]/Cref
			bcs[i + 3N] = Dx(q(t, 0.0)[i]) ~ qₓ₀[j]
		end
	bcs
end

# ╔═╡ 7c9d4d0b-e4ef-43e8-a36c-44ff8b90d0d4
@named pdesysregS = PDESystem(eqs=eqs_same, bcs=bcsreg, domain=domains, ivs=[t, x], dvs=dvs, ps=ps, defaults=pars);

# ╔═╡ 82ee92cf-f38e-4c02-a743-51b8a014cfd8
@time reg42S = discretize_problem(pdesysregS, dx=0.025, grid_align=edge_align);

# ╔═╡ 7a7be574-484e-4e36-bf23-ce7c83b0ab1b
end1 = findfirst(x->x>=BVstop, sol[t]*prob42S.ps[τ]*prob42S.ps[Q]*1e3)

# ╔═╡ 4f5e7e54-ccf9-4067-8af4-fce16fa592f0
let
	prob = prob42S
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
	solx = prob.problem_type.discretespace.axies[x]
	tstop = findfirst(x -> x >= BVstop, BV) - 1 # [idx]; 	# tstop |> println
	BVmax = 70; 	Nt = length(BV)

	# pH
	pHmax = 8
	idxmin = 3
	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
	RC = sol(:, solx, dv=q(t,x)[1])*100
	
	p1 = plot(xlims=(0, BVmax), xlabel="", grid=:off,
			  ylims=(0, pHmax), ylabel="pH", yticks=0:2:pHmax)
	plot!(BV[idxmin:end], pH[idxmin:end]; 
		  lab="", c=col[1], lw=3)

	# initial pH
	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
	prob.ps[qmax]*(1-prob.ps[ε])*prob.ps[V]*1e3
	hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
	annotate!(2, 2*pH0, ("initial pH: $(round(pH0;digits=2))", 9, :left))

	# first detection of pH change
	# pHchange = (round(10*(pH0 - floor(pH0));digits=1) + 1.0)/10 + 1.0
	hline!([pH[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	vline!([BV[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	annotate!(0.95*BV[tstop], 0.95*pHmax, 
			  ("pH $(round(pH[tstop];digits=2)) at $(round(Int, BV[tstop])) L", 9, :right))
	scatter!([BV[tstop]], [pH[tstop]]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

	# Breakthrough
	Cmax = 500
	p2 = plot(xlims=(0, BVmax), ylims=(0, Cmax), grid=:off,
			  ylabel="Concentration [mg/L]")
	plot!(BV, Na; lab=lab[2], col=col[2], lw=3)

	hline!([Cₓ₀[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=1, ls=:dot)

	# breakthrough detection
	NaBT = Na[tstop]
	NaBTint = integrate(BV[1:tstop], Na[1:tstop])

	vline!([BV[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([NaBT]; lab="", c=cp[N+2], lw=3, ls=:dashdot)
	annotate!(0.95*BV[tstop], 0.5Cmax, 
			  ("BT at $(round(Int,BV[tstop])) L: $(round(NaBTint;digits=2)) mg ($(round(NaBT;digits=2)) mg/L)", 9, :right))
	scatter!([BV[tstop]], [Na[tstop]]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

	# Resin capacity
	cap = zeros(Nt, N)
	[cap[:,i] = integrate(solx, sol(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	
	p3 = plot(xlims=(0, BVmax), xlabel="Treated volume [L]", grid=:off,
			  ylims=(0, 100), ylabel="Resin capacity, %")

	areaplot!(BV, cap; lab=permutedims(lab), c=permutedims(col), alpha=0.5, 
			  legendposition=:bottomleft)

	plot!(BV, cap[:,1]; lab="", lw=3, c=col[1])
	vline!([BV[end1]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([cap[end1]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	scatter!([BV[end1]], [cap[end1]];
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)	
	annotate!(1.05*BV[tstop], 2.5cap[tstop], 
			  ("capacity at $(round(Int,BV[tstop])) L: $(round(cap[tstop];digits=2))%", 9, :left))
	# annotate!((BV[end1]<=3 ? 0.2 : -0.2)+BVreg[end2], (BVreg[end2]<=3 ? -10 : -18)+cap[end2], ("resin capacity at $(round(BVreg[end2];digits=1))L: $(round(cap[end2];digits=1))%", 9, (BVreg[end2]<=3 ? :left : :right )))

	# # Resin capacity
	# cap = [integrate(sol_x, RC[j, :]) for j in 1:Nt]
	# p3 = plot(xlims=(0, BVmax), xlabel="Volume [L]", ylims=(0, 100),
	# 		  ylabel="Resin capacity, %")
	# plot!(BV, cap; lab="", c=:black, lw=3)
	# vline!([BV[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	# hline!([cap[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	# annotate!(1.05*BV[tstop], 2.5cap[tstop], 
	# 		  ("resin capacity at $(round(Int,BV[tstop])) L: $(round(cap[tstop];digits=2))%", 9, :left))
	# scatter!([BV[tstop]], [cap[tstop]]; 
	# 		 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				  dpi = 300, size=(680, 600), 
				  suptitle="Fixed-time breakthrough stop control scenario")
	figure |> x->aside(x; v_offset=600)
	
	
	# # Breakthrough curves
	# Cmax = maximum([sol(0:tstop, 1, dv=C(t,x)[i])*prob.ps[Cref]/Cₓ₀[c[i]] for i in 2:N])
	# for i in 2:N
	# 	plot!(BV, sol(:, 1, dv=C(t,x)[i])*prob.ps[Cref]/Cₓ₀[c[i]]; 
	# 		  lab=lab[i], c=cp[i-1], lw=3, 
	# 		  # ylims=(0, ceil(maximum(10Cmax))/10),
	# 		  ylims=(0, 5)
	# 		 )
	# end

	# hline!([1.]; lab="", c=:black, ls=:dot, lw=0.75)
	# vline!([BVstop]; lab="", c=:black, ls=:dot, lw=2)

	# # Accumulated ions
	# Mac = zeros(N)
	# for i in 2:N
	# 	Mac[i] = integrate(BV[1:tstop], sol(0:sol[t][tstop], 1, dv=C(t,x)[i])*
	# 		prob.ps[Cref]*Mi[c[i]]/Zi[c[i]]*1) # [mg]
	# end
	
	# p2 = plot()
	# bar!(lab[2:N], Mac[2:N]; 
	# 	 lab="", c=col[2:N], title="Total ion breakthrough", 
	# 	 ylims=(0, ceil(maximum(Mac))), ylabel="Accumulated ion mass [g]")
	# for i in 2:N 
	# 	annotate!(lab[i], Mac[i]+ceil(maximum(Mac))/15, 
	# 			  ("$(round(Mac[i];digits=1))", 9))
	# end

	# Breakthough concentrations
	# BTstop = 30:5:40
	# n = length(BTstop); BTmin = minimum(BTstop); BTmax = maximum(BTstop)
	# BTidx = [findfirst(x->x>=i, BV) for i in BTstop] # [idx]
	# BT = [sol(j, 1, dv=C(t,x)[i])*prob.ps[Cref]*Mi[c[i]]*1e3/Zi[c[i]] for j in BTidx, i in 1:N]
	
	# p3 = plot()
	# BTper = 50
	# BTlim = 10ceil(maximum(BT/10))
	# comp = 2:3 
	# m = length(comp) ; mmin = minimum(comp) ; mmax = maximum(comp)
	# bw = (BTmax - BTmin)/n

	# plot!(BV, [sol(j, 1, dv=C(t,x)[i])*prob.ps[Cref]*Mi[c[i]]*1e3/Zi[c[i]] for j in 0:sol[t][end], i in comp]; 
	# 	  lab=permutedims(lab[comp]), c=permutedims(col[comp]), lw=3,
	# 	  title="Effluent concentration", grid=:off, legendposition=:topleft,
	# 	  xlabel="Volume [L]", xlims=(0, BTmax),  #xticks=BTstop,
	# 	  ylims=(0, BTlim), #yticks=0:(min(50, BTlim/5)):BTlim,
	# 	  ylabel="Concentration [mg/L]")
	# hline!([Cₓ₀mass[:Na]/10]; lab="", c=:black, ls=:dot, lw=0.75)
	# BTNa = findfirst(x->x>=0.1, sol(:, 1, dv=C(t,x)[2]))
	# vline!([BV[BTNa]]; lab="", c=:black, ls=:dot, lw=2)
	# annotate!(0.98*BV[BTNa], 0.9BTlim, 
	# 		  ("10% breakthrough: $(round(BV[BTNa]; digits=1)) L", 9, :right))
	
	# groupedbar!(BTstop, BT[:, comp]; barwidths=bw, bar_position=:dodge,
	# 			lab=permutedims(lab[comp]), c=permutedims(col[comp]), title="Effluent concentration", grid=:off,
	# 			xlabel="Volume [L]", xticks=BTstop,
	# 			ylims=(0, BTlim), yticks=0:(min(50, BTlim/5)):BTlim,
	# 			ylabel="Concentration [mg/L]")
	
	# for i in comp 
	# 	annotate!(BTstop .+ bw/m*(i .- mmin .- (mmax-mmin)/2), 
	# 			  BT[:, i] .+ 20, ceil.(Int, BT[:, i]); annotationfontsize=9)
	# end

	# Resin capacity/efficiency
	# p4 = plot()
	# solx = prob42S.problem_type.discretespace.axies[x]
	# plot!(BV, 
	# 	  [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]* 100 ;
	# 	  lab="", c=:grey, lw=3, grid=:off, title="Resin capacity", 
	# 	  xlabel="Volume [L]", xlims=(0, BTmax),
	# 	  ylabel="Resin capacity, %", ylims=(0, 100), )
	# capacity = 1 .- [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]
	# plot!(BV, 100capacity; lab="", c=:darkred, lw=3, ls=:dot)
	# vline!([BV[BTNa]]; lab="", c=:black, ls=:dot, lw=2)
	# annotate!(0.98*BV[BTNa], 95, 
	# 		  ("capacity at 10% breakthrough: $(round(100capacity[BTNa];digits=1))%", 9, :right))

	# l = @layout [a{0.25h}; b; c; d]
	
	# p5 = plot(p1, p2, p3, p4, layout=l, dpi=300, size=(680, 800),
	# 	 suptitle="1. Fixed-time stop breakthrough control scenario")
	# p5 |> x->aside(x; v_offset=600)
end

# ╔═╡ e2cb876a-716c-47ce-944b-cb365ca99348
sol[t][end1-1]*prob42S.ps[τ]*prob42S.ps[Q]*1e3

# ╔═╡ 7d787166-df1d-4ba0-83a5-036665ea061c
BVstop/1e3/prob42S.ps[Q]/prob42S.ps[τ]

# ╔═╡ 788ae965-4a3f-48a3-a048-881531dc4f19
sol[t][end1-1]

# ╔═╡ 335a82ca-8cb2-478d-88c5-f69f7e14fce0
solreg, reg = let
	prob = prob42S; 	reg = reg42S; 	
	
	Cₓ₀reg = Dict([i => zero_eq for i in c]); Cₓ₀reg[:H] = 1.4
	reg = updateIC!(sol, reg, sol[t][end1]; reverse=false)
	
	update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
			qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]),
	)
	@time solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
	solreg, reg
end;

# ╔═╡ 13f5e12c-9d95-4ff0-89d8-cf5b84928a56
let
	reg = reg42S
	# solreg = solreg1
	slope = 0.04 # equivalent to ~5% elution Na
	tmin = 10
	tload, index = pHdetection(solreg; slope=slope, tmin=tmin, Cref=reg.ps[Cref], flip=true)
	BV = solreg[t][tmin:end]
	pH = -log10.(solreg[C(t,x)[1]][tmin:end,end] * reg.ps[Cref])
	dpH = diff(pH)./diff(BV)
	Na = solreg[C(t,x)[2]][tmin:end,end] * reg.ps[Cref] 
	
	p1 = plot(BV, pH; lab="pH", lw=2, c=:grey, ls=:dash, ylims=(-0.5, 2.5))
	# hline!([pH[index]]; lab=:none, ls=:dash, c=:black, lw=1)
	scatter!([tload],[pH[index-tmin+1]]; marker=:cross, c=:red, lab="pH = $(round(pH[index];digits=1))")
	plot!(BV, Na; lab="Na", c=cp[1], lw=1)
	vline!([tload]; lab=:none, ls=:dash, c=:black)
	scatter!([tload],[Na[index-tmin+1]]; marker=:cross, c=:black, 
			 lab="$(round(Na[index]*100;digits=1))%")
	hline!([Na[index-tmin+1]]; c=:blue, ls=:dashdot, lab=:none)
	
	
	p2 = plot(BV[2:end], dpH; lab=:none, ylims=(-0.50, 0.0))
	hline!([-slope]; lab="slope = $slope", ls=:dashdot, c=:red)
	vline!([tload]; lab=:none, ls=:dash, c=:black)
	scatter!([tload],[-slope]; marker=:cross, c=:black, lab=:none) 			# slope
	scatter!([tload],[dpH[index-tmin+1]]; marker=:cross, c=:red, lab=:none) # actual
	plot(p1, p2, layout=(2,1), dpi=300, size=(680, 400), xlims=(0, 50))
end

# ╔═╡ b03cc95b-c458-4700-92d9-d48458577a14
# let
# 	reg = updateIC!(sol, reg42S, sol[t][end1]; reverse=false)
# 	solx = reg42S.problem_type.discretespace.axies[x]; Np = length(solx)
# 	integrate(solx, sol(sol[t][end1-1], solx, dv=q(t,x)[1])) |> println
# 	integrate(solx, solreg(0, solx, dv=q(t,x)[1])) |> println
	
# 	plot(solx, sol(sol[t][end1-1], solx, dv=q(t,x)[1]); lw=2, lab="interpolant")
# 	plot!(sol[x], sol[q(t,x)[1]][end1, :]; lw=2, lab="original solution")
# 	scatter!(solx, sol(sol[t][end1-1], solx, dv=q(t,x)[1]); lab="actual values")
# 	# scatter!(sol[x][end:-1:1], solreg[q(t,x)[1]][1,:]; lab="regeneration")
# 	scatter!(solx, reg.u0[@. 2Np*(N-1)+2*(1:Np)-1]; lab="initial values")
# 	scatter!(solx[end:-1:1], solreg(0, solx, dv=q(t,x)[1]); lab="solved")
# end

# ╔═╡ 3fc5ad89-84b5-40b0-aab9-9ab817a5ffdb
md"""
### 1. Fixed-time cycles
"""

# ╔═╡ a470658e-444f-4cc2-81a8-0c0d7a7e03df
md"""
This scenario will assume different fixed times of operation for loading and regeneration.

- Loading: 30-40 L of treated effluent
- Regeneration: 0.1–4 L of HCl 5% regenerant
"""

# ╔═╡ 3b85787c-c67f-4589-ae0e-59401f7b2751
md"""
$(@bind BVregstop confirm(Slider(0.1:0.1:4, default=1, show_value=true))) Regenerant volume: $BVregstop L
"""

# ╔═╡ e135d615-b494-473d-bb39-7175670683b0
end2 = findfirst(x->x>=BVregstop, solreg[t]*reg42S.ps[τ]*reg42S.ps[Q]*1e3)

# ╔═╡ 71fd718e-7d75-4655-a8d4-fdce5290b557
solreg[t][end2]*reg42S.ps[τ]*reg42S.ps[Q]*1e3

# ╔═╡ 18ff9503-091f-4afa-9dd9-4c8df64051e7
let
	prob = prob42S;  	reg = reg42S; 	
	solx = reg.problem_type.discretespace.axies[x]

	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
	BVreg = solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3 # [L]
	pH = -log10.(solreg(:, 1, dv=C(t,x)[1])*reg.ps[Cref])
	BC = [solreg(:, 1, dv=C(t,x)[i])*reg.ps[Cref]*Mi[c[i]]/Zi[c[i]] for i in 1:N]
	BVmax = 5; 	Nt = length(BVreg); 
	pHmax = 8;  BVmin = 10

	findfirst(x -> x > 0.1e-3, BC[2]) |> println

	# pH
	p1 = plot(xlims=(0, BVmax), ylims=(-0.5, pHmax), ylabel="pH", grid=:off)
	plot!(BVreg[BVmin:end], pH[BVmin:end]; lab="", c=col[1], lw=3)
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)

	# pH regenerant
	hline!([-log10(Cₓ₀reg[:H])]; lab="", c=:black, ls=:dot, lw=1)
	annotate!(0.98BVmax, 2, ("final pH: $(round(-log10(Cₓ₀reg[:H]);digits=1))", 9, :right))

	# Elution curve
	Cmax = 10
	p2 = plot(xlims=(0, BVmax), 
			  ylims=(-0.05, Cmax), yticks=0:2:Cmax, ylabel="Concentration [g/L]")
	plot!(BVreg, [BC[i] for i in 2:N]; 
		  lab=permutedims(lab[2:N]), c=permutedims(col[2:N]), lw=3, grid=:off)
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	annotate!(BVreg[end2] + 0.3, 20, 
			  ("BT Na: $(round(BC[2][end2,end];sigdigits=1)) g/L", 9, :left))

	# Resin capacity
	cap = zeros(Nt, N)
	[cap[:,i] = integrate(solx, solreg(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	
	p3 = plot(xlims=(0, BVmax), xlabel="Regenerant volume [L]", grid=:off,
			  ylims=(0, 100), ylabel="Resin capacity, %")

	areaplot!(BVreg, cap; lab=permutedims(lab), c=permutedims(col), alpha=0.5, 
			  legendposition=:bottomright)

	plot!(BVreg, cap[:,1]; lab="", lw=3, c=col[1])
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([cap[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	scatter!([BVreg[end2]], [cap[end2]];
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	annotate!((BVreg[end2]<=3 ? 0.2 : -0.2)+BVreg[end2], (BVreg[end2]<=3 ? -10 : -18)+cap[end2], ("resin capacity at $(round(BVreg[end2];digits=1))L: $(round(cap[end2];digits=1))%", 9, (BVreg[end2]<=3 ? :left : :right )))

	figure = plot!(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				   suptitle="Variable regeneration stop, loading stop at $BVstop L",
				   dpi=300, size=(680, 600))

	figure |> aside
end

# ╔═╡ 8ba53a75-e3d1-4149-a41e-e19fac5baaea
# sol2, prob2 = let

# 	prob2 = updateIC!(solreg, prob42S, solreg[t][end2]; reverse=false)
	
# 	update!(pars, prob2, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
# 			qmax => p0[1], K => vcat(1.0, p0[2:N]), 
# 	)
# 	@time sol = solve(prob2, Rodas5P(), saveat=1, progress=true, progress_steps=1)
# 	sol, prob2
# end;

# ╔═╡ ecd28db0-92db-4817-9307-9fad80dafa91
# let
# 	prob2 = updateIC!(solreg, reg42S, solreg[t][end2]; reverse=false)
# 	solx = prob2.problem_type.discretespace.axies[x]; Np = length(solx)
# 	integrate(solx, sol(solreg[t][end2], solx, dv=q(t,x)[1])) |> println
# 	integrate(solx, sol2(0, solx, dv=q(t,x)[1])) |> println

# 	plot(ylims=(0, 1.02))
# 	plot!(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]); lw=2, lab="interpolant")
# 	plot!(sol[x], solreg[q(t,x)[1]][end2, :]; lw=2, lab="original solution")
# 	scatter!(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]); lab="actual values")
# 	scatter!(sol[x][end:-1:1], sol2[q(t,x)[1]][1,:]; lab="loading")
# 	scatter!(solx, prob2.u0[@. 2Np*(N-1)+2*(1:Np)-1]; lab="initial values")
# 	scatter!(solx[end:-1:1], sol2(0, solx, dv=q(t,x)[1]); lab="solved")
# end

# ╔═╡ a96a1387-8365-481a-ad90-10a867f7a30e
# let
# 	prob = prob2; 	sol = sol2
# 	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
# 	solx = prob.problem_type.discretespace.axies[x]
# 	tstop = findfirst(x -> x >= BVstop, BV) - 1 # [idx]; 	# tstop |> println
# 	BVmax = 70; 	Nt = length(BV)

# 	# pH
# 	pHmax = 8
# 	idxmin = 3
# 	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
# 	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
# 	RC = sol(:, solx, dv=q(t,x)[1])*100
	
# 	p1 = plot(xlims=(0, BVmax), xlabel="", ylims=(0, pHmax), ylabel="pH")
# 	plot!(BV[idxmin:end], pH[idxmin:end]; 
# 		  lab="", c=col[1], lw=3)

# 	# initial pH
# 	pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
# 	prob.ps[qmax]*(1-prob.ps[ε])*prob.ps[V]*1e3
# 	hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
# 	annotate!(2, 1.5*pH0, ("initial pH: $(round(pH0;digits=2))", 9, :left))

# 	# first detection of pH change
# 	# pHchange = (round(10*(pH0 - floor(pH0));digits=1) + 1.0)/10 + 1.0
# 	hline!([pH[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	vline!([BV[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	annotate!(0.95*BV[tstop], 0.95*pHmax, 
# 			  ("pH $(round(pH[tstop];digits=2)) at $(round(Int, BV[tstop])) L", 9, :right))
# 	scatter!([BV[tstop]], [pH[tstop]]; 
# 			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

# 	# Breakthrough
# 	Cmax = 500
# 	p2 = plot(xlims=(0, BVmax), ylims=(0, Cmax),
# 			  ylabel="Concentration [mg/L]")
# 	plot!(BV, Na; lab=lab[2], col=col[2], lw=3)

# 	hline!([Cₓ₀[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=1, ls=:dot)

# 	# breakthrough detection
# 	NaBT = Na[tstop]
# 	NaBTint = integrate(BV[1:tstop], Na[1:tstop])

# 	vline!([BV[tstop]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	hline!([NaBT]; lab="", c=cp[N+2], lw=3, ls=:dashdot)
# 	annotate!(0.95*BV[tstop], 0.5Cmax, 
# 			  ("BT at $(round(Int,BV[tstop])) L: $(round(NaBTint;digits=2)) mg ($(round(NaBT;digits=2)) mg/L)", 9, :right))
# 	scatter!([BV[tstop]], [Na[tstop]]; 
# 			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

# 	# Resin capacity
# 	cap = zeros(Nt, N)
# 	[cap[:,i] = integrate(solx, sol(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	
# 	p3 = plot(xlims=(0, BVmax), xlabel="Regenerant volume [L]", grid=:off,
# 			  ylims=(0, 100), ylabel="Resin capacity, %")

# 	areaplot!(BV, cap; lab=permutedims(lab), c=permutedims(col), alpha=0.5, 
# 			  legendposition=:bottomleft)

# 	plot!(BV, cap[:,1]; lab="", lw=3, c=col[1])
# 	vline!([BV[end1]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	hline!([cap[end1]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	scatter!([BV[end1]], [cap[end1]];
# 			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)	
# 	annotate!(1.05*BV[tstop], 2.5cap[tstop], 
# 			  ("capacity at $(round(Int,BV[tstop])) L: $(round(cap[tstop];digits=2))%", 9, :left))

# 	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
# 				  dpi = 300, size=(680, 600), 
# 				  suptitle="Fixed-time breakthrough stop control scenario")
# 	figure |> aside
# end

# ╔═╡ 78bff468-778b-45f6-9758-287adb4a477e
# solreg2, reg2 = let
# 	prob = prob2
	
# 	reg2 = updateIC!(sol2, prob, sol2[t][end1]; reverse=false)
	
# 	update!(pars, reg2, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
# 			qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]),
# 	)
# 	@time solreg2 = solve(reg2, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
# 	solreg2, reg2
# end;

# ╔═╡ 44690f9e-76b0-47d9-8dc3-401ab191344e
# let
# 	reg2 = updateIC!(sol2, reg42S, sol2[t][end1]; reverse=false)
# 	solx = prob2.problem_type.discretespace.axies[x]; Np = length(solx)
# 	integrate(solx, sol2(sol2[t][end1], solx, dv=q(t,x)[1])) |> println
# 	integrate(solx, solreg2(0, solx, dv=q(t,x)[1])) |> println

# 	plot(ylims=(0, 1.02))
# 	plot!(solx, sol2(sol2[t][end1], solx, dv=q(t,x)[1]); lw=2, lab="interpolant")
# 	plot!(sol[x], sol2[q(t,x)[1]][end1, :]; lw=2, lab="original solution")
# 	scatter!(solx, sol2(sol2[t][end1], solx, dv=q(t,x)[1]); lab="actual values")
# 	scatter!(sol[x][end:-1:1], solreg2[q(t,x)[1]][1,:]; lab="regeneration")
# 	scatter!(solx, reg2.u0[@. 2Np*(N-1)+2*(1:Np)-1]; lab="initial values")
# 	scatter!(solx[end:-1:1], solreg2(0, solx, dv=q(t,x)[1]); lab="solved")
# end

# ╔═╡ 16187793-3049-4be1-a65b-1f42b147f02b
# let
# 	prob = prob2;  	reg = reg2; 	sol = sol2;  	solreg = solreg2
# 	solx = reg.problem_type.discretespace.axies[x]

# 	update!(pars, reg2, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
# 			qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]),
# 	)

# 	BV = sol[t]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]
# 	BVreg = solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3 # [L]
# 	pH = -log10.(solreg(:, 1, dv=C(t,x)[1])*reg.ps[Cref])
# 	BC = [solreg(:, 1, dv=C(t,x)[i])*reg.ps[Cref]*Mi[:Na] for i in 1:N]
# 	BVmax = 5; 	Nt = length(BVreg); 
# 	pHmax = 8;  BVmin = 10

# 	# pH
# 	p1 = plot(xlims=(0, BVmax), ylims=(-0.5, pHmax), ylabel="pH", grid=:off)
# 	plot!(BVreg[BVmin:end], pH[BVmin:end]; lab="", c=col[1], lw=3)
# 	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)

# 	# pH regenerant
# 	hline!([-log10(Cₓ₀reg[:H])]; lab="", c=:black, ls=:dot, lw=1)
# 	annotate!(0.98BVmax, 2, ("final pH: $(round(-log10(Cₓ₀reg[:H]);digits=1))", 9, :right))

# 	# Elution curve
# 	Cmax = 10.05
# 	p2 = plot(xlims=(0, BVmax), ylims=(-0.05, Cmax), ylabel="Concentration [g/L]")
# 	plot!(BVreg, [BC[i] for i in 2:N]; 
# 		  lab=permutedims(lab[2:N]), c=permutedims(col[2:N]), lw=3, grid=:off)
# 	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)

# 	# Resin capacity
# 	cap = zeros(Nt, N)
# 	[cap[:,i] = integrate(solx, solreg(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	
# 	p3 = plot(xlims=(0, BVmax), xlabel="Regenerant volume [L]", grid=:off,
# 			  ylims=(0, 100), ylabel="Resin capacity, %")

# 	areaplot!(BVreg, cap; lab=permutedims(lab), c=permutedims(col), alpha=0.5, 
# 			  legendposition=:bottomright)

# 	plot!(BVreg, cap[:,1]; lab="", lw=3, c=col[1])
# 	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	hline!([cap[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
# 	scatter!([BVreg[end2]], [cap[end2]];
# 			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
# 	annotate!((BVreg[end2]<=3 ? 0.2 : -0.2)+BVreg[end2], (BVreg[end2]<=3 ? -10 : -18)+cap[end2], ("resin capacity at $(round(BVreg[end2];digits=1))L: $(round(cap[end2];digits=1))%", 9, (BVreg[end2]<=3 ? :left : :right )))

# 	figure = plot!(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
# 				   suptitle="Fixed-time regeneration stop control scenario",
# 				   dpi=300, size=(680, 600))

# 	figure |> aside
# end

# ╔═╡ 208e50ff-84b1-4a50-ab1c-85b619b447fd
md"""
##### Fixed time
"""

# ╔═╡ 04e38cc6-5d7b-41df-8855-d1b82fc768b9
Nc = 5

# ╔═╡ 4d9d5d5f-e955-49b2-afe3-89f5acd7ce30
function fixed_cycles(; Nc=2, BVstop=20, BVregstop=1)

	prob = prob42S
	reg = reg42S
	solx = prob.problem_type.discretespace.axies[x]

	load = []
	regn = []
	temp = []

	for n in 1:Nc
		
		# Loading
	
		# initial cycle
		if n > 1
			solreg = temp
			end2 = findfirst(x->x>=BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)
			prob = updateIC!(solreg, prob, solreg[t][end2])
		end
	
		update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
				qmax => p0[1], K => vcat(1.0, p0[2:N]))
	
		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		end1 = findfirst(x->x>=BVstop, sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)

		use = integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
		push!(load, use)
	
		
		# Regeneration
		
		reg = updateIC!(sol, reg, sol[t][end1])
		
		update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
				qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
		
		solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
		temp = solreg

		end2 = findfirst(x->x>=BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)

		cap = integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
		push!(regn, cap)
	end
		
	return load, regn, Nc, BVstop, BVregstop
end

# ╔═╡ 136fc997-88c6-42d9-8116-aa63dfe29b78
fixed_time_cycles = fixed_cycles(; Nc=10, BVstop=30, BVregstop=1);

# ╔═╡ 8ba71079-5d31-4dda-963a-22a481b68495
# ╠═╡ disabled = true
#=╠═╡
cap = let
	load, regn = fixed_time_cycles
	solx = prob42S.problem_type.discretespace.axies[x]
	Nt = length(regn[1][t])
	cap = zeros(Nt, N, Nc)
	for n in 1:Nc
		solreg = regn[n]
		[cap[:, i, n] = integrate(solx, solreg(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	end
	cap
end;
  ╠═╡ =#

# ╔═╡ b2735d20-d71f-414a-aa94-31c92b09ff6e
# ╠═╡ disabled = true
#=╠═╡
use = let
	load, regn = fixed_time_cycles
	solx = prob42S.problem_type.discretespace.axies[x]
	Nt = length(load[1][t])
	use = zeros(Nt, N, Nc)
	for n in 1:Nc
		sol = load[n]
		[use[:, i, n] = integrate(solx, sol(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	end
	# plot(100 .- use[end1, 1, :])
	use
end;
  ╠═╡ =#

# ╔═╡ 72e5e318-eb42-465b-80ad-9ead7b888ba1
let
	load, regn, Nc, BVstop, BVregstop = fixed_time_cycles
	
	# resin capacity remaining
	p1 = plot(xticks=1:Nc, ylims=(60, 100), topmargin=3mm, grid=:off,
			  ylabel="Used resin capacity, %",
			  title="(a) Loading efficiency with premature stop at $(BVstop) L")

	bar!(1:Nc, 100 .- load; lab="", c=:lightblue)
	for n in 1:Nc
		annotate!(n, 102 - load[n], ("$(round(100-load[n];digits=1))", 9))
	end
	

	# resin capacity restored
	p2 = plot(xticks=1:Nc, ylims=(60, 100), grid=:off,
			  ylabel="Restored resin capacity, %",
			  title="(b) Regeneration efficiency, $(round(BVregstop)) L spent of HCl 5%")
	
	bar!(1:Nc, regn; lab="", c=:wheat)
	for n in 1:Nc
		annotate!(n, regn[n] + 2., ("$(round(regn[n];digits=1))", 9))
	end

	figure = plot(p1, p2, layout=(2,1), dpi=300, size=(680, 600),
				  xlabel="No. of cycles loading-regeneration",
				  suptitle="Performance with fixed-time stop criteria")
	figure |> aside
end

# ╔═╡ 6225b970-86c9-484d-89b3-3cf385f7af35
md"""
#### 1.1 Variable loading time
"""

# ╔═╡ 68c844a4-4f7e-4424-8679-5102cade24b9
BVstop_variable = 30:2:36

# ╔═╡ 5e436a4b-61dd-4234-ad56-51a11e665903
function variable_loading(BVstop; Nc=2, BVregstop=2)

	prob = prob42S
	reg = reg42S
	solx = prob42S.problem_type.discretespace.axies[x]

	Ns = length(BVstop)

	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

		for n in 1:Nc
			
			# Loading
		
			# initial cycle
			if n == 1
				prob = prob42S
			elseif n > 1
				solreg = temp
				end2 = time2[m, n-1]
				prob = updateIC!(solreg, prob, solreg[t][end2])
			end
		
			update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
					qmax => p0[1], K => vcat(1.0, p0[2:N]))
		
			sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)
	
			# loading time
			end1 = findfirst(x->x>=BVstop[m], sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)
			time1[m, n] = end1
	
			# resin utilization
			use = integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
			load[m, n] = use
	
			# breakthrough
			BT[m, n] = sol(sol[t][end1],1,dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
			
			# Regeneration
			reg = updateIC!(sol, reg, sol[t][end1])
			
			update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
					qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
			
			solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
			temp = solreg
	
			# regeneration time
			end2 = findfirst(x->x>=BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)
			time2[m, n] = end2
	
			# regenerated capacity
			cap = integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
			regn[m, n] = cap
		end
	
	end
		
	return load, regn, time1, time2, Nc, BT, BVregstop
end

# ╔═╡ 7a550bb5-fec5-4918-bc26-1f3dd4e6118e
cycles_variable_loading = variable_loading(BVstop_variable; Nc=5, BVregstop=2);

# ╔═╡ 547bb551-2604-4f8d-8e8b-35aba51e7a45
# cap_variable = let
# 	load, regn, time1, time2 = fixed_time_cycles_variable
# 	solx = prob42S.problem_type.discretespace.axies[x]
# 	Nt = length(regn[1][t])
# 	Ns = length(BVstop_variable)
# 	cap = zeros(Nt, N, Ns)
# 	for m in 1:Ns
# 		solreg = regn[end]
# 		[cap[:, i, m] = integrate(solx, solreg(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
# 	end
# 	cap
# end;

# ╔═╡ 65e470f2-91ef-4da6-a58f-41e10b011e5e
# use_variable = let
# 	load, regn, time1, time2 = fixed_time_cycles_variable
# 	solx = prob42S.problem_type.discretespace.axies[x]
# 	Nt = length(load[1][t])
# 	Ns = length(BVstop_variable)
# 	use = zeros(Nt, N, Ns)
# 	for m in 1:Ns
# 		sol = load[end]
# 		[use[:, i, m] = integrate(solx, sol(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
# 	end
# 	# plot(100 .- use[end1, 1, :])
# 	use
# end;

# ╔═╡ 57ab6988-0602-43a0-b9a1-bb88dc5846f8
md"""
!!! note
	The loading efficiency increases with the treated volume, but also the breakthrough for a fixed volume of spent regenerant. This is because, for the same amount of regeneration, the restored capacity decreases for an increasing volume treated since more ions are removed.
"""

# ╔═╡ 620403e0-8526-45d8-848d-235754c8b82b
let
	load, regn, time1, time2, Nc, BT, BVregstop = cycles_variable_loading
	BVstop = BVstop_variable
	Ns = length(BVstop)
	bw = (BVstop[end]-BVstop[1])/Ns
	gap = 3
	fs = 8
	lab = ["Cycle 1" "Cycle $Nc"]
		
	# resin capacity remaining
	p1 = plot(xticks=BVstop, ylims=(50, 100), grid=:off, topmargin=3mm,
			  xlabel="Treated volume during loading [L]",
			  ylabel="Used capacity, %", 
			  legendposition=:topleft, legendcolumns=1, legendfontsize=7,
			  title="(a) Loading efficiency: used capacity vs. treated volume")

	groupedbar!(BVstop, 100 .- load[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=:lightblue, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, 100 .- load[m, 1] .+ gap, 
				  ("$(round(100 - load[m, 1];digits=1))", fs, :center))
		annotate!(BVstop[m] + bw/4, 100 .- load[m, Nc] .+ gap, 
				  ("$(round(100 - load[m, Nc];digits=1))", fs, :center))
	end

	# Na breakthrough
	Cmax = 60
	p2 = plot(xticks=BVstop, ylims=(0, Cmax), grid=:off, 
			  xlabel="Treated volume during loading [L]",
			  ylabel="Concentration [mg/L]",
			  legendposition=:topleft, legendcolumns=1, legendfontsize=7,
			  title="(b) Breakthrough concentration of Na at loading stop")
	# bar!(BVstop, BT[1:Nc:end]; lab="", c=col[2])
	groupedbar!(BVstop, BT[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=col[2], alpha=[0.3 1])
	pc = 0.08
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, BT[m, 1] + pc*Cmax, 
				  ("$(round(BT[m, 1];sigdigits=2))", fs, :center))
		annotate!(BVstop[m] + bw/4, BT[m, Nc] + pc*Cmax, 
				  ("$(round(BT[m, Nc];sigdigits=2))", fs, :center))
	end
	
	# resin capacity restored
	p3 = plot(xticks=BVstop, ylims=(85, 100), grid=:off, 
			  xlabel="Treated volume during loading [L]",
			  ylabel="Restored capacity, %",
			  legendposition=:topright, legendcolumns=1, legendfontsize=7,
			  title="(c) Regeneration efficiency with $(round(Int,BVregstop)) L regenerant (HCl 5%)", topmargin=0mm)

	# bar!(BVstop, regn[1:Nc:end]; lab="", c=:wheat)
	groupedbar!(BVstop, regn[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	gap = 2
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, regn[m, 1] .+ gap, 
				  ("$(round(regn[m, 1];digits=1))", fs, :center))
		annotate!(BVstop[m] + bw/4, regn[m, Nc] .+ gap, 
				  ("$(round(regn[m, Nc];digits=1))", fs, :center))
	end

	figure = plot(p1, p2, p3, layout=(3,1), dpi=300, size=(680, 800),
				  suptitle="Variable loading volume, fixed regenerant volume")
	figure |> aside
end

# ╔═╡ 8786b6c6-77b0-4eed-adbc-4790be9ae284
md"""
#### 1.2 Variable regenerant volume
"""

# ╔═╡ 82d7e380-e1bf-43da-afdd-a32ffefa2af9
BVregstop_variable = 1.5:0.5:3.0

# ╔═╡ 12c5278e-3c18-4138-8ed4-1fa5a19a459f
function fixed_cycles_variable_reg(BVregstop; Nc=2, BVstop=35)

	prob = prob42S
	reg = reg42S
	solx = prob42S.problem_type.discretespace.axies[x]

	Ns = length(BVregstop)

	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

		for n in 1:Nc
			
			# Loading
		
			# initial cycle
			if n == 1
				prob = prob42S
			elseif n > 1
				solreg = temp
				end2 = time2[m, n-1]
				prob = updateIC!(solreg, prob, solreg[t][end2])
			end
			
			update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
					qmax => p0[1], K => vcat(1.0, p0[2:N]))
		
			sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)
	
			# loading time
			time1[m, n] = 
				end1 = findfirst(x->x>=BVstop, sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)
	
			# resin utilization
			load[m, n] = 
				integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
	
			# breakthrough
			BT[m, n] = sol(sol[t][end1],1,dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
			
			# Regeneration
			reg = updateIC!(sol, reg, sol[t][end1])

			update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
					qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
			
			solreg = temp = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
	
			# regeneration time
			time2[m, n] = end2 = 
				findfirst(x->x>=BVregstop[m],solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)
	
			# regenerated capacity
			regn[m, n] = 
				integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
		end

	end
		
	return load, regn, time1, time2, Nc, BT, BVstop
end

# ╔═╡ 96808966-f7a1-473b-8601-0ee6cb649730
variable_reg = fixed_cycles_variable_reg(BVregstop_variable; Nc=5, BVstop=30);

# ╔═╡ d8473938-0a1c-4b46-9020-75306e06fa3c
let
	load, regn, time1, time2, Nc, BT, BV = variable_reg
	BVstop = BVregstop_variable
	Ns = length(BVstop)
	fs = 8
	bw = (BVstop[end]-BVstop[1])/Ns
	gap = 3
	pc = 0.08
	lab = ["Cycle 1" "Cycle $Nc"]
		
	# resin capacity remaining
	p1 = plot(xticks=BVstop, ylims=(60, 100), grid=:off, topmargin=3mm,
			  ylabel="Used capacity, %", 
			  title="(a) Loading efficiency after $BV L treated effluent")

	groupedbar!(BVstop, 100 .- load[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=:lightblue, alpha=[0.3 1],
			    legendposition=:topright, legendcolumns=1, legendfontsize=7)
	
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, 100 .- load[m, 1] .+ gap, 
				  ("$(round(100 - load[m, 1];digits=1))", fs, :center))
		annotate!(BVstop[m] + bw/4, 100 .- load[m, Nc] .+ gap, 
				  ("$(round(100 - load[m, Nc];digits=1))", fs, :center))
	end

	# Na breakthrough
	Cmax = 2
	p2 = plot(xticks=BVstop, ylims=(0, Cmax), grid=:off,
			  ylabel="Concentration [mg/L]",
			  legendposition=:topright, legendcolumns=1, legendfontsize=7,
			  title="(b) Na breakthrough after $BV L treated effluent")

	groupedbar!(BVstop, BT[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=col[2], alpha=[0.3 1])
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, BT[m, 1] + pc*Cmax, 
				  ("$(round(BT[m, 1];sigdigits=2))", fs, :center))
		annotate!(BVstop[m] + bw/4, BT[m, Nc] + pc*Cmax, 
				  ("$(round(BT[m, Nc];sigdigits=2))", fs, :center))
	end

	# resin capacity restored
	p3 = plot(xticks=BVstop, ylims=(85, 100), grid=:off, topmargin=3mm,
			  ylabel="Restored capacity, %",
			  legendposition=:topleft, legendcolumns=1, legendfontsize=7,
			  title="(c) Regeneration efficiency vs. regenerant volume (HCl 5%)")

	groupedbar!(BVstop, regn[:, 1:(Nc-1):end];
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(BVstop[m] - bw/4, regn[m, 1] .+ 0.4gap, 
				  ("$(round(regn[m, 1];digits=1))", fs, :center))
		annotate!(BVstop[m] + bw/4, regn[m, Nc] .+ 0.4gap, 
				  ("$(round(regn[m, Nc];digits=1))", fs, :center))
	end

	figure = plot(p1, p2, p3, layout=(3,1), dpi=300, size=(680, 800),
				  xlabel="Spent regenerant volume [L]",
				  suptitle="Fixed loading time, variable regenerant volume")
	figure |> aside
end

# ╔═╡ 0243f891-343c-4dcc-a407-e09bf2fb3516
md"""
### 2. pH detection
"""

# ╔═╡ 2e4e1e1d-ef12-41dd-af3b-0ebaecd3b7ad
md"""
#### 2.1 Fixed regenerant volume
"""

# ╔═╡ 4cb7c562-b581-4ea4-bf64-925ab9819597
function pH_detection_fixed_regenerant(pHchange; Nc=2, BVregstop=2)

	prob = prob42S
	reg = reg42S
	solx = prob42S.problem_type.discretespace.axies[x]

	Ns = length(pHchange)

	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	BV = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

		for n in 1:Nc
			
			# Loading
		
			# initial cycle
			if n == 1
				prob = prob42S
			elseif n > 1
				solreg = temp
				end2 = time2[m, n-1]
				prob = updateIC!(solreg, prob, solreg[t][end2])
			end
			
			update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
					qmax => p0[1], K => vcat(1.0, p0[2:N]))
		
			sol = solve(prob, Rodas5P(), saveat=1, progress=false)

			# pH
			idxmin = 3
			pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
			pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
			end1 = findfirst(x->x>=pH0 + pHchange[m], pH[idxmin:end]) + idxmin - 1
			time1[m, n] = end1
	
			# resin utilization
			use = integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
			load[m, n] = use
	
			# breakthrough
			BT[m, n] = sol(sol[t][end1],1,dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3

			# volume treated
			BV[m, n] = sol[t][end1]*prob.ps[τ]*prob.ps[Q]*1e3
			
			# Regeneration
			if n == 1
				reg = reg42S
			elseif n > 1
			 	end1 = time1[m, n-1]
				reg = updateIC!(sol, reg, sol[t][end1])
			end

			update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
					qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
			
			solreg = solve(reg, Rodas5P(), saveat=0.01, progress=false)
			temp = solreg
	
			# regeneration time
			end2 = findfirst(x->x>=BVregstop,solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)
			time2[m, n] = end2
	
			# regenerated capacity
			cap = integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
			regn[m, n] = cap
		end

	end
		
	return load, regn, time1, time2, Nc, BT, BV, BVregstop
end

# ╔═╡ 230ed788-9526-4182-969b-0d01c20992fe
pHdetect = [0.01, 0.05, 0.1]

# ╔═╡ 8966be55-99cb-4951-ba88-315a03aa7db0
cycles_pH_detection = pH_detection_fixed_regenerant(pHdetect; Nc=3, BVregstop=2)

# ╔═╡ 827ab335-bb6c-46f0-8f37-c27bde773ef6
let
	load, regn, time1, time2, Nc, BT, BV, BVregstop = cycles_pH_detection
	BVstop = pHdetect
	Ns = length(BVstop)
	fs = 8
	bw = 0.8
	gap = 1.5
	pc = 0.08
	lab = ["1 cycle" "$Nc cycles"]
		
	# resin capacity remaining
	p1 = plot(xticks=BVstop, ylims=(80, 100), grid=:off, topmargin=0mm,
			  ylabel="Used capacity, %", 
			  title="(a) Loading efficiency, stop with pH change detection")

	groupedbar!(0:2, 100 .- load[:, 1:(Nc-1):end]; xticks=(0:2, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:lightblue, alpha=[0.3 1],
			    legendposition=:topleft, legendcolumns=1, legendfontsize=7)
	
	for m in 1:Ns
		annotate!(m-1 - bw/4, 100 .- load[m, 1] .+ gap, 
				  ("$(round(100 - load[m, 1];digits=1))", fs, :center))
		annotate!(m-1 + bw/4, 100 .- load[m, Nc] .+ gap, 
				  ("$(round(100 - load[m, Nc];digits=1))", fs, :center))
	end

	# Na breakthrough
	Cmax = 100
	p2 = plot(xticks=BVstop, ylims=(0, Cmax), grid=:off,
			  ylabel="Concentration [mg/L]",
			  legendposition=:outerright, legendcolumns=1, legend=false,
			  title="(c) Na breakthrough, stop with pH change detection")

	groupedbar!(0:2, BT[:, 1:(Nc-1):end]; xticks=(0:2, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=col[2], alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m-1 - bw/4, BT[m, 1] + pc*Cmax, 
				  ("$(round(BT[m, 1];digits=1))", fs, :center))
		annotate!(m-1 + bw/4, BT[m, Nc] + pc*Cmax, 
				  ("$(round(BT[m, Nc];digits=1))", fs, :center))
	end

	# volume treated
	p3 = plot(xticks=BVstop, ylims=(0, 50), grid=:off, topmargin=0mm,
			  ylabel="Volume treated [L]", legend=false,
			  title="(b) Effluent volume treated, variable stop moment")

	gap = 4
	groupedbar!(0:2, BV[:, 1:(Nc-1):end]; xticks=(0:2, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m-1 - bw/4, BV[m, 1] .+ gap, 
				  ("$(round(BV[m, 1];digits=1))", fs, :center))
		annotate!(m-1 + bw/4, BV[m, Nc] .+ gap, 
				  ("$(round(BV[m, Nc];digits=1))", fs, :center))
	end

	figure = plot(p1, p3, p2, layout=(3,1), dpi=300, size=(680, 800),
				  xlabel="pH detection [pH units]",
				  suptitle="Variable pH detection, fixed regenerant volume")
	figure |> aside
end

# ╔═╡ 5d49d4f7-9647-494c-9db7-c4f5b9ea7908
md"""
#### 2.2 pH detection for both phases
"""

# ╔═╡ a1725264-ff20-42fd-b4ee-66c955c2fac4
pH_change_both = [0.005, 0.01, 0.05, 0.1]

# ╔═╡ a70b158f-9b41-42ea-b600-f4c513968357
function pH_detection_both(pHchange; Nc=2)

	prob = prob42S
	reg = reg42S
	solx = prob42S.problem_type.discretespace.axies[x]

	Ns = length(pHchange)

	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	BV = zeros(Ns, Nc)
	BVreg = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

		for n in 1:Nc
			
			# Loading
		
			# initial cycle
			if n == 1
				prob = prob42S
			elseif n > 1
				solreg = temp
				end2 = time2[m, n-1]
				prob = updateIC!(solreg, prob, solreg[t][end2])
			end
			
			update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
					qmax => p0[1], K => vcat(1.0, p0[2:N]))
		
			sol = solve(prob, Rodas5P(), saveat=1, progress=false)

			# pH
			idxmin = 3
			pH0 = -log10(sum([Cₓ₀[c[i]] for i in 2:N]))
			pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
			end1 = findfirst(x->x>=pH0 + pHchange[m], pH[idxmin:end]) + idxmin - 1
			time1[m, n] = end1
	
			# resin utilization
			use = integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
			load[m, n] = use
	
			# breakthrough
			BT[m, n] = sol(sol[t][end1],1,dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3

			# volume treated
			BV[m, n] = sol[t][end1]*prob.ps[τ]*prob.ps[Q]*1e3
			
			# Regeneration
			if n == 1
				reg = reg42S
			elseif n > 1
			 	end1 = time1[m, n-1]
				reg = updateIC!(sol, reg, sol[t][end1])
			end

			update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
					qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
			
			solreg = solve(reg, Rodas5P(), saveat=0.01, progress=false)
			temp = solreg
	
			# regeneration time
			pH0 = -log10(Cₓ₀reg[:H])
			pH = -log10.(solreg(:, 1, dv=C(t,x)[1])*reg.ps[Cref])
			end2 = findfirst(x->-x>=-(pH0+pHchange[m]), pH[idxmin:end])+idxmin-1
			BVreg[m, n] = solreg[t][end2]*reg.ps[τ]*reg.ps[Q]*1e3
			time2[m, n] = end2
	
			# regenerated capacity
			cap = integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
			regn[m, n] = cap
		end

	end
		
	return load, regn, time1, time2, Nc, BT, BV, BVreg
end

# ╔═╡ b8337a8b-26d6-4bcb-a710-899d2f1f3799
cycles_pH_detection_both = pH_detection_both(pH_change_both; Nc=3)

# ╔═╡ 8c9d581c-e9ef-43fb-b8ce-c3f5bfa50f90
let
	load, regn, time1, time2, Nc, BT, BV, BVreg = cycles_pH_detection_both
	BVstop = pH_change_both
	Ns = length(BVstop)
	fs = 8
	bw = 0.8
	gap = 1.5
	pc = 0.08
	lab = ["1 cycle" "$Nc cycles"]

	# volume treated
	p1 = plot(xticks=BVstop, ylims=(0, 50), grid=:off, topmargin=0mm,
			  ylabel="Volume [L]", legendposition=:outertop, legendcolumns=2,
			  title="(a) Effluent volume treated")

	gap = 4
	groupedbar!(1:Ns, BV[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BV[m, 1] .+ gap, 
				  ("$(round(BV[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BV[m, Nc] .+ gap, 
				  ("$(round(BV[m, Nc];digits=1))", fs, :center))
	end

	# Na breakthrough
	Cmax = 100
	p2 = plot(xticks=BVstop, ylims=(0, Cmax), grid=:off,
			  ylabel="Concentration [mg/L]",
			  legendposition=:outerright, legendcolumns=1, legend=false,
			  title="(b) Na breakthrough concentration")

	groupedbar!(1:Ns, BT[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=col[2], alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BT[m, 1] + pc*Cmax, 
				  ("$(round(BT[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BT[m, Nc] + pc*Cmax, 
				  ("$(round(BT[m, Nc];digits=1))", fs, :center))
	end

	# regenerant volume spent
	p3 = plot(xticks=BVstop, ylims=(0, 4), grid=:off, topmargin=0mm,
			  ylabel="Volume [L]", legend=false,
			  title="(c) Regenerant volume spent")

	gap = 0.4
	groupedbar!(1:Ns, BVreg[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BVreg[m, 1] .+ gap, 
				  ("$(round(BVreg[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BVreg[m, Nc] .+ gap, 
				  ("$(round(BVreg[m, Nc];digits=1))", fs, :center))
	end
	

	figure = plot(p1, p2, p3, layout=(3,1), dpi=300, size=(680, 800),
				  xlabel="pH detection [pH units]",
				  suptitle="Loading/regeneration stop with pH detection")
	figure |> aside
end

# ╔═╡ e11cef45-c956-4cdb-9568-005b1c0e8cf2
md"""
### 3. Na detection
"""

# ╔═╡ e760a66a-d418-4b13-af96-37232574b035
md"""
With a sodium analyzer at the outlet, both loading and regeneration can be stopped after detection of a certain threshold: over the threshold for breakthrough concentration, and under the threshold for the eluted concentration.
"""

# ╔═╡ 225d9033-75fd-4c18-a35d-ed9e6306d2c3
function Na_detection_both(Na_detect; Nc=2)

	prob = prob42S
	reg = reg42S
	solx = prob42S.problem_type.discretespace.axies[x]

	Ns = length(Na_detect)

	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	BV = zeros(Ns, Nc)
	BVreg = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

		for n in 1:Nc
			
			# Loading
		
			# initial cycle
			if n == 1
				prob = prob42S
			elseif n > 1
				solreg = temp
				end2 = time2[m, n-1]
				prob = updateIC!(solreg, prob, solreg[t][end2])
			end
			
			update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
					qmax => p0[1], K => vcat(1.0, p0[2:N]))
		
			sol = solve(prob, Rodas5P(), saveat=1, progress=false)

			# Na
			Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3
			idxmin = findfirst(x -> x > Na_detect[m], Na)
			end1 = findfirst(x -> x >= Na_detect[m], Na[idxmin:end])+idxmin-1 
			time1[m, n] = end1
	
			# resin utilization
			use = integrate(solx, sol(sol[t][end1], solx, dv=q(t,x)[1]))*100 
			load[m, n] = use
	
			# breakthrough
			BT[m, n] = sol(sol[t][end1],1,dv=C(t,x)[2])*prob.ps[Cref]*Mi[:Na]*1e3

			# volume treated
			BV[m, n] = sol[t][end1]*prob.ps[τ]*prob.ps[Q]*1e3
			
			# Regeneration
			if n == 1
				reg = reg42S
			elseif n > 1
			 	end1 = time1[m, n-1]
				reg = updateIC!(sol, reg, sol[t][end1])
			end

			update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
					qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
			
			solreg = solve(reg, Rodas5P(), saveat=0.01, progress=false)
			temp = solreg
	
			# regeneration time
			idxmin = 200
			Na = solreg(:, 1, dv=C(t,x)[2])*reg.ps[Cref]*Mi[:Na]*1e6 # [g/L]
			end2 = findfirst(x-> -x > -Na_detect[m]*1e-3, Na[idxmin:end])+idxmin-1
			BVreg[m, n] = solreg[t][end2]*reg.ps[τ]*reg.ps[Q]*1e3
			time2[m, n] = end2
	
			# regenerated capacity
			cap = integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))*100
			regn[m, n] = cap
		end

	end
	
	return load, regn, time1, time2, Nc, BT, BV, BVreg
end

# ╔═╡ 05fcd099-a261-419c-9d68-635ebecd2413
Na_detect = [0.02, 0.1, 1, 10, 100]

# ╔═╡ 55ad8e2c-4e30-4a70-8ec2-6d4424fc070b
cycles_Na_detection_both = Na_detection_both(Na_detect; Nc=3)

# ╔═╡ 8a61660d-0ee0-41a3-9515-2afab8de9133
let
	load, regn, time1, time2, Nc, BT, BV, BVreg = cycles_Na_detection_both
	BVstop = Na_detect
	BVstop = ["20 ppb", "100 ppb", "1 ppm", "10 ppm", "100 ppm"]
	Ns = length(BVstop)
	fs = 8
	bw = 0.8
	gap = 1.5
	pc = 0.08
	lab = ["1 cycle" "$Nc cycles"]

	# volume treated
	p1 = plot(xticks=BVstop, ylims=(0, 50), grid=:off, topmargin=0mm,
			  ylabel="Volume [L]", legendposition=:outertop, legendcolumns=2,
			  title="(a) Effluent volume treated")

	gap = 4
	groupedbar!(1:Ns, BV[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BV[m, 1] .+ gap, 
				  ("$(round(BV[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BV[m, Nc] .+ gap, 
				  ("$(round(BV[m, Nc];digits=1))", fs, :center))
	end

	# Na breakthrough
	Cmax = 100
	p2 = plot(xticks=BVstop, ylims=(0, Cmax), grid=:off,
			  ylabel="Concentration [mg/L]",
			  legendposition=:outerright, legendcolumns=1, legend=false,
			  title="(b) Na breakthrough concentration")

	groupedbar!(1:Ns, BT[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=col[2], alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BT[m, 1] + pc*Cmax, 
				  ("$(round(BT[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BT[m, Nc] + pc*Cmax, 
				  ("$(round(BT[m, Nc];digits=1))", fs, :center))
	end

	# regenerant volume spent
	p3 = plot(xticks=BVstop, ylims=(0, 3), grid=:off, topmargin=0mm,
			  ylabel="Volume [L]", legend=false,
			  title="(c) Regenerant volume spent")

	gap = 0.4
	groupedbar!(1:Ns, BVreg[:, 1:(Nc-1):end]; xticks=(1:Ns, BVstop),
			    lab=lab, barwidths=bw, linealpha=1,	c=:wheat, alpha=[0.3 1])
	for m in 1:Ns
		annotate!(m - bw/4, BVreg[m, 1] .+ gap, 
				  ("$(round(BVreg[m, 1];digits=1))", fs, :center))
		annotate!(m + bw/4, BVreg[m, Nc] .+ gap, 
				  ("$(round(BVreg[m, Nc];digits=1))", fs, :center))
	end
	

	figure = plot(p1, p3, layout=(2, 1), dpi=300, size=(680, 600),
				  xlabel="Na breakthrough threshold",
				  suptitle="Loading/regeneration stop with Na detection")
	figure |> aside
end

# ╔═╡ ec9597bb-26db-4c1e-9896-ea221a4f7f9b
md"""
### 4. Model-based prediction criteria
"""

# ╔═╡ 96ea5cd6-3f76-4916-8698-4f37f0a3eec2
md"""
Table with comparison of all criteria, two by two, and a few cycles
"""

# ╔═╡ d7b1b4f5-876d-4472-b69e-ceeec2c29d78
md"""
#### 4.1 Resin utilization
"""

# ╔═╡ 3e1e978f-3588-4747-a6fe-ed165c9a99ea
function fixed_capacity(; Nc=2, RCstop=90, BVregstop=1)

	prob = prob42S
	reg = reg42S
	solx = prob.problem_type.discretespace.axies[x]

	Ns = length(RCstop)
	
	load = zeros(Ns, Nc)
	regn = zeros(Ns, Nc)
	BT = zeros(Ns, Nc)
	time1 = zeros(Int, Ns, Nc)
	time2 = zeros(Int, Ns, Nc)
	temp = []

	for m in 1:Ns

	for n in 1:Nc
		
		# Loading
	
		# initial cycle
		if n == 1
			prob = prob42S
		elseif n > 1
			solreg = temp
			end2 = time2[m, n-1]
			prob = updateIC!(solreg, prob, solreg[t][end2])
		end
	
		update!(pars, prob, Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
				qmax => p0[1], K => vcat(1.0, p0[2:N]))
	
		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		# remaining resin capacity
		RC = integrate(solx, sol(:, solx, dv=q(t,x)[1]); dims=1)
		
		time1[m, n] = end1 = findfirst(x -> x >= RCstop[m], 1 .- RC)

		load[m, n] = sol[t][end1]*prob.ps[τ]*prob.ps[Q]*1e3 # [L]

		# breakthrough concentration
		BT[m, n] = sol(sol[t][end1], 1, dv=C(t,x)[2]) 		
	
		
		# Regeneration
		reg = updateIC!(sol, reg, sol[t][end1])
		
		update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
				qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]))
		
		solreg = temp =
			solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)

		time2[m, n] = end2 = 
			findfirst(x-> x >= BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)

		# restored capacity
		regn[m, n] = 
			integrate(solx, solreg(solreg[t][end2], solx, dv=q(t,x)[1]))
	end

	end
		
	return load, regn, time1, time2, BT, Nc, RCstop, BVregstop
end

# ╔═╡ 84d5df10-0cf3-4f1a-adf0-4e4d1a1b13cc
cycle_fixed_capacity = fixed_capacity(Nc=5, RCstop=0.6:0.05:0.75, BVregstop=1)

# ╔═╡ 2b7429de-c798-425a-bf37-8bc8bbaa67b3
md"""
!!! note
	In the first cycle, the resin capacity is fully restored and the loading efficiency is maximum: the treated volume decreases with each cycle, since the capacity is not fully restored.

	A higher restored capacity leads to a larger volume of treated effluent during loading, but also to higher breakthrough. Consequently, to reach the same restored capacity, a larger volume is needed in consecutive cycles.
"""

# ╔═╡ fc2358d9-b5cf-401a-ab0f-a68b08cb4359
let
	load, regn, time1, time2, BT, Nc, RCstop, BVregstop = cycle_fixed_capacity

	Ns = length(RCstop)
	fs = 8
	bw = 4
	gap = 3
	pc = 0.08
	lab = ["Cycle 1" "Cycle $Nc"]
	
	# resin capacity remaining
	p1 = plot(xticks=1:Nc, ylims=(0, 50), yticks=0:10:50, topmargin=3mm, 
			  legendposition=:topleft, grid=:off, legendfontsize=7,
			  xlabel="Used resin capacity, %",
			  ylabel="Treated volume [L]",
			  title="(a) Loading efficiency for a given % used resin capacity")

	groupedbar!(100RCstop, load[:, 1:(Nc-1):end]; xticks=100RCstop,
				lab=lab, c=:lightblue, alpha=[0.3 1], linealpha=1, barwidths=bw)
	
	for m in 1:Ns
		annotate!(100RCstop[m] - bw/4, load[m, 1] .+ gap, 
				  ("$(round(load[m, 1];sigdigits=3))", fs, :center))
		annotate!(100RCstop[m] + bw/4, load[m, Nc] .+ gap, 
				  ("$(round(load[m, Nc];sigdigits=3))", fs, :center))
	end
	

	# breakthrough
	p2 = plot(xticks=1:Nc, ylims=(0, 70), topmargin=3mm, grid=:off,
			  legendposition=:topleft, legendfontsize=7,
			  xlabel="Used resin capacity, %",
			  ylabel="Concentration [mg/L]",
			  title="(b) Effluent quality: Na breakthrough concentration, mg/L")

	groupedbar!(100RCstop, BT[:, 1:(Nc-1):end]*Mi[:Na]*1e3; xticks=100RCstop,
				lab=lab, c=cp[1], alpha=[0.3 1], linealpha=1)
	
	for m in 1:Ns
		annotate!(100RCstop[m] - bw/4, BT[m, 1]*Mi[:Na]*1e3 .+ 1.5gap, 
				  ("$(round(BT[m, 1]*Mi[:Na]*1e3;digits=1))", fs, :center))
		annotate!(100RCstop[m] + bw/4, BT[m, Nc]*Mi[:Na]*1e3 .+ 1.5gap, 
				  ("$(round(BT[m, Nc]*Mi[:Na]*1e3;digits=1))", fs, :center))
	end
	

	# resin capacity restored
	p3 = plot(xticks=1:Nc, ylims=(60, 100), grid=:off,
			  legendposition=:topright, legendcolumns=1, legendfontsize=7,
			  xlabel="Used resin capacity during loading, %",
			  ylabel="Restored capacity, %",
			  title="(c) Regeneration efficiency with $BVregstop L spent of HCl 5%")

	groupedbar!(100RCstop, 100regn[:, 1:(Nc-1):end]; xticks=100RCstop,
				lab=lab, c=:wheat, alpha=[0.3 1], linealpha=1)
	
	for m in 1:Ns
		annotate!(100RCstop[m] - bw/4, 100regn[m, 1] .+ gap, 
				  ("$(round(100regn[m, 1];digits=1))", fs, :center))
		annotate!(100RCstop[m] + bw/4, 100regn[m, Nc] .+ gap, 
				  ("$(round(100regn[m, Nc];digits=1))", fs, :center))
	end

	figure = plot(p1, p2, p3, layout=(3,1), dpi=300, size=(680, 800),
				  suptitle="Model-based stop criteria: resin utilization")
	figure |> aside
end

# ╔═╡ bf92682d-4e9e-459c-83e2-3bcfd3729360
md"""
#### 4.2 Resin regeneration
"""

# ╔═╡ d24ddf3f-f5f6-4d6a-bbc7-e2a1ac30ba95


# ╔═╡ deda98c1-ca29-40cd-9f7e-b9b2b5852afd
md"""
# Model-based control quality checks
"""

# ╔═╡ 4d53a2a4-20db-4147-98e7-fdc35d310b42
md"""
- Model prediction of breakthrough based on 95% capacity e.g.

- Uncertainty: compare with model-generated data and modified parameters (e.g. K, qmax)
"""

# ╔═╡ 1d4a3d18-71ca-464e-ab06-b70a0d89140c


# ╔═╡ 60c61a84-3bba-42e4-9aa9-b491f11f11ca
md"""
# Varying influent & model errors
"""

# ╔═╡ 76dd8dd9-f193-4e2e-a782-7c6f4059dc9b
σ = Dict([i => 0. for i in c])

# ╔═╡ 57c724f7-4272-44b0-82d2-edd8548627c2
σ₀ = 10 # [mg/L] so that 95% of samples are within (-10, +10) mg/L

# ╔═╡ e9bea4f3-6d60-4b79-a184-b0ba57ae9f35
[σ[i] = σ₀/(Mi[i]*1e3)*Zi[i] for i in c[2:N]]; σ[:H] = 1e-7; σ

# ╔═╡ c8e3f872-817a-4c6c-89d0-0707aa6f7e2d
begin
	Cₓ₀[:H] = 10.0^(-6); 	Cₓ₀mass[:Na] = 100.0; 	Cₓ₀mass[:K] = 100.0;
							Cₓ₀mass[:Ca] = 100.0; 	Cₓ₀mass[:Mg] = 100.0;
	[Cₓ₀[i] = Cₓ₀mass[i]/(Mi[i]*1e3)*Zi[i] for i in c[2:5]] # [eq/L]
	Cₓ₀
end

# ╔═╡ 9e6eeeb1-dc27-4fa5-b5ea-ebb759c37915
function noise(t)
	# Random.seed!(1234);

	# w = LinearInterpolation(0:1:tmax, randn(length(0:1:tmax)))
	
	return randn()
end	

# ╔═╡ e2ea5ed4-db1b-40ec-a93d-8c9ea4f9828b
noise(0.1)

# ╔═╡ 60c14f0d-3f6d-42d9-b737-9a3a9cc8aca5
@register_symbolic noise(t)

# ╔═╡ a31d2322-a171-4343-b121-48d3bba29689
function step(t)
	if t > 200
		return 1
	else
		return 0
	end
end

# ╔═╡ a63ec162-12e7-42a4-9344-6576de0b50de
@register_symbolic step(t)

# ╔═╡ 61e3201e-788f-4846-9ef7-b320f4e67aa1
bcs_var = let
	Random.seed!(5);
	bcs = Array{Symbolics.Equation}(undef, 2N + 2N) # 2 ICs + 2x2 BCs
		for (i, j) in zip(1:N, c)
			# Initial conditions
			bcs[i] 		= C(0.0, x)[i] ~ C₀[j] 
			bcs[i + N] 	= q(0.0, x)[i] ~ q₀[j] 
			# Boundary conditions - Dirichlet for Cₓ₀, Neumann for the rest
			bcs[i + 2N] = C(t, 0.0)[i] ~ (Cx0[i] + σ[j]*sin(0.5randn()*t))/Cref
				# σ₀/(Mi[:Na]*1e3)*Zi[:Na] * noise(t) / Cref 
			bcs[i + 3N] = Dx(q(t, 0.0)[i]) ~ qx0[i]
		end
	bcs
end

# ╔═╡ 33b03f13-636b-4038-8c75-42b777ad50e5
@named pdesys_var = PDESystem(eqs=eqs_same, bcs=bcs_var, domain=domains, 
							   ivs=[t, x], dvs=dvs, ps=ps, defaults=pars);

# ╔═╡ 4e310059-4f33-45ea-9d60-da1d7824cbd1
@time prob42V = discretize_problem(pdesys_var, dx=0.025, grid_align=center_align);

# ╔═╡ c9da859c-d33f-42d8-bf57-09f5d3fc7c71
md"""
## Loading

"""

# ╔═╡ 3a1cd15b-d101-4e92-8c3b-2339ec7d0cc5
md"""
!!! note 
	With sinusoidal noise in the influent, there is little difference in the observed breakthrough and resin utilization, since on average the adsorbed ion amounts are the same.
"""

# ╔═╡ 4c139460-3116-49c1-8726-40d1cb78bd93
function MonteCarlo(p)

	prob = prob42S

	update!(pars, prob, Q => 120/60/1e6, Cx0 => [Cₓ₀[i] for i in c], 
			qmax => p[1], K => vcat(1.0, p[2:N]); verbose=false)
	
	return sol = solve(prob, Rodas5P(), saveat=1.)
end

# ╔═╡ 2c857135-9ab3-42c5-80c2-3edd5bfd155f
MC, CI, samples = let
	Np = length(p0)
	samples = 300
	CI = [0.34, 5.44, 7.56, 8.83, 10.7]/100 .* p0 
	
	MC = zeros(Np, samples)
	for i in 1:Np
		MC[i, :] = QuasiMonteCarlo.sample(samples, 0., 1., 
										 truncated(Normal(p0[i], CI[i]/2), 
												   lower=0.))
	end
	MC, CI, samples
end

# ╔═╡ 6112953c-9566-4f08-b480-f9adcf96eb4a
solv = let
	Random.seed!(seed_model);

	prob = prob42V
	p = p0 .+ CI .* randn(length(p0))
	update!(pars, prob,	Q => 120/60*1e-6, Cx0 => [Cₓ₀[i] for i in c], 
			qmax => p[1], K => vcat(1.0, p[2:N]), 
	)
	@time sol = 
		solve(prob, Rodas5P(), saveat=0.1, progress=true, progress_steps=1)
end;

# ╔═╡ 91e0f021-db31-4f7a-8cc1-fea9615d681d
let
	α = 0.05
	p = [plot() for i in 2:N]
	for i in 2:N
		p[i-1] = histogram(MC[i,:]; normalize=:probability,
						 c=cp[i-1], lab="", title=L"$\alpha_%$(i)$")

		# calibrated
		vline!([p0[i]]; lab="", c=:black, ls=:solid, lw=3)
		vline!([p0[i]+CI[i]]; lab="", c=:black, lw=3)
		vline!([p0[i]-CI[i]]; lab="", c=:black, lw=3)

		# sampled
		vline!([mean(MC[i, :])]; lab="", c=:red, ls=:dot, lw=2)
		vline!([quantile(MC[i, :], 1-α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
		vline!([quantile(MC[i, :], α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
	end

	# resin capacity
	i = 1
	p5 = plot()
	histogram!(MC[i,:]; normalize=:probability,
					 c=:grey, lab="", title=L"$q_{max}$")
	
	vline!([p0[i]]; lab="", c=:black, ls=:solid, lw=3)
	vline!([p0[i]+CI[i]]; lab="", c=:black, lw=3)
	vline!([p0[i]-CI[i]]; lab="", c=:black, lw=3)
	
	vline!([mean(MC[i, :])]; lab="", c=:red, ls=:dot, lw=2)
	vline!([quantile(MC[i, :], 1-α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
	vline!([quantile(MC[i, :], α/2)]; lab="", c=:red, ls=:dashdot, lw=2)

	plot(p..., p5, layout=(3,2), dpi=300, size=(680, 600), suptitle="Sampling")
	plot!() |> aside
end

# ╔═╡ 97005938-587d-4a58-83c0-0258cb9c8b74
CMC = let
	Ns = samples
	CMC = zeros(length(sol[t]), Ns, N)

	for n in 1:Ns
		p = MC[:, n, :]
		sol = MonteCarlo(p)
		[CMC[:, n, i] = sol(:, 1, dv=C(t,x)[i]) for i in 1:N]
	end
	CMC
end

# ╔═╡ 7f3eb011-a1ae-4b24-a4e7-ab08043f80db
# let
# 	Ns = size(CMC, 2) 
# 	α = 0.05
# 	plot(dpi=300, xlims=(0, 800), grid=:off)

# 	# Monte Carlo simulations
# 	plot!(sol[t], [CMC[:, :, i] for i in 2:N]; lab="", lw=0.5)

# 	# Average
# 	for i in 2:N
# 	plot!(sol[t], mean(CMC[:,:,i], dims=2); 
# 		  lab="Mean ".*lab[i], lw=2, c=col[i])

# 	plot!(sol[t], mean(CMC[:,:,i], dims=2); 
# 		  lab="", lw=2, c=col[i], fillcolor=col[i], fillalpha=0.7,
# 		  fillrange=([quantile(CMC[j,:,i], α/2) for j in 1:length(sol[t])], 
# 					 [quantile(CMC[j,:,i], 1-α/2) for j in 1:length(sol[t])]))
# 	end
	
# 	# plot!(sol[t], [quantile(CMC[i,:], 1-α/2) for i in 1:length(sol[t])];
# 	# 	  lab="", c=:black, ls=:dashdot)
# 	# plot!(sol[t], [quantile(CMC[i,:], α/2) for i in 1:length(sol[t])];
# 	# 	  lab="", c=:black, ls=:dashdot)
	
# 	plot!() |> aside
# end

# ╔═╡ 4353e933-4f0a-4165-a0ab-6e4d93f52f2f
@bindname seed_model confirm(NumberField(0:100, default=16), label="Update")

# ╔═╡ 8142e92b-39d3-4dca-8fe0-588d0e70aba8
@bindname seed_test NumberField(0:100, default=16)

# ╔═╡ 068dfd04-3923-4b39-a952-1bf9206a74cb
let
	# seed_model = 15
	Random.seed!(seed_test);

	names = vcat("qm", ["α"*string(i) for i in 2:5])
	Np = length(names)
	δ = CI./p0*100

	"Calibrated model parameters:" |> println
	join(@sprintf("%s: %.2f ± %.2f (%.2f%%)\n", names[i], p0[i], CI[i], δ[i]) for i in 1:Np) |> println

	p = p0 .+ CI .* randn(Np)
		
	"10% random variability:" |> println
	join(@sprintf("%s: %.2f ± %.2f (%.2f%%)\n", names[i], p[i], p[i]-p0[i], 100*(p[i]/p0[i]-1)) for i in 1:Np) |> println
end

# ╔═╡ 3b04b7e9-a9ae-4069-8edf-fc374e9e696e
let
	α = 0.05
	BV = sol[t]*prob42S.ps[τ]*prob42S.ps[Q]*1e3 # [L]
	BVv = solv[t]*prob42V.ps[τ]*prob42V.ps[Q]*1e3 # [L]
	Vmax = 80 # [L]
	Vres = prob42V.ps[V]*(1-prob42V.ps[ε])*1e3 # [L]

	qe = [prob42V.ps[qmax]*Mi[c[i]]/Zi[c[i]]*prob42V.ps[K][i]*Cₓ₀[c[i]]/
			sum([prob42V.ps[K][j]*Cₓ₀[c[j]] for j in 1:N]) for i in 1:N] # [g/L]

	# pH
	p1 = plot(xlims=(0, Vmax), ylims=(0, 8), legendposition=:topleft, 
			  xlabel="Treated volume [L]", ylabel="pH", grid=:off)

	plot!(BVv[3:end], 
		  -log10.(prob42V.ps[Cref]*solv(solv[t][3:end], 1, dv=C(t,x)[1]));
		  lw=4, c=col[1], lab="", ls=:dot)

	plot!(BV[3:end], -log10.(mean(CMC[3:end,:,1], dims=2)*prob42S.ps[Cref]); 
		  lab="pH", lw=2, c=col[1], fillcolor=col[1], fillalpha=0.3,
		  fillrange=([quantile(-log10.(CMC[j,:,1]*
			  prob42S.ps[Cref]), α/2) for j in 3:length(sol[t])], 
					 [quantile(-log10.(CMC[j,:,1]*
			  prob42S.ps[Cref]), 1-α/2) for j in 3:length(sol[t])]))

	# Breakthrough curves
	p2 = plot(xlims=(0, Vmax), ylims=(0, 6.0), legendposition=:topright, 
			  xlabel="Treated volume [L]", ylabel="Concentration, C/C₀",
			  title="Variable influent concentration: ±$σ₀% C/C₀")
	
	# influent
	hline!([1]; lab="", c=:black, lw=0.5, ls=:dot, grid=:off)
	plot!(BVv, [solv(:, 0., dv=C(t,x)[i])*prob42V.ps[Cref]/Cₓ₀[c[i]] for i in 2:N]; lab="", c=permutedims(col[2:N]), lw=0.5)
	
	# variable
	plot!(BVv, [solv(solv[t], 1, dv=C(t,x)[i])*prob42V.ps[Cref]/Cₓ₀[c[i]] for i in 2:N]; c=permutedims(col[2:N]), lab="", lw=2.5, ls=:solid)

	plot!(BVv, [solv(solv[t], 0., dv=C(t,x)[i])* 
		  prob42V.ps[Cref]/Cₓ₀[c[i]] for i in 2:N]; 
		  lab="", c=permutedims(col[2:N]), lw=1.5,
		  inset = bbox(0.17, 0.15, 0.3, 0.3), subplot=2, 
		  xlims=(0, 5), ylims=(0.75, 1.25), yticks=0.8:0.2:1.2,
		  title="", titlefontsize=9)
	hline!([1]; lab="", c=:black, lw=0.5, ls=:dot, subplot=2)

	plot!([], []; inset=bbox(0.4, 0.775, 12/80, 12/500), subplot=3, box=:box, 
		  grid=:off, ticks=false, backgroundinside=false, legend=false, 
		  bordercolor=:grey)

	plot!([], []; inset=bbox(0.0975, 0.65, 5/80, 0.4/5), subplot=4, box=:box, 
		  grid=:off, ticks=false, backgroundinside=false, legend=false, 
		  bordercolor=:grey)

	# model
	for i in 2:N
	plot!(BV, mean(CMC[:,:,i], dims=2)*prob42S.ps[Cref]/Cₓ₀[c[i]]; 
		  lab=lab[i], lw=2, c=col[i], alpha=0.5, ls=:dashdot,
		  fillcolor=col[i], fillalpha=0.3,
		  fillrange=([quantile(CMC[j,:,i]*
			  prob42S.ps[Cref]/Cₓ₀[c[i]], α/2) for j in 1:length(sol[t])], 
					 [quantile(CMC[j,:,i]*
			  prob42S.ps[Cref]/Cₓ₀[c[i]], 1-α/2) for j in 1:length(sol[t])]))
	end

	plot!(twinx(), [], []; lab="", xlims=(0, Vmax), ylims=(0, 600), grid=:off, 
		  ylabel="Concentration [mg/L]")


	# Predicted breakthrough
	Cmax = 0.06
	p5 = plot(xlims=(30, 42), xticks=30:2:42, 
			  ylims=(-0.0002, 0.06), yticks=0:0.01:0.06, grid=:off,
			  xlabel="Treated volume [L]", ylabel="Concentration, C/C₀",
			  title="5% breakthrough time: mean prediction + 95% CI")

	p = vcat(prob42V.ps[qmax], prob42V.ps[K][2:N])
	Δ = round.(CI;digits=2)
	δ = round.(100CI./p;sigdigits=2)
	# annotate!(30.5, 0.75Cmax, ("Model parameters:", 9, :left))
	# annotate!(30.5, 0.65Cmax, ("qmax: $(p0[1]) ± $(Δ[1]) ($(δ[1])%) eq/L", 9, :left))
	# annotate!(30.5, 0.55Cmax, ("α(Na): $(p0[2]) ± $(Δ[2]) ($(δ[2])%)", 9, :left))
	# annotate!(30.5, 0.45Cmax, ("α(K):   $(p0[3]) ± $(Δ[3]) ($(δ[3])%)", 9, :left))
	# annotate!(30.5, 0.35Cmax, ("α(Mg): $(p0[4]) ± $(Δ[4]) ($(δ[4])%)", 9, :left))
	# annotate!(30.5, 0.25Cmax, ("α(Ca): $(p0[5]) ± $(Δ[5]) ($(δ[5])%)", 9, :left))

	idxP = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.5) for j in 1:length(sol[t])])
	idxC1 = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.975) for j in 1:length(sol[t])])
	idxC2 = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.025) for j in 1:length(sol[t])])
	idxM = findfirst(x -> x >= 0.05, solv(:, 1, dv=C(t,x)[2]))

	annotate!(30.5, 0.95Cmax, ("Prediction: $(round(BV[idxP];digits=1)) ± $(round(BV[idxP]-BV[idxC1];digits=1)) L", 9, :left))
	annotate!(30.5, 0.75Cmax, ("Validation: $(round(BVv[idxM];digits=1)) L", 9, :left))

	hline!([0.05]; lab="", c=cp[N+2], ls=:dot, lw=3)
	vline!([BVv[idxM]]; lab="", c=cp[N+2], ls=:dot, lw=3)

	# variable
	plot!(BVv, [solv(solv[t], 1, dv=C(t,x)[i])*
		prob42V.ps[Cref]/Cₓ₀[c[i]] for i in 2:N]; 
		  c=permutedims(col[2:N]), lab="", lw=4, ls=:solid)
	# constant
	# plot!(BV, [sol(sol[t], 1, dv=C(t,x)[i])*
	# 	  prob42S.ps[Cref]/Cₓ₀[c[i]] for i in 2:N]; 
	# 	  c=permutedims(col[2:N]), lab="", lw=3, ls=:dot)

	# Uncertainty
	for i in 2:N
	plot!(BV, [quantile(CMC[j,:,i], 0.5) for j in 1:length(sol[t])]*prob42V.ps[Cref]/Cₓ₀[c[i]]; 
		  lab=lab[i], lw=2, c=col[i], ls=:dashdot, alpha=0.5,
		  fillcolor=col[i], fillalpha=0.2,
		  fillrange=([quantile(CMC[j,:,i]*
			  prob42V.ps[Cref]/Cₓ₀[c[i]], α/2) for j in 1:length(sol[t])], 
					 [quantile(CMC[j,:,i]*
			prob42V.ps[Cref]/Cₓ₀[c[i]], 1-α/2) for j in 1:length(sol[t])]))
	# plot!(BV, [quantile(CMC[j,:,i]*
	# 	Mi[c[i]]/Zi[c[i]], 1-α/2) for j in 1:length(sol[t])];
	# 	  lab="", c=col[i], ls=:dashdot)
	# plot!(BV, [quantile(CMC[j,:,i]*
	# 	Mi[c[i]]/Zi[c[i]], α/2) for j in 1:length(sol[t])];
	# 	  lab="", c=col[i], ls=:dashdot)
	end

	# % breakthrough
	plot!(twinx(), [], []; lab="", xlims=(30, 40), xticks=30:2:40, grid=:off,
		  ylims=(0, 6), yticks=0:1:6, ylabel="Concentration [mg/L]")

	
	# Accumulated in effluent
	# p3 = plot(xlims=(0, Vmax), ylims=(0, 0.1Vmax),
	# 		  xlabel="Treated volume [L]", ylabel="Effluent [mg]")

	# plot!(BV, BV*Cₓ₀mass[:Na]*1e-3; lab="", lw=1, c=:black, ls=:dashdot)
	# # plot!(BVv[2:end], [integrate(BVv[1:j], solv(solv[t][1:j], 0, dv=C(t,x)[2])) for j in 2:length(BV)]*prob42V.ps[Cref]*Mi[:Na]; c=:black, lw=2)
	
	# plot!(BV, [BV*Cₓ₀mass[c[i]]*1e-3 .- qe[i]*Vres for i in 2:N]; 
	# 	  lab="", ls=:dot, c=permutedims(col[2:N]))

	# plot!(BV[2:end], [[integrate(BV[1:j], sol(sol[t][1:j], 1, dv=C(t,x)[i])) for j in 2:length(BV)]*prob42S.ps[Cref]*Mi[c[i]]/Zi[c[i]] for i in 2:N];
	# 	  lab=permutedims(lab[2:N]), c=permutedims(col[2:N]), lw=2)
	
	# plot!(BVv[2:end], [[integrate(BVv[1:j], solv(solv[t][1:j], 1, dv=C(t,x)[i])) for j in 2:length(BVv)]*prob42V.ps[Cref]*Mi[c[i]]/Zi[c[i]] for i in 2:N];
	# 	  lab="", c=permutedims(col[2:N]), lw=2, ls=:dot)
	

	# Accumulated in resin
	# p4 = plot(xlims=(0, Vmax), ylims=(0, 0.1Vmax),
	# 		  xlabel="Treated volume [L]", ylabel="Resin [mg]")

	# plot!(BV, BV*Cₓ₀mass[:Na]*1e-3; lab="", lw=1, c=:black, ls=:dashdot)

	# hline!([qe[i]*Vres for i in 1:N]; lab="", ls=:dot, lw=1.5, c=col[N:-1:2])
	
	# plot!(BV, [integrate(sol[x][2:end], sol(:, sol[x][2:end], dv=q(t,x)[i]); dims=1)* prob42S.ps[Cref]*Mi[c[i]]*1e3/Zi[c[i]] for i in 2:N] * Vres;
	# 	  lab="", c=permutedims(col[2:N]), lw=2)

	# plot!(BVv, [integrate(solv[x], solv(:, solv[x], dv=q(t,x)[i]); dims=1)* prob42V.ps[Cref]*Mi[c[i]]*1e3/Zi[c[i]] for i in 2:N] * Vres;
	# 	  lab="", c=permutedims(col[2:N]), lw=2, ls=:dot)
	

	figure = plot(p1, p2, p5,
				  layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				  suptitle="Prediction: <10% error in model parameters",
				  dpi=300, size=(680, 800))
	figure |> x-> aside(x; v_offset=500)
end

# ╔═╡ 3da033f0-5b5c-4b73-b485-5102eb8fe366
md"""
## Regeneration
"""

# ╔═╡ acedb039-3ff6-4d2a-9da3-16492f6e9c30
md"""
Variable elution:
"""

# ╔═╡ 0d2737e3-4878-46f1-bf1c-4a361e455eff
idxM = findfirst(x -> x >= 0.05, solv(:, 1, dv=C(t,x)[2]))

# ╔═╡ ad25ec42-9ba5-440e-8cd4-54521f5e88cf
solregv, regv = let
	Random.seed!(seed_model);

	prob = prob42V; 	reg = reg42S

	Cₓ₀reg = Dict([i => zero_eq for i in c]); Cₓ₀reg[:H] = 1.4
	reg = updateIC!(solv, reg, solv[t][idxM]; reverse=false)
	
	p = p0 .+ CI .* randn(length(p0))
	
	update!(pars, reg,	Q => 25/60*1e-6, # * 0.1 * sin(randn()*t), 
			Cx0 => [Cₓ₀reg[i] for i in c], 
			qmax => p[1], K => vcat(1.0, p[2:N]), Cref => prob.ps[Cref])
	
	@time solreg = 
		solve(reg, Rodas5P(), saveat=0.1, progress=true, progress_steps=1)
	solreg, reg
end;

# ╔═╡ fb263e2f-7dba-42dc-bb6c-5cbac8779fc7
md"""
Model with errors:
"""

# ╔═╡ 38c93997-7c9e-4ed7-a05c-95f4b86c083c
idxP = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.5) for j in 1:length(sol[t])])

# ╔═╡ 58071477-2e5f-4690-b16f-ddf407d208d7
idxC1 = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.975) for j in 1:length(sol[t])])

# ╔═╡ d1328678-5a0f-4b27-9778-680f974322d0
idxC2 = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.025) for j in 1:length(sol[t])])

# ╔═╡ 90ada351-91f6-4de4-a79b-52927a82fe6c
end5 = findfirst(x -> x >= 0.05, [quantile(CMC[j,:,2], 0.5) for j in 1:length(sol[t])])

# ╔═╡ bd158374-97ef-4819-8606-24a4414984f8
CMC[end5, :, 2] |> x-> quantile(x, 0.025)

# ╔═╡ b69d5017-0842-47a0-9e21-dcb4c173abe7
sol[t][end5]*prob42S.ps[τ]*prob42S.ps[Q]*1e3

# ╔═╡ 6948180c-0638-4014-a6f4-3ebceb291b54
solregm, regm = let
	Random.seed!(seed_model);

	prob = prob42S; 	reg = reg42S
	
	Cₓ₀reg = Dict([i => zero_eq for i in c]); Cₓ₀reg[:H] = 1.4
	reg = updateIC!(sol, reg, sol[t][end5]; reverse=false)
	
	update!(pars, reg, Q => 25/60*1e-6, Cx0 => [Cₓ₀reg[i] for i in c], 
			qmax => p0[1], Cref => prob.ps[Cref], K => vcat(1., p0[2:N]),
	)
	@time solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
	solreg, reg
end;

# ╔═╡ 5a565526-013a-4ad5-ac51-bb581f4c67c6
function MonteCarloReg(p)

	reg = regm

	Cₓ₀reg = Dict([i => zero_eq for i in c]); Cₓ₀reg[:H] = 1.4
	update!(pars, reg, Q => 25/60/1e6, Cx0 => [Cₓ₀reg[i] for i in c], 
			qmax => p[1], K => vcat(1.0, p[2:N]); verbose=false)
	
	return sol = solve(reg, Rodas5P(), saveat=0.01)
end

# ╔═╡ d51bad9f-09c7-45f3-b580-f75419058fd7
MCR, samplereg = let
	Np = length(p0)
	samples = 50
	CI = [0.34, 5.44, 7.56, 8.83, 10.7]/100 .* p0 
	
	MC = zeros(Np, samples)
	for i in 1:Np
		MC[i, :] = QuasiMonteCarlo.sample(samples, 0., 1., 
										 truncated(Normal(p0[i], CI[i]/2), 
												   lower=0.))
	end
	MC, samples
end

# ╔═╡ 16900190-d8fb-4d6b-967b-3b1af9e89153
let
	MC = MCR
	α = 0.05
	p = [plot() for i in 2:N]
	for i in 2:N
		p[i-1] = histogram(MC[i,:]; normalize=:probability,
						 c=cp[i-1], lab="", title=L"$\alpha_%$(i)$")

		# calibrated
		vline!([p0[i]]; lab="", c=:black, ls=:solid, lw=3)
		vline!([p0[i]+CI[i]]; lab="", c=:black, lw=3)
		vline!([p0[i]-CI[i]]; lab="", c=:black, lw=3)

		# sampled
		vline!([mean(MC[i, :])]; lab="", c=:red, ls=:dot, lw=2)
		vline!([quantile(MC[i, :], 1-α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
		vline!([quantile(MC[i, :], α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
	end

	# resin capacity
	i = 1
	p5 = plot()
	histogram!(MC[i,:]; normalize=:probability,
					 c=:grey, lab="", title=L"$q_{max}$")
	
	vline!([p0[i]]; lab="", c=:black, ls=:solid, lw=3)
	vline!([p0[i]+CI[i]]; lab="", c=:black, lw=3)
	vline!([p0[i]-CI[i]]; lab="", c=:black, lw=3)
	
	vline!([mean(MC[i, :])]; lab="", c=:red, ls=:dot, lw=2)
	vline!([quantile(MC[i, :], 1-α/2)]; lab="", c=:red, ls=:dashdot, lw=2)
	vline!([quantile(MC[i, :], α/2)]; lab="", c=:red, ls=:dashdot, lw=2)

	plot(p..., p5, layout=(3,2), dpi=300, size=(680, 600), suptitle="Sampling")
	plot!() |> aside
end

# ╔═╡ c18275a4-3ca0-4179-ae88-ed49391be2a8
CMCreg, CMCregM, CMCregC1, CMCregC2 = let
	α = 0.05
	Ns = samplereg
	Nt = length(solregm[t])
	CMC = zeros(length(solregm[t]), Ns, N)

	for n in 1:Ns
		p = MCR[:, n, :]
		solregm = MonteCarloReg(p)
		[CMC[:, n, i] = solregm(:, 1, dv=C(t,x)[i]) for i in 1:N]
	end

	CMCreg = CMC
	CMCregM = [quantile(CMC[j,:,i], 0.5) for j in 1:Nt, i in 1:N]
	CMCregC1 = [quantile(CMC[j,:,i], 1-α/2) for j in 1:Nt, i in 1:N]
	CMCregC2 = [quantile(CMC[j,:,i], α) for j in 1:Nt, i in 1:N]

	CMCreg, CMCregM, CMCregC1, CMCregC2
end

# ╔═╡ 5b449a34-071d-487b-a4e8-836f1c3afbdf
# let
# 	CMC = CMCreg
# 	Ns = size(CMC, 2) 
# 	α = 0.05

# 	BVreg = solregm[t]*regm.ps[τ]*regm.ps[Q]*1e3
	
# 	plot(dpi=300, xlims=(0, 3), grid=:off)

# 	# Monte Carlo simulations
# 	# plot!(BVreg, [CMC[:, :, i] for i in 2:N]; lab="", lw=0.5)

# 	# Average
# 	for i in 2:N
# 	plot!(BVreg, CMCregM[:, i]; 
# 		  lab="Mean ".*lab[i], lw=2, c=col[i])

# 	plot!(BVreg, CMCregM[:, i]; 
# 		  lab="", lw=2, c=col[i], fillcolor=col[i], fillalpha=0.7,
# 		  fillrange=(CMCregC1[:, i], CMCregC2[:, i]))
# 	end
	
# 	plot!() |> aside
# end

# ╔═╡ 94b01121-2ba6-4859-a7e2-df82d424b01f
let
	prob = prob42S;  	reg = reg42S; 	
	solx = regm.problem_type.discretespace.axies[x]

	BV = solregv[t]*reg.ps[τ]*reg.ps[Q]*1e3 # [L]
	BVreg = solregm[t]*reg.ps[τ]*reg.ps[Q]*1e3 # [L]
	pH = -log10.(solregm(solregm[t], 1, dv=C(t,x)[1])*prob.ps[Cref])
	
	BVmax = 2; 	Nt = length(BVreg); 
	pHmax = 8;  BVmin = 3
	α = 0.05

	# (a) pH
	p1 = plot(xlims=(0, BVmax), xticks=0:0.5:BVmax,
			  ylims=(-0.5, pHmax), ylabel="pH", grid=:off)

	# measured
	plot!(BV[BVmin:end], 
		  -log10.(solregv(solregv[t][BVmin:end], 1, dv=C(t,x)[1])*prob.ps[Cref]);
		  lab="", c=:grey, ls=:solid, lw=3)
	# model	
	plot!(BVreg[BVmin:end], pH[BVmin:end]; 
		  lab="", c=col[1], lw=2, ls=:dot,
		  fillcolor=col[1], fillalpha=0.4,
			  fillrange=(-log10.(CMCregC1[:, 1]*regm.ps[Cref]), 
						 -log10.(CMCregC2[:, 1]*regm.ps[Cref])))
	
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([-log10(Cₓ₀reg[:H])]; lab="", c=:black, ls=:dot, lw=1)
	annotate!(0.98BVmax, 2, ("final pH: $(round(-log10(Cₓ₀reg[:H]);digits=1))", 9, :right))

	
	# (b) Elution curve
	Cmax = 10
	p2 = plot(xlims=(0, BVmax), xticks=0:0.5:BVmax, grid=:off,
			  ylims=(-0.05, Cmax), ylabel="Concentration [g/L]")

	# Measured
	plot!(BV, [prob.ps[Cref]*Mi[c[i]]/Zi[c[i]]*
			   solregv(solregv[t], 1, dv=C(t,x)[i]) for i in 2:N]; 
		  lab="", c=permutedims(col[2:N]), lw=3, ls=:solid)

	# Model
	for i in 2:N
		plot!(BVreg, CMCregM[:, i]*regm.ps[Cref]*Mi[c[i]]/Zi[c[i]]; 
			  lab=lab[i], lw=2, c=col[i], ls=:dashdot, alpha=0.5, grid=:off,
			  fillcolor=col[i], fillalpha=0.4,
			  fillrange=(CMCregC1[:, i]*regm.ps[Cref]*Mi[c[i]]/Zi[c[i]], 
						 CMCregC2[:, i]*regm.ps[Cref]*Mi[c[i]]/Zi[c[i]]))
	end
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	# annotate!(BVreg[end2] + 0.3, 5, 
	# 		  ("BT Na: $(round(BC[2][end2,end];sigdigits=1)) g/L", 9, :left))


	# (c) Resin capacity
	cap = zeros(Nt, N)
	[cap[:,i] = integrate(solx, solregm(:, solx, dv=q(t,x)[i]); dims=1) * 100 for i in 1:N]
	
	p3 = plot(xlims=(0, BVmax), xticks=0:0.5:BVmax, 
			  xlabel="Regenerant volume [L]", grid=:off,
			  ylims=(0, 100), ylabel="Resin capacity, %")

	areaplot!(BVreg, cap; lab=permutedims(lab), c=permutedims(col), alpha=0.5, 
			  legendposition=:bottomright)

	plot!(BVreg, cap[:,1]; lab="", lw=3, c=col[1])
	vline!([BVreg[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	hline!([cap[end2]]; lab="", c=cp[N+2], ls=:dot, lw=3)
	scatter!([BVreg[end2]], [cap[end2]];
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	annotate!((BVreg[end2]<=3 ? 0.2 : -0.2)+BVreg[end2], (BVreg[end2]<=3 ? -10 : -18)+cap[end2], ("resin capacity at $(round(BVreg[end2];digits=1))L: $(round(cap[end2];digits=1))%", 9, (BVreg[end2]<=3 ? :left : :right )))

	figure = plot!(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				   suptitle="",
				   dpi=300, size=(680, 600))

	figure |> aside
end

# ╔═╡ 620542bd-441b-46bd-94e7-1ce10496a725
md"""
# Probabilistic model
"""

# ╔═╡ 5d4dc962-dfe9-4833-96b0-21c3ff9c9540
ypred = zeros(length(t4), N)

# ╔═╡ 67c19e62-f714-4c72-9ad0-8c59f226607d
@model function pprob(prob, data_t, data_y)
	Random.seed!(1234)

	Nt = length(data_t)
	yobs = zeros(Nt, N)
	ypred = zeros(Nt, N)
	

	# Prior distributions for parameters
	qm ~ truncated(Normal(p0[1], CI[1]), lower=0)

	α = Vector{Real}(undef, N)	
	α[1] = 1.0
	[α[i] ~ truncated(Normal(p0[i], CI[i]), lower=0) for i in 2:N]

	σ ~ truncated(Normal(0.0, 0.05), lower = 1e-8)

	# Solve model with sampled parameters
	update!(pars, prob, Q => 120/60/1e6, Cx0 => [Cₓ₀[i] for i in c], 
			qmax => qm, K => vcat(1.0, α[2:N]); verbose=false)
	
	sol = solve(prob, Rodas5P(), saveat=data_t*60/prob.ps[τ], 
				progress=true, progress_steps=1)

	if sol.retcode != ReturnCode.Success
        Turing.@addlogprob! -Inf
        return nothing
	end

	# Store data
	for i in 2:N
		ypred[:, i] = sol(:, 1, dv=C(t,x)[i])*
			prob.ps[Cref]*Mi[c[i]]*1e3/Zi[c[i]]
	end

	# Posterior
	for i in eachindex(Nt), j in 2:N
        yobs[i, j] ~ Normal(ypred[i, j], σ)
    end
end

# ╔═╡ b50e4e4e-1d8c-458e-9791-1981329a2994
model = pprob(prob42S, t4, y4[:, 2:N])

# ╔═╡ e1f0257a-10e8-4b6d-ab77-ae99a7da4767


# ╔═╡ 7e688179-e2bf-49fc-ada3-2d7ea06fa3c9
chain = sample(model, NUTS(), MCMCSerial(), 300, 3; progress=true)

# ╔═╡ a260c65a-b482-4be8-80e1-621a333de233


# ╔═╡ c978ee21-2747-47fe-b02e-2e604acc9eed
md"""
# Regeneration with RO concentrate
"""

# ╔═╡ e74dde90-96a3-4562-81e7-cb72a5bd2720


# ╔═╡ cd5eeff7-0829-4a28-837b-d444a65c7238
md"""
# Appendix
"""

# ╔═╡ e7506e41-70eb-4a9d-aa68-3250d7d95714
html"""
    <style>
    main {
        max-width: 1400px;
        margin-left: 30px;
        margin-right: 120px;
		padding-left: max(160px, 12%);
    	padding-right: max(10px, 30%);
    }
	pluto-output aside img {
    	max-width: 230%;
    	height: auto;
    	object-fit: contain;
	}
	</style>
    """

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
DomainSets = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
GlobalSensitivity = "af5da776-676b-467e-8baf-acd8249e4f0f"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MethodOfLines = "94925ecb-adb7-4558-8ed8-f975c56a0bf4"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
NumericalIntegration = "e7bfaba1-d571-5449-8927-abc22e82249b"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QuasiMonteCarlo = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Turing = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"

[compat]
Distributions = "~0.25.123"
DomainSets = "~0.7.16"
FiniteDiff = "~2.29.0"
GlobalSensitivity = "~2.8.0"
Interpolations = "~0.16.2"
LaTeXStrings = "~1.4.0"
MethodOfLines = "~0.11.9"
ModelingToolkit = "~10.25.0"
NumericalIntegration = "~0.3.4"
Optim = "~2.0.0"
OrdinaryDiffEq = "~6.102.1"
Plots = "~1.41.4"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.79"
QuasiMonteCarlo = "~0.3.5"
StatsPlots = "~0.15.8"
Turing = "~0.43.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "445b4b8a43ee5981540d76e4ae43482c25e7db3a"

[[deps.ADTypes]]
git-tree-sha1 = "f7304359109c768cf32dc5fa2d371565bb63b68a"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.21.0"
weakdeps = ["ChainRulesCore", "ConstructionBase", "EnzymeCore"]

    [deps.ADTypes.extensions]
    ADTypesChainRulesCoreExt = "ChainRulesCore"
    ADTypesConstructionBaseExt = "ConstructionBase"
    ADTypesEnzymeCoreExt = "EnzymeCore"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Dates", "Distributed", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "UUIDs"]
git-tree-sha1 = "8ac6182431567907e0d5170bcac6dd48fa541f78"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "5.15.1"

    [deps.AbstractMCMC.extensions]
    AbstractMCMCOnlineStatsExt = "OnlineStats"
    AbstractMCMCTensorBoardLoggerExt = "TensorBoardLogger"

    [deps.AbstractMCMC.weakdeps]
    OnlineStats = "a15396b6-48d5-5d58-9928-6d29437db91e"
    TensorBoardLogger = "899adc3e-224a-11e9-021f-63837185c80f"

[[deps.AbstractPPL]]
deps = ["AbstractMCMC", "Accessors", "BangBang", "DensityInterface", "JSON", "LinearAlgebra", "MacroTools", "OrderedCollections", "Random", "StatsBase"]
git-tree-sha1 = "e7be2de9646c1f61332de9f1e32c7dedf1e00831"
uuid = "7a57a42e-76ec-4ea3-a279-07e840d6d9cf"
version = "0.14.2"
weakdeps = ["Distributions"]

    [deps.AbstractPPL.extensions]
    AbstractPPLDistributionsExt = ["Distributions", "LinearAlgebra"]

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "MacroTools"]
git-tree-sha1 = "856ecd7cebb68e5fc87abecd2326ad59f0f911f3"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.43"

    [deps.Accessors.extensions]
    AxisKeysExt = "AxisKeys"
    IntervalSetsExt = "IntervalSets"
    LinearAlgebraExt = "LinearAlgebra"
    StaticArraysExt = "StaticArrays"
    StructArraysExt = "StructArrays"
    TestExt = "Test"
    UnitfulExt = "Unitful"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "7e35fca2bdfba44d797c53dfe63a51fabf39bfc0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.4.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AdvancedHMC]]
deps = ["AbstractMCMC", "ArgCheck", "DocStringExtensions", "LinearAlgebra", "LogDensityProblems", "LogDensityProblemsAD", "ProgressMeter", "Random", "Setfield", "Statistics", "StatsBase", "StatsFuns"]
git-tree-sha1 = "15e1bafe97eb0eeb77c8dae63cdcfa133986b254"
uuid = "0bf59076-c3b1-5ca4-86bd-e02cd72cde3d"
version = "0.8.3"

    [deps.AdvancedHMC.extensions]
    AdvancedHMCADTypesExt = "ADTypes"
    AdvancedHMCCUDAExt = "CUDA"
    AdvancedHMCComponentArraysExt = "ComponentArrays"
    AdvancedHMCMCMCChainsExt = "MCMCChains"
    AdvancedHMCOrdinaryDiffEqExt = "OrdinaryDiffEq"

    [deps.AdvancedHMC.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ComponentArrays = "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"

[[deps.AdvancedMH]]
deps = ["AbstractMCMC", "Distributions", "DocStringExtensions", "FillArrays", "LinearAlgebra", "LogDensityProblems", "Random", "Requires"]
git-tree-sha1 = "62ddbccf0ce5c26f8ef3cebe4bedef6b1599d616"
uuid = "5b7e9947-ddc0-4b3f-9b55-0d8042f74170"
version = "0.8.10"
weakdeps = ["DiffResults", "ForwardDiff", "MCMCChains", "StructArrays"]

    [deps.AdvancedMH.extensions]
    AdvancedMHForwardDiffExt = ["DiffResults", "ForwardDiff"]
    AdvancedMHMCMCChainsExt = "MCMCChains"
    AdvancedMHStructArraysExt = "StructArrays"

[[deps.AdvancedPS]]
deps = ["AbstractMCMC", "Distributions", "Random", "Random123", "Requires", "SSMProblems", "StatsFuns"]
git-tree-sha1 = "d92dd3fb4cc2748860ae8d5dd1d324cf0715a53b"
uuid = "576499cb-2369-40b2-a588-c64705576edc"
version = "0.7.2"
weakdeps = ["Libtask"]

    [deps.AdvancedPS.extensions]
    AdvancedPSLibtaskExt = "Libtask"

[[deps.AdvancedVI]]
deps = ["ADTypes", "Accessors", "ChainRulesCore", "DiffResults", "DifferentiationInterface", "Distributions", "DocStringExtensions", "FillArrays", "Functors", "LinearAlgebra", "LogDensityProblems", "Optimisers", "ProgressMeter", "Random", "StatsBase"]
git-tree-sha1 = "d69d7d9e1756fff9dd5d3fd26add46ee5ac62be4"
uuid = "b5ca4192-6429-45e5-a2d9-87aec30a685c"
version = "0.6.2"

    [deps.AdvancedVI.extensions]
    AdvancedVIBijectorsExt = ["Bijectors", "Optimisers"]
    AdvancedVIEnzymeExt = ["Enzyme", "ChainRulesCore"]
    AdvancedVIMooncakeExt = ["Mooncake", "ChainRulesCore"]
    AdvancedVIReverseDiffExt = ["ReverseDiff", "ChainRulesCore"]

    [deps.AdvancedVI.weakdeps]
    Bijectors = "76274a88-744f-5084-9051-94815aaf08c4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "f9e9a66c9b7be1ad7372bbd9b062d9230c30c5ce"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.5.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "d57bd3762d308bded22c3b82d033bff85f6195c6"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.4.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "libblastrampoline_jll"]
git-tree-sha1 = "7f54761502ff149a9d492e4acefe9805898e29b3"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.2+0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "d81ae5489e13bc03567d4fbbb06c546a5e53c857"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.22.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = ["CUDSS", "CUDA"]
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceMetalExt = "Metal"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "e0b47732a192dd59b9d079a06d04235e2f833963"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.12.2"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "4126b08903b777c88edf1754288144a0492c05ad"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.8"

[[deps.BangBang]]
deps = ["Accessors", "ConstructionBase", "InitialValues", "LinearAlgebra"]
git-tree-sha1 = "a49f9342fc60c2a2aaa4e0934f06755464fcf438"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.4.6"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTablesExt = "Tables"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BigCombinatorics]]
git-tree-sha1 = "a7d7d6194b3582455be7dd9495561626b3803cc5"
uuid = "7b33fef7-ef1d-5d54-bb67-cf96d4c8a166"
version = "0.3.6"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

[[deps.Bijectors]]
deps = ["AbstractPPL", "ArgCheck", "ChainRulesCore", "ChangesOfVariables", "DifferentiationInterface", "Distributions", "DocStringExtensions", "EnzymeCore", "FillArrays", "Functors", "InverseFunctions", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "MappedArrays", "Random", "Reexport", "Roots", "SparseArrays", "Statistics", "Test"]
git-tree-sha1 = "aa6d7edaf6ec7c432eb308020472c22624a63a39"
uuid = "76274a88-744f-5084-9051-94815aaf08c4"
version = "0.15.22"

    [deps.Bijectors.extensions]
    BijectorsDistributionsADExt = "DistributionsAD"
    BijectorsForwardDiffExt = "ForwardDiff"
    BijectorsLazyArraysExt = "LazyArrays"
    BijectorsMooncakeExt = "Mooncake"
    BijectorsReverseDiffChainRulesExt = ["ChainRules", "ReverseDiff"]
    BijectorsReverseDiffExt = "ReverseDiff"

    [deps.Bijectors.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    DistributionsAD = "ced4e74d-a319-5a8a-b0ac-84af2272839c"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "f21cfd4950cb9f0587d5067e69405ad2acd27b87"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.6"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "0f606a9894e2bcda541ceb82a91a13c5d450ed97"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "1.9.3"

    [deps.BlockArrays.extensions]
    BlockArraysAdaptExt = "Adapt"
    BlockArraysBandedMatricesExt = "BandedMatrices"

    [deps.BlockArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"

[[deps.BracketingNonlinearSolve]]
deps = ["CommonSolve", "ConcreteStructs", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "750f782fcc7e09283be7d8a7aa687a95e4911b60"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.6.2"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.BracketingNonlinearSolve.extensions]
    BracketingNonlinearSolveChainRulesCoreExt = ["ChainRulesCore", "ForwardDiff"]
    BracketingNonlinearSolveForwardDiffExt = "ForwardDiff"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "PrecompileTools", "Preferences", "Static"]
git-tree-sha1 = "f3a21d7fc84ba618a779d1ed2fcca2e682865bab"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.2.7"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.ChainRules]]
deps = ["Adapt", "ChainRulesCore", "Compat", "Distributed", "GPUArraysCore", "IrrationalConstants", "LinearAlgebra", "Random", "RealDot", "SparseArrays", "SparseInverseSubset", "Statistics", "StructArrays", "SuiteSparse"]
git-tree-sha1 = "3b704353e517a957323bd3ac70fa7b669b5f48d4"
uuid = "082447d4-558c-5d27-93f4-14fc19e9eca2"
version = "1.72.6"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e4c6a16e77171a5f5e25e9646617ab1c276c5607"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Chairmarks]]
deps = ["Printf", "Random"]
git-tree-sha1 = "9a49491e67e7a4d6f885c43d00bb101e6e5a434b"
uuid = "0ca39b1e-fe0b-4e98-acfc-b1656634c4de"
version = "1.3.1"
weakdeps = ["Statistics"]

    [deps.Chairmarks.extensions]
    StatisticsChairmarksExt = ["Statistics"]

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "3aa4bf1532aa2e14e0374c4fd72bed9a9d0d0f6c"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.10"
weakdeps = ["InverseFunctions", "Test"]

    [deps.ChangesOfVariables.extensions]
    ChangesOfVariablesInverseFunctionsExt = "InverseFunctions"
    ChangesOfVariablesTestExt = "Test"

[[deps.CloseOpenIntervals]]
deps = ["Static", "StaticArrayInterface"]
git-tree-sha1 = "05ba0d07cd4fd8b7a39541e31a7b0254704ea581"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.13"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"
weakdeps = ["StyledStrings"]

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["PrecompileTools"]
git-tree-sha1 = "351d6f4eaf273b753001b2de4dffb8279b100769"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.9.1"

[[deps.CommonSolve]]
git-tree-sha1 = "78ea4ddbcf9c241827e7035c3a03e2e456711470"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.6"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.CommonWorldInvalidations]]
git-tree-sha1 = "ae52d1c52048455e85a387fbee9be553ec2b68d0"
uuid = "f70d9fcc-98c5-4d4a-abd7-e4cdeebd8ca8"
version = "1.0.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "9d8a54ce4b17aa5bdce0ea5c34bc5e7c340d16ad"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.18.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ComplexityMeasures]]
deps = ["Combinatorics", "DelayEmbeddings", "Distances", "FFTW", "LinearAlgebra", "Neighborhood", "QuadGK", "Random", "Reexport", "SparseArrays", "SpecialFunctions", "StateSpaceSets", "StaticArrays", "Statistics", "StatsBase", "StatsFuns", "Wavelets"]
git-tree-sha1 = "df895f2a7fc136d2719eb6934f66d964b2a594af"
uuid = "ab4b797d-85ee-42ba-b621-05d793b346a2"
version = "3.8.5"

[[deps.CompositeTypes]]
git-tree-sha1 = "bce26c3dab336582805503bed209faab1c279768"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.4"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcreteStructs]]
git-tree-sha1 = "f749037478283d372048690eb3b5f92a79432b34"
uuid = "2569d6c7-a4a2-43d3-a901-331e8e4be471"
version = "0.2.3"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Copulas]]
deps = ["BigCombinatorics", "Combinatorics", "Distributions", "ForwardDiff", "HCubature", "InteractiveUtils", "LambertW", "LinearAlgebra", "LogExpFunctions", "MvNormalCDF", "PolyLog", "QuadGK", "Random", "Roots", "SpecialFunctions", "StatsBase", "StatsFuns", "WilliamsonTransforms"]
git-tree-sha1 = "b7ec542876a22344c0708f6b5be95df1fb7d9c4b"
uuid = "ae264745-0b69-425e-9d9d-cf662c5eec93"
version = "0.1.32"
weakdeps = ["Plots", "RecipesBase"]

    [deps.Copulas.extensions]
    CopulasPlotsExt = ["Plots", "RecipesBase"]

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"
weakdeps = ["OffsetArrays"]

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelayEmbeddings]]
deps = ["Distances", "Distributions", "LinearAlgebra", "Neighborhood", "Random", "Reexport", "StateSpaceSets", "Statistics", "StatsBase"]
git-tree-sha1 = "a37b1ad1e45b107e5176c8856114a15ee7fe9c6e"
uuid = "5732040d-69e3-5649-938a-b6b4f237613f"
version = "2.9.0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "BracketingNonlinearSolve", "ConcreteStructs", "DocStringExtensions", "EnzymeCore", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "0514bf55835444420ce81f8e32c5e2369d4af456"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.199.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseForwardDiffExt = ["ForwardDiff"]
    DiffEqBaseGTPSAExt = "GTPSA"
    DiffEqBaseGeneralizedGeneratedExt = "GeneralizedGenerated"
    DiffEqBaseMPIExt = "MPI"
    DiffEqBaseMeasurementsExt = "Measurements"
    DiffEqBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    DiffEqBaseMooncakeExt = "Mooncake"
    DiffEqBaseReverseDiffExt = "ReverseDiff"
    DiffEqBaseSparseArraysExt = "SparseArrays"
    DiffEqBaseTrackerExt = "Tracker"
    DiffEqBaseUnitfulExt = "Unitful"

    [deps.DiffEqBase.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    GeneralizedGenerated = "6b9d7cbe-bcb9-11e9-073f-15a7a543e2eb"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.DiffEqCallbacks]]
deps = ["ConcreteStructs", "DataStructures", "DiffEqBase", "DifferentiationInterface", "LinearAlgebra", "Markdown", "PrecompileTools", "RecipesBase", "RecursiveArrayTools", "SciMLBase", "StaticArraysCore"]
git-tree-sha1 = "f17b863c2d5d496363fe36c8d8535cc6a33c9952"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.12.0"
weakdeps = ["Functors"]

    [deps.DiffEqCallbacks.extensions]
    DiffEqCallbacksFunctorsExt = "Functors"

[[deps.DiffEqNoiseProcess]]
deps = ["CommonSolve", "DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "PoissonRandom", "QuadGK", "Random", "Random123", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "76fbb6985dda2aba32c97148540ad2c043e741e3"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.26.0"

    [deps.DiffEqNoiseProcess.extensions]
    DiffEqNoiseProcessOptimExt = "Optim"
    DiffEqNoiseProcessReverseDiffExt = "ReverseDiff"

    [deps.DiffEqNoiseProcess.weakdeps]
    Optim = "429524aa-4258-5aef-a3af-852621145aeb"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.DifferentiationInterface]]
deps = ["ADTypes", "LinearAlgebra"]
git-tree-sha1 = "7ae99144ea44715402c6c882bfef2adbeadbc4ce"
uuid = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
version = "0.7.16"

    [deps.DifferentiationInterface.extensions]
    DifferentiationInterfaceChainRulesCoreExt = "ChainRulesCore"
    DifferentiationInterfaceDiffractorExt = "Diffractor"
    DifferentiationInterfaceEnzymeExt = ["EnzymeCore", "Enzyme"]
    DifferentiationInterfaceFastDifferentiationExt = "FastDifferentiation"
    DifferentiationInterfaceFiniteDiffExt = "FiniteDiff"
    DifferentiationInterfaceFiniteDifferencesExt = "FiniteDifferences"
    DifferentiationInterfaceForwardDiffExt = ["ForwardDiff", "DiffResults"]
    DifferentiationInterfaceGPUArraysCoreExt = "GPUArraysCore"
    DifferentiationInterfaceGTPSAExt = "GTPSA"
    DifferentiationInterfaceMooncakeExt = "Mooncake"
    DifferentiationInterfacePolyesterForwardDiffExt = ["PolyesterForwardDiff", "ForwardDiff", "DiffResults"]
    DifferentiationInterfaceReverseDiffExt = ["ReverseDiff", "DiffResults"]
    DifferentiationInterfaceSparseArraysExt = "SparseArrays"
    DifferentiationInterfaceSparseConnectivityTracerExt = "SparseConnectivityTracer"
    DifferentiationInterfaceSparseMatrixColoringsExt = "SparseMatrixColorings"
    DifferentiationInterfaceStaticArraysExt = "StaticArrays"
    DifferentiationInterfaceSymbolicsExt = "Symbolics"
    DifferentiationInterfaceTrackerExt = "Tracker"
    DifferentiationInterfaceZygoteExt = ["Zygote", "ForwardDiff"]

    [deps.DifferentiationInterface.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DiffResults = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
    Diffractor = "9f5e2b26-1114-432f-b630-d3fe2085c51c"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastDifferentiation = "eb9bf01b-bf85-4b60-bf87-ee5de06c00be"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    GTPSA = "b27dd330-f138-47c5-815b-40db9dd9b6e8"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.DispatchDoctor]]
deps = ["MacroTools", "Preferences"]
git-tree-sha1 = "42cd00edaac86f941815fe557c1d01e11913e07c"
uuid = "8d63f2c5-f18a-4cf2-ba9d-b3f60fc568c8"
version = "0.4.28"
weakdeps = ["ChainRulesCore", "EnzymeCore"]

    [deps.DispatchDoctor.extensions]
    DispatchDoctorChainRulesCoreExt = "ChainRulesCore"
    DispatchDoctorEnzymeCoreExt = "EnzymeCore"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "fbcc7610f6d8348428f722ecbe0e6cfe22e672c6"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.123"
weakdeps = ["ChainRulesCore", "DensityInterface", "Test"]

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "c249d86e97a7e8398ce2068dce4c078a1c3464de"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.7.16"

    [deps.DomainSets.extensions]
    DomainSetsMakieExt = "Makie"
    DomainSetsRandomExt = "Random"

    [deps.DomainSets.weakdeps]
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DynamicPPL]]
deps = ["ADTypes", "AbstractMCMC", "AbstractPPL", "Accessors", "BangBang", "Bijectors", "Chairmarks", "Compat", "ConstructionBase", "DifferentiationInterface", "Distributions", "DocStringExtensions", "FillArrays", "InteractiveUtils", "LinearAlgebra", "LogDensityProblems", "MacroTools", "OrderedCollections", "PrecompileTools", "Printf", "Random", "Statistics", "Test"]
git-tree-sha1 = "14ded43f8a62568b17efc62d9f04491a2d4d8c9b"
uuid = "366bfd00-2699-11ea-058f-f148b4cae6d8"
version = "0.40.14"

    [deps.DynamicPPL.extensions]
    DynamicPPLEnzymeCoreExt = ["EnzymeCore"]
    DynamicPPLForwardDiffExt = ["ForwardDiff"]
    DynamicPPLMCMCChainsExt = ["MCMCChains"]
    DynamicPPLMarginalLogDensitiesExt = ["MarginalLogDensities"]
    DynamicPPLMooncakeExt = ["Mooncake", "DifferentiationInterface"]
    DynamicPPLReverseDiffExt = ["ReverseDiff"]

    [deps.DynamicPPL.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
    MarginalLogDensities = "f0c3360a-fb8d-11e9-1194-5521fd7ee392"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.DynamicPolynomials]]
deps = ["Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "Test"]
git-tree-sha1 = "3f50fa86c968fc1a9e006c07b6bc40ccbb1b704d"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.4"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "PrecompileTools", "TestItems", "Tricks"]
git-tree-sha1 = "faa92764e98ab0a1991b06c9b5bb3f1b70b4bf45"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.11.0"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EllipticalSliceSampling]]
deps = ["AbstractMCMC", "ArrayInterface", "Distributions", "Random", "Statistics"]
git-tree-sha1 = "e611b7fdfbfb5b18d5e98776c30daede41b44542"
uuid = "cad2338a-1db2-11e9-3401-43bc07c9ede2"
version = "2.0.0"

[[deps.EnumX]]
git-tree-sha1 = "7bebc8aad6ee6217c78c5ddcf7ed289d65d0263e"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.6"

[[deps.EnzymeCore]]
git-tree-sha1 = "990991b8aa76d17693a98e3a915ac7aa49f08d1a"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.18"
weakdeps = ["Adapt", "ChainRulesCore"]

    [deps.EnzymeCore.extensions]
    AdaptExt = "Adapt"
    EnzymeCoreChainRulesCoreExt = "ChainRulesCore"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterface", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "PrecompileTools", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "cc294ead6a85e975a8519dd4a0a6cb294eeb18d1"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.30.0"
weakdeps = ["StaticArrays"]

    [deps.ExponentialUtilities.extensions]
    ExponentialUtilitiesStaticArraysExt = "StaticArrays"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.ExproniconLite]]
git-tree-sha1 = "c13f0b150373771b0fdc1713c97860f8df12e6c2"
uuid = "55351af7-c7e9-48d6-89ff-24e801d99491"
version = "0.10.14"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FastBroadcast]]
deps = ["ArrayInterface", "LinearAlgebra", "Polyester", "Static", "StaticArrayInterface", "StrideArraysCore"]
git-tree-sha1 = "ab1b34570bcdf272899062e1a56285a53ecaae08"
uuid = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
version = "0.3.5"

[[deps.FastClosures]]
git-tree-sha1 = "acebe244d53ee1b461970f8910c235b259e772ef"
uuid = "9aa1b823-49e4-5ca5-8b0f-3971ec8bab6a"
version = "0.3.2"

[[deps.FastGaussQuadrature]]
deps = ["LinearAlgebra", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "0044e9f5e49a57e88205e8f30ab73928b05fe5b6"
uuid = "442a2c76-b920-505d-bb47-c5924d526838"
version = "1.1.0"

[[deps.FastPower]]
git-tree-sha1 = "862831f78c7a48681a074ecc9aac09f2de563f71"
uuid = "a4df4552-cc26-4903-aec0-212e50a0e84b"
version = "1.3.1"

    [deps.FastPower.extensions]
    FastPowerEnzymeExt = "Enzyme"
    FastPowerForwardDiffExt = "ForwardDiff"
    FastPowerMeasurementsExt = "Measurements"
    FastPowerMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    FastPowerMooncakeExt = "Mooncake"
    FastPowerReverseDiffExt = "ReverseDiff"
    FastPowerTrackerExt = "Tracker"

    [deps.FastPower.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2f979084d1e13948a3352cf64a25df6bd3b4dca3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.16.0"
weakdeps = ["PDMats", "SparseArrays", "StaticArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStaticArraysExt = "StaticArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FindFirstFunctions]]
deps = ["PrecompileTools"]
git-tree-sha1 = "27b495de668ccea58de6b06d6d13181396598ea0"
uuid = "64ca27bc-2ba2-4a57-88aa-44e436879224"
version = "1.8.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "9340ca07ca27093ff68418b7558ca37b05f8aeb1"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.29.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "eef4c86803f47dcb61e9b8790ecaa96956fdd8ae"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.2"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Functors]]
deps = ["Compat", "ConstructionBase", "LinearAlgebra", "Random"]
git-tree-sha1 = "60a0339f28a233601cb74468032b5c302d5067de"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.5.2"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ee0585b62671ce88e48d3409733230b401c9775c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.22"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "7dd7173f7129a1b6f84e0f03e0890cd1189b0659"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.22+0"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "a694e2a57394e409f7a11ee0977362a9fafcb8c7"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.6"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "6b4d2dc81736fe3980ff0e8879a9fc7c33c44ddf"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.2+0"

[[deps.Glob]]
git-tree-sha1 = "83cb0092e2792b9e3a865b6655e88f5b862607e2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.4.0"

[[deps.GlobalSensitivity]]
deps = ["Combinatorics", "ComplexityMeasures", "Copulas", "Distributions", "FFTW", "ForwardDiff", "KernelDensity", "LinearAlgebra", "Parameters", "QuasiMonteCarlo", "Random", "RecursiveArrayTools", "Statistics", "StatsBase", "ThreadsX", "Trapz"]
git-tree-sha1 = "e9f673523e5b6f731010ea4ba93d0d7f9db78ed6"
uuid = "af5da776-676b-467e-8baf-acd8249e4f0f"
version = "2.8.0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "031d63d09bd3e6e319df66bb466f5c3e8d147bee"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.13.4"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HCubature]]
deps = ["Combinatorics", "DataStructures", "LinearAlgebra", "QuadGK", "StaticArrays"]
git-tree-sha1 = "8ee627fb73ecba0b5254158b04d4745611b404a1"
uuid = "19dc6840-f33b-545b-b366-655c7e3ffd49"
version = "1.8.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "68c173f4f449de5b438ee67ed0c9c748dc31a2ec"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.28"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImplicitDiscreteSolve]]
deps = ["DiffEqBase", "OrdinaryDiffEqCore", "Reexport", "SciMLBase", "SimpleNonlinearSolve", "SymbolicIndexingInterface", "UnPack"]
git-tree-sha1 = "3e9ef0da0cabc23fc74e24cb233e184023f3b3ce"
uuid = "3263718b-31ed-49cf-8a0f-35a466e8af96"
version = "1.2.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "4c1acff2dc6b6967e7e750633c50bc3b8d83e617"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "65d505fa4c0d7072990d659ef3fc086eb6da8208"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.16.2"
weakdeps = ["ForwardDiff", "Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsForwardDiffExt = "ForwardDiff"
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "d966f85b3b7a8e49d034d27a189e9a4874b4391a"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.13"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "a779299d77cd080bf77b97535acecd73e1c5e5cb"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.17"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.Jieko]]
deps = ["ExproniconLite"]
git-tree-sha1 = "2f05ed29618da60c06a87e9c033982d4f71d0b6c"
uuid = "ae98c720-c025-4a4a-838c-29b094483192"
version = "0.2.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.JuliaFormatter]]
deps = ["CommonMark", "Glob", "JuliaSyntax", "PrecompileTools", "TOML"]
git-tree-sha1 = "eeb150b63d8bfdca04a0c367a8b95d3863fdc892"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "2.2.1"

[[deps.JuliaSyntax]]
git-tree-sha1 = "937da4713526b96ac9a178e2035019d3b78ead4a"
uuid = "70703baa-626e-46a2-a12c-08ffd08c73b4"
version = "0.4.10"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "Setfield", "StaticArrays", "SymbolicIndexingInterface"]
git-tree-sha1 = "c19bb0be5fe979b4476ab5531393b4fb0752ace5"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.21.1"

    [deps.JumpProcesses.extensions]
    JumpProcessesKernelAbstractionsExt = ["Adapt", "KernelAbstractions"]

    [deps.JumpProcesses.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTA", "Interpolations", "StatsBase"]
git-tree-sha1 = "4260cfc991b8885bf747801fb60dd4503250e478"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.11"

[[deps.Krylov]]
deps = ["LinearAlgebra", "Printf", "SparseArrays"]
git-tree-sha1 = "125d65fe5042faf078383312dd060adf11d90802"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.LambertW]]
deps = ["IrrationalConstants"]
git-tree-sha1 = "c6cd4fe119e21ebb1ca8adfaf32c989465c8d920"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "1.0.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LatticeRules]]
deps = ["Random"]
git-tree-sha1 = "7f5b02258a3ca0221a6a9710b0a0a2e8fb4957fe"
uuid = "73f95e8e-ec14-4e6a-8b18-0d2e271c4e55"
version = "0.0.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "95ba48564903b43b2462318aa243ee79d81135ff"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3acf07f130a76f87c041cfb2ff7d7284ca67b072"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.2+0"

[[deps.Libtask]]
deps = ["MistyClosures", "Test"]
git-tree-sha1 = "8a6ffa41b9c56aeb88c2962b12f35dd18ea2623f"
uuid = "6f1fad26-d15e-5dc8-ae53-837a1d7b8c9f"
version = "0.9.17"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2a7a12fc0a4e7fb773450d17975322aa77142106"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.2+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "9f7253c0574b4b585c8909232adb890930da980a"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.6"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "738bdcacfef25b3a9e4a39c28613717a6b23751e"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.6.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "1fdeafa25801ec2e0caa88fbe0afd7c41d3d087a"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.57.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveBLISExt = ["blis_jll", "LAPACK_jll"]
    LinearSolveBandedMatricesExt = "BandedMatrices"
    LinearSolveBlockDiagonalsExt = "BlockDiagonals"
    LinearSolveCUDAExt = "CUDA"
    LinearSolveCUDSSExt = "CUDSS"
    LinearSolveCUSOLVERRFExt = ["CUSOLVERRF", "SparseArrays"]
    LinearSolveCliqueTreesExt = ["CliqueTrees", "SparseArrays"]
    LinearSolveEnzymeExt = ["EnzymeCore", "SparseArrays"]
    LinearSolveFastAlmostBandedMatricesExt = "FastAlmostBandedMatrices"
    LinearSolveFastLapackInterfaceExt = "FastLapackInterface"
    LinearSolveForwardDiffExt = "ForwardDiff"
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockDiagonals = "0a1fb500-61f7-11e9-3c65-f5ef3456f9f0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    CUSOLVERRF = "a8cc9031-bad2-4722-94f5-40deabb4245c"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    FastAlmostBandedMatrices = "9d29842c-ecb8-4973-b1e9-a27b1157504e"
    FastLapackInterface = "29a986be-02c6-4525-aec4-84b980013641"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "d9625f27ded4ad726ceca7819394a4cc77ed25b3"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.2.0"

[[deps.LogDensityProblemsAD]]
deps = ["DocStringExtensions", "LogDensityProblems"]
git-tree-sha1 = "7b83f3ad0a8105f79a067cafbfd124827bb398d0"
uuid = "996a588d-648d-4e1f-a8f0-a84b347e47b1"
version = "1.13.1"

    [deps.LogDensityProblemsAD.extensions]
    LogDensityProblemsADADTypesExt = "ADTypes"
    LogDensityProblemsADDifferentiationInterfaceExt = ["ADTypes", "DifferentiationInterface"]
    LogDensityProblemsADEnzymeExt = "Enzyme"
    LogDensityProblemsADFiniteDifferencesExt = "FiniteDifferences"
    LogDensityProblemsADForwardDiffBenchmarkToolsExt = ["BenchmarkTools", "ForwardDiff"]
    LogDensityProblemsADForwardDiffExt = "ForwardDiff"
    LogDensityProblemsADReverseDiffExt = "ReverseDiff"
    LogDensityProblemsADTrackerExt = "Tracker"
    LogDensityProblemsADZygoteExt = "Zygote"

    [deps.LogDensityProblemsAD.weakdeps]
    ADTypes = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
    BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDifferences = "26cc04aa-876d-5657-8c51-4c34ba976000"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"
weakdeps = ["ChainRulesCore", "ChangesOfVariables", "InverseFunctions"]

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "DataAPI", "Dates", "Distributions", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "060d6bc7cf60e621dfd056ed2c1a2db1e68db0fe"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "7.7.0"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "DataStructures", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "2f464b68e84673727b4e4216a6254fba7da5cf4e"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.3.17"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MLJModelInterface]]
deps = ["InteractiveUtils", "REPL", "Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "c275fae2e693206b4527dd9d2382aa15359ef3ed"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.12.1"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.MappedArrays]]
git-tree-sha1 = "0ee4497a4e80dbd29c058fcee6493f5219556f40"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.3"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MaybeInplace]]
deps = ["ArrayInterface", "LinearAlgebra", "MacroTools"]
git-tree-sha1 = "54e2fdc38130c05b42be423e90da3bade29b74bd"
uuid = "bb5d69b7-63fc-4a16-80bd-7e42200c7bdb"
version = "0.1.4"
weakdeps = ["SparseArrays"]

    [deps.MaybeInplace.extensions]
    MaybeInplaceSparseArraysExt = "SparseArrays"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.MethodOfLines]]
deps = ["Combinatorics", "DiffEqBase", "DomainSets", "IfElse", "Interpolations", "Latexify", "LinearAlgebra", "ModelingToolkit", "OrdinaryDiffEq", "PDEBase", "PrecompileTools", "RuntimeGeneratedFunctions", "SciMLBase", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "TermInterface"]
git-tree-sha1 = "ec5b791da14e103fe0c0d1f27ee44d7a4f8e0293"
uuid = "94925ecb-adb7-4558-8ed8-f975c56a0bf4"
version = "0.11.9"

[[deps.MicroCollections]]
deps = ["Accessors", "BangBang", "InitialValues"]
git-tree-sha1 = "44d32db644e84c75dab479f1bc15ee76a1a3618f"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.2.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.MistyClosures]]
git-tree-sha1 = "d1a692e293c2a0dc8fda79c04cad60582f3d4de3"
uuid = "dbe65cb8-6be2-42dd-bbc5-4196aaced4f4"
version = "2.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.ModelingToolkit]]
deps = ["ADTypes", "AbstractTrees", "ArrayInterface", "BlockArrays", "ChainRulesCore", "Combinatorics", "CommonSolve", "Compat", "ConstructionBase", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DiffEqNoiseProcess", "DiffRules", "DifferentiationInterface", "Distributed", "Distributions", "DocStringExtensions", "DomainSets", "DynamicQuantities", "EnumX", "ExprTools", "FindFirstFunctions", "ForwardDiff", "FunctionWrappers", "FunctionWrappersWrappers", "Graphs", "ImplicitDiscreteSolve", "InteractiveUtils", "JuliaFormatter", "JumpProcesses", "Latexify", "Libdl", "LinearAlgebra", "MLStyle", "Moshi", "NaNMath", "OffsetArrays", "OrderedCollections", "OrdinaryDiffEqCore", "PrecompileTools", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SCCNonlinearSolve", "SciMLBase", "SciMLPublic", "SciMLStructures", "Serialization", "Setfield", "SimpleNonlinearSolve", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics", "URIs", "UnPack", "Unitful"]
git-tree-sha1 = "57cab9a2801c08138b50fa6bb95622de330e9274"
uuid = "961ee093-0014-501f-94e3-6117800e7a78"
version = "10.25.0"

    [deps.ModelingToolkit.extensions]
    MTKBifurcationKitExt = "BifurcationKit"
    MTKCasADiDynamicOptExt = "CasADi"
    MTKDeepDiffsExt = "DeepDiffs"
    MTKFMIExt = "FMI"
    MTKInfiniteOptExt = "InfiniteOpt"
    MTKLabelledArraysExt = "LabelledArrays"
    MTKPyomoDynamicOptExt = "Pyomo"

    [deps.ModelingToolkit.weakdeps]
    BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665"
    CasADi = "c49709b8-5c63-11e9-2fb2-69db5844192f"
    DeepDiffs = "ab62b9b5-e342-54a8-a765-a90f495de1a6"
    FMI = "14a09403-18e3-468f-ad8a-74f8dda2d9ac"
    InfiniteOpt = "20393b10-9daf-11e9-18c9-8db751c92c57"
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    Pyomo = "0e8e1daf-01b5-4eba-a626-3897743a3816"

[[deps.Moshi]]
deps = ["ExproniconLite", "Jieko"]
git-tree-sha1 = "53f817d3e84537d84545e0ad749e483412dd6b2a"
uuid = "2e0e35c7-a2e4-4343-998d-7ef72827ed2d"
version = "0.3.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MuladdMacro]]
git-tree-sha1 = "cac9cc5499c25554cba55cd3c30543cff5ca4fab"
uuid = "46d2c3a1-f734-5fdb-9937-b9b9aeba4221"
version = "0.2.4"

[[deps.MultivariatePolynomials]]
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "d38b8653b1cdfac5a7da3b819c0a8d6024f9a18c"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.13"
weakdeps = ["ChainRulesCore"]

    [deps.MultivariatePolynomials.extensions]
    MultivariatePolynomialsChainRulesCoreExt = "ChainRulesCore"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "7c3ff68a904d0f7404e5d2f7f5bc667934d8d616"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "22df8573f8e7c593ac205455ca088989d0a2c7a0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.6.7"

[[deps.MvNormalCDF]]
deps = ["Distributions", "FillArrays", "LinearAlgebra", "Primes", "Random", "StatsBase"]
git-tree-sha1 = "87ce0f8162ba3a9f1b6410d2180be4f04c0baac1"
uuid = "37188c8d-bc69-4638-b057-733e744175ec"
version = "0.3.2"

[[deps.NLSolversBase]]
deps = ["ADTypes", "DifferentiationInterface", "FiniteDiff", "LinearAlgebra"]
git-tree-sha1 = "b3f76b463c7998473062992b246045e6961a074e"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "8.0.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "e2c3bba08dd6dedfe17a17889131b885b8c082f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.27"

[[deps.Neighborhood]]
deps = ["Distances", "NearestNeighbors", "Random", "Test"]
git-tree-sha1 = "fdea60ca30d724e76cc3b3d90d7f9d29d3d5cab5"
uuid = "645ca80c-8b79-4109-87ea-e1f58159d116"
version = "0.2.4"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SciMLLogging", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "da41d5cd90c1118f5a8914392cf561263a1d09cb"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.15.0"

    [deps.NonlinearSolve.extensions]
    NonlinearSolveFastLevenbergMarquardtExt = "FastLevenbergMarquardt"
    NonlinearSolveFixedPointAccelerationExt = "FixedPointAcceleration"
    NonlinearSolveLeastSquaresOptimExt = "LeastSquaresOptim"
    NonlinearSolveMINPACKExt = "MINPACK"
    NonlinearSolveNLSolversExt = "NLSolvers"
    NonlinearSolveNLsolveExt = ["NLsolve", "LineSearches"]
    NonlinearSolvePETScExt = ["PETSc", "MPI", "SparseArrays"]
    NonlinearSolveSIAMFANLEquationsExt = "SIAMFANLEquations"
    NonlinearSolveSpeedMappingExt = "SpeedMapping"
    NonlinearSolveSundialsExt = "Sundials"

    [deps.NonlinearSolve.weakdeps]
    FastLevenbergMarquardt = "7a0df574-e128-4d35-8cbd-3d84502bf7ce"
    FixedPointAcceleration = "817d07cb-a79a-5c30-9a31-890123675176"
    LeastSquaresOptim = "0fc2ff8b-aaa3-5acd-a817-1944a5e08891"
    LineSearches = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
    MINPACK = "4854310b-de5a-5eb6-a2a5-c1dee2bd17f9"
    MPI = "da04e1cc-30fd-572f-bb4f-1f8673147195"
    NLSolvers = "337daf1e-9722-11e9-073e-8b9effe078ba"
    NLsolve = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    SIAMFANLEquations = "084e46ad-d928-497d-ad5e-07fa361a48c4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SpeedMapping = "f1835b91-879b-4a3f-a438-e4baacf14412"
    Sundials = "c3572dad-4567-51f8-b174-8c6c989267f4"

[[deps.NonlinearSolveBase]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "Compat", "ConcreteStructs", "DifferentiationInterface", "EnzymeCore", "FastClosures", "LinearAlgebra", "Markdown", "MaybeInplace", "Preferences", "Printf", "RecursiveArrayTools", "SciMLBase", "SciMLJacobianOperators", "SciMLLogging", "SciMLOperators", "SciMLStructures", "Setfield", "StaticArraysCore", "SymbolicIndexingInterface", "TimerOutputs"]
git-tree-sha1 = "9323e063ce69851b28e5c96539b2e2024f9a7e4c"
uuid = "be0214bd-f91f-a760-ac4e-3421ce2b2da0"
version = "2.11.2"

    [deps.NonlinearSolveBase.extensions]
    NonlinearSolveBaseBandedMatricesExt = "BandedMatrices"
    NonlinearSolveBaseChainRulesCoreExt = "ChainRulesCore"
    NonlinearSolveBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    NonlinearSolveBaseForwardDiffExt = "ForwardDiff"
    NonlinearSolveBaseLineSearchExt = "LineSearch"
    NonlinearSolveBaseLinearSolveExt = "LinearSolve"
    NonlinearSolveBaseMooncakeExt = "Mooncake"
    NonlinearSolveBaseReverseDiffExt = "ReverseDiff"
    NonlinearSolveBaseSparseArraysExt = "SparseArrays"
    NonlinearSolveBaseSparseMatrixColoringsExt = "SparseMatrixColorings"
    NonlinearSolveBaseTrackerExt = "Tracker"

    [deps.NonlinearSolveBase.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    LineSearch = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
    LinearSolve = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SparseMatrixColorings = "0a514795-09f3-496d-8182-132a7b665d35"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.NonlinearSolveFirstOrder]]
deps = ["ADTypes", "ArrayInterface", "CommonSolve", "ConcreteStructs", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLJacobianOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "df31d105d8e7254447256a44606f2a7e98b61aba"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "1.11.1"

[[deps.NonlinearSolveQuasiNewton]]
deps = ["ArrayInterface", "CommonSolve", "ConcreteStructs", "LinearAlgebra", "LinearSolve", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "SciMLOperators", "StaticArraysCore"]
git-tree-sha1 = "ade27e8e9566b6cec63ee62f6a6650a11cf9a2eb"
uuid = "9a2c21bd-3a47-402d-9113-8faf9a0ee114"
version = "1.12.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveQuasiNewton.extensions]
    NonlinearSolveQuasiNewtonForwardDiffExt = "ForwardDiff"

[[deps.NonlinearSolveSpectralMethods]]
deps = ["CommonSolve", "ConcreteStructs", "LineSearch", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "eafd027b5cd768f19bb5de76c0e908a9065ddd36"
uuid = "26075421-4e9a-44e1-8bd1-420ed7ad02b2"
version = "1.6.0"
weakdeps = ["ForwardDiff"]

    [deps.NonlinearSolveSpectralMethods.extensions]
    NonlinearSolveSpectralMethodsForwardDiffExt = "ForwardDiff"

[[deps.NumericalIntegration]]
deps = ["Interpolations", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9c6c03622c410853ac378cf17f71016d23a823a7"
uuid = "e7bfaba1-d571-5449-8927-abc22e82249b"
version = "0.3.4"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "117432e406b5c023f665fa73dc26e79ec3630151"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.17.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c9cbeda6aceffc52d8a0017e71db27c7a7c0beaf"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.5+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["ADTypes", "EnumX", "FillArrays", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "PositiveFactorizations", "Printf", "SparseArrays", "Statistics"]
git-tree-sha1 = "e4f98846b70ef55e111ac8c40add135256c0cc47"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "2.0.0"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Optimisers]]
deps = ["ChainRulesCore", "ConstructionBase", "Functors", "LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "36b5d2b9dd06290cd65fcf5bdbc3a551ed133af5"
uuid = "3bd65402-5787-11e9-1adc-39752487f4e2"
version = "0.4.7"

    [deps.Optimisers.extensions]
    OptimisersAdaptExt = ["Adapt"]
    OptimisersEnzymeCoreExt = "EnzymeCore"
    OptimisersReactantExt = "Reactant"

    [deps.Optimisers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Reactant = "3c362404-f566-11ee-1572-e11a4b42c853"

[[deps.Optimization]]
deps = ["ADTypes", "ArrayInterface", "ConsoleProgressMonitor", "DocStringExtensions", "LinearAlgebra", "Logging", "LoggingExtras", "OptimizationBase", "Printf", "Reexport", "SciMLBase", "SparseArrays", "TerminalLoggers"]
git-tree-sha1 = "5366b739224064a03399be1a21a18e4553631181"
uuid = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
version = "5.5.1"

[[deps.OptimizationBase]]
deps = ["ADTypes", "ArrayInterface", "DifferentiationInterface", "DocStringExtensions", "FastClosures", "LinearAlgebra", "PDMats", "PrecompileTools", "Reexport", "SciMLBase", "SciMLLogging", "SparseArrays", "SparseConnectivityTracer", "SparseMatrixColorings", "SymbolicIndexingInterface"]
git-tree-sha1 = "642ea067d729f1f18795a434550b40de0a4fbbc8"
uuid = "bca83a33-5cc9-4baa-983d-23429ab6bcbb"
version = "5.1.1"

    [deps.OptimizationBase.extensions]
    OptimizationChainRulesCoreExt = "ChainRulesCore"
    OptimizationEnzymeExt = ["ChainRulesCore", "Enzyme"]
    OptimizationFiniteDiffExt = "FiniteDiff"
    OptimizationForwardDiffExt = "ForwardDiff"
    OptimizationMLDataDevicesExt = "MLDataDevices"
    OptimizationMLUtilsExt = "MLUtils"
    OptimizationMooncakeExt = "Mooncake"
    OptimizationReverseDiffExt = "ReverseDiff"
    OptimizationSymbolicAnalysisExt = "SymbolicAnalysis"
    OptimizationZygoteExt = "Zygote"

    [deps.OptimizationBase.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FiniteDiff = "6a86dc24-6348-571c-b903-95158fe2bd41"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLDataDevices = "7e8f7934-dd98-4c1a-8fe8-92b47a384d40"
    MLUtils = "f1d291b0-491e-4a28-83b9-f70985020b54"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SymbolicAnalysis = "4297ee4d-0239-47d8-ba5d-195ecdf594fe"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.OptimizationOptimJL]]
deps = ["Optim", "OptimizationBase", "Reexport", "SciMLBase", "SparseArrays"]
git-tree-sha1 = "3e6aa79dd37700e32f95e9b312f0742343cd9437"
uuid = "36348300-93cb-4f02-beb5-3c3902f8871e"
version = "0.4.13"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.OrdinaryDiffEq]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "CommonSolve", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "ExponentialUtilities", "FastBroadcast", "FastClosures", "FillArrays", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "InteractiveUtils", "LineSearches", "LinearAlgebra", "LinearSolve", "Logging", "MacroTools", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqAdamsBashforthMoulton", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqDefault", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqExplicitRK", "OrdinaryDiffEqExponentialRK", "OrdinaryDiffEqExtrapolation", "OrdinaryDiffEqFIRK", "OrdinaryDiffEqFeagin", "OrdinaryDiffEqFunctionMap", "OrdinaryDiffEqHighOrderRK", "OrdinaryDiffEqIMEXMultistep", "OrdinaryDiffEqLinear", "OrdinaryDiffEqLowOrderRK", "OrdinaryDiffEqLowStorageRK", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqNordsieck", "OrdinaryDiffEqPDIRK", "OrdinaryDiffEqPRK", "OrdinaryDiffEqQPRK", "OrdinaryDiffEqRKN", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqSDIRK", "OrdinaryDiffEqSSPRK", "OrdinaryDiffEqStabilizedIRK", "OrdinaryDiffEqStabilizedRK", "OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "Polyester", "PreallocationTools", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "SimpleUnPack", "SparseArrays", "Static", "StaticArrayInterface", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "89cd4e81d7a668f8858fba6779212f41a0360260"
uuid = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
version = "6.102.1"

[[deps.OrdinaryDiffEqAdamsBashforthMoulton]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqLowOrderRK", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "09aae1486c767caa6bce9de892455cbdf5a6fbc8"
uuid = "89bda076-bce5-4f1c-845f-551c83cdda9a"
version = "1.5.0"

[[deps.OrdinaryDiffEqBDF]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqSDIRK", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "TruncatedStacktraces"]
git-tree-sha1 = "ce8db53fd1e4e41c020fd53961e7314f75e4c21c"
uuid = "6ad6398a-0878-4a85-9266-38940aa047c8"
version = "1.10.1"

[[deps.OrdinaryDiffEqCore]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "DataStructures", "DiffEqBase", "DocStringExtensions", "EnumX", "FastBroadcast", "FastClosures", "FastPower", "FillArrays", "FunctionWrappersWrappers", "InteractiveUtils", "LinearAlgebra", "Logging", "MacroTools", "MuladdMacro", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleUnPack", "Static", "StaticArrayInterface", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "4b68f9ca0cfa68cb9ee544df96391d47ca0e62a9"
uuid = "bbf590c4-e513-4bbe-9b18-05decba2e5d8"
version = "1.36.0"

    [deps.OrdinaryDiffEqCore.extensions]
    OrdinaryDiffEqCoreEnzymeCoreExt = "EnzymeCore"
    OrdinaryDiffEqCoreMooncakeExt = "Mooncake"

    [deps.OrdinaryDiffEqCore.weakdeps]
    EnzymeCore = "f151be2c-9106-41f4-ab19-57ee4f262869"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"

[[deps.OrdinaryDiffEqDefault]]
deps = ["ADTypes", "DiffEqBase", "EnumX", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqBDF", "OrdinaryDiffEqCore", "OrdinaryDiffEqRosenbrock", "OrdinaryDiffEqTsit5", "OrdinaryDiffEqVerner", "PrecompileTools", "Preferences", "Reexport", "SciMLBase"]
git-tree-sha1 = "7d5ddeee97e1bdcc848f1397cbc3d03bd57f33e7"
uuid = "50262376-6c5a-4cf5-baba-aaf4f84d72d7"
version = "1.8.0"

[[deps.OrdinaryDiffEqDifferentiation]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "FunctionWrappersWrappers", "LinearAlgebra", "LinearSolve", "OrdinaryDiffEqCore", "SciMLBase", "SciMLOperators", "SparseMatrixColorings", "StaticArrayInterface", "StaticArrays"]
git-tree-sha1 = "320b5f3e4e61ca0ad863c63c803f69973ba6efce"
uuid = "4302a76b-040a-498a-8c04-15b101fed76b"
version = "1.16.1"
weakdeps = ["SparseArrays"]

    [deps.OrdinaryDiffEqDifferentiation.extensions]
    OrdinaryDiffEqDifferentiationSparseArraysExt = "SparseArrays"

[[deps.OrdinaryDiffEqExplicitRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "4c0633f587395d7aaec0679dc649eb03fcc74e73"
uuid = "9286f039-9fbf-40e8-bf65-aa933bdc4db0"
version = "1.4.0"

[[deps.OrdinaryDiffEqExponentialRK]]
deps = ["ADTypes", "DiffEqBase", "ExponentialUtilities", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "3b81416ff11e55ea0ae7b449efc818256d9d450b"
uuid = "e0540318-69ee-4070-8777-9e2de6de23de"
version = "1.8.0"

[[deps.OrdinaryDiffEqExtrapolation]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastPower", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "9e1b11cf448a2c1bca640103c1c848a20aa2f967"
uuid = "becaefa8-8ca2-5cf9-886d-c06f3d2bd2c4"
version = "1.9.0"

[[deps.OrdinaryDiffEqFIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "FastGaussQuadrature", "FastPower", "LinearAlgebra", "LinearSolve", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "b968d66de3de5ffcf18544bc202ca792bad20710"
uuid = "5960d6e9-dd7a-4743-88e7-cf307b64f125"
version = "1.16.0"

[[deps.OrdinaryDiffEqFeagin]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "815b54211201ec42b8829e0275ab3c9632d16cbe"
uuid = "101fe9f7-ebb6-4678-b671-3a81e7194747"
version = "1.4.0"

[[deps.OrdinaryDiffEqFunctionMap]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "fe750e4b8c1b1b9e1c1319ff2e052e83ad57b3ac"
uuid = "d3585ca7-f5d3-4ba6-8057-292ed1abd90f"
version = "1.5.0"

[[deps.OrdinaryDiffEqHighOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "42096f72136078fa02804515f1748ddeb1f0d47d"
uuid = "d28bc4f8-55e1-4f49-af69-84c1a99f0f58"
version = "1.5.0"

[[deps.OrdinaryDiffEqIMEXMultistep]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Reexport", "SciMLBase"]
git-tree-sha1 = "a5dcd75959dada0005b1707a5ca9359faa1734ba"
uuid = "9f002381-b378-40b7-97a6-27a27c83f129"
version = "1.7.0"

[[deps.OrdinaryDiffEqLinear]]
deps = ["DiffEqBase", "ExponentialUtilities", "LinearAlgebra", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "925fc0136e8128fd19abf126e9358ec1f997390f"
uuid = "521117fe-8c41-49f8-b3b6-30780b3f0fb5"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowOrderRK]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "3cc4987c8e4725276b55a52e08b56ded4862917e"
uuid = "1344f307-1e59-4825-a18e-ace9aa3fa4c6"
version = "1.6.0"

[[deps.OrdinaryDiffEqLowStorageRK]]
deps = ["Adapt", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "e6bd0a7fb6643a57b06a90415608a81aaf7bd772"
uuid = "b0944070-b475-4768-8dec-fb6eb410534d"
version = "1.7.0"

[[deps.OrdinaryDiffEqNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "DiffEqBase", "FastBroadcast", "FastClosures", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MuladdMacro", "NonlinearSolve", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "PreallocationTools", "RecursiveArrayTools", "SciMLBase", "SciMLOperators", "SciMLStructures", "SimpleNonlinearSolve", "StaticArrays"]
git-tree-sha1 = "f59c1c07cfa674c1d3f5dd386c4274d9bc2be221"
uuid = "127b3ac7-2247-4354-8eb6-78cf4e7c58e8"
version = "1.15.0"

[[deps.OrdinaryDiffEqNordsieck]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqTsit5", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "c90aa7fa0d725472c4098096adf6a08266c2f682"
uuid = "c9986a66-5c92-4813-8696-a7ec84c806c8"
version = "1.4.0"

[[deps.OrdinaryDiffEqPDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "Polyester", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "9d599d2eafdf74ab26ea6bf3feb28183a2ade143"
uuid = "5dd0a6cf-3d4b-4314-aa06-06d4e299bc89"
version = "1.6.0"

[[deps.OrdinaryDiffEqPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "Reexport", "SciMLBase"]
git-tree-sha1 = "8e35132689133255be6d63df4190b5fc97b6cf2b"
uuid = "5b33eab2-c0f1-4480-b2c3-94bc1e80bda1"
version = "1.4.0"

[[deps.OrdinaryDiffEqQPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "63fb643a956b27cd0e33a3c6d910c3c118082e0f"
uuid = "04162be5-8125-4266-98ed-640baecc6514"
version = "1.4.0"

[[deps.OrdinaryDiffEqRKN]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "a31c41f9dbea7c7179c6e544c25c7e144d63868c"
uuid = "af6ede74-add8-4cfd-b1df-9a4dbb109d7a"
version = "1.5.0"

[[deps.OrdinaryDiffEqRosenbrock]]
deps = ["ADTypes", "DiffEqBase", "DifferentiationInterface", "FastBroadcast", "FiniteDiff", "ForwardDiff", "LinearAlgebra", "LinearSolve", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static"]
git-tree-sha1 = "f34bc2f58656843596d09a4c4de8c20724ebc2f1"
uuid = "43230ef6-c299-4910-a778-202eb28ce4ce"
version = "1.18.1"

[[deps.OrdinaryDiffEqSDIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "LinearAlgebra", "MacroTools", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "RecursiveArrayTools", "Reexport", "SciMLBase", "TruncatedStacktraces"]
git-tree-sha1 = "20caa72c004414435fb5769fadb711e96ed5bcd4"
uuid = "2d112036-d095-4a1e-ab9a-08536f3ecdbf"
version = "1.7.0"

[[deps.OrdinaryDiffEqSSPRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "StaticArrays"]
git-tree-sha1 = "3bce87977264916bd92455754ab336faec68bf8a"
uuid = "669c94d9-1f4b-4b64-b377-1aa079aa2388"
version = "1.7.0"

[[deps.OrdinaryDiffEqStabilizedIRK]]
deps = ["ADTypes", "DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "OrdinaryDiffEqDifferentiation", "OrdinaryDiffEqNonlinearSolve", "OrdinaryDiffEqStabilizedRK", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "75abe7462f4b0b2a2463bb512c8a5458bbd39185"
uuid = "e3e12d00-db14-5390-b879-ac3dd2ef6296"
version = "1.6.0"

[[deps.OrdinaryDiffEqStabilizedRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays"]
git-tree-sha1 = "7e94d3d1b3528b4bcf9e0248198ee0a2fd65a697"
uuid = "358294b1-0aab-51c3-aafe-ad5ab194a2ad"
version = "1.4.0"

[[deps.OrdinaryDiffEqSymplecticRK]]
deps = ["DiffEqBase", "FastBroadcast", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "RecursiveArrayTools", "Reexport", "SciMLBase"]
git-tree-sha1 = "e8dd5ab225287947016dc144a5ded1fb83885638"
uuid = "fa646aed-7ef9-47eb-84c4-9443fc8cbfa8"
version = "1.7.0"

[[deps.OrdinaryDiffEqTsit5]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "778c7d379265f17f40dbe9aaa6f6a2a08bc7fa3e"
uuid = "b1df2697-797e-41e3-8120-5422d3b24e4a"
version = "1.5.0"

[[deps.OrdinaryDiffEqVerner]]
deps = ["DiffEqBase", "FastBroadcast", "LinearAlgebra", "MuladdMacro", "OrdinaryDiffEqCore", "Polyester", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "Static", "TruncatedStacktraces"]
git-tree-sha1 = "185578fa7c38119d4318326f9375f1cba0f0ce53"
uuid = "79d7bb75-1356-48c1-b8c0-6832512096c2"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDEBase]]
deps = ["DomainSets", "ModelingToolkit", "SciMLBase", "SymbolicIndexingInterface", "SymbolicUtils", "Symbolics"]
git-tree-sha1 = "dc3f4a0efd478a391dbffd987c7e6b469c874678"
uuid = "a7812802-0625-4b9e-961c-d332478797e5"
version = "0.1.19"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e4cff168707d441cd6bf3ff7e4832bdf34278e4a"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.37"
weakdeps = ["StatsBase"]

    [deps.PDMats.extensions]
    StatsBaseExt = "StatsBase"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "063ef757a1e0e15af77bbe92be92da672793fd4e"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PoissonRandom]]
deps = ["LogExpFunctions", "Random"]
git-tree-sha1 = "67afbcbe9e184d6729a92a022147ed4cf972ca7b"
uuid = "e409e4f3-bfea-5376-8464-e040bb5c01ab"
version = "0.4.7"

[[deps.PolyLog]]
git-tree-sha1 = "4c1d0f01e381fe73fce9b31a96ea76b49762242f"
uuid = "85e3b03c-9856-11eb-0374-4dc1f8670e7f"
version = "2.6.0"
weakdeps = ["ChainRulesCore", "ForwardDiff"]

    [deps.PolyLog.extensions]
    PolyLogChainRulesExt = "ChainRulesCore"
    PolyLogForwardDiffExt = "ForwardDiff"

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "6f7cd22a802094d239824c57d94c8e2d0f7cfc7d"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.18"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "972089912ba299fba87671b025cd0da74f5f54f7"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.0"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "c05b4c6325262152483a1ecb6c69846d2e01727b"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.34"

    [deps.PreallocationTools.extensions]
    PreallocationToolsForwardDiffExt = "ForwardDiff"
    PreallocationToolsReverseDiffExt = "ReverseDiff"
    PreallocationToolsSparseConnectivityTracerExt = "SparseConnectivityTracer"

    [deps.PreallocationTools.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseConnectivityTracer = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "REPL", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "624de6279ab7d94fc9f672f0068107eb6619732c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "3.3.2"

    [deps.PrettyTables.extensions]
    PrettyTablesTypstryExt = "Typstry"

    [deps.PrettyTables.weakdeps]
    Typstry = "f0ed7684-a786-439e-b1e3-3b82803b501e"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "f0803bc1171e455a04124affa9c21bba5ac4db32"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.6"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "fbb92c6c56b34e1a2c4c36058f68f332bec840e7"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "34f7e5d2861083ec7596af8b8c092531facf2192"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.8.2+2"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "da7adf145cce0d44e892626e647f9dcbe9cb3e10"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.8.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "9eca9fc3fe515d619ce004c83c31ffd3f85c7ccf"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.8.2+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "8f528b0851b5b7025032818eb5abbeb8a736f853"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.8.2+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.QuasiMonteCarlo]]
deps = ["Accessors", "ConcreteStructs", "LatticeRules", "LinearAlgebra", "PrecompileTools", "Primes", "Random", "Sobol"]
git-tree-sha1 = "017a6731a1754173013b98f458adb631e368a9d2"
uuid = "8a4e6c94-4038-4cdc-81c3-7e6ffdb2a71b"
version = "0.3.5"
weakdeps = ["Distributions"]

    [deps.QuasiMonteCarlo.extensions]
    QuasiMonteCarloDistributionsExt = "Distributions"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Random123]]
deps = ["Random", "RandomNumbers"]
git-tree-sha1 = "dbe5fd0b334694e905cb9fda73cd8554333c46e2"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.7.1"

[[deps.RandomNumbers]]
deps = ["Random"]
git-tree-sha1 = "c6ec94d2aaba1ab2ff983052cf6a606ca5985902"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.6.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "LinearAlgebra", "PrecompileTools", "RecipesBase", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "81750fa442cc5cb3ab1ba89b1cdf8eda543e6905"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.47.0"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsFastBroadcastExt = "FastBroadcast"
    RecursiveArrayToolsForwardDiffExt = "ForwardDiff"
    RecursiveArrayToolsKernelAbstractionsExt = "KernelAbstractions"
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    RecursiveArrayToolsReverseDiffExt = ["ReverseDiff", "Zygote"]
    RecursiveArrayToolsSparseArraysExt = ["SparseArrays"]
    RecursiveArrayToolsStatisticsExt = "Statistics"
    RecursiveArrayToolsStructArraysExt = "StructArrays"
    RecursiveArrayToolsTablesExt = ["Tables"]
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    FastBroadcast = "7034ab61-46d4-4ed7-9d0f-46aef9175898"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "02d31ad62838181c1a3a5fd23a1ce5914a643601"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.3"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.ResettableStacks]]
deps = ["StaticArrays"]
git-tree-sha1 = "31c086583c92ab32d82ebef0d09fbcd6dd2c54a7"
uuid = "ae5879a3-cd67-5da8-be7f-38c6eb64a37b"
version = "1.2.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "5b3d50eb374cea306873b371d3f8d3915a018f0b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.9.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.Roots]]
deps = ["Accessors", "CommonSolve", "Printf"]
git-tree-sha1 = "8a433b1ede5e9be9a7ba5b1cc6698daa8d718f1d"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.2.10"

    [deps.Roots.extensions]
    RootsChainRulesCoreExt = "ChainRulesCore"
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"
    RootsUnitfulExt = "Unitful"

    [deps.Roots.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "2f609ec2295c452685d3142bc4df202686e555d2"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.16"

[[deps.SCCNonlinearSolve]]
deps = ["CommonSolve", "PrecompileTools", "Reexport", "SciMLBase", "SymbolicIndexingInterface"]
git-tree-sha1 = "733200f6b48ee91261cceefdc9237815dc782ee3"
uuid = "9dfe8606-65a1-4bb3-9748-cb89d1561431"
version = "1.9.0"
weakdeps = ["ChainRulesCore"]

    [deps.SCCNonlinearSolve.extensions]
    SCCNonlinearSolveChainRulesCoreExt = "ChainRulesCore"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SSMProblems]]
deps = ["AbstractMCMC", "Distributions", "Random"]
git-tree-sha1 = "cbf723e4c486375cf91236db53a7beefe8291951"
uuid = "26aad666-b158-4e64-9d35-0e672562fa48"
version = "0.6.1"

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "1b007ef25d4c18da29f956c2cad0a45f1c98f0f3"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.136.0"

    [deps.SciMLBase.extensions]
    SciMLBaseChainRulesCoreExt = "ChainRulesCore"
    SciMLBaseDifferentiationInterfaceExt = "DifferentiationInterface"
    SciMLBaseDistributionsExt = "Distributions"
    SciMLBaseEnzymeExt = "Enzyme"
    SciMLBaseForwardDiffExt = "ForwardDiff"
    SciMLBaseMLStyleExt = "MLStyle"
    SciMLBaseMakieExt = "Makie"
    SciMLBaseMeasurementsExt = "Measurements"
    SciMLBaseMonteCarloMeasurementsExt = "MonteCarloMeasurements"
    SciMLBaseMooncakeExt = "Mooncake"
    SciMLBasePartialFunctionsExt = "PartialFunctions"
    SciMLBasePyCallExt = "PyCall"
    SciMLBasePythonCallExt = "PythonCall"
    SciMLBaseRCallExt = "RCall"
    SciMLBaseReverseDiffExt = "ReverseDiff"
    SciMLBaseTrackerExt = "Tracker"
    SciMLBaseZygoteExt = ["Zygote", "ChainRulesCore"]

    [deps.SciMLBase.weakdeps]
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DifferentiationInterface = "a0c0ee7d-e4b9-4e03-894e-1c5f64a51d63"
    Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    MLStyle = "d8e11817-5142-5d16-987a-aa16d5891078"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    MonteCarloMeasurements = "0987c9cc-fe09-11e8-30f0-b96dd679fdca"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PartialFunctions = "570af359-4316-4cb7-8c74-252c00c2016b"
    PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
    PythonCall = "6099a3de-0909-46bc-b1f4-468b9a2dfc0d"
    RCall = "6f49c342-dc21-5d91-9882-a32aef131414"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SciMLJacobianOperators]]
deps = ["ADTypes", "ArrayInterface", "ConcreteStructs", "ConstructionBase", "DifferentiationInterface", "FastClosures", "LinearAlgebra", "SciMLBase", "SciMLOperators"]
git-tree-sha1 = "e96d5e96debf7f80a50d0b976a13dea556ccfd3a"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.12"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "7eebb9985e35b123e12025a3a2ad020cd6059f71"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.8.0"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "794c760e6aafe9f40dcd7dd30526ea33f0adc8b7"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.15.1"
weakdeps = ["SparseArrays", "StaticArraysCore"]

    [deps.SciMLOperators.extensions]
    SciMLOperatorsSparseArraysExt = "SparseArrays"
    SciMLOperatorsStaticArraysCoreExt = "StaticArraysCore"

[[deps.SciMLPublic]]
git-tree-sha1 = "0ba076dbdce87ba230fff48ca9bca62e1f345c9b"
uuid = "431bcebd-1456-4ced-9d72-93c2757fff0b"
version = "1.0.1"

[[deps.SciMLStructures]]
deps = ["ArrayInterface", "PrecompileTools"]
git-tree-sha1 = "607f6867d0b0553e98fc7f725c9f9f13b4d01a32"
uuid = "53ae85a6-f571-4167-b2af-e1d143709226"
version = "1.10.0"

[[deps.ScientificTypesBase]]
deps = ["InteractiveUtils"]
git-tree-sha1 = "e785eaa35a0f5518a388f9010e66fda64ea95ede"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.1.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ebe7e59b37c400f694f52b58c93d26201387da70"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.9"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleNonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "MaybeInplace", "NonlinearSolveBase", "PrecompileTools", "Reexport", "SciMLBase", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "315da09948861edbc6d18e066c08903487bb580d"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.10.0"

    [deps.SimpleNonlinearSolve.extensions]
    SimpleNonlinearSolveChainRulesCoreExt = "ChainRulesCore"
    SimpleNonlinearSolveReverseDiffExt = "ReverseDiff"
    SimpleNonlinearSolveTrackerExt = "Tracker"

    [deps.SimpleNonlinearSolve.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "be8eeac05ec97d379347584fa9fe2f5f76795bcb"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.5"

[[deps.SimpleUnPack]]
git-tree-sha1 = "58e6353e72cde29b90a69527e56df1b5c3d8c437"
uuid = "ce78b400-467f-4804-87d8-8f486da07d0a"
version = "1.1.0"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SparseConnectivityTracer]]
deps = ["ADTypes", "DocStringExtensions", "FillArrays", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "590b72143436e443888124aaf4026a636049e3f5"
uuid = "9f842d2f-2579-4b1d-911e-f412cf18a3f5"
version = "1.2.1"

    [deps.SparseConnectivityTracer.extensions]
    SparseConnectivityTracerChainRulesCoreExt = "ChainRulesCore"
    SparseConnectivityTracerLogExpFunctionsExt = "LogExpFunctions"
    SparseConnectivityTracerNNlibExt = "NNlib"
    SparseConnectivityTracerNaNMathExt = "NaNMath"
    SparseConnectivityTracerSpecialFunctionsExt = "SpecialFunctions"

    [deps.SparseConnectivityTracer.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    LogExpFunctions = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    NaNMath = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.SparseInverseSubset]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "52962839426b75b3021296f7df242e40ecfc0852"
uuid = "dc90abb0-5640-4711-901d-7e5b23a2fada"
version = "0.1.2"

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "6ed48d9a3b22417c765dc273ae3e1e4de035e7c8"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.23"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = "CUDA"
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f2685b435df2613e25fc10ad8c26dddb8640f547"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.6.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StateSpaceSets]]
deps = ["Distances", "LinearAlgebra", "Neighborhood", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "3cb1c581c396d534d86525cf9cb3a1e54f906f83"
uuid = "40b095a5-5852-4c12-98c7-d43bf788e795"
version = "2.5.3"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Static"]
git-tree-sha1 = "96381d50f1ce85f2663584c8e886a6ca97e60554"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.8.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eee1b9ad8b29ef0d936e3ec9838c7ec089620308"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.16"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "89f86d9376acd18a1a4fbef66a56335a3a7633b8"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.5.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "91f091a8716a6bb38417a6e6f274602a19aaa685"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.5.2"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "88cf3587711d9ad0a55722d339a013c4c56c5bbc"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.8"

[[deps.StrideArraysCore]]
deps = ["ArrayInterface", "CloseOpenIntervals", "IfElse", "LayoutPointers", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface", "ThreadingUtilities"]
git-tree-sha1 = "83151ba8065a73f53ca2ae98bc7274d817aa30f2"
uuid = "7792a7ef-975c-4747-a70f-980b88e8d1da"
version = "0.5.8"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "d05693d339e37d6ab134c5ab53c29fce5ee5d7d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.4"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "a2c37d815bf00575332b7bd0389f771cb7987214"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.7.2"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = ["GPUArraysCore", "KernelAbstractions"]
    StructArraysLinearAlgebraExt = "LinearAlgebra"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.SymbolicIndexingInterface]]
deps = ["Accessors", "ArrayInterface", "RuntimeGeneratedFunctions", "StaticArraysCore"]
git-tree-sha1 = "94c58884e013efff548002e8dc2fdd1cb74dfce5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.3.46"
weakdeps = ["PrettyTables"]

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

[[deps.SymbolicLimits]]
deps = ["SymbolicUtils"]
git-tree-sha1 = "f75c7deb7e11eea72d2c1ea31b24070b713ba061"
uuid = "19f23fe9-fdab-4a78-91af-e7b7767979c3"
version = "0.2.3"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "ArrayInterface", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "ExproniconLite", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicIndexingInterface", "TaskLocalValues", "TermInterface", "TimerOutputs", "Unityper"]
git-tree-sha1 = "a85b4262a55dbd1af39bb6facf621d79ca6a322d"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "3.32.0"

    [deps.SymbolicUtils.extensions]
    SymbolicUtilsLabelledArraysExt = "LabelledArrays"
    SymbolicUtilsReverseDiffExt = "ReverseDiff"

    [deps.SymbolicUtils.weakdeps]
    LabelledArrays = "2ee39098-c373-598a-b85f-a56591580800"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.Symbolics]]
deps = ["ADTypes", "ArrayInterface", "Bijections", "CommonWorldInvalidations", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "DynamicPolynomials", "LaTeXStrings", "Latexify", "Libdl", "LinearAlgebra", "LogExpFunctions", "MacroTools", "Markdown", "NaNMath", "OffsetArrays", "PrecompileTools", "Primes", "RecipesBase", "Reexport", "RuntimeGeneratedFunctions", "SciMLBase", "SciMLPublic", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArraysCore", "SymbolicIndexingInterface", "SymbolicLimits", "SymbolicUtils", "TermInterface"]
git-tree-sha1 = "7e53ce26142025d12565217082830f2f8caa5275"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "6.58.0"

    [deps.Symbolics.extensions]
    SymbolicsD3TreesExt = "D3Trees"
    SymbolicsForwardDiffExt = "ForwardDiff"
    SymbolicsGroebnerExt = "Groebner"
    SymbolicsLuxExt = "Lux"
    SymbolicsNemoExt = "Nemo"
    SymbolicsPreallocationToolsExt = ["PreallocationTools", "ForwardDiff"]
    SymbolicsSymPyExt = "SymPy"
    SymbolicsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Symbolics.weakdeps]
    D3Trees = "e3df1716-f71e-5df9-9e2d-98e193103c45"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    Groebner = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
    Lux = "b2108857-7c20-44ae-9111-449ecde12c47"
    Nemo = "2edaba10-b0f1-5616-af89-8c11ac63239a"
    PreallocationTools = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "f2c1efbc8f3a609aadf318094f8fc5204bdaf344"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TaskLocalValues]]
git-tree-sha1 = "67e469338d9ce74fc578f7db1736a74d93a49eb8"
uuid = "ed4db957-447d-4319-bfb6-7fa9ae7ecf34"
version = "0.1.3"

[[deps.TaylorDiff]]
deps = ["ChainRules", "ChainRulesCore", "MacroTools", "SymbolicUtils", "Symbolics"]
git-tree-sha1 = "8f0e85143bcf1101c046c4ea88804ee6ef400cc7"
uuid = "b36ab563-344f-407b-a36a-4f200bebf99c"
version = "0.3.5"

    [deps.TaylorDiff.extensions]
    TaylorDiffNNlibExt = ["NNlib"]
    TaylorDiffSFExt = ["SpecialFunctions"]
    TaylorDiffZygoteExt = ["Zygote"]

    [deps.TaylorDiff.weakdeps]
    NNlib = "872c559c-99b0-510c-b3b7-b6c96a88d5cd"
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TestItems]]
git-tree-sha1 = "42fd9023fef18b9b78c8343a4e2f3813ffbcefcb"
uuid = "1c621080-faea-4a02-84b6-bbd5e436b8fe"
version = "1.0.0"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "d969183d3d244b6c33796b5ed01ab97328f2db85"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.5"

[[deps.ThreadsX]]
deps = ["Accessors", "ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "SplittablesBase", "Transducers"]
git-tree-sha1 = "70bd8244f4834d46c3d68bd09e7792d8f571ef04"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.12"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "3748bd928e68c7c346b52125cf41fff0de6937d0"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.29"

    [deps.TimerOutputs.extensions]
    FlameGraphsExt = "FlameGraphs"

    [deps.TimerOutputs.weakdeps]
    FlameGraphs = "08572546-2f56-4bcf-ba4e-bab62c3a3f89"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Transducers]]
deps = ["Accessors", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "SplittablesBase", "Tables"]
git-tree-sha1 = "4aa1fdf6c1da74661f6f5d3edfd96648321dade9"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.85"

    [deps.Transducers.extensions]
    TransducersAdaptExt = "Adapt"
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

[[deps.Trapz]]
git-tree-sha1 = "79eb0ed763084a3e7de81fe1838379ac6a23b6a0"
uuid = "592b5752-818d-11e9-1e9a-2b8ca4a44cd1"
version = "2.0.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

[[deps.Turing]]
deps = ["ADTypes", "AbstractMCMC", "AbstractPPL", "Accessors", "AdvancedHMC", "AdvancedMH", "AdvancedPS", "AdvancedVI", "BangBang", "Bijectors", "Compat", "DataStructures", "DifferentiationInterface", "Distributions", "DocStringExtensions", "DynamicPPL", "EllipticalSliceSampling", "ForwardDiff", "Libtask", "LinearAlgebra", "LogDensityProblems", "MCMCChains", "Optimization", "OptimizationOptimJL", "OrderedCollections", "Printf", "Random", "Reexport", "SciMLBase", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "05e4dd6cfcaaa38e423c054a82aa3a43abd9a5b3"
uuid = "fce5fe82-541a-59a6-adf8-730c64b5f9a0"
version = "0.43.2"

    [deps.Turing.extensions]
    TuringDynamicHMCExt = "DynamicHMC"

    [deps.Turing.weakdeps]
    DynamicHMC = "bbc10e6e-7c05-544b-b16e-64fede858acb"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "57e1b2c9de4bd6f40ecb9de4ac1797b81970d008"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.28.0"
weakdeps = ["ConstructionBase", "ForwardDiff", "InverseFunctions", "LaTeXStrings", "Latexify", "NaNMath", "Printf"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    LatexifyExt = ["Latexify", "LaTeXStrings"]
    NaNMathExt = "NaNMath"
    PrintfExt = "Printf"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "25008b734a03736c41e2a7dc314ecb95bd6bbdb0"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.6"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wavelets]]
deps = ["DSP", "FFTW", "LinearAlgebra", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "d0ec97a100abbe47a5e9a02361841da49cce6029"
uuid = "29a6e085-ba6d-5f35-a997-948ac2efa89a"
version = "0.10.1"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "e9aeb174f95385de31e70bd15fa066a505ea82b9"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.7"

[[deps.WilliamsonTransforms]]
deps = ["Distributions", "Roots", "TaylorDiff"]
git-tree-sha1 = "87eacb5bed90fe0aa4773ba63d9de9d84c28dc11"
uuid = "48feb556-9bdd-43a2-8e10-96100ec25e22"
version = "0.2.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─b8d2511d-0865-4d5a-999e-bf40b2314dd8
# ╠═eee0fd64-96b6-4753-b262-9c9c4bfd81a1
# ╠═051f3f5c-a493-4e22-b77c-45222096cdd1
# ╠═9b1b7dbe-5fbc-4426-9d91-664a67938d34
# ╠═794e1837-5b56-42a1-af76-c9b566523ec4
# ╠═c21a3c39-d908-4444-9ebe-63f4569a0b31
# ╠═2862b00a-a238-4448-a7e5-e5f121e4b0d8
# ╠═14ddcaa9-2d26-40ea-a1cd-a8f1c144d272
# ╠═63850d62-693a-4651-851e-3e4665a950d1
# ╠═3844c32a-1212-421e-a25c-82a859aebb10
# ╠═00dd236d-444d-4227-9068-8d1cdb26f0dd
# ╠═d833b79d-ca76-43a9-b3fc-a3e5e6a989dd
# ╠═34266a3a-004e-495c-8eac-3df81d7f2232
# ╟─581c7094-6e13-4118-a3a7-6d4f9a2abff4
# ╟─8b071f14-cfcf-4ba6-a65d-85962b6dedc2
# ╟─57dec26e-3e14-45ec-83b3-34567fc16414
# ╠═3354bf0d-ff7e-4f47-9b07-c9dd7995ac52
# ╠═9a763df8-dfd8-42b3-9519-e2401f08bd7e
# ╠═00503082-7ea5-486a-b121-5e130c2fde4e
# ╠═791bdb1e-cea6-4f94-8bce-eb50fc38ddfb
# ╟─78e93282-dcd5-4517-98f7-86059c31a551
# ╠═aa4be22a-5dc8-4592-b4fa-b8fbe53447d7
# ╠═cfed0866-324d-4f74-9aae-83966de26969
# ╠═0668d48e-bad2-4d1d-88a5-d0c80e6bfc06
# ╠═565221c8-15a9-4482-bcf4-885009acbe41
# ╟─9ce4958a-914a-412a-a283-2f93d1912b9f
# ╠═79f86646-ecf8-47cd-972c-08b1dccd4268
# ╠═a03fc089-0c67-4e52-a283-b8553c074f53
# ╠═5faf1dfa-4b4c-44db-8159-b8629a1346d6
# ╠═ad35a9f3-f94b-4bb4-aae9-56f67a50f51a
# ╠═19d1c246-925b-4373-90ab-2f415416a5b6
# ╟─ed5ff484-b0c2-42dc-b91e-ae5cc19cca38
# ╠═75b26c32-0c4f-4583-8635-acdcfb8ddcfc
# ╠═5ef51740-4097-47be-b0ed-559bbeceac42
# ╠═aa0d5203-09c4-40fc-9792-20a56ebe1d4a
# ╠═81d0e150-ae98-4998-bc1f-c58e2c58e40b
# ╠═148119ed-49c8-4106-9c41-60788d3ff7d6
# ╠═215cf20b-c2b5-4e80-9584-9b9510d1df7a
# ╠═03178e82-d523-486c-ac5e-49ed20871dca
# ╠═1966ab11-75b6-4855-906e-16af2cd77f8d
# ╠═21954982-4d00-410a-a1e4-ba247a90d744
# ╠═6890ad0b-a50a-48de-95eb-d294bfffd6de
# ╠═1afa3007-1c4d-420d-a564-97954b703f66
# ╠═eb29fafd-b527-4ec1-a76e-292df208e5aa
# ╠═bbe1895a-cbce-4a0a-a7e6-147b2361b499
# ╠═fa5e1464-735d-4999-9893-0e50f4360f6f
# ╠═9accd87d-68dc-4b73-8eda-193255022ada
# ╠═7ee93e74-65c5-48c9-b456-6186bc81c71c
# ╠═3e07fe4c-131e-4d34-bbdd-836972049a16
# ╠═8dc31ff5-525e-40ac-8e50-2906c89eefff
# ╠═481221d6-960d-4f37-920c-7c8031e0c6b8
# ╠═f3fb5bb6-4500-4174-9c4a-636817625ce3
# ╠═9da4312a-00eb-400e-ba3d-39f03a25a9e2
# ╠═45679967-8e1b-4f3c-af73-41939bf99956
# ╠═0ee7fee1-5a45-45b1-ba8d-944a31414ad5
# ╠═d3c12be8-bfd3-4312-adbc-50dad911b61b
# ╠═cffe894d-d84e-4673-8c94-6d739410841a
# ╠═ddd71bc0-a30a-4114-8d5b-5400c796a343
# ╠═d9878823-74d3-4a41-a043-89c9f772d4d6
# ╠═cd46bc3f-334f-4dc5-859a-b581de1ebdcc
# ╠═66f76310-6861-45c7-8fb0-3bc5879d9f45
# ╠═9e8d5262-e0e4-40ef-a8be-5c38cc650747
# ╠═71ea145c-23fa-409f-a235-354110004312
# ╠═6a455ed0-8147-49a2-a729-93a1749e0f47
# ╠═215e8b8d-b113-4349-8522-8ebf395e54cd
# ╠═e9d33fe7-3bb0-4ac7-909c-a4aef8e8e72f
# ╠═65caa2ef-2e28-49f8-8e02-55a144abca8b
# ╠═2e6b1b04-294a-4395-8e69-0348f4442a5d
# ╠═0c894ba4-0753-4327-854a-cf4093a9020e
# ╠═b7e24c21-0d85-4c24-8000-4ba700cbd87e
# ╠═b4b71052-9609-4fb1-a541-e9e21919f275
# ╠═886f1479-c180-4f0f-bd22-1c1632382c27
# ╠═b33073d7-7676-4184-a426-0cb19bd75c7d
# ╠═398f3af2-fe1f-4f87-976b-5e9b2f041d58
# ╠═3b2bd4d5-1733-472e-ab2e-bf756488209a
# ╠═e6862433-a14d-41b1-8f0f-970c47f6f14f
# ╠═c94b8d18-d346-412a-ab1d-5cc7768b9d88
# ╟─8f918e3e-edf7-4e77-80db-0248a18b4983
# ╠═344e091a-fff8-47c0-8bfc-85c33c8bfcb5
# ╠═d2de844f-4124-4cf4-80c7-52300ded90d3
# ╟─1abc008b-4a60-44fc-b624-52c4840aebbb
# ╟─3fe71419-b2e4-4705-b5fd-655c8e8c69df
# ╠═8975d13c-bf27-4a78-812f-21122ff45aab
# ╠═01039316-bdd0-4b08-a75a-ad0d63b781b2
# ╠═cd0d9c44-9299-4c9d-a6d2-2f197d21c201
# ╟─77e04104-139c-4be3-862c-5d2655c7bdca
# ╠═303dab3b-b49e-4cc2-9c32-22d34c23a9a0
# ╠═6d9c39d7-b725-4722-bec7-684173c8e3d5
# ╠═72b1485f-2d7f-4a06-adce-4e05c37aced7
# ╠═e6ec0b6c-ae86-406f-8a18-b629ab39dda6
# ╠═ca9830b5-0b56-4ec7-a2a0-4841fd7f30c7
# ╟─cb90e8da-6d8a-49b7-9677-2dbbcf0745e2
# ╟─b3fadf2c-4a39-4dcc-9f6a-77a38584c000
# ╠═0eeaec2e-28ab-4826-973d-8a650d31373a
# ╠═96891cfa-bd75-43e3-8a42-41b019d2fbc8
# ╠═e320c41b-a06c-4953-b258-de6ffa983a6a
# ╠═34c81e22-e726-4183-b10a-f9ea21e0acc1
# ╠═a0660c1b-5497-43bb-8da4-d32b6c2cb759
# ╟─90ea266e-c8c8-4848-9b9e-58e867c71b23
# ╠═fd909090-e1be-42ee-ba88-45c37793dfbd
# ╠═e5238e11-583b-492c-9cc9-79d6c48583d0
# ╟─b2786e24-3119-42b9-9b09-eddb2c14c80d
# ╠═fc389dab-2d23-4923-845a-d895399e57da
# ╠═bb99d2a3-dc1e-4418-add5-47bebca9124a
# ╟─05ff089a-6d32-482b-85d2-c985feea65fd
# ╠═60499827-399c-4053-aa56-731f20012ad5
# ╟─aa96d243-3dd3-4a0d-bb57-c9494cb35b8e
# ╠═e7bcaaac-64c6-456b-b52b-f6006bf05f37
# ╠═e7d5f938-b446-4e71-93ae-ec77a55f18ed
# ╟─659597cf-5258-4997-80de-1a866988a3b4
# ╠═4bc7ce04-9f4e-4779-88fc-823d3ffb6ecb
# ╠═f98fc0ba-72c7-44a5-9b15-0b9b30e6b13e
# ╠═ef6206f7-305e-4199-b023-97fa55107dc4
# ╠═37698f81-a6a2-4843-91b2-09eb8ef902e9
# ╠═754f07f9-9fed-430e-804b-d40a292dded4
# ╠═37f89836-caf0-466a-abcf-d14041cc559c
# ╠═88a4f5ba-9857-45fc-87bd-9da63b1da4af
# ╟─48360ff2-9923-41a4-bac1-448c7f5a6131
# ╠═0d4433db-9a67-4a5a-bd40-4dbfa7bf0f43
# ╠═40ce2635-e918-455c-889f-4ab43a2dfa15
# ╠═a8de4251-3319-429f-aace-c18adccfb9ff
# ╠═0bf5f346-4046-433d-bbcc-7fb04ac0ef0d
# ╠═6542cd88-7816-4b78-a6ab-8e76967d98ff
# ╠═090782ac-3e17-453b-953d-911615b6a45d
# ╟─a8bb2a58-08f3-4495-a22e-f4fed33a936e
# ╠═c8b83cee-caa9-4087-8487-aa93a39c0e9b
# ╠═b50e0b99-434c-415d-9639-afdf5e8726d3
# ╠═b9aadf8f-5a77-46e5-986a-bb4c2a314335
# ╟─410d04f4-8158-41a1-aadc-59277cf260dd
# ╟─460c29de-39bf-4dcc-afb0-be20d117860a
# ╟─f174d441-2934-478e-beeb-3b2e50c796ca
# ╟─ad0e7db9-e293-4aa6-8e33-683a9ed05440
# ╠═df52ae4c-97d4-48f7-9a5d-ef281b26577d
# ╟─1fcd177b-dee2-4e56-8817-fa57bc4eea0c
# ╟─569fa407-e152-4c5b-aa8c-dd696ac2326d
# ╠═ce14146a-f9dc-425b-bc5d-dc628f33dc56
# ╟─44f56d11-a5a1-444e-826b-a86515ce78d9
# ╟─3a359f7e-919e-41aa-94e6-c1004fbe0a62
# ╟─e2973d09-defd-42cc-93f4-e55e38a3e05a
# ╟─67b555d4-edac-4731-b167-06e2265f1dcc
# ╟─13f5e12c-9d95-4ff0-89d8-cf5b84928a56
# ╟─509a586c-0be9-4a11-aec3-ac3e71cf29d7
# ╟─f596e2e6-81a1-4fd6-9a2e-d9ce5b03ec3a
# ╟─7e1c7fa6-e329-417d-95e2-ba2896a05cd2
# ╟─b0b839a1-1bd4-43cc-9f43-325189278f5f
# ╟─ae7704b0-ad6a-41b0-9f05-67bababe40d6
# ╟─47bbd2ab-0829-4136-84de-f98f7b3408a7
# ╟─9f1329e5-04b0-4395-983a-8dad2cf4bef1
# ╟─d7f657c3-e540-4e2d-933a-5a3534957015
# ╟─33e1d1ad-85f4-48ad-933b-a85f2203c300
# ╟─67c6d1c1-ba68-4966-a3cb-f23d38a9049a
# ╟─df88c1dc-ad21-4ef4-b41e-416ee3543a4a
# ╟─9ce03269-5734-40cf-a87f-7b9896d93bdd
# ╟─a9853f0c-6e44-42b2-96bb-24fc22164df0
# ╟─8eb84797-2805-4f94-911d-fedf2ca10fc6
# ╟─4e58cf1e-623a-4296-b3fa-a7dd4375c481
# ╟─4b8afbf2-6bc7-43b4-b5f1-0bd3e7854698
# ╟─f3b30039-70c0-4f8e-85f2-d94bb9a160f3
# ╟─54b90a93-992c-42bf-8672-61905b874c47
# ╟─4242c3d6-20e4-4a27-9998-6498dafaffc0
# ╟─dd1a1903-872f-42cb-b138-ebc044c40c14
# ╠═a0f7e2f3-2693-4241-a4c7-9af5617e04dd
# ╠═bff5109f-88ad-4269-8ce5-4a9810f4053e
# ╠═f1dbaf46-580f-4511-b012-6efb7f582230
# ╟─c07242d9-5015-4cce-a2d7-24d34200b6d1
# ╟─86754b3c-a9bf-4138-a85a-aae4a3e7c78d
# ╟─220c7e50-022e-45fc-a79f-d44596e89b52
# ╟─1464fa94-e1c6-4626-89fe-e28d10155d16
# ╠═4f5e7e54-ccf9-4067-8af4-fce16fa592f0
# ╟─e10178ba-d810-408c-a65c-0d206280e541
# ╟─8f74a1f1-7fa3-4672-ae77-51c904fe2cc2
# ╟─cc953d89-ebf3-4812-a9b0-f9e73e40f611
# ╠═d4b5e8f6-f288-47de-b05f-377181ce3326
# ╟─8b6f1c83-cc1d-466f-94e0-9f04ba0e7cbb
# ╟─e8c5ecfd-78a8-451a-94da-d26a12dbbe3a
# ╟─b14e753a-7c41-48d5-88c3-5f2dac364e47
# ╟─60125c5d-4997-43fd-a48e-f06373ec54b0
# ╟─f08b716a-cbbe-4c1d-b9b0-03e415a9d1ed
# ╠═2bf1ffa1-d687-4cca-8072-98908b0e4a4c
# ╟─d37bdba9-ac5b-4023-8aa2-8ca70300549d
# ╟─584c3a98-9bd7-4d46-9b6c-a098694e835d
# ╠═3fb3da00-ee78-4339-a458-b6f947f6900f
# ╠═d8ed6fd2-0ae6-4b8a-889b-87ce71580edb
# ╟─82c98042-7e70-4243-b8b4-30143322a077
# ╠═b2ca0f67-05b7-4a94-a679-5a649d002fc7
# ╠═97b01366-9ff8-40c3-b6dd-fee9510eb58f
# ╠═3d573f46-8def-4a6d-83b7-aa577996b079
# ╠═aaeeca5e-08f6-479d-a300-02c4132752df
# ╟─c5dabff0-fe48-42c9-b292-d3880fd2df1f
# ╟─66466e20-66d0-47dc-81fc-f2123d9bce09
# ╠═c8df5818-a945-4f40-a5db-db866b7397f5
# ╠═39f3f4a9-5513-4973-80b1-f07c152f7d28
# ╠═eaa12a3d-23e6-41d6-9e9c-913b9b2c0b82
# ╠═7c9d4d0b-e4ef-43e8-a36c-44ff8b90d0d4
# ╠═82ee92cf-f38e-4c02-a743-51b8a014cfd8
# ╠═7a7be574-484e-4e36-bf23-ce7c83b0ab1b
# ╠═e2cb876a-716c-47ce-944b-cb365ca99348
# ╠═7d787166-df1d-4ba0-83a5-036665ea061c
# ╠═788ae965-4a3f-48a3-a048-881531dc4f19
# ╠═335a82ca-8cb2-478d-88c5-f69f7e14fce0
# ╠═e135d615-b494-473d-bb39-7175670683b0
# ╠═71fd718e-7d75-4655-a8d4-fdce5290b557
# ╟─b03cc95b-c458-4700-92d9-d48458577a14
# ╠═5519a8d7-4abe-449d-83f4-d50ecce71494
# ╟─3fc5ad89-84b5-40b0-aab9-9ab817a5ffdb
# ╟─a470658e-444f-4cc2-81a8-0c0d7a7e03df
# ╟─3b85787c-c67f-4589-ae0e-59401f7b2751
# ╠═18ff9503-091f-4afa-9dd9-4c8df64051e7
# ╠═8ba53a75-e3d1-4149-a41e-e19fac5baaea
# ╟─ecd28db0-92db-4817-9307-9fad80dafa91
# ╟─a96a1387-8365-481a-ad90-10a867f7a30e
# ╟─78bff468-778b-45f6-9758-287adb4a477e
# ╟─44690f9e-76b0-47d9-8dc3-401ab191344e
# ╟─16187793-3049-4be1-a65b-1f42b147f02b
# ╟─208e50ff-84b1-4a50-ab1c-85b619b447fd
# ╠═04e38cc6-5d7b-41df-8855-d1b82fc768b9
# ╠═4d9d5d5f-e955-49b2-afe3-89f5acd7ce30
# ╠═136fc997-88c6-42d9-8116-aa63dfe29b78
# ╟─8ba71079-5d31-4dda-963a-22a481b68495
# ╟─b2735d20-d71f-414a-aa94-31c92b09ff6e
# ╠═72e5e318-eb42-465b-80ad-9ead7b888ba1
# ╟─6225b970-86c9-484d-89b3-3cf385f7af35
# ╠═68c844a4-4f7e-4424-8679-5102cade24b9
# ╟─5e436a4b-61dd-4234-ad56-51a11e665903
# ╠═7a550bb5-fec5-4918-bc26-1f3dd4e6118e
# ╟─547bb551-2604-4f8d-8e8b-35aba51e7a45
# ╟─65e470f2-91ef-4da6-a58f-41e10b011e5e
# ╟─57ab6988-0602-43a0-b9a1-bb88dc5846f8
# ╠═620403e0-8526-45d8-848d-235754c8b82b
# ╟─8786b6c6-77b0-4eed-adbc-4790be9ae284
# ╠═82d7e380-e1bf-43da-afdd-a32ffefa2af9
# ╟─12c5278e-3c18-4138-8ed4-1fa5a19a459f
# ╠═96808966-f7a1-473b-8601-0ee6cb649730
# ╠═d8473938-0a1c-4b46-9020-75306e06fa3c
# ╟─0243f891-343c-4dcc-a407-e09bf2fb3516
# ╟─2e4e1e1d-ef12-41dd-af3b-0ebaecd3b7ad
# ╟─4cb7c562-b581-4ea4-bf64-925ab9819597
# ╠═230ed788-9526-4182-969b-0d01c20992fe
# ╠═8966be55-99cb-4951-ba88-315a03aa7db0
# ╠═827ab335-bb6c-46f0-8f37-c27bde773ef6
# ╟─5d49d4f7-9647-494c-9db7-c4f5b9ea7908
# ╠═a1725264-ff20-42fd-b4ee-66c955c2fac4
# ╟─a70b158f-9b41-42ea-b600-f4c513968357
# ╠═b8337a8b-26d6-4bcb-a710-899d2f1f3799
# ╠═8c9d581c-e9ef-43fb-b8ce-c3f5bfa50f90
# ╟─e11cef45-c956-4cdb-9568-005b1c0e8cf2
# ╟─e760a66a-d418-4b13-af96-37232574b035
# ╟─225d9033-75fd-4c18-a35d-ed9e6306d2c3
# ╠═05fcd099-a261-419c-9d68-635ebecd2413
# ╠═55ad8e2c-4e30-4a70-8ec2-6d4424fc070b
# ╠═8a61660d-0ee0-41a3-9515-2afab8de9133
# ╟─ec9597bb-26db-4c1e-9896-ea221a4f7f9b
# ╟─96ea5cd6-3f76-4916-8698-4f37f0a3eec2
# ╟─d7b1b4f5-876d-4472-b69e-ceeec2c29d78
# ╠═3e1e978f-3588-4747-a6fe-ed165c9a99ea
# ╠═84d5df10-0cf3-4f1a-adf0-4e4d1a1b13cc
# ╟─2b7429de-c798-425a-bf37-8bc8bbaa67b3
# ╠═fc2358d9-b5cf-401a-ab0f-a68b08cb4359
# ╟─bf92682d-4e9e-459c-83e2-3bcfd3729360
# ╠═d24ddf3f-f5f6-4d6a-bbc7-e2a1ac30ba95
# ╟─deda98c1-ca29-40cd-9f7e-b9b2b5852afd
# ╟─4d53a2a4-20db-4147-98e7-fdc35d310b42
# ╠═1d4a3d18-71ca-464e-ab06-b70a0d89140c
# ╟─60c61a84-3bba-42e4-9aa9-b491f11f11ca
# ╠═76dd8dd9-f193-4e2e-a782-7c6f4059dc9b
# ╠═57c724f7-4272-44b0-82d2-edd8548627c2
# ╠═e9bea4f3-6d60-4b79-a184-b0ba57ae9f35
# ╠═c8e3f872-817a-4c6c-89d0-0707aa6f7e2d
# ╠═9e6eeeb1-dc27-4fa5-b5ea-ebb759c37915
# ╠═e2ea5ed4-db1b-40ec-a93d-8c9ea4f9828b
# ╠═60c14f0d-3f6d-42d9-b737-9a3a9cc8aca5
# ╠═a31d2322-a171-4343-b121-48d3bba29689
# ╠═a63ec162-12e7-42a4-9344-6576de0b50de
# ╠═61e3201e-788f-4846-9ef7-b320f4e67aa1
# ╠═33b03f13-636b-4038-8c75-42b777ad50e5
# ╠═4e310059-4f33-45ea-9d60-da1d7824cbd1
# ╠═6112953c-9566-4f08-b480-f9adcf96eb4a
# ╟─c9da859c-d33f-42d8-bf57-09f5d3fc7c71
# ╟─3a1cd15b-d101-4e92-8c3b-2339ec7d0cc5
# ╠═9da02a78-6118-47f6-aa32-0157e2cad312
# ╠═4c139460-3116-49c1-8726-40d1cb78bd93
# ╠═2c857135-9ab3-42c5-80c2-3edd5bfd155f
# ╠═91e0f021-db31-4f7a-8cc1-fea9615d681d
# ╠═97005938-587d-4a58-83c0-0258cb9c8b74
# ╟─7f3eb011-a1ae-4b24-a4e7-ab08043f80db
# ╠═4353e933-4f0a-4165-a0ab-6e4d93f52f2f
# ╠═8142e92b-39d3-4dca-8fe0-588d0e70aba8
# ╠═068dfd04-3923-4b39-a952-1bf9206a74cb
# ╠═3b04b7e9-a9ae-4069-8edf-fc374e9e696e
# ╟─3da033f0-5b5c-4b73-b485-5102eb8fe366
# ╟─acedb039-3ff6-4d2a-9da3-16492f6e9c30
# ╠═0d2737e3-4878-46f1-bf1c-4a361e455eff
# ╠═ad25ec42-9ba5-440e-8cd4-54521f5e88cf
# ╟─fb263e2f-7dba-42dc-bb6c-5cbac8779fc7
# ╠═38c93997-7c9e-4ed7-a05c-95f4b86c083c
# ╠═58071477-2e5f-4690-b16f-ddf407d208d7
# ╠═d1328678-5a0f-4b27-9778-680f974322d0
# ╠═bd158374-97ef-4819-8606-24a4414984f8
# ╠═b69d5017-0842-47a0-9e21-dcb4c173abe7
# ╠═90ada351-91f6-4de4-a79b-52927a82fe6c
# ╠═6948180c-0638-4014-a6f4-3ebceb291b54
# ╠═5a565526-013a-4ad5-ac51-bb581f4c67c6
# ╠═d51bad9f-09c7-45f3-b580-f75419058fd7
# ╟─16900190-d8fb-4d6b-967b-3b1af9e89153
# ╠═c18275a4-3ca0-4179-ae88-ed49391be2a8
# ╠═5b449a34-071d-487b-a4e8-836f1c3afbdf
# ╠═94b01121-2ba6-4859-a7e2-df82d424b01f
# ╟─620542bd-441b-46bd-94e7-1ce10496a725
# ╠═55e3366b-e1e8-44b2-ad02-44146ce0d839
# ╠═5d4dc962-dfe9-4833-96b0-21c3ff9c9540
# ╠═67c19e62-f714-4c72-9ad0-8c59f226607d
# ╠═b50e4e4e-1d8c-458e-9791-1981329a2994
# ╠═e1f0257a-10e8-4b6d-ab77-ae99a7da4767
# ╠═7e688179-e2bf-49fc-ada3-2d7ea06fa3c9
# ╠═a260c65a-b482-4be8-80e1-621a333de233
# ╟─c978ee21-2747-47fe-b02e-2e604acc9eed
# ╠═e74dde90-96a3-4562-81e7-cb72a5bd2720
# ╟─cd5eeff7-0829-4a28-837b-d444a65c7238
# ╠═e7506e41-70eb-4a9d-aa68-3250d7d95714
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
