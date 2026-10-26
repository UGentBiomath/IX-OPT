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

# ╔═╡ 92f1185b-ddad-40a8-bd02-cc09ca969519
using PlutoUI; TableOfContents()

# ╔═╡ 85d74ede-e1ac-4ea7-a514-98899159ef45
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq

# ╔═╡ 3c0e40f3-01bc-4503-9577-4f51df700191
using NumericalIntegration, Statistics, Interpolations

# ╔═╡ 40997d48-71fd-4387-b503-7b35aad78e54
using Plots, LaTeXStrings, Plots.PlotMeasures, StatsPlots

# ╔═╡ e755c941-e82c-4b4a-9b07-8b12506de465
using PlutoLinks

# ╔═╡ 7f353f6c-a053-4919-8bdf-5c4cb6ed0fc1
using Random

# ╔═╡ f5023c00-491f-46c7-8a60-39baee9534c5
md"""
## Libraries
"""

# ╔═╡ 1706c1f9-717a-43af-beea-de483934a3a8
md"""
## Modules
"""

# ╔═╡ ff12d5a8-2781-4abf-b10c-62e3863dde7a
IX = @ingredients("IX_mcp_module.jl")

# ╔═╡ be304ff8-b720-4eeb-99f1-5fa6c9a41c2a
Plot = @ingredients("IX_mcp_plotting.jl")

# ╔═╡ 4f44aab8-916b-11f1-1bd9-9bec891531eb
md"""
# IX model – 5 components
"""

# ╔═╡ e4950454-6422-4315-9476-1e98a13f321a
md"""
## Average
"""

# ╔═╡ ca4aa7f8-d9dc-47c0-93d0-36491f633fd3
md"""
### Loading
"""

# ╔═╡ bf420b21-d616-4e9c-8c7d-19243ec09550
@unpack t, x = IX.ivs

# ╔═╡ 983d8d90-a177-4d58-af6e-b1821e7ac46f
@unpack c, N = IX.components(N=5)

# ╔═╡ 714300ab-09c1-4c6d-98da-3bef18957a7d
@unpack C, q = IX.dvs(N=N)

# ╔═╡ 9f9b153f-7cd7-438c-a0a7-d02ab0210aec
model = IX.sys_average(N=5)

# ╔═╡ 2e6ddf43-dbd2-4189-a4a5-ac4931ac9b6c
dx = 1/40

# ╔═╡ 67ddea74-84a1-49d4-ba3c-9c8fa611ceb6
grid_align = edge_align

# ╔═╡ 30d31108-5adc-4433-b4ce-dea1fdae8329
@time prob = IX.MOLdiscretize(model, dx=dx, grid_align=grid_align);

# ╔═╡ 0b5a79d8-0ef2-49b3-82ac-c6d5d9afa701
@time sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1);

# ╔═╡ b9e25013-e140-4542-987b-e6ddb604aed5
update = IX.update

# ╔═╡ 0ede36ed-f4de-4b78-a55e-78fd279688ba
@unpack params, L, D, A, V, Q, ε, u, τ, Dg, Cref, qref, qmax, rₚ, Cₓ₀, Z, M, Dz, Pe, Dp, R, K, C₀, q₀ = IX.pars(N=N)

# ╔═╡ 08d67be8-9471-44dc-8ca7-339a2ea5728d
params[Dg]

# ╔═╡ 03350557-e499-4ff1-8e7a-e93f3f992160
Mi = IX.Mi

# ╔═╡ 8ffc05e1-e92b-44e7-9b70-ba1a8ebcfa49
Zi = IX.Zi

# ╔═╡ 94919636-7f55-4f2e-b9bc-1e19bb6065cc
ps = collect(keys(params))

# ╔═╡ e6e97ce5-7be1-48c8-a6b4-8834be453e87
md"""
### Regeneration
"""

# ╔═╡ c3ad2eda-0770-4bff-94ef-3d67a6dc2716
bcsreg = IX.bcs_average_reg(N=N, Cₓ₀=IX.Cₓ₀reg, C₀=IX.C₀, q₀=IX.q₀load)

# ╔═╡ 357aac49-2645-41d4-b71b-99b26e19f050
regmodel = IX.sys_average_reg(N=N, bcs=bcsreg, domain=(t∈(0, 50.), x∈(0, 1.)));

# ╔═╡ df70d886-1cdc-400b-baf4-5a8e0b782e76
@time regtemp = IX.MOLdiscretize(regmodel, dx=dx; grid_align=grid_align);

# ╔═╡ 7d3f91be-e4c1-4721-9a53-41d8a3fe15d6
md"""
## Full-conservation
"""

# ╔═╡ 9a369bd3-e629-4bdb-9d0c-45fa8cc7722e
full_model = IX.sys_full(N=5)

# ╔═╡ 29337edf-bf3b-4f7e-b58f-e0b604d9b377
@time probF = IX.MOLdiscretize(full_model, dx=dx, grid_align=grid_align);

# ╔═╡ 68370136-628b-4221-8c94-d2ab29b2502c
@time solF = solve(probF, Rodas5P(), saveat=1, progress=true, progress_steps=1);

# ╔═╡ b004c772-2b2a-4845-a6a6-ac78155ed0e6
md"""
### Regeneration
"""

# ╔═╡ a5548d15-f6fe-4e76-8d3e-1a4344d75bd0
regmodelF = IX.sys_full_reg(N=N, bcs=bcsreg, domain=(t∈(0, 50.), x∈(0, 1.)));

# ╔═╡ 27157081-7a3c-4aa4-bbd2-c783283dd11d
@time regtempF = IX.MOLdiscretize(regmodelF, dx=dx; grid_align=grid_align);

# ╔═╡ 8dc820b8-ca06-49f8-bf59-3de899f19c7d
md"""
# Model selection
"""

# ╔═╡ 828a8999-9f52-4e60-b7a4-f9a00d53776f
md"""
Effective-mobility vs. full charge-conservation model
"""

# ╔═╡ ac4fc18d-1a18-44c6-9417-397d9901c7a2
pars = IX.pars(N=N)

# ╔═╡ 7dc086b5-45bb-4bdb-b6a3-78b507419d9d
let
	@unpack Dp = IX.pars(N=N)
	Dpi = prob.ps[Dp]
	
	Dpt = sum([Dpi[i]*sol(:, 1, dv=q(t,x)[i]) for i in 1:N])./sum([sol(:, 1, dv=q(t,x)[j]) for j in 1:N])

	plot(Dpt)
end

# ╔═╡ 3af2d2af-289e-4314-87b5-0e6050add6ff
prob.ps[Dp]

# ╔═╡ 3b6e72c2-1907-4ca5-8b3b-718d6784d7f8
# let
# 	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]
# 	color = vcat(:grey, 1:4)
# 	label = ["H⁺", "Na⁺", "K⁺", "Mg²⁺", "Ca²⁺"]
	
# 	# Average
# 	p1 = plot(xlims=(0, 200), xlabel="Treated volume [BV]", 
# 			  ylims=(-0.5, 100.5), ylabel=L"% exchange sites, $q/q_{max}$", 
# 			  legendposition=:topleft)

# 	for i in 1:N
# 		plot!(BV, 100sol(:, 1, dv=q(t,x)[i]); lw=2.5, c=color[i], lab=label[i])
# 	end

# 	plot!(BV, 100sum([sol(:, 1, dv=q(t,x)[i]) for i in 1:N]); 
# 		  ls=:dashdot, lw=3, c=:black, lab=L"\sum \bar q_i")

# 	p2 = plot(xlims=(0, 200), xlabel="Treated volume [BV]", 
# 			  ylims=(1e-12, 1e-9), ylabel=L"Ion mobility in resin, $D_p$ [m²s⁻¹]",
# 			  legendposition=:inline)

# 	@unpack Dp = IX.pars(N=N)
# 	Dpi = prob.ps[Dp]
	
# 	Dpt = sum([Dpi[i]*sol(:, 1, dv=q(t,x)[i]) for i in 1:N])./sum([sol(:, 1, dv=q(t,x)[j]) for j in 1:N])

# 	plot!(BV, Dpt; lw=2.5, c=:black, lab="", yscale=:log10)
# 	for i in 1:N
# 		hline!([Dpi[i]]; lab=label[i], lw=2, ls=:dot, c=color[i])
# 	end

# 	plot(p1, p2, layout=(2,1))

# 	plot!(dpi=300, size=(680, 500), grid=:off, widen=false)
# end

# ╔═╡ fb566e4a-9b99-475c-b2cf-0f8345e3a319
md"""
# Discretization analysis
"""

# ╔═╡ 2a9a124e-4422-4352-9b8e-b97e2033b4c7
md"""
See notebook...
"""

# ╔═╡ 4ccb0a7c-18cb-43ab-a395-0d0cd047a720
md"""
# Operational scenarios
"""

# ╔═╡ 9ad29232-660e-4b4a-9539-458fcf26afb0


# ╔═╡ 904456f6-c25b-47b7-8a44-c3d11bae309f
md"""
# Operating capacity
"""

# ╔═╡ 72f12e0a-f0fa-4820-a0dc-7e02a3aad68d
let
	plot(dpi=300, size=(680, 400), widen=false, box=:axes, grid=:off, margin=5mm)

	BT = findfirst(x -> x >= 0.2, sol(:, 1, dv=C(t,x)[5]))

	plot!(sol[t], sol(:, 1, dv=C(t,x)[5]); 
		  lw=2.5, lab="", c=:black,
		  xlims=(0, 500), xticks=([sol[t][BT]], []), ylims=(-0.05, 1.2), yticks=[],
		  xlabel="Volume of effluent treated",
		  ylabel="Effluent concentration")
	hline!([sol(550, 1, dv=C(t,x)[5])]; lab="", lw=1, ls=:dot, c=:black)
	hline!([0.2]; lw=1, lab="", ls=:dash, c=:black)

	vline!([sol[t][BT]]; lab="", lw=1, c=:black, ls=:dash)
	scatter!([sol[t][BT]], [sol(sol[t][BT], 1, dv=C(t,x)[5])];
			 lab="", c=:black, marker=:xcross, markerstrokewidth=3, markersize=6)

	annotate!(100, 0.05, ("Leakage concentration", 9, :left))
	annotate!(100, 0.25, ("Breakthrough concentration", 9, :left))
	annotate!(100, 1.08, ("Influent concentration", 9, :left))
	

	savefig("breakthrough.png")
	plot!()
end

# ╔═╡ 8926d4c6-6400-4a2d-bccc-864fd4ffa379
md"""
# Optimization scenarios
"""

# ╔═╡ 20680a91-5093-449a-a15a-2bf2bbd940f1
md"""
## Stopping criteria
"""

# ╔═╡ 4952a6b7-8efd-4e00-ace2-421e84952117
md"""
### 1. Fixed-time operation
"""

# ╔═╡ b3a8e288-b396-4846-842c-e5348f06125d
col = vcat(:grey, 1:N-1)

# ╔═╡ 07fc17b6-94fc-4b95-b8dc-8c69eb9554b9
md"""
## Fixed loading-regeneration cycles
"""

# ╔═╡ f8084726-27c7-45bd-b376-ce0f61879130
md"""
### Breakthrough, resin capacity
"""

# ╔═╡ d4c76776-194c-439b-873f-1091f1bdceba
md"""
Breakthrough for Na⁺ ions:
$(@bindname NaBT Select([0.05, 0.1, 0.15, 0.20], default=0.1)) = $(100*NaBT)% inlet Na concentration
"""

# ╔═╡ 5ff8098e-42a0-4f26-8b90-336b7ffc1110
BTidxreg = findfirst(x -> x >= NaBT, solF(:, 1, dv=C(t,x)[2]))

# ╔═╡ c6fe9b00-6bcc-4be0-a34c-b883da867028
regprobF = remake(regtempF; u0=IX.update_u0(probF, solF, solF[t][BTidxreg]; reverse=true));

# ╔═╡ af271f05-a813-4d81-a83c-55ec1902a542
@time solregF = solve(regprobF, Rodas5P(), saveat=0.01);

# ╔═╡ 7844f7e5-c5c2-4150-b688-890aec27f59f
md"""
Default HCl regeneration volume:
$(@bindname Vreg Slider(0:0.1:3, default=1, show_value=true)) = $Vreg L
"""

# ╔═╡ 46bfb780-21f4-4522-bd84-6aef8808b25e
BTidx = findfirst(x -> x >= NaBT, sol(:, 1, dv=C(t,x)[2]))

# ╔═╡ 49bce273-ec27-432b-8c30-2e8ae4f84467
sol[t][BTidx]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]

# ╔═╡ 902410c1-edb7-4e23-a460-c969169cced5
sol[t][BTidx]*prob.ps[τ]*prob.ps[Q]*1e3

# ╔═╡ f6d8b3d3-4c8a-4a66-8a28-40dff51207e4
regprob = remake(regtemp; u0=IX.update_u0(prob, sol, sol[t][BTidx]; reverse=true));

# ╔═╡ 3c1b1d7a-af93-429b-8077-9cb91531e2c7
@time solreg = solve(regprob, Rodas5P(), saveat=0.01);

# ╔═╡ e8c5f1cc-46c7-46a0-9aed-fee9ba6d71c4
let
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]
	color = vcat(:grey, 1:4)
	label = ["H⁺", "Na⁺", "K⁺", "Mg²⁺", "Ca²⁺"]

	# Breakthrough
	
	# Average
	p1 = plot(xlims=(0, 200), xlabel="Treated effluent volume [BV]", 
			  ylims=(-0.01, 4.7), ylabel=L"Normalized ion concentration, $C/C_0$",
			  title="(a) Breakthrough curves")

	for i in 1:N
		plot!(BV, sol(:, 1, dv=C(t,x)[i]); lw=3, c=color[i], lab=label[i])
	end

	# Full-conservation
	for i in 1:N
		plot!(BV, solF(:, 1, dv=C(t,x)[i]); lw=2, ls=:dot, c=color[i], lab="")
	end

	# Elution curves

	BVreg = solreg[t]*regprob.ps[τ]*regprob.ps[Q]/regprob.ps[V]

	p2 = plot(xlims=(0, 10), ylims=(-0.01, 330), 
			  xlabel="Spent regenerant volume [BV]", 
			  title="(b) Elution curves", legendposition=:bottomright)

	# Average
	for i in 1:N
		plot!(BVreg, solreg(:, 1, dv=C(t,x)[i]); lw=3, c=color[i], lab=label[i])
	end

	# Full-conservation
	for i in 1:N
		plot!(BVreg, solregF(:, 1, dv=C(t,x)[i]); lw=2, ls=:dot, c=color[i], lab="")
	end
	
	plot(p1, p2)
	
	plot!(dpi=300, size=(680, 400), grid=:off, widen=false)
	savefig("model_selection.png")
	plot!()
end

# ╔═╡ c37ed861-46fa-4c87-b8ee-e1c57523c824
RCidx = findfirst(x -> x >= Vreg, solreg[t]*regprob.ps[τ]*regprob.ps[Q]*1e3)

# ╔═╡ 3f98f4aa-2dab-4b8b-9993-970f1ac0e814
md"""
Treated volume: $(@bindname BVstopfixed Slider(20:1:80, default=25, show_value=true)) $BVstopfixed L
"""

# ╔═╡ 191e5a51-1261-4b0b-8bf0-3d9be58dae1b
md"""
Spent regenerant volume: $(@bindname BVregstopfixed Slider(0:0.1:3, default=1, show_value=true)) $BVregstopfixed L
"""

# ╔═╡ 692bf75f-8dc6-4949-8928-788ec55eb8d6
let
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]
	BVreg = solreg[t]*regprob.ps[τ]*regprob.ps[Q]/regprob.ps[V]
	BT = Plot.breakthrough(prob, sol; comp=2)
	cap = 100*Plot.capacity(prob, sol; comp=1)
	c = :black
	ls = :dot
	# Breakthrough
	p1 = Plot.breakthrough_curves(prob, sol; 
								  show_pH=true,
								  comp=2:N, lw=3, time="BV",
								  breakthrough_line=NaBT, 
								  legendposition=:topleft,
								  widen=false,
								  title="(a) Loading")
	p2 = Plot.distribution(prob, sol; time="BV", Vmax=200, xticks=0:50:200)

	vline!(p1.subplots[1], [BV[BTidx]]; lab="", c=c, lw=2, ls=ls)
	vline!(p2.subplots[1], [BV[BTidx]]; lab="", c=c, lw=2, ls=ls)
	vline!(p1.subplots[2], [BV[BTidx]]; lab="", c=c, lw=2, ls=ls)

	# vline!(p1.subplots[1], [BVstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	# vline!(p1.subplots[2], [BVstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	# vline!(p2.subplots[1], [BVstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	
	hline!(p1.subplots[1], [Plot.pH(prob, sol)[BTidx]]; lab="", c=c, lw=2, ls=ls)
	hline!(p1.subplots[2], [BT[BTidx]]; lab="", c=c, lw=2, ls=ls)
	hline!(p2.subplots[1], [100Plot.capacity(prob, sol; comp=1)[BTidx]]; lab="", c=c, lw=2, ls=ls)
	
	scatter!(p1.subplots[2], [BV[BTidx]], [BT[BTidx]]; 
			 lab="", c=:black, m=:x, markersize=6, markerstrokewidth=2.5)
	scatter!(p1.subplots[1], [BV[BTidx]], [Plot.pH(prob, sol)[BTidx]]; 
			 lab="", c=:black, m=:x, markersize=6, markerstrokewidth=2.5)
	scatter!(p2.subplots[1], [BV[BTidx]], [cap[BTidx]]; 
			 lab="", c=:black, m=:x, markersize=6, markerstrokewidth=2.5)

	l = @layout [ a{0.6h}; b]
	L = plot(p1, p2, layout=l, xlims=(0, 200))
	plot!(size=(680, 700))

	# Regeneration
	BT = Plot.breakthrough(regprob, solreg; comp=2)
	cap = 100*Plot.capacity(regprob, solreg; comp=1)
	
	p3 = Plot.elution_curves(regprob, solreg; 
								  show_pH=true,
								  comp=2:N, lw=3, time="BV",
								  breakthrough_line=0.1, 
							 		xlims=(0, 10), xticks=0:2:10,
							 		# xlabel="Treated volume [BV]",
								  widen=false, title="(b) Regeneration")
	p4 = Plot.distribution(regprob, solreg; 
						   time="BV",
						   elution=true, xlims=(0, 10), xticks=0:2:10, 
						   legendposition=:bottomright)
	
	vline!(p3.subplots[1], [BVreg[RCidx]]; lab="", c=c, lw=2, ls=ls)
	vline!(p3.subplots[2], [BVreg[RCidx]]; lab="", c=c, lw=2, ls=ls)
	vline!(p4.subplots[1], [BVreg[RCidx]]; lab="", c=c, lw=2, ls=ls)

	# vline!(p3.subplots[1], [BVregstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	# vline!(p3.subplots[2], [BVregstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	# vline!(p4.subplots[1], [BVregstopfixed]; lab="", c=c, lw=3, ls=:dashdot)
	
	hline!(p3.subplots[1], [Plot.pH(regprob, solreg)[RCidx]]; lab="", c=c, lw=2, ls=ls)
	hline!(p3.subplots[2], [BT[RCidx]]; lab="", c=c, lw=2, ls=ls)
	hline!(p4.subplots[1], [cap[RCidx]]; lab="", c=c, lw=2, ls=ls)
	
	scatter!(p3.subplots[1], [BVreg[RCidx]], [Plot.pH(regprob, solreg)[RCidx]]; 
			 lab="", c=c, m=:x, markersize=6, markerstrokewidth=2.5)
	scatter!(p3.subplots[2], [BVreg[RCidx]], [BT[RCidx]]; 
			 lab="", c=c, m=:x, markersize=6, markerstrokewidth=2.5)
	scatter!(p4.subplots[1], [BVreg[RCidx]], [cap[RCidx]]; 
			 lab="", c=c, m=:x, markersize=6, markerstrokewidth=2.5)
	
	plot(p3, xlims=(0, 9))
	R = plot(p3, p4, layout=l)

	plot(L, R, layout=(1, 2), size=(680, 550))
	savefig("fixed_curves.png")
	plot!()
end

# ╔═╡ ad5c8b41-82c8-4034-9e24-726597aca74f
md"""
### Consecutive running cycles
"""

# ╔═╡ 1687efa2-b296-47cc-8b7d-4036ec34dd04
md"""
In this criterion, a fixed duration of both the loading and regeneration cycles is used.
"""

# ╔═╡ de199fc0-62d5-4bb5-9272-41e024197544
md"""
$(@bindname Nc NumberField(1:20, default=10)) $Nc cycles of loading-regeneration
"""

# ╔═╡ 41f8477b-77a3-49d6-bb5f-6b474959c0cf
md"""
$(@bindname BVstop confirm(Slider(0:1:60, default=25, show_value=true))) $BVstop L of treated volume
"""

# ╔═╡ be1cb9c5-5c01-43e2-bb93-e67b9f0937d4
md"""
$(@bindname BVregstop confirm(Slider(0:0.1:3, default=0.8, show_value=true))) $BVregstop L of spent regenerant
"""

# ╔═╡ 60e89f62-b975-4f0b-bed4-2d98e7aa0c08
p0 = vcat(prob.ps[qmax], [prob.ps[K][i] for i in 2:N])

# ╔═╡ 83c7674e-e663-4fad-a074-947449fb0272
function fixed_load(; prob=prob)

	update(params, 
		   Q => 120/60*1e-6, 
		   Cₓ₀ => vcat(1e-6, [params[Cₓ₀][i] for i in c[2:3]]), 
		   qmax => p0[1], K => vcat(1.0, p0[2:N]), 
	)
	@time sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

	return sol
end

# ╔═╡ 497590c0-93d4-4921-8603-0fba2d40aa1f
function fixed_cycles(; Nc=1, BVstop=20, BVregstop=1)

	loadprob = prob

	if prob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x]
	elseif prob.problem_type.disc.grid_align == MethodOfLines.EdgeAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x][2:end]
	end

	load = []; regn = []; temp = []; tempreg = []; loadtemp = []
	loads = []; regns = []; sols = []; regs = []; ends1 = []; ends2 = []

	for n in 1:Nc
		
		# Loading
		loadparams = IX.update(params, Q => 120/60*1e-6, 
								  Cₓ₀ => [IX.Cₓ₀load[i] for i in c],
								  qmax => p0[1], K => vcat(1.0, p0[2:N]))

		# Initial cycle
		if n == 1
			prob = remake(loadprob, p=loadparams)
		elseif n > 1
			solreg = temp; reg = tempreg
			end2 = findfirst(x -> x >= BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)
		
			prob = remake(loadprob, p=loadparams, 
						  u0=IX.update_u0(reg, solreg, solreg[t][end2];reverse=true))
		end

		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		end1 = findfirst(x -> x >= BVstop, sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)

		use = 100*integrate(grid, sol(sol[t][end1], grid, dv=q(t,x)[1])) 
		# use = 100 * Plot.capacity(prob, sol; comp=1)[end1]
		push!(load, use)
		push!(loads, prob)
		push!(sols, sol)
		push!(ends1, end1)
	
		
		# Regeneration
		if regprob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x]
		elseif regprob.problem_type.disc.grid_align ==MethodOfLines.EdgeAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x][2:end]
		end
		regparams = IX.update(params, Q => 25/60*1e-6, 
							  Cₓ₀ => [IX.Cₓ₀reg[i] for i in c], 
							  Cref => prob.ps[Cref], 
							  qmax => p0[1], K => vcat(1., p0[2:N]))

		reg = remake(regprob, p=regparams, 
					 u0=IX.update_u0(prob, sol, sol[t][end1]; reverse=true))
		
		solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
		temp = solreg; tempreg = reg
		grid = solreg[x]

		end2 = findfirst(x -> x >= BVregstop, solreg[t]*reg.ps[τ]*reg.ps[Q]*1e3)

		# cap = 100 * Plot.capacity(reg, sol; comp=1, time=solreg[t][1:end2])
		cap = 100*integrate(grid, solreg(solreg[t][end2], grid, dv=q(t,x)[1]))
		push!(regn, cap)
		push!(regns, reg)
		push!(regs, solreg)
		push!(ends2, end2)
	end
		
	return (; load, regn, Nc, BVstop, BVregstop, loads, regns, sols, regs, ends1, 
			ends2)
end

# ╔═╡ 0566f850-5083-4e40-9ef4-263f6a103c64
@time fixed_time_cycles = fixed_cycles(Nc=Nc, BVstop=BVstop, BVregstop=BVregstop)

# ╔═╡ fcd9a255-e395-4bb3-bf27-c2f335566d26
let
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = fixed_time_cycles

	plot()
	plot!(sols[1][t][1:ends1[1]], 
		  sols[1](sols[1][t][1:ends1[1]], 1, dv=C(t,x)[2])*loads[1].ps[Cref]*Mi[:Na]*1e3;
		  lab="Cycle 1", lw=3, widen=false)
	tacc = 0.
	for i in 2:Nc
		tacc += sols[i][t][ends1[i-1]]
		plot!(sols[i][t][1:ends1[i]] .+ tacc, 
			  sols[i](sols[i][t][1:ends1[i]], 1, dv=C(t,x)[2]) * 	
				loads[i].ps[Cref]*Mi[:Na]*1e3;
			  lab="Cycle $i", legendposition=:inline, lw=3, c=1)
	end
	hline!([sols[end](sols[end][t][ends1[end]], 1, dv=C(t,x)[2]) * 1.025 *
				loads[end].ps[Cref]*Mi[:Na]*1e3]; 
		   lab="", ls=:dash, c=:grey, lw=2, alpha=0.5)
	plot!(dpi=300, size=(680, 400), xlabel="Treated volume [L]", grid=:off,
		  ylabel="Breakthrough concentration [mg/L]", ylims=(0, 40))
	savefig("breakthrough_cycles.png")
	plot!()
end

# ╔═╡ e4d6e820-5bc0-4212-b6ad-22f1c74b4ab7
fixed_time_cycles.load[end] - 100

# ╔═╡ fa8382c1-4f18-45dd-b432-4066c2d9ea23
fixed_time_cycles.regn[end]

# ╔═╡ 63256b4d-6487-4d23-895c-8c7143e5bf2b
let
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = fixed_time_cycles

	# effluent quality after breakthrough
	Ymax = 20
	p1 = plot(xticks=1:Nc, ylims=(0, Ymax), topmargin=3mm, grid=:off,
			  ylabel="Concentration, %",
			  title="(a) Effluent quality after breakthrough, % inlet Na⁺ concentration")
	
	Ci = 100 * [sols[i][C(t,x)[2]][ends1[i], end] for i in 1:Nc]
	
	bar!(1:Nc, Ci; lab="", c=:lightblue, fillalpha=1)
	for n in 1:Nc
		annotate!(n, 0.08Ymax + Ci[n], ("$(round(Ci[n];digits=1))", 8))
	end	
	
	# resin capacity remaining after loading
	Ymax = 30
	p2 = plot(xticks=1:Nc, ylims=(0, Ymax), topmargin=3mm, grid=:off,
			  ylabel="Remaining capacity, %", framestyle=:box, 
			  title="(a) Remaining resin capacity before regeneration, $(round(Int,BVstop/prob.ps[V]*1e-3)) BV treated")

	bar!(1:Nc, load; lab="Resin capacity, H⁺ form", c=:grey, fillalpha=0.7, 
		 legendposition=:topleft, xforegroundcoloraxis=:white)
	for n in 1:Nc
		# annotate!(n, 0.05Ymax + load[n], ("$(round(load[n];digits=1))", 8))
	end

	Ymax = 50
	plot!(twinx(), 1:Nc, Ci; lab="Na⁺ breakthrough concentration", c=1, lw=2, 
		  m=:circle, framestyle=:axes, legendposition=:topright,
		  ylims=(0, Ymax), yticks=0:10:Ymax, ylabel="% Na⁺ inlet concentration")
	

	# resin capacity remaining after regeneration
	Ymin = 0
	Ymax = 100
	p3 = plot(xticks=1:Nc, ylims=(Ymin, Ymax), yticks=Ymin:20:Ymax, 
			  topmargin=3mm, grid=:off,
			  xlabel="No. of cycles loading-regeneration",
			  ylabel="Occupied capacity, %",
			  title="(b) Remaining resin capacity after regeneration, $(round(Int,BVregstop/prob.ps[V]*1e-3)) BV HCl 5% spent")
	
	# bar!(1:Nc, regn; lab="", c=:darkseagreen, fillalpha=1)
	for n in 1:Nc
		# annotate!(n, 0.08(Ymax-Ymin) + regn[n], ("$(round(regn[n];digits=1))", 8))
	end

	Ymax = 12
	Ci = 100 * [Plot.capacity(regns[i], regs[i]; comp=j)[ends2[i]] for i in 1:Nc, j in N:-1:1] |> stack
	color = permutedims(reverse(Plot.color[1:N]))
	label = permutedims(reverse(Plot.label[1:N]))
	# plot!(twinx(), 1:Nc, Ci; lw=2, c=:black, m=:square,
	# 	  ylims=(0, Ymax), yticks=0:4:Ymax, ylabel="% resin capacity",
	# 	  lab="Accumulated Mg⁺ and Ca⁺ ions", legendposition=:topright)
	groupedbar!(1:Nc, Ci; 
				bar_position=:stack, c=color, lab="", fillalpha=0.7,
			    xforegroundcoloraxis=:white, legendcolumns=N)
	[plot!([], []; seriestype=:bar, c=Plot.color[j], lab=Plot.label[j],
		    fillalpha=0.7, legendposition=:outertop) for j in 1:N]

	figure = plot(p2, p3, layout=(2,1), dpi=300, size=(780, 650), 
				  plot_titlevspan=0.1, legendfontsize=10,
				  # xlabel="No. of cycles loading-regeneration", 
				  # suptitle="Column performance with fixed-time stop criteria control"
				 )
	# figure |> aside
	savefig("fixed_cycles.png")
	plot!()
end

# ╔═╡ 6eb6d99b-65e2-47aa-9095-480452bd2954
let 
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = fixed_time_cycles

	@show 100-load[end]
	@show regn[end]
	@show BVstop / loads[end].ps[V]*1e-3
	@show BVregstop / regns[end].ps[V]*1e-3
end

# ╔═╡ a7143b46-8b83-4f16-9b09-9aa75b87f0b7
md"""
## Adaptive control of cycle length
"""

# ╔═╡ cccc18c5-3ad0-4506-bac1-b52fc02aa332
md"""
Based on measurable quantities, such as:
- pH
- Na breakthrough
"""

# ╔═╡ f2aebf87-7197-4c39-bc68-90e6b0cb64ac
md"""
### pH measurement
"""

# ╔═╡ 05ee4bd3-1f29-410c-bfda-0e768f8a47fe
md"""
Assuming a resolution of 0.01 pH units and a sampling frequency of 1 Hz (1 sample per second), electrochemical probes used in industry can detect pH changes at a rate of 0.01 pH units s$^{-1}$.

For a measurement accuracy of ±0.05 pH units, and for an average response time of 3 s (~0.3 Hz), a more realistic lower bound for pH change could be ~0.015 pH units s$^{-1}$.

Based on these assumptions and a slope of 0.01, the breakthrough detection reaches 23 mg/L Na, or equivalently 23 ppm or 0.001 M. This is much higher than the best sodium analyzers.

Another important issue is the daily calibration of the pH probes to avoid sensor drift.
"""

# ╔═╡ 9387f387-3901-4de3-ba62-85d29dbf94fb
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
	# figure |> aside
end

# ╔═╡ 41409955-ed92-4bb3-ae31-875d4699a54f
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

# ╔═╡ 787f5fe0-fc33-4c64-b722-6cc20292968b
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

# ╔═╡ a4025273-fbe3-4f72-9f1e-dbd658483a5b
md"""
### pH control
"""

# ╔═╡ 3291abb6-ca4d-47b7-9b84-b2c9a422a4aa
md"""
pH change: $(@bind pHchange Select([0.01, 0.05, 0.1], default=0.1)) units
"""

# ╔═╡ 56bcb6a3-a53f-427c-a365-f4f528c0f754
md"""
Filter window: $(@bindname w1 Slider(1:15, default=3, show_value=true))
"""

# ╔═╡ b5b4e4a1-43f1-41e5-b442-b57a61ec763b
let
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V] # Bed volumes [L] 
	Vmax = 120
	Nt = length(BV)
	tmax = 500
	ΔpH = pHchange
	σ = pHchange/2
	# w1 = 3
	Random.seed!(1234);
	
	pHmax = 8; 	dpHmin = -0.0008; 	dpHmax = 0.2001
	idxmin = 3
	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	dpH = diff(pH)./diff(BV)
		
	Na = sol(:, 1, dv=C(t,x)[2])#*prob.ps[Cref]#*Mi[:Na]*1e3 # [mg/L]
	K = sol(:, 1, dv=C(t,x)[3])#*prob.ps[Cref]#*Mi[:K]*1e3 # [mg/L]
	Mg = sol(:, 1, dv=C(t,x)[4])#*prob.ps[Cref]#*Mi[:Mg]/Zi[:Mg]*1e3 # [mg/L]
	Naint = [integrate(BV[1:j], Na[1:j]) for j in 2:Nt] # [mg]

	# (a) pH
	p1 = plot(xlims=(0, Vmax), xlabel="Treated volume [BV]", ylims=(0, pHmax), 
			  ylabel="pH", grid=:off)
	plot!(BV[idxmin:end], pH[idxmin:end]; 
		  lab="", c=:grey, lw=3, legendposition=:topright)	
	# plot!([], [], lab="ΔpH/Δt", c=col[1], ls=:dashdot, lw=3, grid=:off)

	# initial pH
	pH0 = -log10(sum([IX.:Cₓ₀load[c[i]] for i in 2:N])) # pH = 1.7
	# hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
	annotate!(1, 0.5*pH0, ("initial pH: $(round(pH0;sigdigits=3))", 9, :left))

	# noisy measurement
	noisy_pH = pH[idxmin:end] .+ σ.*randn(length(pH[idxmin:end]))
	
	# first detection of pH change
	idxpH = findfirst(x -> x >= pH0+ΔpH, filter(noisy_pH, w1)) + idxmin - 1
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2.5)
	hline!([pH0+ΔpH]; lab="", c=N+1, ls=:dot, lw=2.5)
	scatter!([BV[idxpH]], [pH0+ΔpH]; 
			 lab="", m=:xcross, markersize=5, c=:black, markerstrokewidth=2)
	annotate!(0.955*BV[idxpH], 0.5*pH0, ("pH: $(round(pH0+ΔpH;sigdigits=3)) at $(round(Int,BV[idxpH])) BV ($pHchange change)", 9, :right))	

	# Resin capacity remaining and restored
	@show Plot.capacity(prob, sol; comp=1)[idxpH] - 1
	regprob = remake(regtemp; u0=IX.update_u0(prob, sol, sol[t][idxpH]; reverse=true));
	solreg = solve(regprob, Rodas5P(), saveat=0.01);
	RCidx = findfirst(x -> x >= Vreg, solreg[t]*regprob.ps[τ]*regprob.ps[Q]*1e3)
	@show Plot.capacity(regprob, solreg; comp=1)[RCidx]


	# inset
	plot!(BV, pH; inset = (1, bbox(0.075, 0., 0.55+ΔpH, 0.6)), subplot=2,
			 xlims=(70, 95), xticks=60:5:100, ylims=(1.6, 1.9), yticks=1:0.1:2, 
			 lab="model", c=:grey, lw=3, box=:box, grid=:off, legendfonthalign=:center,
		 	 legendposition=:outertop, legendcolumns=-1, legendfontsize=8)
	# noisy measurement
	plot!(BV[idxmin:end], noisy_pH; 
		  subplot=2, c=:grey, ls=:dot, lw=2, lab="pH probe")
	# moving average filter
	plot!(BV[idxmin:end], filter(noisy_pH, w1); 
		  subplot=2, lw=2, c=:darkred, lab="moving average")
	hline!([pH0+ΔpH]; lab="", c=N+1, ls=:dot, lw=2.5, subplot=2)
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2.5, subplot=2)
	# scatter!([BV[idxpH]], [pH0+ΔpH]; subplot=2,
	# 		 lab="", m=:xcross, markersize=9, c=:black, markerstrokewidth=2)

	# pH slope
	# plot!(twinx(), BV[idxmin:end-1], dpH[idxmin:end];
	# 	  lab="", c=col[1], ls=:dashdot, lw=3, grid=:off,
	# 	  xlims=(0, tmax), ylims=(dpHmin, dpHmax), ylabel="ΔpH/Δt")

	# (b) Breakthrough
	Cmax = 4
	p2 = plot(xlims=(0, Vmax), xlabel="Treated volume [BV]", ylims=(-0.02, Cmax),
			  ylabel=L"Na⁺ concentration, $C/C_0$", grid=:off, legendposition=:topleft)
	plot!(BV[1:idxpH], Na[1:idxpH]; lab="Na⁺", c=1, lw=3)
	plot!(BV[idxpH:end], Na[idxpH:end]; lab="", c=1, lw=3, ls=:dot, alpha=0.9)
	plot!(BV[1:idxpH], Na[1:idxpH]; fillrange=0, c=1, alpha=0.3, lab="")
	plot!(BV[1:end], K[1:end]; lab="K⁺", c=2, lw=3, ls=:dot, alpha=0.9)

	# breakthrough detection
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2.5)
	hline!([Na[idxpH]]; lab="", c=N+1, lw=2.5, ls=:dot)
	scatter!([BV[idxpH]], [Na[idxpH]]; 
			 lab="", m=:xcross, markersize=5, c=:black, markerstrokewidth=2)
	annotate!(0.98*BV[idxpH], 0.5 + Na[idxpH],
			  ("Breakthrough of Na⁺ ions after $(round(Int,BV[idxpH])) BV: $(round(Na[idxpH]*prob.ps[Cref]*Mi[:Na]*1e3;digits=1)) mg/L ", 9, :right))
	# (total: $(round(Int, Naint[idxpH-1])*prob.ps[Cref]*Mi[:Na]*1e3) mg)
	
	figure = plot(p1, p2, layout=grid(2, 1, heights=[0.5, 0.5]), 
				  dpi = 300, size=(680, 600), 
				  suptitle="Cycle stopping with pH change threshold detection")
	# figure |> x-> aside(x; v_offset=0)
	savefig("pH_control.png")
	plot!()
end

# ╔═╡ 731eae26-849a-4119-ab8a-5b4e8131b5c0
function pH_stop_fixed_reg_cycles(; Nc=1, BVstop=20, BVregstop=3)

	loadprob = prob

	if prob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x]
	elseif prob.problem_type.disc.grid_align == MethodOfLines.EdgeAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x][2:end]
	end

	load = []; regn = []; temp = []; tempreg = []; loadtemp = []
	loads = []; regns = []; sols = []; regs = []; ends1 = []; ends2 = []

	for n in 1:Nc
		
		# Loading
		loadparams = IX.update(params, Q => 120/60*1e-6, 
								  Cₓ₀ => [IX.Cₓ₀load[i] for i in c],
								  qmax => p0[1], K => vcat(1.0, p0[2:N]))

		# Initial cycle
		if n == 1
			prob = remake(loadprob, p=loadparams)
		elseif n > 1
			solreg = temp; reg = tempreg
			end2 = findfirst(x->x>=BVregstop,solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])
		
			prob = remake(loadprob, p=loadparams, 
						  u0=IX.update_u0(reg, solreg, solreg[t][end2];reverse=true))
		end

		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		# pH stop
		ΔpH = pHchange
		σ = pHchange/2
		idxmin = 3
		pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
		pH0 = -log10(sum([IX.:Cₓ₀load[c[i]] for i in 2:N])) # pH = 1.7
		noisy_pH = pH[idxmin:end] .+ σ.*randn(length(pH[idxmin:end]))
		end1 = findfirst(x -> x >= pH0+ΔpH, filter(noisy_pH, w1)) + idxmin - 1
		# end1 = findfirst(x -> x >= BVstop, sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)

		use = 100*integrate(grid, sol(sol[t][end1], grid, dv=q(t,x)[1])) 
		# use = 100 * Plot.capacity(prob, sol; comp=1)[end1]
		push!(load, use)
		push!(loads, prob)
		push!(sols, sol)
		push!(ends1, end1)
	
		
		# Regeneration
		if regprob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x]
		elseif regprob.problem_type.disc.grid_align ==MethodOfLines.EdgeAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x][2:end]
		end
		regparams = IX.update(params, Q => 25/60*1e-6, 
							  Cₓ₀ => [IX.Cₓ₀reg[i] for i in c], 
							  Cref => prob.ps[Cref], 
							  qmax => p0[1], K => vcat(1., p0[2:N]))

		reg = remake(regprob, p=regparams, 
					 u0=IX.update_u0(prob, sol, sol[t][end1]; reverse=true))
		
		solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
		temp = solreg; tempreg = reg
		grid = solreg[x]

		end2 = findfirst(x -> x >= BVregstop,solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])

		# cap = 100 * Plot.capacity(reg, sol; comp=1, time=solreg[t][1:end2])
		cap = 100*integrate(grid, solreg(solreg[t][end2], grid, dv=q(t,x)[1]))
		push!(regn, cap)
		push!(regns, reg)
		push!(regs, solreg)
		push!(ends2, end2)
	end
		
	return (; load, regn, Nc, BVstop, BVregstop, loads, regns, sols, regs, ends1, 
			ends2)
end

# ╔═╡ 850c2a32-d4c5-4dc3-88e5-4fcb31cedb18
@time pH_cycles = pH_stop_fixed_reg_cycles(Nc=5, BVstop=BVstop, BVregstop=3.)

# ╔═╡ 2743c754-8fd9-4d35-b82e-f577399dca16
let 
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = pH_cycles

	@show 100-load[end]
	@show regn[end]
	@show sols[end][t][ends1[end]]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]
	@show BVregstop #/ regns[end].ps[V]*1e-3
	@show sols[end](:,1, dv=C(t,x)[2])[ends1[end]]*100
end

# ╔═╡ f877762e-7ab5-4b3a-9065-1dc728d54b4d
md"""
### pH slope control
"""

# ╔═╡ bb8013c2-8299-4a6b-acc8-0b8d9a5d98b4
md"""
pH slope resolution (units s⁻¹): 
$(@bind dpHslope confirm(Slider(0.001:0.0005:0.02, default=0.01, show_value=true)))
"""

# ╔═╡ d16741e3-5138-4c71-8121-b786556a89e2
md"""
Filter window (min): $(@bind w1slope Slider(1:15, default=3, show_value=true))
"""

# ╔═╡ 4001ddf9-2cd7-4704-bfe3-ab24b6c12d62
md"""
Slope filter window (min): $(@bind w2slope Slider(1:30, default=15, show_value=true))
"""

# ╔═╡ 9c7a047a-e940-4599-9694-04d133ed2749
let
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V] # Bed volumes [L]
	Nt = length(BV)
	Vmax = 120
	tmax = 500
	ΔpH = pHchange
	σ = pHchange/2
	w1 = w1slope
	w2 = w2slope
	Random.seed!(1234);
	SF = 1.1
	
	pHmax = 8
	imin = 3
	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	dpH = diff(pH)./diff(BV/60/(prob.ps[Q]*1e3)) # [min⁻¹]
		
	Na = sol(:, 1, dv=C(t,x)[2])#*prob.ps[Cref]*Mi[:Na]*1e3 # [mg/L]
	K = sol(:, 1, dv=C(t,x)[3])#*prob.ps[Cref]*Mi[:K]*1e3 # [mg/L]
	Mg = sol(:, 1, dv=C(t,x)[4])#*prob.ps[Cref]*Mi[:Mg]/Zi[:Mg]*1e3 # [mg/L]
	Naint = [integrate(BV[1:j], Na[1:j]) for j in 2:Nt] # [mg]

	# (a) pH
	p1 = plot(xlims=(0, Vmax), xlabel="Treated volume [BV]", legendposition=:topright,
			  ylims=(0, pHmax), ylabel="pH", grid=:off, title="(a) pH change threshold: ΔpH = $ΔpH, $w1-min moving average")
	plot!(BV[imin:end], pH[imin:end]; lab="", c=:grey, lw=3)	

	# initial pH
	pH0 = -log10(sum([IX.:Cₓ₀load[c[i]] for i in 2:N])) # pH = 1.7
	annotate!(1, 0.5*pH0, ("initial pH: $(round(pH0;digits=2))", 9, :left))

	# noisy pH measurement (upsampled)
	# pHint = LinearInterpolation(BV, pH)
	# BVup = 0:(3/60):BV[end]
	# pHup = pHint(BVup)
	# noisy_pH_up = pHup .+ σ.*randn(length(pHup))

	# noisy pH measurement
	noisy_pH = pH[imin:end] .+ σ.*randn(length(pH[imin:end]))
	noisy_dpH = filter(diff(noisy_pH)./diff(BV[imin:end]/60/(prob.ps[Q]*1e3)), w2)
	plot!(BV[imin:end], noisy_pH; c=:darkblue, ls=:dot, lw=1, lab="")
	
	# first detection of pH change
	idxpH = findfirst(x -> x >= pH0+ΔpH, filter(noisy_pH, w1)) + imin - 1
	@show BV[idxpH]/prob.ps[V]*1e-3
	@show Na[idxpH]*prob.ps[Cref]*Mi[:Na]*1e3
	
	hline!([filter(noisy_pH, w1)[idxpH-imin+1]]; 
		   lab="", c=N+1, ls=:dot, lw=3)
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=3)
	scatter!([BV[idxpH]], [pH0+ΔpH]; 
			 lab="", m=:xcross, markersize=6, c=:black, markerstrokewidth=2)
	annotate!(0.955*BV[idxpH], 0.5*pH0, ("pH: $(round(filter(noisy_pH, w1)[idxpH-imin+1];digits=2)) after $(round(Int,BV[idxpH])) BV", 9, :right))	

	# inset
	plot!(BV[imin:end], pH[imin:end]; 
		  inset = (1, bbox(0.075, 0.0, 0.55+ΔpH, 0.6)), subplot=2, 
		  xlims=(60, 95), xticks=60:5:100, ylims=(1.6, 1.9), yticks=1.6:0.1:1.9, 
		  lab="model", c=:grey, lw=3, box=:box, grid=:off, ylabel="", 
		  legendposition=:outertop, legendcolumns=3, legendfontsize=8)
	
	# noisy measurement
	plot!(BV[imin:end], noisy_pH; 
		  subplot=2, c=:blue, ls=:dot, lw=1.5, lab="pH probe")
	# moving average filter
	plot!(BV[imin:end], filter(noisy_pH, w1); 
		  subplot=2, lw=2.5, c=:darkred, lab="moving average")

	hline!([pH0+ΔpH]; subplot=2, lab="", ls=:dot, c=:black, lw=1)
	hline!([filter(noisy_pH, w1)[idxpH-imin+1]]; 
		   lab="", c=N+1, ls=:dot, lw=2.5, subplot=2)
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2.5, subplot=2)

	# up-sample
	# plot!(BVup, noisy_pH_up; subplot=2, lab="", c=:orange, lw=2)
	# plot!(BVup, filter(noisy_pH_up, 3w); 
	# 	  subplot=2, lw=2, c=:darkred, lab="0.3 Hz")

	
	# (b) pH slope
	dpHmin = -0.05; 	dpHmax = 0.15
	p2 = plot(yticks=dpHmin:0.05:dpHmax, xlabel="Treated volume [BV]",
			  xlims=(0, Vmax), ylims=(dpHmin, dpHmax), ylabel="ΔpH/Δt [min⁻¹]", 
			  grid=:off, minorgrid=:off, minorticks=1, legend=:bottomright,
			  title="(b) pH slope threshold: ΔpH/Δt = $dpHslope min⁻¹, $w2 min,  $(round(Int, 100*(SF-1)))% SF", leftmargin=3mm)
	
	plot!(BV[imin+1:end], dpH[imin:end]; lab="", c=:green, lw=3)
	# hline!([0.01]; lab="", ls=:dot, lw=0.75, c=:black)
	plot!(BV[imin+1:end], noisy_dpH; c=:darkblue, ls=:dot, lw=1, lab="")

	# pH slope detection
	idxpH2 = findfirst(x -> x >= SF*dpHslope, noisy_dpH[3imin:end]) + 3imin + 2

	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2)
	vline!([BV[idxpH2]]; c=N+2, lw=3, lab="", ls=:dot)
	hline!([dpH[idxpH2]]; c=N+2, lw=3, lab="", ls=:dot)
	scatter!([BV[idxpH2]], [dpH[idxpH2]]; 
			 lab="", m=:star, markersize=9, c=:black, markerstrokewidth=2)
	annotate!(0.95*BV[idxpH2], -1.8*dpH[idxpH2], ("ΔpH/Δt: $(round(SF*dpHslope;digits=3)) after $(round(Int,BV[idxpH2])) BV", 9, :right))	
	
	# inset
	plot!(BV[2:end], dpH; inset = (1, bbox(0.075, 0.0, 0.65, 0.55)), 
		  subplot=2, xlims=(60, 95), ylims=(-0.02, 0.02),yticks=-0.02:0.01:0.02,
		  lab="model", c=:green, lw=3, box=:box, grid=:off, ylabel="",
		  legendposition=:outertop, legendcolumns=2, legendfontsize=8)
	
	hline!([0.01]; subplot=2, lab="", ls=:dot, lw=1, c=:black)
	plot!(BV[imin+1:end], noisy_dpH; 
		  subplot=2, c=:darkblue, ls=:solid, lw=2, lab="moving average")
	vline!([BV[idxpH2]]; c=N+2, lw=2.5, lab="", ls=:dot, subplot=2)
	hline!([noisy_dpH[idxpH2-imin]]; 
		   c=N+2, lw=2.5, lab="", ls=:dot, subplot=2)
	
	# plot!(BVup[1:end-1], noisy_dpH_up; 
		  # subplot=2, c=:orange, ls=:solid, lw=2, lab="pH probe")
	# noisy_dpH_up = diff(moving_average(noisy_pH_up, 2w))./diff(BVup)
	# moving average filter
	# plot!(BV[idxmin:end], moving_average(noisy_pH, w); 
	# 	  subplot=2, lw=2, c=:darkred, lab="mov. average"

	
	# (c) Breakthrough
	Cmax = 4
	p3 = plot(xlims=(0, Vmax), xlabel="Treated volume [BV]", ylims=(-0.01, Cmax),
			  ylabel=L"Concentration, $C/C_0$", grid=:off, 
			  legendposition=:topleft)
	
	plot!(BV[1:idxpH], Na[1:idxpH]; lab="Na⁺", c=1, lw=3)
	plot!(BV[idxpH:end], K[idxpH:end]; lab="K⁺", c=2, lw=3, ls=:dot, alpha=0.9)
	plot!(BV[idxpH:end], Na[idxpH:end]; lab="", c=1, lw=3, alpha=0.7)
	
	hline!([IX.:Cₓ₀load[:Na]*Mi[:Na]*1e3]; lab="", c=:black, lw=0.75, ls=:dot)
	
	plot!(BV[1:idxpH], Na[1:idxpH]; fillrange=0, c=1, alpha=0.2, lab="")
	plot!(BV[1:idxpH2], Na[1:idxpH2]; fillrange=0, c=1, alpha=0.5, lab="")
	
	# plot!(BV[1:end], Mg[1:end]; lab=lab[4], c=col[4], lw=3, ls=:dot, alpha=0.9)
	

	# breakthrough detection
	vline!([BV[idxpH]]; lab="", c=N+1, ls=:dot, lw=2)
	# hline!([Na[idxpH]]; lab="", c=cp[N+1], lw=2, ls=:dot)
	annotate!(0.95*BV[idxpH2], Na[idxpH2],
			  ("Na⁺ breakthrough after $(round(Int,BV[idxpH2])) BV: $(round(Int,Na[idxpH2]*prob.ps[Cref]*Mi[:Na]*1e3)) mg/L", 9, :right))

	vline!([BV[idxpH2]]; lab="", c=N+2, ls=:dot, lw=3)
	hline!([1]; lab="", c=:black, lw=1, ls=:dot)
	# scatter!([BV[idxpH2]], [Na[idxpH2]];
	# 		 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	# annotate!(0.95*BV[idxpH2], 0.6Cmax, 
	# 		  ("BT at $(round(BV[idxpH2];digits=1)) L: $(round(Int, Na[idxpH2])) mg ($(round(Naint[idxpH2-1];digits=1)) mg/L)", 9, :right))
	
	
	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.4, 0.4, 0.2]), 
				  dpi = 300, size=(740, 950), #xlabel="Time [min]",
				  suptitle="pH change threshold vs. pH slope change cycle stopping")
	# figure |> x-> aside(x; v_offset=0)
	savefig("pH_slope_control.png")
	plot!()
end

# ╔═╡ c6ea2e44-3dc0-456e-bca9-b947468810c3
md"""
### Na concentration
"""

# ╔═╡ 6c6b75d9-d4e6-4a81-bcf7-58f1dcd8595e
md"""
For very low thresholds, the difference in pH measured at breakthrough is undetectable.

Depending on the resolution and update time of the Na analyzer, the detection can be swift.

The lowest range of detection is 10$^{-6}$ M or ~0.02 ppm Na (23 ppb). The major issue of ion selective electrodes is the interference with other ions such as H+, K+ and Li+. If Na is the first ion to break through the column, and in absence of Li ions, sodium analyzers are an alternative to pH sensors.
"""

# ╔═╡ 7a82894f-a68c-4b21-9434-0e1a463729c8
md""" Threshold for Na breakthrough: 
$(@bindname Na_threshold Select([0.02 => "20 ppb", 0.1 => "100 ppb", 1 => "1 ppm", 10 => "10 ppm", 100 => "100 ppm"]))
"""

# ╔═╡ 3fe68f2c-4db5-403a-ba68-b5ac501305f3
let
	BV = sol[t]*prob.ps[τ]*prob.ps[Q]/prob.ps[V] 
	BVmax = 140
	Nt = length(BV)

	pH = -log10.(sol(:, 1, dv=C(t,x)[1])*prob.ps[Cref])
	Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref] # [eq/L]

	# breakthrough detection in ppb!
	maxNa = Na_threshold

	idxmin = 3
	idxNa = findfirst(x -> x >= maxNa, Na[idxmin:Nt]*Mi[:Na]*1e3) + idxmin - 1
	NaBT = Na[idxNa]*Mi[:Na]*1e3
	NaBTint = integrate(BV[1:idxNa], Na[1:idxNa]*Mi[:Na]*1e3)

	# Resin capacity remaining and restored
	@show Plot.capacity(prob, sol; comp=1)[idxNa] - 1
	regprob = remake(regtemp; u0=IX.update_u0(prob, sol, sol[t][idxNa]; reverse=true));
	solreg = solve(regprob, Rodas5P(), saveat=0.01);
	RCidx = findfirst(x -> x >= Vreg, solreg[t]*regprob.ps[τ]*regprob.ps[Q]*1e3)
	@show Plot.capacity(regprob, solreg; comp=1)[RCidx]

	# pH
	pHmax = 8
	
	p1 = plot(xlims=(0, BVmax), xlabel="", ylims=(0, pHmax), ylabel="pH")
	plot!(BV[idxmin:idxNa], pH[idxmin:idxNa]; lab="", c=:grey, lw=3, grid=:off)
	plot!(BV[idxNa:end], pH[idxNa:end]; lab="", c=:grey, lw=3, ls=:dot, alpha=0.7)

	# initial pH
	pH0 = -log10(sum([IX.:Cₓ₀load[c[i]] for i in 2:N]))
	prob.ps[qmax]*(1-prob.ps[ε])*prob.ps[V]*1e3
	hline!([pH0]; lab="", c=:black, ls=:dot, lw=0.75)
	annotate!(2, 1.5*pH0, ("initial pH: $(round(pH0;digits=1))", 9, :left))

	vline!([BV[idxNa]]; lab="", c=:black, lw=2.5, ls=:dashdot)
	hline!([pH[idxNa]]; lab="", c=:black, lw=2.5, ls=:dashdot)
	scatter!([BV[idxNa]], [pH[idxNa]]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)
	annotate!([0.95BV[idxNa]], [0.95*pHmax], ("pH after $(round(Int,BV[idxNa])) BV: $(round(pH[idxNa];sigdigits=3))", 9, :right))
	

	# Breakthrough
	Cmax = 4
	p2 = plot(xlims=(0, BVmax), xlabel="Treated volume [BV]", ylims=(0, Cmax), 
			  grid=:off, ylabel=L"Concentration, $C/C_0$")
	plot!(BV[1:idxNa], [prob.ps[Cref]*Mi[c[i]]/Zi[c[i]]*1e3*1e-2*
		  sol(:, 1, dv=(C(t,x)[i]))[1:idxNa] for i in 2:N];
		  lab=permutedims(Plot.label[2:N]), c=permutedims(1:N-1), lw=3)
	plot!(BV[idxNa:end], [prob.ps[Cref]*Mi[c[i]]/Zi[c[i]]*1e3*1e-2*
		  sol(:, 1, dv=(C(t,x)[i]))[idxNa:end] for i in 2:N];
		  lab="", c=permutedims(1:N-1), lw=3, ls=:dot, alpha=0.7)

	hline!([IX.:Cₓ₀load[:Na]*Mi[:Na]*1e1]; lab="", c=:black, lw=1, ls=:dot)
	
	vline!([BV[idxNa]]; lab="", c=:black, ls=:dashdot, lw=2.5)
	hline!([NaBT/100]; lab="", c=:black, ls=:dashdot, lw=2.5)
	scatter!([BV[idxNa]], [Na[idxNa]*Mi[:Na]*1e1]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=2)

	annotate!(0.95*BV[idxNa], 0.5 + NaBT*1e-2, 
			  ("Na⁺ breakthrough after $(round(Int,BV[idxNa])) BV: $(round(Int,NaBT)) mg/L", 9, :right))

	# Resin capacity
	p3 = Plot.distribution(prob, sol; time="BV")

	# solx = prob.problem_type.discretespace.axies[x]
	# plot!(BV, 
	# 	  [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]* 100 ;
	# 	  lab="", c=:grey, lw=3, grid=:off, title="Resin capacity", 
	# 	  xlabel="Treated volume [L]", xlims=(0, BVmax),
	# 	  ylabel="Resin capacity, %", ylims=(-0.5, 101), )
	# capacity = 1 .- [integrate(solx, sol(j, solx, dv=q(t,x)[1])) for j in sol[t]]
	# plot!(BV, 100capacity; lab="", c=:darkred, lw=3, ls=:dot)
	vline!([BV[idxNa]]; lab="", c=:black, ls=:dashdot, lw=2.5)
	# annotate!(0.98*BV[idxNa], 95, 
	# 		  ("saturation at breakthrough: $(round(100capacity[idxNa];digits=1))%", 9, :right))
	
	figure = plot(p1, p2, p3, layout=grid(3, 1, heights=[0.2, 0.4, 0.4]), 
				  dpi = 300, size=(680, 700), xlims=(0, BVmax), xticks=0:20:200,
				  suptitle="Cycle stopping with Na⁺ breakthrough threshold")
	# figure |> x -> aside(x; v_offset=0)
	savefig("Na_control.png")
	plot!()
end

# ╔═╡ c2ecf3a2-4cc2-4d1c-a729-6f412de22b47
function Na_stop_fixed_reg_cycles(; Nc=1, BVstop=20, BVregstop=3)

	loadprob = prob

	if prob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x]
	elseif prob.problem_type.disc.grid_align == MethodOfLines.EdgeAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x][2:end]
	end

	load = []; regn = []; temp = []; tempreg = []; loadtemp = []
	loads = []; regns = []; sols = []; regs = []; ends1 = []; ends2 = []

	for n in 1:Nc
		
		# Loading
		loadparams = IX.update(params, Q => 120/60*1e-6, 
								  Cₓ₀ => [IX.Cₓ₀load[i] for i in c],
								  qmax => p0[1], K => vcat(1.0, p0[2:N]))

		# Initial cycle
		if n == 1
			prob = remake(loadprob, p=loadparams)
		elseif n > 1
			solreg = temp; reg = tempreg
			end2 = findfirst(x->x>=BVregstop,solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])
		
			prob = remake(loadprob, p=loadparams, 
						  u0=IX.update_u0(reg, solreg, solreg[t][end2];reverse=true))
		end

		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		# Na stop
		Na = sol(:, 1, dv=C(t,x)[2])*prob.ps[Cref] # [eq/L]
		maxNa = Na_threshold
		idxmin = 3
		end1 = findfirst(x -> x >= maxNa, Na[idxmin:end]*Mi[:Na]*1e3) + idxmin - 1
		# end1 = findfirst(x -> x >= BVstop, sol[t]*prob.ps[τ]*prob.ps[Q]*1e3)

		use = 100*integrate(grid, sol(sol[t][end1], grid, dv=q(t,x)[1])) 
		# use = 100 * Plot.capacity(prob, sol; comp=1)[end1]
		push!(load, use)
		push!(loads, prob)
		push!(sols, sol)
		push!(ends1, end1)
	
		
		# Regeneration
		if regprob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x]
		elseif regprob.problem_type.disc.grid_align ==MethodOfLines.EdgeAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x][2:end]
		end
		regparams = IX.update(params, Q => 25/60*1e-6, 
							  Cₓ₀ => [IX.Cₓ₀reg[i] for i in c], 
							  Cref => prob.ps[Cref], 
							  qmax => p0[1], K => vcat(1., p0[2:N]))

		reg = remake(regprob, p=regparams, 
					 u0=IX.update_u0(prob, sol, sol[t][end1]; reverse=true))
		
		solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
		temp = solreg; tempreg = reg
		grid = solreg[x]

		end2 = findfirst(x -> x >= BVregstop,solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])

		# cap = 100 * Plot.capacity(reg, sol; comp=1, time=solreg[t][1:end2])
		cap = 100*integrate(grid, solreg(solreg[t][end2], grid, dv=q(t,x)[1]))
		push!(regn, cap)
		push!(regns, reg)
		push!(regs, solreg)
		push!(ends2, end2)
	end
		
	return (; load, regn, Nc, BVstop, BVregstop, loads, regns, sols, regs, ends1, 
			ends2)
end

# ╔═╡ db6d781f-8db9-4e95-b965-2aa748fe7d97
@time Na_cycles = Na_stop_fixed_reg_cycles(Nc=5, BVstop=BVstop, BVregstop=3.)

# ╔═╡ f19c2d96-3fc4-4011-86f3-044f02b26c43
let 
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = Na_cycles

	@show 100-load[end]
	@show regn[end]
	@show sols[end][t][ends1[end]]*prob.ps[τ]*prob.ps[Q]/prob.ps[V]
	@show BVregstop #/ regns[end].ps[V]*1e-3
	@show sols[end](:,1, dv=C(t,x)[2])[ends1[end]]*100
end

# ╔═╡ e208347b-5aef-4181-82c0-ebc5b82ec134
25/prob.ps[V]/60/1000

# ╔═╡ 5243b1f8-a698-45a8-a650-b5fb5039b20f
md"""
## Model-based control strategies
"""

# ╔═╡ 7ea4c7fb-f496-426f-b3ec-fad0302e5f88
md"""
Output calculations:
- Resin capacity utilisation
- Breakthrough and elution curves 
- Breakthrough times and duration of each cycle
- Likelihood of breakthrough, updated every 10 minutes
"""

# ╔═╡ 216a4a97-26ab-4e9f-aa43-55f9e8c59aa3
md"""
We recreate the model for the regeneration process with a variable concentration of regenerant Cₓ₀:
"""

# ╔═╡ 08158174-6453-4e72-ad0d-4110d7275379
# Boundary conditions
bcsreg_var = let
	@unpack Cₓ₀ = IX.pars(N=N)
	Dx = IX.Dx
	bcs = Array{Symbolics.Equation}(undef, 2N + 4N) # 2N ICs + N BCs
	for (i, j) in zip(1:N, c)
		# Initial conditions
		bcs[i] 		= C(0.0, x)[i] ~ IX.C₀[j] 
		bcs[i + N] 	= q(0.0, x)[i] ~ IX.q₀load[j] 
		# Boundary conditions - Dirichlet for Cₓ₀, Neumann for the rest
		bcs[i + 2N] = C(t, 0.0)[i] ~ Cₓ₀[i]/Cref
		bcs[i + 3N] = Dx(C(t, 1.0)[i]) ~ 0.
		bcs[i + 4N] = Dx(q(t, 0.0)[i]) ~ 0.
		bcs[i + 5N] = Dx(q(t, 1.0)[i]) ~ 0.
	end
	bcs
end

# ╔═╡ 21aa8bc8-8dcd-4fad-8cc8-e3a49c020577
regmodel2 = IX.sys_average_reg(N=N, bcs=bcsreg_var, domain=(t∈(0, 50.), x∈(0, 1.)));

# ╔═╡ a70be775-8c77-46a2-835d-0eaa378bb4ae
@time regtemp2 = IX.MOLdiscretize(regmodel2, dx=dx; grid_align=grid_align);

# ╔═╡ 142086b2-4d1d-486e-b1ec-9df4810346ef
regparams = IX.update(params, 
					  Cₓ₀ => vcat(1.4, [IX.zero_eq for i in 2:N]),
					  Cref => prob.ps[Cref],
					  Q => 25/60/1e6,
					  );

# ╔═╡ 9c28e7f0-b1be-436a-8293-6b139d203132
regprob2 = remake(regtemp2; 
				  p=regparams,
				  u0=IX.update_u0(prob, sol, sol[t][end]; reverse=true));

# ╔═╡ db8fbe59-5aa7-4dba-918b-189084c93969
@time solreg2 = solve(regprob2, Rodas5P(), saveat=0.01);

# ╔═╡ 84926400-3c00-4024-b4bd-d252d5a911cb
md"""
With 10% of the initial Na concentration, the regeneration efficiency of HCl will be given by the restored resin capacity by using a certain volume of HCl at different concentrations.
"""

# ╔═╡ f6af7beb-aa6e-4318-b975-05191f1afa1c
Plot.distribution(regprob, solreg; elution=true, xlims=(0, 6), xticks=0:1:6, time="BV")

# ╔═╡ 75bbed03-a241-4879-923e-f0eb167361ef
Plot.distribution(regprob2, solreg2; elution=true, xlims=(0, 6), xticks=0:1:6, time="BV")

# ╔═╡ 43465bde-3989-4833-93a2-669ad7f52612
md"""
### Regeneration vs. concentration
"""

# ╔═╡ 5ee56023-c336-4bc6-a015-419ac7805ece
md"""
We calculate the regeneration of capacity by using 1 L of regenerant in exhausted resins.
"""

# ╔═╡ 9bc995d7-44a1-44e8-a5d9-4534e5f92914
Ci_values = 0.002:0.04:(100/36.5)

# ╔═╡ 568c355c-aca0-4c24-8951-f7c10ae9abee
@time cap, eff = let
	cap = zeros(length(Ci_values))
	eff = zeros(length(Ci_values))
	
	for (i, Ci) in zip(1:length(Ci_values), Ci_values)
	
		# Adapt regenerant concentration dosage
		regparams = IX.update(params, 
						  Cₓ₀ => vcat(Ci, [IX.zero_eq for i in 2:N]),
						  Cref => prob.ps[Cref],
						  Q => 25/60/1e6,
						  );
		# Remake problem
		regprob = remake(regprob2; 
					  p=regparams,
					  u0=IX.update_u0(prob, sol, sol[t][end]; reverse=true));
		# Solve problem
		solreg = solve(regprob, Rodas5P(), saveat=0.01);

		# Extract regeneration efficiency
		RCidx = findfirst(x -> x >= 1., solreg[t]*regprob.ps[τ]*regprob.ps[Q]*1e3)
		
		cap[i] = Plot.capacity(regprob, solreg; time=solreg[t], comp=1)[RCidx]

		RC = cap[i]*prob.ps[qmax]*prob.ps[V]*(1-prob.ps[ε])*1e3 # restored capacity 

		eff[i] = RC./Ci # restored capacity [eq] / regenerant used [eq]
	end

	cap, eff
end

# ╔═╡ a5d7e48c-af63-46f5-bdae-943b3b3f7b71
let
	plot(dpi=300, size=(680, 420), widen=false, grid=:off, margin=3mm)

	# Regenerated capacity fraction [%]
	plot!(Ci_values*36.5, 100*cap; 
		  lw=3, c=:grey, lab="capacity with HCl as regenerant", legendposition=:best,
		  xlims=(0, 100), xlabel="Regenerant concentration [g/L]", xticks=0:25:100,
		  ylabel="Regeneration efficiency [%]\n Regenerated column capacity [%]", 
		  ylims=(0, 100),
		  # title="Regeneration efficiency vs. regenerant concentration",
		 )
	plot!(Ci_values*49.0, 100*cap; 
		  lw=3, c=:grey, lab="capacity with H₂SO₄ as regenerant", 
		  alpha=0.8, ls=:dash)

	# H₂SO₄ 
	C50 = findfirst(x -> x >= 50, Ci_values*36.5) - 1
	hline!([100*cap[C50]]; lab="", c=:black, ls=:dash)
	vline!([50]; lab="", c=:black, ls=:dash)
	scatter!([50], [100*cap[C50]]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=3)
	annotate!(1.52*50, 0.85*100*cap[C50], 
			 ("$(round(Int, 100*cap[C50]))% capacity: ~$(round(50/36.5;digits=1)) eq HCl 5% (50 g/L)  \n\n vs. ~$(round(50/36.5*49/36.5;digits=1)) eq of H₂SO₄ 50 g/L required", 9, :center))

	# Efficiency
	plot!(Ci_values*36.5, 100eff;
		  lw=3, c=2, lab="regeneration efficiency", ls=:dashdotdot, alpha=0.8)
	plot!(Ci_values*36.5, 100*cap.*eff; 
		  c=:black, ls=:dot, lw=2, lab="capacity × efficiency", alpha=0.7)
	hline!([100*eff[C50]]; lab="", c=2, ls=:dash, alpha=0.8)
	scatter!([50], [100*eff[C50]]; 
			 lab="", c=2, m=:xcross, markersize=6, markerstrokewidth=3, alpha=0.8)
	annotate!(1.55*50, 1.02*40, 
			  ("~$(round(Int, 100eff[36]))% efficiency: $(round(eff[36];digits=1)) eq H⁺ restored\n\n per 1 eq HCl 5% used", 9, :center))
	# annotate!(25, 20, ("capacity × efficiency", 9, :center))
	
	savefig("efficiency.png")
	plot!()
end

# ╔═╡ f98541fe-ae35-4f95-af2d-e4930cef2138
md"""
### Regeneration efficiency
"""

# ╔═╡ 0af62d26-4453-43ab-824b-d85ca7c43b30
md"""
Regeneration efficiency decreases with increasing regenerant concentration because each extra unit of regenerant only gives a small increase in regenerated capacity.
"""

# ╔═╡ 9a7dd642-8000-45d0-81af-ba168ed0860f
md"""
### Efficiency vs. contact time
"""

# ╔═╡ a3ac46d3-f429-4b6a-9404-39d40d846186
Ci_values2 = [10, 20, 36.5, 50] ./ 36.5  # 0.1:0.4:(50/36.5)

# ╔═╡ d6a6c468-a9e6-4ded-a6ac-8847c3f799c5
100 .- 3Ci_values2 .+ 1500 .* Ci_values2 .|> x -> round(Int, x)

# ╔═╡ 8613c163-850e-4af0-8aa8-20f2bd7d715d
@time cap2, eff2 = let
	cap = zeros(length(solreg2[t]), length(Ci_values2))
	eff = zeros(length(solreg2[t]), length(Ci_values2))
	
	for (i, Ci) in zip(1:length(Ci_values2), Ci_values2)
	
		# Adapt regenerant concentration dosage
		regparams = IX.update(params, 
						  Cₓ₀ => vcat(Ci, [IX.zero_eq for i in 2:N]),
						  Cref => prob.ps[Cref],
						  Q => 25/60/1e6,
						  );
		# Remake problem
		regprob = remake(regprob2; 
					  p=regparams,
					  u0=IX.update_u0(prob, sol, sol[t][end]; reverse=true));
		# Solve problem
		solreg = solve(regprob, Rodas5P(), saveat=0.01);

		# Extract regeneration efficiency
		# RCidx = findfirst(x -> x >= 1., solreg[t]*regprob.ps[τ]*regprob.ps[Q]*1e3)
		
		cap[:, i] = Plot.capacity(regprob, solreg; time=solreg[t], comp=1)

		RC = cap[:, i]*prob.ps[qmax]*prob.ps[V]*(1-prob.ps[ε])*1e3 

		eff[:, i] = RC./(Ci*solreg[t]*regprob2.ps[τ]*regprob2.ps[Q]*1e3)
	end

	cap, eff
end;

# ╔═╡ 94358d87-d685-435a-98f0-734b0d486b73
md"""
For a 40 min contact time at approx. 25 mL/min flow rate, the regeneration efficiency is about 80 % for a completely exhausted SAC resin with 1 L of 50 g/L HCl. For lower concentrations of HCl, the contact time needed to reach the same regeneration efficiency is much longer.
"""

# ╔═╡ c93af4de-9aba-4360-973f-4de6215aa56e
let
	BV = solreg2[t]*regprob2.ps[τ]*regprob2.ps[Q]/regprob2.ps[V] # bed volumes
	plot(dpi=300, size=(680, 420), widen=false, margin=3mm, box=:axes, grid=:off)

	idx = round.(Int, 1600 .- 120 * Ci_values2)
	
	
	for (i, Ci) in zip(1:length(Ci_values2), Ci_values2)
		if i == 3
			plot!(BV[1:idx[i]], 100*cap2[1:idx[i], i]; 
			  lw=3, c=:grey, alpha=2Ci, lab="", ls=:dash)
			annotate!(1.05*BV[idx[i]], 100cap2[idx[i], i], 
					  ("1.0 eq/L (37 g/L HCl or 50 g/L H₂SO₄)", 8, :left))
		# H2SO4
		elseif i == 4
			plot!(BV[1:idx[i]], 100*cap2[1:idx[i], i]; 
			  lw=3, c=:grey, alpha=2Ci, lab="regenerated capacity")
			annotate!(1.05*BV[idx[i]], 1.02*100cap2[idx[i], i], 
					  ("1.4 eq/L (50 g/L HCl or 69 g/L H₂SO₄)", 8, :left))
		else
		# HCl
		plot!(BV[1:idx[i]], 100*cap2[1:idx[i], i]; 
			  lw=3, c=:grey, alpha=2Ci, lab="",
			  xlims=(0, 9), xticks=0:3:9, 
			  xlabel="Bed volumes of spent regenerant [BV]", 
			  ylabel="Regeneration efficiency [%]\n Regenerated column capacity [%]",
			  ylims=(0, 100), legendposition=(0.75, 0.35)
			  # title="Regeneration efficiency vs. spent regenerant volume"
			 )
		annotate!(1.05*BV[idx[i]], 100cap2[idx[i], i], 
				("$(round(Ci;digits=1)) eq/L ($(round(Int,Ci*36.5)) g/L HCl)", 8, :left))
		end
	end
	hline!([89]; lab="", c=:black, ls=:dash)
	vline!([3]; lab="", c=:black, ls=:dash)
	scatter!([3], [89]; 
			 lab="", c=:black, m=:xcross, markersize=6, markerstrokewidth=3)
	idx = findfirst(x -> x >= 3, BV)
	annotate!(2.8, 93, 
			  ("~$(round(Int, 100cap2[idx, end]))% capacity restored", 9, :right))
	annotate!(3.2, 7, 
			  ("3 BV HCl 5% ≈ 40 min contact time, ~$(round(Int, 100eff2[idx, end]))% efficiency", 9, :left))

	# Efficiency 
	idxmin = 2
	# plot!(BV[idxmin:end], [mean(100eff2[idxmin:i, end]) for i in idxmin:length(BV)]; 
	# 	  lab="efficiency, 50 g/L HCl (mean)", 
	# 	  c=2, lw=3, ls=:dashdotdot, alpha=0.8)
	idxmin = 7
	plot!(BV[idxmin:end], 100eff2[idxmin:end, end];
		  c=2, lw=3, ls=:dashdotdot, alpha=0.8, lab="efficiency, 50 g/L HCl")
	plot!(BV, 100eff2[:, end].*cap2[:, end]; 
		  lab="capacity × efficiency", ls=:dot, lw=2, c=:black, alpha=0.7)

	# Efficiency for non-exhausted columns 
	# plot!(BV, 100Plot.capacity(regprob, solreg; time=solreg[t], comp=1))
	
	savefig("efficiency2.png")
	plot!()
end

# ╔═╡ 12f6540f-4de3-40bb-aff2-524705a6f680
md"""
### Calculation of the operating capacity
"""

# ╔═╡ da593e8b-ff28-4af7-af39-03d20ffe12bd
let
	@unpack load, regn, loads, regns, sols, regs, Nc, BVstop, BVregstop, ends1, ends2 = fixed_time_cycles

	plot(dpi=300, size=(680, 400), widen=false, margin=3mm)

	alpha = 0.4:0.1:(Nc./ 10 + 0.4) |> collect; #@show alpha

	for i in 1:Nc
		if i in vcat(1, 2:5, 10)
			plot!(sols[i][t][1:ends1[i]]*loads[i].ps[τ]*loads[i].ps[Q]/loads[i].ps[V], 
				  100*sols[i](sols[i][t][1:ends1[i]], 1, dv=C(t,x)[2]);
				  lab="", c=1, alpha=alpha[i], lw=4)
			annotate!(1.03*sols[i][t][ends1[i]]*loads[i].ps[τ]*loads[i].ps[Q]/loads[i].ps[V],
			100*sols[i](sols[i][t][ends1[i]], 1, dv=C(t,x)[2]),
			("Cycle $i", 9, :left))
		else
			plot!(sols[i][t][1:ends1[i]]*loads[i].ps[τ]*loads[i].ps[Q]/loads[i].ps[V], 
			  100*sols[i](sols[i][t][1:ends1[i]], 1, dv=C(t,x)[2]);
			  lab="", c=1, alpha=1, lw=4)
		end
	end
	hline!([5]; lab="", c=:black, ls=:dash, lw=1.5, grid=:off)
	BTidx = findfirst(x -> x >= 0.05, sols[end](:, 1, dv=C(t,x)[2]))
	vline!([sols[end][t][BTidx]*loads[end].ps[τ]*loads[end].ps[Q]/loads[end].ps[V]]; 
		   lab="", c=:black, ls=:dash, lw=1.5)
	scatter!([sols[end][t][BTidx]*loads[end].ps[τ]*loads[end].ps[Q]/loads[end].ps[V]], [5]; 
			 lab="", m=:xcross, markersize=7, c=:black, markerstrokewidth=4)
	annotate!(0.95 * sols[end][t][BTidx]*loads[end].ps[τ]*loads[end].ps[Q]/loads[end].ps[V], 
			  1.2 * 5, ("Total volume at 5% Na⁺ breakthrough: $(round(sols[end][t][BTidx]*loads[end].ps[τ]*loads[end].ps[Q]/loads[end].ps[V];digits=1)) BV", 9, :right))
	plot!(xlabel="Bed volumes of treated effluent [BV]", 
		  ylabel="Na⁺ breakthrough concentration [mg/L]", 
		  xlims=(0, 90), xticks=0:15:105, ylims=(0, 15))
		  # title="Column tests for operating capacity calculation"
	savefig("capacity.png")
	plot!()
end

# ╔═╡ 40b92aab-5ecf-482e-9d31-7b26cd84500f
md"""
## Trade-offs for optimization
"""

# ╔═╡ 79a969ca-5965-4820-aa3c-e2b51dd38c61
md"""
### Effluent quality vs. regenerant volume
"""

# ╔═╡ 1b47c432-dfb7-45ad-bffc-f548a9c4bf4b
function Na_Vreg_cycles(; Nc=1, Na_load=0.1, Vreg=3.)

	loadprob = prob

	if prob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x]
	elseif prob.problem_type.disc.grid_align == MethodOfLines.EdgeAlignedGrid()
		grid = prob.problem_type.discretespace.grid[x][2:end]
	end

	temp = []; tempreg = []; loadtemp = []
	BVstop = []; BVregstop = []

	cap_load = Vector{Any}(undef, Nc)
	loads = Vector{Any}(undef, Nc)
	sols = Vector{Any}(undef, Nc)
	ends1 = Vector{Any}(undef, Nc)

	cap_reg = Vector{Any}(undef, Nc)
	regns = Vector{Any}(undef, Nc)
	regs = Vector{Any}(undef, Nc)
	ends2 = Vector{Any}(undef, Nc)
	

	for n in 1:Nc
		
		# Loading
		loadparams = IX.update(params, Q => 120/60*1e-6, 
								  Cₓ₀ => [IX.Cₓ₀load[i] for i in c],
								  qmax => p0[1], K => vcat(1.0, p0[2:N]))

		# Initial cycle
		if n == 1
			prob = remake(loadprob, p=loadparams)
		elseif n > 1
			solreg = temp; reg = tempreg
			end2 = findfirst(x -> x >= Vreg, 
							 solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])
		
			prob = remake(loadprob, p=loadparams, 
						  u0=IX.update_u0(reg, solreg, solreg[t][end2];reverse=true))
		end

		sol = solve(prob, Rodas5P(), saveat=1, progress=true, progress_steps=1)

		end1 = findfirst(x -> x >= Na_load, sol(:, 1, dv=C(t,x)[2]))

		cap_load[n] = 100*integrate(grid, sol(sol[t][end1], grid, dv=q(t,x)[1])) 
		# use = 100 * Plot.capacity(prob, sol; comp=1)[end1]
		loads[n] = prob
		sols[n] = sol
		ends1[n] = end1
		# push!(BVstop, sol[t][end1]*prob.ps[τ]*prob.ps[Q]/prob.ps[V])
		
		# Regeneration
		if regprob.problem_type.disc.grid_align == MethodOfLines.CenterAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x]
		elseif regprob.problem_type.disc.grid_align ==MethodOfLines.EdgeAlignedGrid()
			grid = regprob.problem_type.discretespace.grid[x][2:end]
		end
		regparams = IX.update(params, Q => 25/60*1e-6, 
							  Cₓ₀ => [IX.Cₓ₀reg[i] for i in c], 
							  Cref => prob.ps[Cref], 
							  qmax => p0[1], K => vcat(1., p0[2:N]))

		reg = remake(regprob, p=regparams, 
					 u0=IX.update_u0(prob, sol, sol[t][end1]; reverse=true))
		
		solreg = solve(reg, Rodas5P(), saveat=0.01, progress=true, progress_steps=1)
		temp = solreg; tempreg = reg
		grid = solreg[x]

		end2 = findfirst(x -> x >= Vreg, solreg[t]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])

		# cap = 100 * Plot.capacity(reg, sol; comp=1, time=solreg[t][1:end2])
		cap_reg[n] = 100*integrate(grid, solreg(solreg[t][end2], grid, dv=q(t,x)[1]))
		regns[n] = reg
		regs[n] = solreg
		ends2[n] = end2
		# push!(BVregstop, solreg[t][end2]*reg.ps[τ]*reg.ps[Q]/reg.ps[V])
	end
		
	return (; Nc, cap_load, loads, sols, cap_reg, regns, regs, ends1, ends2, 
			BVstop, BVregstop)
end

# ╔═╡ c8303901-9e32-4ba0-85a9-a230322d6830
@time Naload2, Vreg2, VloadNa2, VregnNa2 = let
	Nc = 5
	Naload = [0.02, 0.1, 0.5] #0.01:0.05:0.26
	Vreg = 0:1:5 # [0., 1., 2.0, 3.0, 4.0]  
	
	Vload = zeros(length(Naload), length(Vreg), 2)
	Vregn = zeros(length(Naload), length(Vreg), 2)
	

	for (i, Na) in zip(1:length(Naload), Naload)
		for (j, Vr) in zip(1:length(Vreg), Vreg)
			Na_Vreg_cycles_test = Na_Vreg_cycles(Nc=Nc, Na_load=Na, Vreg=Vr)
		
			@unpack loads, regns, sols, regs, ends1, ends2 = Na_Vreg_cycles_test
	
			# Treated volume
			Vload[i, j, :] = [sols[k][t][ends1[k]] *
				loads[k].ps[τ]*loads[k].ps[Q]/loads[k].ps[V] for k in (1, Nc)]
	
			# Regenerant volume
			Vregn[i, j, :] = [regs[k][t][ends2[k]] * 
				regns[k].ps[τ]*regns[k].ps[Q]/regns[k].ps[V] for k in (1, Nc)]
	
			# Resin capacity
			# RC[i, j, :] = [Plot.capacity(regns[k], regs[k]; comp=1, time=regs[k][t][1:ends2[k]])[end] for k in (1, Nc)]
		end
	end

	Naload, Vreg, Vload, Vregn
end

# ╔═╡ 3839332c-384d-44b1-b48c-349f08bbc9a9
let
	Na = Naload2; Vreg = Vreg2; V = VloadNa2; Vregn = VregnNa2

	label = ["$(round(Int, 100*Na[j])) mg/L" for j in 1:length(Na)]

	alpha = 0.5 .+ (Na .- minimum(Na)) ./ (maximum(Na) .- minimum(Na))

	Vline = Vector{Any}(undef, length(Na)) 
	[Vline[i] = CubicSplineInterpolation(Vreg, V[i, :, end], extrapolation_bc=Line()) for i in 1:length(Na)]

	Vreg = 0:0.01:5.2
	
	plot(dpi=300, size=(680, 400), widen=false, grid=:off, margin=3mm)

	for i in 1:length(Na)
		plot!(Vreg, Vline[i](Vreg); 
			  lw=4, c=1, alpha=alpha[i], lab="",
			  xlabel="Spent regenerant volume [BV]", xlims=(0, 6),
			  ylabel="Treated effluent volume [BV]", ylims=(0, 90))

		annotate!(1.05*Vreg[end], Vline[i](Vreg[end]), 
				 (label[i], 9, :left))
	end

	Vreg = 0:0.01:5.2

	plot!(Vreg[2:end], diff(Vline[3](Vreg))./diff(Vreg).*Vreg[2:end]; 
		  lw=3, ls=:dash, c=2, alpha=0.8, lab="")
	annotate!(2., 40, (L"marginal resource efficiency: $\eta_{chem} = V_{\textrm{regen}}\frac{\Delta V_{\textrm{treated}}}{\Delta V_{\textrm{regen}}}$", 10, :left))

	plot!(Vreg, Vline[1](Vreg)./(Vreg .+ 1.5);
		  lw=3, ls=:dashdotdot, c=3, alpha=0.8, lab="")
	annotate!(3.5, 20, (L"cycle productivity: $\frac{V_{\textrm{treated}}}{V_{\textrm{regen}} + V_{\textrm{waste}}}$", 10, :left))

	# plot!(title="Treated vs. spent regenerant volume & Na⁺ threshold")
	savefig("pareto.png")
	plot!()
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
MethodOfLines = "94925ecb-adb7-4558-8ed8-f975c56a0bf4"
ModelingToolkit = "961ee093-0014-501f-94e3-6117800e7a78"
NumericalIntegration = "e7bfaba1-d571-5449-8927-abc22e82249b"
OrdinaryDiffEq = "1dea7af3-3e70-54e6-95c3-0bf5283fa5ed"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
Interpolations = "~0.16.2"
LaTeXStrings = "~1.4.0"
MethodOfLines = "~0.11.9"
ModelingToolkit = "~10.25.0"
NumericalIntegration = "~0.3.4"
OrdinaryDiffEq = "~6.102.1"
Plots = "~1.41.4"
PlutoLinks = "~0.1.8"
PlutoUI = "~0.7.79"
StatsPlots = "~0.15.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.7"
manifest_format = "2.0"
project_hash = "d093b4017dbaa72d2d17d6a27925e3aa2b6f12c6"

[[deps.ADTypes]]
git-tree-sha1 = "bbc22a9a08a0ef6460041086d8a7b27940ed4ffd"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "1.22.0"
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
git-tree-sha1 = "2eeb2c9bef11013efc6f8f97f32ee59b146b09fb"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.44"

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
git-tree-sha1 = "0761717147821d696c9470a7a86364b2fbd22fd8"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.5.2"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

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
git-tree-sha1 = "54f895554d05c83e3dd59f6a396671dae8999573"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.24.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceAMDGPUExt = "AMDGPU"
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
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Bijections]]
git-tree-sha1 = "a2d308fcd4c2fb90e943cf9cd2fbfa9c32b69733"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.2.2"

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
git-tree-sha1 = "7ad7171d693ae5552ac43862e7f6b61df4471c2b"
uuid = "70df07ce-3d50-431d-a3e7-ca6ddb60ac1e"
version = "1.12.1"
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
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "12177ad6b3cad7fd50c8b3825ce24a99ad61c18f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.26.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

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

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "REPL", "UUIDs"]
git-tree-sha1 = "cfb7a2e89e245a9d5016b70323db412b3a7438d5"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "3.0.2"

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
git-tree-sha1 = "65ea18ada9814f09c5013924c42fe8b53d6ee467"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.10.3"

    [deps.CommonMark.extensions]
    CommonMarkMarkdownASTExt = "MarkdownAST"
    CommonMarkMarkdownExt = "Markdown"

    [deps.CommonMark.weakdeps]
    Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
    MarkdownAST = "d0879d2d-cac2-40c8-9cee-1863dc0c7391"

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

[[deps.Compiler]]
git-tree-sha1 = "382d79bfe72a406294faca39ef0c3cef6e6ce1f1"
uuid = "807dbc54-b67e-4c79-8afb-eafe4df6f2e1"
version = "0.1.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

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
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

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

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

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

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffEqBase]]
deps = ["ArrayInterface", "BracketingNonlinearSolve", "ConcreteStructs", "DocStringExtensions", "FastBroadcast", "FastClosures", "FastPower", "FunctionWrappers", "FunctionWrappersWrappers", "LinearAlgebra", "Logging", "Markdown", "MuladdMacro", "PrecompileTools", "Printf", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLOperators", "SciMLStructures", "Setfield", "Static", "StaticArraysCore", "SymbolicIndexingInterface", "TruncatedStacktraces"]
git-tree-sha1 = "c5fe5125fcba8f98cdc5c7221b6c324883899c07"
uuid = "2b5f629d-d688-5b77-993f-72d75c75574e"
version = "6.213.0"

    [deps.DiffEqBase.extensions]
    DiffEqBaseCUDAExt = "CUDA"
    DiffEqBaseChainRulesCoreExt = "ChainRulesCore"
    DiffEqBaseDynamicQuantitiesExt = "DynamicQuantities"
    DiffEqBaseEnzymeExt = ["ChainRulesCore", "Enzyme"]
    DiffEqBaseFlexUnitsExt = "FlexUnits"
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
    DynamicQuantities = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    FlexUnits = "76e01b6b-c995-4ce6-8559-91e72a3d4e95"
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
git-tree-sha1 = "51f445c25b80c26c83c61520cde214adf8c05227"
uuid = "459566f4-90b8-5000-8ac3-15dfb0a30def"
version = "4.17.0"

    [deps.DiffEqCallbacks.extensions]
    DiffEqCallbacksFunctorsExt = "Functors"

    [deps.DiffEqCallbacks.weakdeps]
    Functors = "d9f16b24-f501-4c13-a1f2-28368ffc5196"

[[deps.DiffEqNoiseProcess]]
deps = ["CommonSolve", "DiffEqBase", "Distributions", "GPUArraysCore", "LinearAlgebra", "Markdown", "PoissonRandom", "QuadGK", "Random", "RecipesBase", "RecursiveArrayTools", "ResettableStacks", "SciMLBase", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "5b5fe6eb9e994ff0de4e84fe149b03c706fd6403"
uuid = "77a26b50-5914-5dd7-bc55-306e6241c503"
version = "5.31.0"

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
git-tree-sha1 = "e421c1938fafab0165b04dc1a9dbe2a26272952c"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.125"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.DynamicPolynomials]]
deps = ["LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Reexport", "StarAlgebras", "Test"]
git-tree-sha1 = "5bfabc3827dfdd164359bad6800c115a81280c00"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.6.6"

[[deps.DynamicQuantities]]
deps = ["DispatchDoctor", "PrecompileTools", "TestItems", "Tricks"]
git-tree-sha1 = "7b098159d57dd03cb83fb5173cfba462d1b5e2de"
uuid = "06fc5a27-2a28-4c7c-a15d-362465fb6821"
version = "1.12.3"

    [deps.DynamicQuantities.extensions]
    DynamicQuantitiesLinearAlgebraExt = "LinearAlgebra"
    DynamicQuantitiesMeasurementsExt = "Measurements"
    DynamicQuantitiesRecursiveArrayToolsExt = "RecursiveArrayTools"
    DynamicQuantitiesSciMLBaseExt = "SciMLBase"
    DynamicQuantitiesScientificTypesExt = "ScientificTypes"
    DynamicQuantitiesUnitfulExt = "Unitful"

    [deps.DynamicQuantities.weakdeps]
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    SciMLBase = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
    ScientificTypes = "321657f4-b219-11e9-178b-2701a2544e81"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.EnumX]]
git-tree-sha1 = "c49898e8438c828577f04b92fc9368c388ac783c"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.7"

[[deps.EnzymeCore]]
git-tree-sha1 = "24bbb6fc8fb87eb71c1f8d00184a60fc22c63903"
uuid = "f151be2c-9106-41f4-ab19-57ee4f262869"
version = "0.8.19"
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
git-tree-sha1 = "9cb7fe11da6adb8683cbacf8aa9b5237941e3a75"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.5+0"

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
git-tree-sha1 = "cac41ca6b2d399adfc95e51240566f8a60a80806"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.1.0+0"

[[deps.FFTA]]
deps = ["AbstractFFTs", "DocStringExtensions", "LinearAlgebra", "MuladdMacro", "Primes", "Random", "Reexport"]
git-tree-sha1 = "65e55303b72f4a567a51b174dd2c47496efeb95a"
uuid = "b86e33f2-c0db-4aa1-a6e0-ab43e668529e"
version = "0.3.1"

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
git-tree-sha1 = "73e879af0e767bd6dfade7c5b09d7b05657a8284"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.30.0"

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
git-tree-sha1 = "cddeab6487248a39dae1a960fff0ac17b2a28888"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "1.3.3"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "9e0fb9e54594c47f278d75063980e43066e26e20"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "83cf05ab16a73219e5f6bd1bdfa9848fa24ac627"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.2.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

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
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Glob]]
git-tree-sha1 = "83cb0092e2792b9e3a865b6655e88f5b862607e2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.4.0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Inflate", "LinearAlgebra", "Random", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "7eb45fe833a5b7c51cf6d89c5a841d5967e44be3"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.14.0"
weakdeps = ["Distributed", "SharedArrays"]

    [deps.Graphs.extensions]
    GraphsSharedArraysExt = "SharedArrays"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

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
git-tree-sha1 = "79d6bd28c8d9bccc2229784f1bd637689b256377"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.14"
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
git-tree-sha1 = "67c6f1f085cb2671c93fe34244c9cccde30f7a26"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.0"

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
git-tree-sha1 = "c0c9b76f3520863909825cbecdef58cd63de705a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.5+0"

[[deps.JuliaFormatter]]
deps = ["CommonMark", "Glob", "JuliaSyntax", "PrecompileTools", "TOML"]
git-tree-sha1 = "ec5e99ad66c51b6fcd3dbab6655b17118a29a588"
uuid = "98e50ef6-434e-11e9-1051-2b60c6c9e899"
version = "2.3.0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "58927c485919bf17ea308d9d82156de1adf4b006"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.10.12"

[[deps.JuliaSyntax]]
git-tree-sha1 = "937da4713526b96ac9a178e2035019d3b78ead4a"
uuid = "70703baa-626e-46a2-a12c-08ffd08c73b4"
version = "0.4.10"

[[deps.JumpProcesses]]
deps = ["ArrayInterface", "DataStructures", "DiffEqBase", "DiffEqCallbacks", "DocStringExtensions", "FunctionWrappers", "Graphs", "LinearAlgebra", "PoissonRandom", "Random", "RecursiveArrayTools", "Reexport", "SciMLBase", "StaticArrays", "SymbolicIndexingInterface"]
git-tree-sha1 = "310fc6ddd46db8946f988b142fc2821f2f6456ff"
uuid = "ccbc3e58-028d-4f4c-8cd5-9ae44345cda5"
version = "9.23.1"

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
git-tree-sha1 = "c4d19f51afc7ba2afbe32031b8f2d21b11c9e26e"
uuid = "ba0b0d4f-ebba-5204-a429-3ac8c609bfb7"
version = "0.10.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "17b94ecafcfa45e8360a4fc9ca6b583b049e4e37"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.1.0+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

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

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static", "StaticArrayInterface"]
git-tree-sha1 = "a9eaadb366f5493a5654e843864c13d8b107548c"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.17"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

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
git-tree-sha1 = "cc3ad4faf30015a3e8094c9b5b7f19e85bdf2386"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.42.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d620582b1f0cbe2c72dd1d5bd195a9ce73370ab1"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.42.0+0"

[[deps.LineSearch]]
deps = ["ADTypes", "CommonSolve", "ConcreteStructs", "FastClosures", "LinearAlgebra", "MaybeInplace", "PrecompileTools", "SciMLBase", "SciMLJacobianOperators", "StaticArraysCore"]
git-tree-sha1 = "25454bc65c3eec4656cbe201c3d9336af49178c7"
uuid = "87fe0de2-c867-4266-b59a-2f0a94fc965b"
version = "0.1.8"
weakdeps = ["LineSearches"]

    [deps.LineSearch.extensions]
    LineSearchLineSearchesExt = "LineSearches"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Printf"]
git-tree-sha1 = "a666999510c794fe1d9dfa629ef33366f11103aa"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.6.1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LinearSolve]]
deps = ["ArrayInterface", "ChainRulesCore", "ConcreteStructs", "DocStringExtensions", "EnumX", "GPUArraysCore", "InteractiveUtils", "Krylov", "Libdl", "LinearAlgebra", "MKL_jll", "Markdown", "OpenBLAS_jll", "PrecompileTools", "Preferences", "RecursiveArrayTools", "Reexport", "SciMLBase", "SciMLLogging", "SciMLOperators", "Setfield", "StaticArraysCore"]
git-tree-sha1 = "15e73fc6ac5ad564842af06fb9f479199e7174fb"
uuid = "7ed4a6bd-45f5-4d41-b270-4a48e9bafcae"
version = "3.63.0"

    [deps.LinearSolve.extensions]
    LinearSolveAMDGPUExt = "AMDGPU"
    LinearSolveAlgebraicMultigridExt = "AlgebraicMultigrid"
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
    LinearSolveGinkgoExt = ["Ginkgo", "SparseArrays"]
    LinearSolveHYPREExt = "HYPRE"
    LinearSolveIterativeSolversExt = "IterativeSolvers"
    LinearSolveKernelAbstractionsExt = "KernelAbstractions"
    LinearSolveKrylovKitExt = "KrylovKit"
    LinearSolveMetalExt = "Metal"
    LinearSolveMooncakeExt = "Mooncake"
    LinearSolvePETScExt = ["PETSc", "SparseArrays"]
    LinearSolveParUExt = ["ParU_jll", "SparseArrays"]
    LinearSolvePardisoExt = ["Pardiso", "SparseArrays"]
    LinearSolveRecursiveFactorizationExt = "RecursiveFactorization"
    LinearSolveSparseArraysExt = "SparseArrays"
    LinearSolveSparspakExt = ["SparseArrays", "Sparspak"]

    [deps.LinearSolve.weakdeps]
    AMDGPU = "21141c5a-9bdb-4563-92ae-f87d6854732e"
    AlgebraicMultigrid = "2169fc97-5a83-5252-b627-83903c6c433c"
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
    Ginkgo = "4c8bd3c9-ead9-4b5e-a625-08f1338ba0ec"
    HYPRE = "b5ffcf37-a2bd-41ab-a3da-4bd9bc8ad771"
    IterativeSolvers = "42fd0dbc-a981-5370-80f2-aaf504508153"
    KernelAbstractions = "63c18a36-062a-441e-b654-da1e3ab1ce7c"
    KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
    LAPACK_jll = "51474c39-65e3-53ba-86ba-03b1b862ec14"
    Metal = "dde4c033-4e86-420c-a63e-0dd931031962"
    Mooncake = "da2b9cff-9c12-43a0-ae48-6db2b0edb7d6"
    PETSc = "ace2c81b-2b5f-4b1e-a30d-d662738edfe0"
    ParU_jll = "9e0b026c-e8ce-559c-a2c4-6a3d5c955bc9"
    Pardiso = "46dd5b70-b6fb-5a00-ae2d-e8fea33afaf2"
    RecursiveFactorization = "f2c3362d-daeb-58d1-803e-2bc74f2840b4"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    Sparspak = "e56a9233-b9d6-4f03-8d0f-1825330902ac"
    blis_jll = "6136c539-28a5-5bf0-87cc-b183200dce32"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.LoweredCodeUtils]]
deps = ["CodeTracking", "Compiler", "JuliaInterpreter"]
git-tree-sha1 = "5d4278f755440f70648d80cc6225f51e78e94094"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.5.1"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

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
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

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
deps = ["DataStructures", "LinearAlgebra", "MutableArithmetics", "StarAlgebras"]
git-tree-sha1 = "4838893d9b035c2f6967c0d533350e1755b58a70"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.5.19"
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
git-tree-sha1 = "7c25249fc13a070f5ba433c50e21e22bb33c6fb0"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.7.1"

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

[[deps.NearestNeighbors]]
deps = ["AbstractTrees", "Distances", "StaticArrays"]
git-tree-sha1 = "e2c3bba08dd6dedfe17a17889131b885b8c082f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.27"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonlinearSolve]]
deps = ["ADTypes", "ArrayInterface", "BracketingNonlinearSolve", "CommonSolve", "ConcreteStructs", "DifferentiationInterface", "FastClosures", "FiniteDiff", "ForwardDiff", "LineSearch", "LinearAlgebra", "LinearSolve", "NonlinearSolveBase", "NonlinearSolveFirstOrder", "NonlinearSolveQuasiNewton", "NonlinearSolveSpectralMethods", "PrecompileTools", "Preferences", "Reexport", "SciMLBase", "SciMLLogging", "SimpleNonlinearSolve", "StaticArraysCore", "SymbolicIndexingInterface"]
git-tree-sha1 = "d27bcf0cebf8786edcc2eaa4455c959e680334e7"
uuid = "8913a72c-1f9b-4ce2-8d82-65094dcecaec"
version = "4.16.0"

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
git-tree-sha1 = "eea7cbe389b168c77df7ff779fb7277019c685c8"
uuid = "5959db7a-ea39-4486-b5fe-2dd0bf03d60d"
version = "2.0.0"

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
git-tree-sha1 = "2ac022577e5eac7da040de17776d51bb770cd895"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

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
git-tree-sha1 = "58e5ed5e386e156bd93e86b305ebd21ac63d2d04"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.1+0"

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

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "844a829c8dc9fd0fe62eced22bc2d0dfd66a3f51"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.1.0"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "aea4eede5ab3ee188906d0cf3bbfa36eb543dccc"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.8"

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

[[deps.Polyester]]
deps = ["ArrayInterface", "BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "ManualMemory", "PolyesterWeave", "Static", "StaticArrayInterface", "StrideArraysCore", "ThreadingUtilities"]
git-tree-sha1 = "16bbc30b5ebea91e9ce1671adc03de2832cff552"
uuid = "f517fe37-dbe3-4b94-8317-1923a5111588"
version = "0.7.19"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "645bed98cd47f72f67316fd42fc47dee771aefcd"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.2.2"

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
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "25cdd1d20cd005b52fc12cb6be3f75faaf59bb9b"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.7"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "5e8e8b0ab68215d7a2b14b9921a946fee794749e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.3"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

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
git-tree-sha1 = "76a9c4d3f7b256dd34074e714ffe0ced8d3d21de"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "3.50.0"

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

[[deps.Revise]]
deps = ["CodeTracking", "FileWatching", "InteractiveUtils", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Preferences", "REPL", "UUIDs"]
git-tree-sha1 = "5f4f629c085b87e71125eec6773f5f872c74a47a"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.14.2"
weakdeps = ["Distributed"]

    [deps.Revise.extensions]
    DistributedExt = "Distributed"

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

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cfcdc949c4660544ab0fdeed169561cb22f835f4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.18"

[[deps.SCCNonlinearSolve]]
deps = ["CommonSolve", "PrecompileTools", "Reexport", "SciMLBase", "SymbolicIndexingInterface"]
git-tree-sha1 = "521a298294a6fd11068df0148b4d57bb31550453"
uuid = "9dfe8606-65a1-4bb3-9748-cb89d1561431"
version = "1.13.0"
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

[[deps.SciMLBase]]
deps = ["ADTypes", "Accessors", "Adapt", "ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Moshi", "PreallocationTools", "PrecompileTools", "Preferences", "Printf", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLLogging", "SciMLOperators", "SciMLPublic", "SciMLStructures", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface"]
git-tree-sha1 = "bdbc7c751d36fdb104d437361950d884d489af85"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "2.153.1"

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
git-tree-sha1 = "7156a5b51cba1bea33a82a036939ead4131f92bc"
uuid = "19f34311-ddf3-4b8b-af20-060888a46c0e"
version = "0.1.13"

[[deps.SciMLLogging]]
deps = ["Logging", "LoggingExtras", "Preferences"]
git-tree-sha1 = "0161be062570af4042cf6f69e3d5d0b0555b6927"
uuid = "a6db7da4-7206-11f0-1eab-35f2a5dbe1d1"
version = "1.9.1"

    [deps.SciMLLogging.extensions]
    SciMLLoggingTracyExt = "Tracy"

    [deps.SciMLLogging.weakdeps]
    Tracy = "e689c965-62c8-4b79-b2c5-8359227902fd"

[[deps.SciMLOperators]]
deps = ["Accessors", "ArrayInterface", "DocStringExtensions", "LinearAlgebra"]
git-tree-sha1 = "0e34162268883db01c04f988895a80d0659071bb"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "1.17.0"
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
git-tree-sha1 = "744c3f0fb186ad28376199c1e72ca39d0c614b5d"
uuid = "727e6d20-b764-4bd8-a329-72de5adea6c7"
version = "2.11.0"

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

[[deps.SparseMatrixColorings]]
deps = ["ADTypes", "DocStringExtensions", "LinearAlgebra", "PrecompileTools", "Random", "SparseArrays"]
git-tree-sha1 = "f63d76c7b7c329cf11badd564fd8ba877b09c3fe"
uuid = "0a514795-09f3-496d-8182-132a7b665d35"
version = "0.4.27"

    [deps.SparseMatrixColorings.extensions]
    SparseMatrixColoringsCUDAExt = ["CUDA", "cuSPARSE"]
    SparseMatrixColoringsCliqueTreesExt = "CliqueTrees"
    SparseMatrixColoringsColorsExt = "Colors"
    SparseMatrixColoringsGPUArraysExt = "GPUArrays"
    SparseMatrixColoringsJuMPExt = ["JuMP", "MathOptInterface"]

    [deps.SparseMatrixColorings.weakdeps]
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CliqueTrees = "60701a23-6482-424a-84db-faee86b9b1f8"
    Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
    GPUArrays = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
    JuMP = "4076af6c-e467-56ae-b986-b466b2749572"
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
    cuSPARSE = "b26da814-b3bc-49ef-b0ee-c816305aa060"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StarAlgebras]]
deps = ["LinearAlgebra", "MutableArithmetics", "SparseArrays"]
git-tree-sha1 = "235b1f9d287bbf34083b3d0829343a7942c0ad1c"
uuid = "0c0c59c1-dc5f-42e9-9a8b-b5dc384a6cd1"
version = "0.3.0"

[[deps.Static]]
deps = ["CommonWorldInvalidations", "IfElse", "PrecompileTools", "SciMLPublic"]
git-tree-sha1 = "49440414711eddc7227724ae6e570c7d5559a086"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "1.3.1"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "SciMLPublic", "Static"]
git-tree-sha1 = "aa1ea41b3d45ac449d10477f65e2b40e3197a0d2"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.9.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "246a8bb2e6667f832eea063c3a56aef96429a3db"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.18"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

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

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "86f5831495301b2a1387476cb30f86af7ab99194"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.8.0"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
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

    [deps.SymbolicIndexingInterface.extensions]
    SymbolicIndexingInterfacePrettyTablesExt = "PrettyTables"

    [deps.SymbolicIndexingInterface.weakdeps]
    PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "d673e0aca9e46a2f63720201f55cc7b3e7169b16"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "2.0.0"

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

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "ea3e54c2bdde39062abf5a9758a23735558705e1"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.4.0"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "248a7031b3da79a127f14e5dc5f417e26f9f6db7"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.1.0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

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
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

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
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

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
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

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
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

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
git-tree-sha1 = "429722587208f02b1cecbddcd20133df2f1ed796"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.47.0+0"

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
git-tree-sha1 = "850b06095ee71f0135d644ffd8a52850699581ed"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.3+0"

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
git-tree-sha1 = "e51150d5ab85cee6fc36726850f0e627ad2e4aba"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.58+0"

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
# ╟─f5023c00-491f-46c7-8a60-39baee9534c5
# ╠═92f1185b-ddad-40a8-bd02-cc09ca969519
# ╠═85d74ede-e1ac-4ea7-a514-98899159ef45
# ╠═3c0e40f3-01bc-4503-9577-4f51df700191
# ╠═40997d48-71fd-4387-b503-7b35aad78e54
# ╟─1706c1f9-717a-43af-beea-de483934a3a8
# ╠═e755c941-e82c-4b4a-9b07-8b12506de465
# ╠═ff12d5a8-2781-4abf-b10c-62e3863dde7a
# ╠═be304ff8-b720-4eeb-99f1-5fa6c9a41c2a
# ╟─4f44aab8-916b-11f1-1bd9-9bec891531eb
# ╟─e4950454-6422-4315-9476-1e98a13f321a
# ╟─ca4aa7f8-d9dc-47c0-93d0-36491f633fd3
# ╠═bf420b21-d616-4e9c-8c7d-19243ec09550
# ╠═983d8d90-a177-4d58-af6e-b1821e7ac46f
# ╠═714300ab-09c1-4c6d-98da-3bef18957a7d
# ╠═9f9b153f-7cd7-438c-a0a7-d02ab0210aec
# ╠═08d67be8-9471-44dc-8ca7-339a2ea5728d
# ╠═2e6ddf43-dbd2-4189-a4a5-ac4931ac9b6c
# ╠═67ddea74-84a1-49d4-ba3c-9c8fa611ceb6
# ╠═30d31108-5adc-4433-b4ce-dea1fdae8329
# ╠═0b5a79d8-0ef2-49b3-82ac-c6d5d9afa701
# ╠═b9e25013-e140-4542-987b-e6ddb604aed5
# ╠═0ede36ed-f4de-4b78-a55e-78fd279688ba
# ╠═03350557-e499-4ff1-8e7a-e93f3f992160
# ╠═8ffc05e1-e92b-44e7-9b70-ba1a8ebcfa49
# ╠═94919636-7f55-4f2e-b9bc-1e19bb6065cc
# ╟─e6e97ce5-7be1-48c8-a6b4-8834be453e87
# ╠═c3ad2eda-0770-4bff-94ef-3d67a6dc2716
# ╠═357aac49-2645-41d4-b71b-99b26e19f050
# ╠═df70d886-1cdc-400b-baf4-5a8e0b782e76
# ╟─7d3f91be-e4c1-4721-9a53-41d8a3fe15d6
# ╠═9a369bd3-e629-4bdb-9d0c-45fa8cc7722e
# ╠═29337edf-bf3b-4f7e-b58f-e0b604d9b377
# ╠═68370136-628b-4221-8c94-d2ab29b2502c
# ╟─b004c772-2b2a-4845-a6a6-ac78155ed0e6
# ╠═5ff8098e-42a0-4f26-8b90-336b7ffc1110
# ╠═a5548d15-f6fe-4e76-8d3e-1a4344d75bd0
# ╠═27157081-7a3c-4aa4-bbd2-c783283dd11d
# ╠═c6fe9b00-6bcc-4be0-a34c-b883da867028
# ╠═af271f05-a813-4d81-a83c-55ec1902a542
# ╟─8dc820b8-ca06-49f8-bf59-3de899f19c7d
# ╟─828a8999-9f52-4e60-b7a4-f9a00d53776f
# ╟─e8c5f1cc-46c7-46a0-9aed-fee9ba6d71c4
# ╠═ac4fc18d-1a18-44c6-9417-397d9901c7a2
# ╠═7dc086b5-45bb-4bdb-b6a3-78b507419d9d
# ╠═3af2d2af-289e-4314-87b5-0e6050add6ff
# ╟─3b6e72c2-1907-4ca5-8b3b-718d6784d7f8
# ╟─fb566e4a-9b99-475c-b2cf-0f8345e3a319
# ╟─2a9a124e-4422-4352-9b8e-b97e2033b4c7
# ╟─4ccb0a7c-18cb-43ab-a395-0d0cd047a720
# ╠═9ad29232-660e-4b4a-9539-458fcf26afb0
# ╟─904456f6-c25b-47b7-8a44-c3d11bae309f
# ╠═72f12e0a-f0fa-4820-a0dc-7e02a3aad68d
# ╟─8926d4c6-6400-4a2d-bccc-864fd4ffa379
# ╟─20680a91-5093-449a-a15a-2bf2bbd940f1
# ╟─4952a6b7-8efd-4e00-ace2-421e84952117
# ╠═b3a8e288-b396-4846-842c-e5348f06125d
# ╠═83c7674e-e663-4fad-a074-947449fb0272
# ╟─07fc17b6-94fc-4b95-b8dc-8c69eb9554b9
# ╟─f8084726-27c7-45bd-b376-ce0f61879130
# ╟─d4c76776-194c-439b-873f-1091f1bdceba
# ╟─7844f7e5-c5c2-4150-b688-890aec27f59f
# ╠═46bfb780-21f4-4522-bd84-6aef8808b25e
# ╠═49bce273-ec27-432b-8c30-2e8ae4f84467
# ╠═902410c1-edb7-4e23-a460-c969169cced5
# ╠═f6d8b3d3-4c8a-4a66-8a28-40dff51207e4
# ╠═3c1b1d7a-af93-429b-8077-9cb91531e2c7
# ╠═c37ed861-46fa-4c87-b8ee-e1c57523c824
# ╟─3f98f4aa-2dab-4b8b-9993-970f1ac0e814
# ╟─191e5a51-1261-4b0b-8bf0-3d9be58dae1b
# ╟─692bf75f-8dc6-4949-8928-788ec55eb8d6
# ╟─ad5c8b41-82c8-4034-9e24-726597aca74f
# ╟─1687efa2-b296-47cc-8b7d-4036ec34dd04
# ╟─497590c0-93d4-4921-8603-0fba2d40aa1f
# ╟─de199fc0-62d5-4bb5-9272-41e024197544
# ╟─41f8477b-77a3-49d6-bb5f-6b474959c0cf
# ╟─be1cb9c5-5c01-43e2-bb93-e67b9f0937d4
# ╠═60e89f62-b975-4f0b-bed4-2d98e7aa0c08
# ╠═0566f850-5083-4e40-9ef4-263f6a103c64
# ╠═fcd9a255-e395-4bb3-bf27-c2f335566d26
# ╠═e4d6e820-5bc0-4212-b6ad-22f1c74b4ab7
# ╠═fa8382c1-4f18-45dd-b432-4066c2d9ea23
# ╟─63256b4d-6487-4d23-895c-8c7143e5bf2b
# ╠═6eb6d99b-65e2-47aa-9095-480452bd2954
# ╟─a7143b46-8b83-4f16-9b09-9aa75b87f0b7
# ╟─cccc18c5-3ad0-4506-bac1-b52fc02aa332
# ╟─f2aebf87-7197-4c39-bc68-90e6b0cb64ac
# ╟─05ee4bd3-1f29-410c-bfda-0e768f8a47fe
# ╠═7f353f6c-a053-4919-8bdf-5c4cb6ed0fc1
# ╠═9387f387-3901-4de3-ba62-85d29dbf94fb
# ╠═41409955-ed92-4bb3-ae31-875d4699a54f
# ╠═787f5fe0-fc33-4c64-b722-6cc20292968b
# ╟─a4025273-fbe3-4f72-9f1e-dbd658483a5b
# ╟─3291abb6-ca4d-47b7-9b84-b2c9a422a4aa
# ╟─56bcb6a3-a53f-427c-a365-f4f528c0f754
# ╟─b5b4e4a1-43f1-41e5-b442-b57a61ec763b
# ╟─731eae26-849a-4119-ab8a-5b4e8131b5c0
# ╠═850c2a32-d4c5-4dc3-88e5-4fcb31cedb18
# ╠═2743c754-8fd9-4d35-b82e-f577399dca16
# ╟─f877762e-7ab5-4b3a-9065-1dc728d54b4d
# ╟─bb8013c2-8299-4a6b-acc8-0b8d9a5d98b4
# ╟─d16741e3-5138-4c71-8121-b786556a89e2
# ╟─4001ddf9-2cd7-4704-bfe3-ab24b6c12d62
# ╟─9c7a047a-e940-4599-9694-04d133ed2749
# ╟─c6ea2e44-3dc0-456e-bca9-b947468810c3
# ╟─6c6b75d9-d4e6-4a81-bcf7-58f1dcd8595e
# ╟─7a82894f-a68c-4b21-9434-0e1a463729c8
# ╠═3fe68f2c-4db5-403a-ba68-b5ac501305f3
# ╟─c2ecf3a2-4cc2-4d1c-a729-6f412de22b47
# ╠═db6d781f-8db9-4e95-b965-2aa748fe7d97
# ╠═f19c2d96-3fc4-4011-86f3-044f02b26c43
# ╠═e208347b-5aef-4181-82c0-ebc5b82ec134
# ╟─5243b1f8-a698-45a8-a650-b5fb5039b20f
# ╟─7ea4c7fb-f496-426f-b3ec-fad0302e5f88
# ╟─216a4a97-26ab-4e9f-aa43-55f9e8c59aa3
# ╠═08158174-6453-4e72-ad0d-4110d7275379
# ╠═21aa8bc8-8dcd-4fad-8cc8-e3a49c020577
# ╠═a70be775-8c77-46a2-835d-0eaa378bb4ae
# ╠═142086b2-4d1d-486e-b1ec-9df4810346ef
# ╠═9c28e7f0-b1be-436a-8293-6b139d203132
# ╠═db8fbe59-5aa7-4dba-918b-189084c93969
# ╟─84926400-3c00-4024-b4bd-d252d5a911cb
# ╠═f6af7beb-aa6e-4318-b975-05191f1afa1c
# ╠═75bbed03-a241-4879-923e-f0eb167361ef
# ╟─43465bde-3989-4833-93a2-669ad7f52612
# ╟─5ee56023-c336-4bc6-a015-419ac7805ece
# ╠═9bc995d7-44a1-44e8-a5d9-4534e5f92914
# ╠═568c355c-aca0-4c24-8951-f7c10ae9abee
# ╟─a5d7e48c-af63-46f5-bdae-943b3b3f7b71
# ╟─f98541fe-ae35-4f95-af2d-e4930cef2138
# ╟─0af62d26-4453-43ab-824b-d85ca7c43b30
# ╟─9a7dd642-8000-45d0-81af-ba168ed0860f
# ╠═a3ac46d3-f429-4b6a-9404-39d40d846186
# ╠═d6a6c468-a9e6-4ded-a6ac-8847c3f799c5
# ╟─8613c163-850e-4af0-8aa8-20f2bd7d715d
# ╟─94358d87-d685-435a-98f0-734b0d486b73
# ╠═c93af4de-9aba-4360-973f-4de6215aa56e
# ╟─12f6540f-4de3-40bb-aff2-524705a6f680
# ╠═da593e8b-ff28-4af7-af39-03d20ffe12bd
# ╟─40b92aab-5ecf-482e-9d31-7b26cd84500f
# ╟─79a969ca-5965-4820-aa3c-e2b51dd38c61
# ╟─1b47c432-dfb7-45ad-bffc-f548a9c4bf4b
# ╠═c8303901-9e32-4ba0-85a9-a230322d6830
# ╠═3839332c-384d-44b1-b48c-349f08bbc9a9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
