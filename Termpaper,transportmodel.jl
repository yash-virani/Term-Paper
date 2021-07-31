using JuMP
using Clp
using DataFrames
using Plots, StatsPlots

### setup of sets and parameters ###
# Timesteps as array of integers
T = [1, 2]

# Plant/Zone set as strings
P = ["Coal","Gas","Wind"]
NONDISP = ["Wind"]
Z = ["z1", "z2"]

mc = Dict(
    "Coal" => 30,
    "Gas" => 80,
    "Wind" => 0)

g_max = Dict(
    "Coal" =>  150,
    "Gas" => 40,
    "Wind" => 100)

# Zone to Plant
z2p = Dict("z1" => ["Coal"], "z2" => ["Gas", "Wind"])
# maps the inverse of n2p as a dict: power plant to node
p2z = Dict("Coal" => "z1", "Gas" => "z2", "Wind" => "z2")
p2z = Dict(p => k for (k,v) in z2p for p in v)

# Nodes indexed by strings, demand timeseries as array > Dict
demand = Dict("z1" => [120, 120], "z2" => [40, 40])
availability = Dict("Wind" => [1, 0.8])
# net transfer capacity (NTC) between Zones
# NTC is symmetric in both directions
ntc = Dict(("z1","z2") => 50, ("z2","z1") => 50)

# % Model definition

# Initialize the JuMP model 
m = Model(Clp.Optimizer);
# Define variables G, EX
@variables m begin 
    G[P, T] >= 0 # posivie variables for generation and exchange
    EX[Z, Z, T] >= 0
end

# Objective Function
@objective(m, Min,
    sum(mc[p] * G[p,t] for p in P, t in T));

# Capacity Constraints
@constraint(m, Gmax[p=P, t=T],
    G[p,t] <= (p in NONDISP ? availability[p][t]*g_max[p] : g_max[p]));

# Exchange is bound by NTCs
@constraint(m, Exchange[z=Z, zz=Z, t=T], 
    EX[z,zz,t] <= ((z,zz) in keys(ntc) ? ntc[z,zz] : 0));

# Energy Balance.
@constraint(m, ElectricityBalance[z=Z, t=T], # the balance is created for each zone
    sum(G[p,t] for p in z2p[z]) # z2p holds which plants are at zone z
    + sum(EX[zz,z,t] - EX[z,zz,t] for zz in Z) # Commercial Exchange as sum of Import - Exports
    == demand[z][t]); # Has su supply all demand. 

optimize!(m); # Run the optimization

# Result analysis
result_G = value.(G)
g = DataFrame(
        (variable="generation_"*p, zone=p2z[p], t=t, value=result_G[p,t])
    for p in P, t in T)
result_EX = value.(EX);

exchange = DataFrame(
    (variable="exchange", zone=z, t=t,
     value = sum(result_EX[zz,z,t] - result_EX[z,zz,t] for zz in Z))
    for z in Z, t in T);

energybalance = vcat(g, exchange);

marginal = dual.(ElectricityBalance);


### Visualize the Results 

colors = [:brown :red :blue :purple]

df_zone_1 = filter(x-> x.zone == "z1", energybalance)
df_zone_2 = filter(x-> x.zone == "z2", energybalance)

x1 = df_zone_1[:,:t]
y1 = df_zone_1[:,:value]
g1 = df_zone_1[:,:variable]

z1 = groupedbar(x1,y1, group=g1, color=colors,
                bar_position=:stack, legend=:bottomleft,
                yaxis="MW", title="Zone 1", grid=false)

marginal_z1 = marginal["z1",:].data
z1_twin = twinx(z1)
scatter!(z1_twin, marginal_z1, color=:black, legend=false, ylim=(-10,50), grid=false)

x2 = df_zone_2[:,:t]
y2 = df_zone_2[:,:value]
g2 = df_zone_2[:,:variable]

z2 = groupedbar(x2,y2, group=g2, color=colors, bar_position=:stack, legend=true,
                title="Zone 2", grid=false)

marginal_z2 = marginal["z2",:].data
z2_twin = twinx(z2)
scatter!(z2_twin, marginal_z2,
    color=:black,
    legend=false,
    ylim=(-10,50),
    grid=false,
    yaxis="Price")

plot(z1, z2, grid=(2,1))
savefig("result_transport.pdf")