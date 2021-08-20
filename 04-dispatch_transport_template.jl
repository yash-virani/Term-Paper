using JuMP
# using Clp
# import Pkg
# Pkg.add("Gurobi")
# Pkg.build("Gurobi")
using Gurobi
using Plots, StatsPlots
using DataFrames, CSV 

using Statistics
#  Preprocessing
### data load ###
datapath = joinpath(@__DIR__, "data")
# datapath = "data"

timeseries_path = joinpath(datapath, "time_series")
static_path = joinpath(datapath, "static")

plants =  CSV.read(joinpath(static_path,"plants.csv"), DataFrame)
ntc_data =  CSV.read(joinpath(static_path,"ntc.csv"), DataFrame)


timeseries = Dict(splitext(files)[1] => CSV.read(joinpath(timeseries_path, files), DataFrame)
    for files in readdir(timeseries_path))

######################################################################
### create sets based on input data ###

### sets ###
T = 1:size(timeseries["demand_2015"], 1) |> collect
# T = 1:4000 |> collect
TECH = unique(plants[:,:id] |> Vector)
DISP = unique(plants[plants[:,:disp] .== 1, :id])
NONDISP = unique(plants[plants[:,:disp] .== 0, :id])
S = plants[plants[:,:storage_capacity] .!= 0, :id]
Z = unique(ntc_data[:,:from_country])
CU = plants[plants[:,:tech] .== "CU", :id] |> Vector
CU_NEG = plants[plants[:,:tech] .== "CU_neg", :id] |> Vector

### parameters ###

### utility functions and mappings ###
zipcols(df::DataFrame, x::Symbol) = df[:,x] |> Vector
zipcols(df::DataFrame, x::Vector) = zip(eachcol(df[:,x])...) |> collect

function dictzip(df::DataFrame, x::Pair)
    dictkeys = zipcols(df, x[1])
    dictvalues = zipcols(df, x[2])
    return zip(dictkeys, dictvalues) |> collect |> Dict
end

coldict(df::DataFrame) = Dict(string(name) => Vector(vec) for (name,vec) in pairs(eachcol(df)))

### Static Parameters ###
mc = dictzip(plants, :id => :mc_el)
g_max = dictzip(plants, :id => :g_max)
d_max = dictzip(plants, :id => :d_max)
eta = dictzip(plants, :id => :eta)
storage_capacity = dictzip(plants, :id => :storage_capacity)

map_country2id = Dict()
for z in Z
    map_country2id[z] = (plants[plants.country .== z, "id"])
end

map_id2country = Dict()
for plant_id in plants[:,:id] 
    map_id2country[plant_id] = unique((plants[plants.id .== plant_id, "country"]))[1]
end

map_id2tech  = Dict()
for plant_id in plants[:,:id] 
    map_id2tech[plant_id] = unique((plants[plants.id .== plant_id, "tech"]))[1]
end
# Hyrdo Inflows
psp_inflow = coldict(timeseries["hydro_psp_inflow_countries"])
PSP_list = plants[plants[:,:tech] .== "hydro_psp", :id]
res_inflow = coldict(timeseries["hydro_res_inflow_countries"])
RES_list = plants[plants[:,:tech] .== "hydro_res", :id]

### time series ###
coldict(df::DataFrame) = Dict(string(name) => Vector(vec) for (name,vec) in pairs(eachcol(df)))

elec_demand = coldict(timeseries["demand_2015"])

# Creating feed in for non dispatchables, !!!!!!!!!! Hydropower is still missing due to lack of Timeseries data!!!!!!!!!!!!!!!
feed_in_wind_off = Dict()
feed_in_wind_on = Dict()
feed_in_pv = Dict()
feed_in_ror = Dict()
feed_in = Dict()

for z in Z
    if ("wind onshore" in  plants[(plants.country .== z),"tech"]) & (z in names(timeseries["wind_on_2015"]))
        feed_in_wind_off[z] = timeseries["wind_on_2015"][!,z] .* (plants[(plants.country .== z) .& (plants.tech .== "wind onshore"), "g_max"])
    end

    if ("wind offshore" in  plants[(plants.country .== z),"tech"]) & (z in names(timeseries["wind_off_2015"]))
        feed_in_wind_on[z] = timeseries["wind_off_2015"][!,z] .* (plants[(plants.country .== z) .& (plants.tech .== "wind offshore"), "g_max"])
    end

    if ("solar" in  plants[(plants.country .== z),"tech"]) & (z in names(timeseries["pv_2015"]))
        feed_in_pv[z] = timeseries["pv_2015"][!,z] .* (plants[(plants.country .== z) .& (plants.tech .== "solar"), "g_max"])
    end
    if ("ror" in  plants[(plants.country .== z),"tech"])
        feed_in_ror[z] = [0.65 for i=1:8760] .* (plants[(plants.country .== z) .& (plants.tech .== "ror"), "g_max"])
        # feed_in_ror[z] = [0.65 for i=1:size(T)[1]] .* (plants[(plants.country .== z) .& (plants.tech .== "ror"), "g_max"])
    end
    feed_in[z] = [0 for i=1:8760]
    # feed_in[z] = [0 for i=1:size(T)[1]]
end

# Summing the different feedins 
for z in Z
    if haskey(feed_in_wind_off, z) feed_in[z] = feed_in[z] .+ feed_in_wind_off[z] end
    if haskey(feed_in_wind_on, z) feed_in[z] = feed_in[z] .+ feed_in_wind_on[z] end
    if haskey(feed_in_pv, z) feed_in[z] = feed_in[z] .+ feed_in_pv[z] end
    if haskey(feed_in_ror, z) feed_in[z] = feed_in[z] .+ feed_in_ror[z] end

end

feed_in_by_z_nondisp= Dict()

for x in keys(feed_in_wind_off)
    feed_in_by_z_nondisp[x,"wind offshore"] = feed_in_wind_off[x]
end

for x in keys(feed_in_wind_on)
    feed_in_by_z_nondisp[x,"wind onshore"] = feed_in_wind_on[x]
end

for x in keys(feed_in_pv)
    feed_in_by_z_nondisp[x,"solar"] = feed_in_pv[x]
end

### sucessor: give the next index of an array, return 1 for the last element###
successor(arr, x) = x == length(arr) ? 1 : x + 1

### ntc ###
ntc = dictzip(ntc_data, [:from_country, :to_country] => :ntc)

# actual model creation
###############################################################################
m = Model(Gurobi.Optimizer)
set_optimizer_attribute(m, "TimeLimit", 600)


@variables m begin
    g_max[disp] >= G[disp=DISP, T] >= 0
    d_max[s] >= D_stor[s=S,T] >= 0
    storage_capacity[s] >= L_stor[s=S,T] >= 0
    EX[z=Z,zz=Z,T] >= 0
    CU[Z,T] >= 0
    BALANCE_P[Z,T] >= 0

end

@objective(m, Min,
    sum(mc[disp] * G[disp,t] for disp in DISP, t in T)
    + sum((CU[z,t] + BALANCE_P[z,t])*1000 for z in Z, t in T)
    );

@constraint(m, ElectricityBalance[z=Z, t=T],
    sum(G[disp,t] for disp in intersect(map_country2id[z],DISP))
    + feed_in[z][t]
    - sum(D_stor[s,t] for s in intersect(map_country2id[z],S))
    - CU[z,t]
    + BALANCE_P[z,t]
    + sum(EX[zz,z,t] - EX[z,zz,t] for zz in Z)
    ==
    elec_demand[z][t])

# NTC upper-bound for exchange between zones (lower bound not need as ntc is declared as greater than 0) 
@constraint(m, Exchange[z=Z, zz=Z, t=T], 
    EX[z,zz,t] <= ((z,zz) in keys(ntc) ? ntc[z,zz] : 0));


@constraint(m, StorageLevel[s=S, t=T],
    L_stor[s, successor(T,t)]
    ==
    L_stor[s, t]
    + sqrt(eta[s])*D_stor[s,t]
    - (1/sqrt(eta[s])) * G[s,t]
    + (s in PSP_list ? psp_inflow[map_id2country[s]][t] : 0)
    + (s in RES_list  ? res_inflow[map_id2country[s]][t] : 0)
    );

optimize!(m)


# post processing
########################################################
generic_names(len::Int) = [Symbol("x$i") for i in 1:len]

function DataFrame(x::JuMP.Containers.DenseAxisArray,
    names::AbstractVector{Symbol}=generic_names(length(x.axes));
    kwargs...)

    if length(names) != length(x.axes)
        throw(ArgumentError("Length of argument 'names' is $(length(names)) and
        does not fit the variable's dimension of $(length(x.axes))"))
    end

    push!(names, :value)

    iter = Iterators.product(x.axes...) |> collect |> vec
    df = DataFrame([(i..., value(x[i...])) for i in iter])
    rename!(df, names)
    return df
end

result_G = DataFrame(G, [:id, :hour])
insertcols!(result_G, 2, :zone => [map_id2country[id] for id in result_G[!,:id]])
insertcols!(result_G, 3, :technology => [map_id2tech[id] for id in result_G[!,:id]])

# Debugging



# Creating new result DF to collect exchanges between zones
result_EX_data = DataFrame(EX, [:from_Z, :to_Z, :hour] )
# Filtering out entries from and to same zone (as thos will by definition allways be 0)
result_EX_data = result_EX_data[result_EX_data.from_Z .!= result_EX_data.to_Z, :]
# dividing the exchange to only import data
res_imp = copy(result_EX_data) # fetching a copy to not alter the orginal DF
insertcols!(res_imp, 3, :technology => "import") #adding column with marking 'import'
res_imp = select!(res_imp, Not(:from_Z)) # selecting only import data
rename!(res_imp, :to_Z => :zone) # renaming the column to uniform 'zone' label

# dividing the exchange to only export data
res_exp = copy(result_EX_data) # fetching a copy to not alter the orginal DF
insertcols!(res_exp, 3, :technology => "export") #adding column with marking 'export'
res_exp = select!(res_exp, Not(:to_Z)) # selecting only export data
res_exp.value = -res_exp.value # As export have to show as negative, *(-1) is applied
rename!(res_exp, :from_Z => :zone) # renaming the column to uniform 'zone' label

# recombining the DFs import and export
result_EX = vcat(res_imp, res_exp)

##########################################################################
# TABLE 1

plants_mc = dictzip(plants, [:id, :country] => :mc_el)

plants_mc_id = dictzip(plants, :id => :mc_el )
#manipulate result_G
result_G_TEST = copy(result_G)
result_G_TEST = select!(result_G_TEST, Not(:technology))


#production (annual, summer, winter) per id, country
id_h_value = dictzip(result_G_TEST, [:id,:hour] => :value)
tech_country = dictzip(plants, :id => :country)
# Summer Timing (Astronomical Season): 21. June[171*24 = 4104] - 23. September[265*24 = 6360]
# Winter Timing (Astronomical Season): 1. Januar - 20. März[78*24 = 1872] ; 22. Dezember - 31. Dezember[10*24 = 240]
T_s = 4104:6360
T_w1 = 0:1872
T_w2 = 8520:8760
T_w = vcat(T_w1,T_w2)

sum_tech = DataFrame(
    (technology = tech,
    country = tech_country[tech],
    value_annual = sum( (tech,t) in keys(id_h_value) ? id_h_value[tech,t] : 0 for t in T),
    value_summer = sum( (tech,t) in keys(id_h_value) ? id_h_value[tech,t] : 0 for t in T_s),
    value_winter = sum( (tech,t) in keys(id_h_value) ? id_h_value[tech,t] : 0 for t in T_w)
    )
    for tech in TECH
)

avg_annual_z = dictzip(sum_tech, [:technology, :country] => :value_annual)
avg_summer_z = dictzip(sum_tech, [:technology, :country] => :value_summer)
avg_winter_z = dictzip(sum_tech, [:technology, :country] => :value_winter)


#total annual production in country z
total_production_z = DataFrame(
    (country = z,
    total_gen_annual = sum((id,z) in keys(avg_annual_z) ? avg_annual_z[id,z] : 0 for id in sum_tech.technology),
    total_gen_summer = sum((id,z) in keys(avg_summer_z) ? avg_summer_z[id,z] : 0 for id in sum_tech.technology),
    total_gen_winter = sum((id,z) in keys(avg_winter_z) ? avg_winter_z[id,z] : 0 for id in sum_tech.technology)
    )
    for z in Z
)

total_annual_z = dictzip(total_production_z, :country => :total_gen_annual)
total_summer_z = dictzip(total_production_z, :country => :total_gen_summer)
total_winter_z = dictzip(total_production_z, :country => :total_gen_winter)


# Durchschnittlichen Marginal cost * producution
# Durchschnittlichen Marginal cost sum(mc_cost*production_technology for Z)/ (production in country Z)
table1 = DataFrame(
    (zone = z,
    avg_price_annual= sum((id,z) in keys(avg_annual_z) ? avg_annual_z[id,z]*plants_mc[id,z] : 0 for id in sum_tech.technology)/total_annual_z[z],
    avg_price_summer= sum((id,z) in keys(avg_summer_z) ? avg_summer_z[id,z]*plants_mc[id,z] : 0 for id in sum_tech.technology)/total_summer_z[z],
    avg_price_winter= sum((id,z) in keys(avg_winter_z) ? avg_winter_z[id,z]*plants_mc[id,z] : 0 for id in sum_tech.technology)/total_winter_z[z],
    )
    for z in Z
)


# Variance for Germany (DE)
find_DE = result_G_TEST[findall(x -> occursin("DE",x), result_G_TEST[!,:zone]),:]
dict_DE = dictzip(find_DE, [:id, :hour] => :value)
for (id,t) in keys(dict_DE)
    dict_DE[id,t] *= plants_mc_id[id]
end


#Standardabweichung der täglichen Kosten in Deutschland 
#--> Eigentlich müssten die nicht prod. TECH rausgenommen weren (verfälscht)
sqrt(var(collect(values(dict_DE))))


##########################################################################
# MC Cost per Zone

id_per_t = dictzip(result_G_TEST, [:zone, :hour, :id] => :value)

findmax(collect(("UK",3064,id) in keys(id_per_t) && id_per_t["UK",3064,id] > 0 ? plants_mc_id[id] : 0 for id in result_G_TEST.id))

mc = Dict()
for z in Z
    for t in T
        test = collect(("UK",3064,id) in keys(id_per_t) && id_per_t["UK",3064,id] > 0 ? id_per_t["UK",3064,id] : 0 for id in result_G_TEST.id)
        mc[z,t] => (findmax(test))
    end
end

# TEST

# countries= DataFrame(
#     (country = z,
#     dataframe_z = (result_G_TEST[findall(x -> occursin(z,x), result_G_TEST[!,:zone]),:])
#     dict_z = dictzip()
#     )
#     for z in Z
# )


# result_G_summer = copy(result_G)
# result_G_summer = filter!(row -> 6360>row.hour>4103 , result_G_TEST)
# mc_annual_s = dictzip(result_G_TEST, [:id, :zone,:hour] => :value)

# T_s = 4104:6360

# T_w1 = 0:1872
# T_w2 = 8520:8760
# T_w = vcat(T_w1,T_w2)


# avg_price_annual_z = DataFrame()
# avg_price_annual_summer_z = DataFrame()
# avg_price_annual_winter_z = DataFrame()
# # ANNUAL
# for z in Z
#     avg_price_annual_z[!,z] = sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in result_G.id, t in T)/8760
# end

# # SUMMER
# for z in Z
#     avg_price_annual_summer_z[!,z] = 
#         sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in result_G.id, t in T_s)/8760
# end

# # WINTER
# for z in Z
#     avg_price_annual_winter_z[!,z] = 
#         sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in result_G.id, t in T_w)/8760
# end

# table1 = DataFrame(
#     (zone = z,
#     avg_price_annual= sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in TECH, t in T)/8760,
#     avg_price_summer= sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in TECH, t in T_s)/8760,
#     avg_price_winter= sum((id,z) in keys(plants_mc) && (id,z,t) in keys(mc_annual) ? plants_mc[id, z]*mc_annual[id,z,t] : 0 for id in TECH, t in T_w)/8760,
#     )
#     for z in Z
# )


# TABLE 2
ntc_time = dictzip(result_EX_data, [:from_Z, :to_Z, :hour] => :value)

table2 = DataFrame(
    (from_zone = z,
    to_zone = zz,
    ntc_max = ntc[z,zz],
    avg_sum_year = sum(ntc_time[z,zz,t] for t in T)/8760,
    utilisation_percent =  ntc[z,zz] != 0 ? (sum(ntc_time[z,zz,t] for t in T)/8760)/ntc[z,zz] : 0,
    utilisation_equal100 = ntc[z,zz] != 0 ? sum(ntc[z,zz] == ntc_time[z,zz,t] ? 1 : 0 for t in T)/8760 : 0,
    utilisation_bigger90 = ntc[z,zz] != 0 ? sum(ntc[z,zz]*0.9 <= ntc_time[z,zz,t] ? 1 : 0 for t in T)/8760 : 0,
    )
    for z in Z, zz in Z if z!=zz
)
CSV.write("/Users/lukasschirren/desktop/FromTo_ntc_utilisation.csv", table2)


# TABLE 3
#TEST
#export_FR = res_exp[findall(x -> occursin("FR", x), res_exp[!, :zone]), :]
#sum_FR = sum(export_FR.value)/8760
#export_FR_ntc_MAX = ntc_data()
#CSV.write("/Users/lukasschirren/desktop/inflow_AT.csv", export_DE)

# max NTC capacities (import/export)
fixed_ntc_in = DataFrame()
fixed_ntc_out = DataFrame()
for z in Z
    fixed_ntc_in[!,z] = (ntc_data[findall(x -> occursin(z,x), ntc_data[!,:to_country]),:]).ntc
    fixed_ntc_out[!,z] = (ntc_data[findall(x -> occursin(z,x), ntc_data[!,:from_country]),:]).ntc
end


# Inflow/Outflow over the whole year
import_z = DataFrame()
export_z = DataFrame()
for z in Z
    import_z[!,z] = (res_imp[findall(x -> occursin(z, x), res_imp[!, :zone]), :]).value
    export_z[!,z] = (res_exp[findall(x -> occursin(z, x), res_exp[!, :zone]), :]).value
end

table3 = DataFrame(
    (zone=z,
    max_import_ntc = sum(fixed_ntc_in[!,z]),
    avg_import = sum(import_z[!,z])/8760,
    max_export_ntc = sum(fixed_ntc_out[!,z]),
    avg_export = sum(export_z[!,z])/8760,
    utilisation_in = (sum(import_z[!,z])/8760)/sum(fixed_ntc_in[!,z]),
    utilisation_out = -(sum(export_z[!,z])/8760)/sum(fixed_ntc_out[!,z])
    )
    for z in Z
)
CSV.write("/Users/lukasschirren/desktop/Country_ntc_import_export.csv", table3)

##########################################################################

gen_by_tech = combine(groupby(result_G, [:zone, :technology, :hour]),
    :value => sum => :value)

result_feedin = DataFrame(
        (id=id,
        zone=map_id2country[id],
        technology=map_id2tech[id],
        hour=t,
        value=haskey(feed_in_by_z_nondisp,[map_id2country[id],map_id2tech[id]]) ?  feed_in_by_z_nondisp[map_id2country[id],map_id2tech[id]][t] : 0
        )
    for id in NONDISP, t in T)

gen_by_tech2 = combine(groupby(result_feedin, [:zone, :technology, :hour]),
    :value => sum => :value)

result_CU = DataFrame(CU, [:zone, :hour])
result_CU[!,:technology] .= "curtailment"

############################### BALANCE_P muss noch im weiteren aufgenommen werden ######################################
result_BALANCE_P = DataFrame(BALANCE_P, [:zone, :hour])
result_BALANCE_P[!,:technology] .= "Balancing_Power"
############################### BALANCE_P muss noch im weiteren aufgenommen werden ######################################

result_D = DataFrame(D_stor, [:id, :hour])
insertcols!(result_D, 2, :zone => [map_id2country[id] for id in result_D[!,:id]])
insertcols!(result_D, 3, :technology => [map_id2tech[id] for id in result_D[!,:id]])

demand_by_storage = combine(groupby(result_D, [:zone, :technology, :hour]),
    :value => sum => :value)
electricity_demand = DataFrame((zone=k, technology="demand", hour=x, value=v[x])
    for (k,v) in elec_demand for x in 1:length(v))

result_generation = vcat(gen_by_tech2, gen_by_tech)
result_demand = vcat(demand_by_storage, result_CU, electricity_demand)

#added color for exchange (same color as it is clearly visible that export is negative and import positive)
colordict = Dict(
    "solar" => :yellow, "wind offshore" => :lightblue, "wind onshore" => :lightblue, "hydro_psp" => :darkblue, "hydro_res" => :darkblue, "ror" => :blue, 
    "ccgt" => :lightgrey, "ccot" => :lightgrey,"generator" => :lightgrey, "steam" => :brown, "geothermal" => :pink, "demand" => :darkgrey,
    "curtailment" => :red, "import" => :purple, "export" => :purple) 

function plot_energybalance(df_gen::DataFrame, df_dem::DataFrame, df_ex::DataFrame, z::AbstractString)
    
    df_generation = filter(x-> x.zone == z, df_gen)
    df_demand = filter(x-> x.zone == z, df_dem)
    # added new DF for import and export
    df_exchange = filter(x-> x.zone == z, df_ex)

    table_gen = unstack(df_generation, :hour, :technology, :value)
    labels = names(table_gen[:,Not(:hour)]) |> permutedims
    colors = [colordict[tech] for tech in labels]
    data_gen = Array(table_gen[:,Not(:hour)])

    p = areaplot(data_gen, label=labels, color=colors, width=0,
        leg=:outertopright, title="Dispatch of zone $z")

    table_dem = unstack(df_demand, :hour, :technology, :value)
    labels = names(table_dem[:,Not(:hour)]) |> permutedims
    colors = [colordict[tech] for tech in labels]
    data_dem = - Array(table_dem[:,Not(:hour)])

    areaplot!(p, data_dem, label=labels, color=colors, width=0)
    hline!(p, [0], color=:black, label="", width=2)

    # added new table, labels and color specification  for plotting
    table_ex = unstack(df_exchange, :hour, :technology, :value, allowduplicates = true)
    labels = names(table_ex[:,Not(:hour)]) |> permutedims
    colors = [colordict[tech] for tech in labels]
    data_ex = Array(table_ex[:,Not(:hour)])

    areaplot!(p, data_ex, label=labels, color=colors, width=0)

    return p
end

result_generation[result_generation[:,:zone] .== "DE",:]


z1 = plot_energybalance(result_generation, result_demand, result_EX, "DE");
# savefig("results_z1.pdf")


df_exchange = filter(x-> x.zone == "z1", result_EX)

# Exploring the exchange results between zones
SUMexp_z1 = sum(result_EX_data[result_EX_data.from_Z .== "z1", :].value);
SUMexp_z2 = sum(result_EX_data[result_EX_data.from_Z .== "z2", :].value);
NETimp_z1 = round(SUMexp_z2 - SUMexp_z1, digits=2);
NETexp_z2 = NETimp_z1;
println("The Zone 1 has a net import of: $NETimp_z1 units, the zone 2 has a net export of $NETexp_z2 units")
println("The direction of commercial Exchange is therefor from zone 2 to zone 1")

averageEX = round((sum(result_EX_data.value)*2)/length(result_EX_data.value), digits=4);
println("The average load of the interconnector is: $averageEX")

fullUtlizationH = length(result_EX_data[result_EX_data.value .==0.5, :].value);
fullUtlizationPercent = round((fullUtlizationH/length(result_EX_data.value))*100,digits= 2);
println("The interconnector is fully utilized in $fullUtlizationH hours, which is $fullUtlizationPercent % of the time")
