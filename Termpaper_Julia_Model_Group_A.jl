using JuMP
using Gurobi
using Plots, StatsPlots
using DataFrames, CSV 
using Statistics

### utility functions ###
function dictzip(df::DataFrame, x::Pair)
    dictkeys = zipcols(df, x[1])
    dictvalues = zipcols(df, x[2])
    return zip(dictkeys, dictvalues) |> collect |> Dict
end

coldict(df::DataFrame) = Dict(string(name) => Vector(vec) for (name,vec) in pairs(eachcol(df)))

generic_names(len::Int) = [Symbol("x$i") for i in 1:len]

zipcols(df::DataFrame, x::Symbol) = df[:,x] |> Vector

zipcols(df::DataFrame, x::Vector) = zip(eachcol(df[:,x])...) |> collect

### sucessor: give the next index of an array, return 1 for the last element###
successor(arr, x) = x == length(arr) ? 1 : x + 1

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

function plot_ldc2(dem, feed, wind_off, wind_on, pv, ror)
    
    resid_load_data = DataFrame(load = dem,
                                feedin = feed,
                                feed_in_wind_off = wind_off,
                                feed_in_wind_on = wind_on,
                                feed_in_pv = pv,
                                feed_in_ror = ror
                                )
   

    df_copy = copy(resid_load_data)
    df_copy[!, :residual_load] = df_copy[!, :load] .- df_copy[!,:feedin]
    df_copy[!, :load_wo_wind] = df_copy[!, :load] .- df_copy[!,:feed_in_wind_on] .- df_copy[!,:feed_in_wind_off]
    df_copy[!, :load_wo_wind_pv] = df_copy[!, :load_wo_wind] .- df_copy[!,:feed_in_pv]
    
    sorted_load = sort(df_copy[!,:load], rev=true)
    sorted_residual = sort(df_copy[!,:residual_load], rev=true)
    sorted_load_wo_wind = sort(df_copy[!,:load_wo_wind], rev=true)
    sorted_load_wo_wind_pv = sort(df_copy[!,:load_wo_wind_pv], rev=true)


    total_load = sum(df_copy[!,:load])
    rl = df_copy[!,:residual_load]
    pos_rl = filter(x-> x >= 0, rl)
    share = sum(rl / total_load)
    share_pos = sum(pos_rl / total_load)
    
    p = plot(
        sorted_load,
        color=:black,
        width=3,
        label="Load",
        xlabel="Number of hours",
        ylabel="MW",
        fillrange=sorted_load_wo_wind,
        fillcolor=:lightblue,
        fillalpha=0.2,
        leg=:outerbottom,
        title="Share of nondispatchable generation of total load: $(round(Int, 100*(1 - share))) % \n load covered bei nondispatchable generation: $(round(Int, 100*(1 - share_pos))) %",
        titlefontsize=6
    )

    plot!(p,
        sorted_load_wo_wind,
        color=:black,
        width=1,
        label="Load with wind infeed",
        xlabel="Number of hours",
        ylabel="MW",
        fillrange=sorted_residual,
        fillcolor=:yellow,
        fillalpha=0.2)

    plot!(p,
        sorted_load_wo_wind_pv,
        color=:black,
        width=1,
        label="Load with solar and wind infeed",
        xlabel="Number of hours",
        ylabel="MW",
        fillrange=sorted_residual,
        fillcolor=:darkblue,
        fillalpha=0.2)

    plot!(p,
        sorted_residual,
        color=:black,
        width=2,
        label="Residual load with solar, wind and ror infeed",
        xlabel="Number of hours",
        ylabel="MW",
        fill=0,
        fillcolor=:red,
        fillalpha=0.2)
    
    hline!([0], width=2, color = :black, label="")
    
    return p
end

function create_line_constrain_map(result_EX_data, ntc)
    line_constrains_map = copy(result_EX_data)
    line_constrains_map = line_constrains_map[line_constrains_map.from_Z .!= line_constrains_map.to_Z, :]
    ntc_max_collection = []
    constrained_coll = []
    counrty_set_coll = []
    for row in eachrow(line_constrains_map)
        append!(ntc_max_collection,ntc[(row.from_Z, row.to_Z)])
        append!(constrained_coll, ntc[(row.from_Z, row.to_Z)] == row.value ? 1 : 0)
        append!(counrty_set_coll, [sort([row.from_Z, row.to_Z])])
    end
    insertcols!(line_constrains_map, 4, :ntc_max => ntc_max_collection )
    insertcols!(line_constrains_map, 5, :constrained => constrained_coll )
    insertcols!(line_constrains_map, 6, :country_sets => counrty_set_coll ) 
    line_constrains_map = line_constrains_map[(line_constrains_map.ntc_max .!= 0), :]

    line_constrains_map = combine(groupby(copy(line_constrains_map), [:country_sets, :hour]), :constrained .=> [maximum])
    line_constrains_map = combine(groupby(copy(line_constrains_map), [:country_sets]), :constrained_maximum .=> [mean])
    insertcols!(line_constrains_map, 1, :from => [x[1] for x in line_constrains_map.country_sets])
    insertcols!(line_constrains_map, 2, :to => [x[2] for x in line_constrains_map.country_sets])
    select!(line_constrains_map, Not(:country_sets))
    return line_constrains_map
end

function run_model(enable_2030::Bool, scen_name::String, ntc_factor,extra_storage_factor, Time_steps::Int)

    results_path  = joinpath(@__DIR__, "Results")
    if !isdir(results_path)
        mkdir(results_path)
    end


    #  Preprocessing
    ### data load ###
    datapath = joinpath(@__DIR__, "data")

    timeseries_path = joinpath(datapath, "time_series")
    static_path = joinpath(datapath, "static")

    plants =  CSV.read(joinpath(static_path,"plants.csv"), DataFrame)

    ntc_data =  CSV.read(joinpath(static_path,"ntc.csv"), DataFrame)
    ntc_data[!,"ntc"] .= ntc_data[!,"ntc"] .* ntc_factor

    timeseries = Dict(splitext(files)[1] => CSV.read(joinpath(timeseries_path, files), DataFrame)
        for files in readdir(timeseries_path))

    ############### switch this on for 2030 renewable capacity values #######################
    if enable_2030
        renewables_2030 = CSV.read(joinpath(static_path,"renewables_2030.csv"), DataFrame)
        for x in eachrow(plants), i in eachrow(renewables_2030)
            if (x.country == i.country) & (x.tech == i.tech) 
                x["g_max"] =  i[:g_max_2030]
            end
        end
    end

    ##################### Fuel reduction according to scenario ######################
    if scen_name == "Mixed_FR_PL_exept"
        fuel_reduction = ["uran", "lignite", "hard coal"]
        reduce_dict = Dict("uran" => 0.2, "lignite" => 0.2, "hard coal" => 0.7)
        for x in eachrow(plants)
            if x.fuel in fuel_reduction
                if (x.country == "PL") & (x.fuel == "hard coal")
                    x["g_max"] = x["g_max"] * 0.8
                elseif (x.country == "PL") & (x.fuel == "lignite")
                    x["g_max"] = x["g_max"] * 0.7
                elseif (x.country == "FR") & (x.fuel == "uran")
                    x["g_max"] =  x["g_max"] *0.7
                else
                    x["g_max"] =  x["g_max"] * reduce_dict[x.fuel]
                end
            end
        end
    end

    if scen_name == "Atomkraft_Nein_Danke"
        for x in eachrow(plants)    
            for x in eachrow(plants)
                if x.fuel == "uran"
                    x["g_max"] = 0
                end
            end
        end
    end

    if scen_name == "Minus_50"
        fuel_reduction = ["uran", "lignite", "hard coal"]
        for x in eachrow(plants)
            if x.fuel in fuel_reduction
                x["g_max"] =  x["g_max"] * 0.5
            end
        end
    end

    ############### Add new extra storage  #######################
    if extra_storage_factor != 0
        mc_el_sto = 20 # sets the marginal costs of the new storage in EUR per MWh
        g_max_fac = 1 # sets the factor to multiply the max generation by (times cumulated generation capacity in country  divided by 100) 
        eta_sto = 0.9 # sets the efficiency of storing in and out of new storage
        d_max_fac = 1 # sets the factor to multiply the max demand by (times cumulated generation capacity in country divided by 100) 
        storage_capacity_fac = 1 #  # sets the factor to multiply the max demand by (times cumulated generation capacity in country) 
        for zone in unique(ntc_data[:,:from_country])
            sum_gen_cap = sum(plants[plants[:,:country] .== zone, :g_max])
            push!(plants,
                [(zone * "_storage_NA") 
                zone 
                "undetermined Storage" 
                "NA" 
                mc_el_sto 
                g_max_fac * (sum_gen_cap/100) * extra_storage_factor
                eta_sto 
                1 
                d_max_fac*(sum_gen_cap/100) * extra_storage_factor
                storage_capacity_fac * sum_gen_cap* extra_storage_factor]) 
        end
    end

    ######################################################################
    ### create sets based on input data ###

    ### sets ###
    # T = 1:size(timeseries["demand_2015"], 1) |> collect
    T = 1:Time_steps |> collect
    TECH = unique(plants[:,:id] |> Vector)
    DISP = unique(plants[plants[:,:disp] .== 1, :id])
    NONDISP =  unique(plants[plants[:,:disp] .== 0, :id])
    S = plants[plants[:,:storage_capacity] .!= 0, :id]
    Z = unique(ntc_data[:,:from_country])
    CU = plants[plants[:,:tech] .== "CU", :id] |> Vector

    ### parameters ###

    ### Static Parameters ###
    mc = dictzip(plants, :id => :mc_el)
    g_max = dictzip(plants, :id => :g_max)
    d_max = dictzip(plants, :id => :d_max)
    eta = dictzip(plants, :id => :eta)
    storage_capacity = dictzip(plants, :id => :storage_capacity)

    ### Mappings of static data #####
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
        
    map_id2mc_el  = Dict()
    for plant_id in plants[:,:id] 
        map_id2mc_el[plant_id] = ((plants[plants.id .== plant_id, "mc_el"]))[1]
    end


    ### time series data processing ###

    ##### Hyrdo Inflows #####
    psp_inflow = coldict(timeseries["hydro_psp_inflow_countries"])
    PSP_list = plants[plants[:,:tech] .== "hydro_psp", :id]
    res_inflow = coldict(timeseries["hydro_res_inflow_countries"])
    RES_list = plants[plants[:,:tech] .== "hydro_res", :id]

    elec_demand = coldict(timeseries["demand_2015"])

    ####### Creating feed in for non dispatchables ######
    feed_in_wind_off = Dict()
    feed_in_wind_on = Dict()
    feed_in_pv = Dict()
    feed_in_ror = Dict()
    feed_in = Dict()

    for z in Z
        if ("wind_onshore" in  plants[(plants.country .== z),"tech"]) & (z in names(timeseries["wind_on_2015"]))
            feed_in_wind_off[z] = timeseries["wind_on_2015"][!,z] .* (plants[(plants.country .== z) .& (plants.tech .== "wind_onshore"), "g_max"])
        end

        if ("wind_offshore" in  plants[(plants.country .== z),"tech"]) & (z in names(timeseries["wind_off_2015"]))
            feed_in_wind_on[z] = timeseries["wind_off_2015"][!,z] .* (plants[(plants.country .== z) .& (plants.tech .== "wind_offshore"), "g_max"])
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
        feed_in_by_z_nondisp[x,"wind_offshore"] = feed_in_wind_off[x]
    end

    for x in keys(feed_in_wind_on)
        feed_in_by_z_nondisp[x,"wind_onshore"] = feed_in_wind_on[x]
    end

    for x in keys(feed_in_pv)
        feed_in_by_z_nondisp[x,"solar"] = feed_in_pv[x]
    end


    ###### NTCs #######
    ntc = dictzip(ntc_data, [:from_country, :to_country] => :ntc)

    ######### actual model creation ##########################

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
        + sum((BALANCE_P[z,t])*1000 for z in Z, t in T)
        + sum((CU[z,t]) for z in Z, t in T)
        + sum(EX[z,zz,t] for z in Z, zz in Z, t in T)
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

# fetching the results from model
result_G = DataFrame(G, [:id, :hour])
insertcols!(result_G, 2, :zone => [map_id2country[id] for id in result_G[!,:id]])
insertcols!(result_G, 3, :technology => [map_id2tech[id] for id in result_G[!,:id]])

result_CU = DataFrame(CU, [:zone, :hour])
result_CU[!,:technology] .= "curtailment"

result_BALANCE_P = DataFrame(BALANCE_P, [:zone, :hour])
result_BALANCE_P[!,:technology] .= "Balancing_Power"

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

    ############# Saving Data ######################
    
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
table2_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"FromTo_ntc_utilisation.csv"

CSV.write(joinpath(results_path,table2_name), table2)


# TABLE 3

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

table3_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"Country_ntc_import_export.csv"
CSV.write(joinpath(results_path, table3_name), table3)

# plotting load duration curve
country = "DE"
dem = elec_demand[country]
feed = feed_in[country]
wind_off = haskey(feed_in_wind_off, country) ? feed_in_wind_off[country] : zeros(length(T))
wind_on = haskey(feed_in_wind_on, country) ?  feed_in_wind_on[country] : zeros(length(T))
pv = haskey(feed_in_pv, country) ?  feed_in_pv[country] : zeros(length(T))
ror = haskey(feed_in_ror, country) ?  feed_in_ror[country] : zeros(length(T))

ldc_plot = plot_ldc2(dem, feed, wind_off, wind_on, pv, ror)

ldc_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"ldc_DE.png"
savefig(ldc_plot,joinpath(results_path,ldc_name))

# fetching mc results from generation of model
result_mc = DataFrame(G, [:id, :hour])
result_mc = result_mc[result_mc.value .!= 0, :]
insertcols!(result_mc, 2, :mc_el => [map_id2mc_el[id] for id in result_mc[!,:id]])
insertcols!(result_mc, 3, :zone => [map_id2country[id] for id in result_mc[!,:id]])
select!(result_mc, Not([:value, :id]))

result_balance_P = DataFrame(BALANCE_P, [:zone, :hour])
result_balance_P = result_balance_P[result_balance_P.value .!= 0, :]

# adding coss of 1000 EUR for the balance power / load shedding 
insertcols!(result_balance_P, 2, :mc_el => 1000)
select!(result_balance_P, Not(:value))
result_mc = vcat(result_mc, result_balance_P)

# saving table with highest marginal cost in each and zone (merit oder)
res_mc_grouped_by_zone_hourly = combine(groupby(copy(result_mc), [:zone, :hour]), :mc_el .=> [mean, maximum])
res_mc_grouped_by_zone_year_mean = combine(groupby(res_mc_grouped_by_zone_hourly, :zone), :mc_el_mean .=> [mean, std])
res_mc_grouped_by_zone_year_max = combine(groupby(res_mc_grouped_by_zone_hourly, :zone), :mc_el_maximum .=> [mean, std])

res_mc_grouped_by_zone_year_mean_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"res_mc_grouped_by_zone_year_mean.csv"
res_mc_grouped_by_zone_year_max_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"res_mc_grouped_by_zone_year_max.csv"
CSV.write(joinpath(results_path,res_mc_grouped_by_zone_year_max_name), res_mc_grouped_by_zone_year_max)
CSV.write(joinpath(results_path,res_mc_grouped_by_zone_year_mean_name), res_mc_grouped_by_zone_year_mean)


# saving the line constrains for the map creation
line_constrains_map = create_line_constrain_map(result_EX_data, ntc)
line_constrains_map_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"line_constrain_map.csv"
CSV.write(joinpath(results_path,line_constrains_map_name), line_constrains_map)

# saving the objective value
obj_val_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"objective_value.csv"
ob_val = DataFrame("Objective value:" => objective_value(m))
CSV.write(joinpath(results_path,obj_val_name), ob_val)
print(ob_val)

# saving the balanceing power / loadshedding inteances
result_balance_P = DataFrame(BALANCE_P, [:zone, :hour])
result_balance_P = result_balance_P[result_balance_P.value .!= 0, :]
balance_P_name = string(enable_2030)*"_"*scen_name*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"balance_power.csv"
CSV.write(joinpath(results_path, balance_P_name), result_balance_P)

end
# set first as true for 2030 renewable values, second: list of fuels to remove, third: ntc faktor to be multiplied, fourth: extra storage to be added (factor)
# setting the basic scenarios for fuel removal
Time_steps = 50
scenario_list = [[false, "Base", 1, 0, Time_steps], #1
                [true, "Base", 1, 0, Time_steps], #2
                [true, "Minus_50", 1 ,0, Time_steps], #3
                [true, "Atomkraft_Nein_Danke", 1, 0, Time_steps], #4
                [true, "Mixed_FR_PL_exept", 1, 0, Time_steps], #5
                ]

# itterating through scearios + ntc and storage changes 
# !!!!!!!!!! takes approx 4-5h !!!!!!!!!!!!!!!!!!
# for ntc_fac in [1,0.8,1.2], sto_fac in [0,0.5,1,3], scen in copy(scenario_list)
#     scen[3] = ntc_fac
#     scen[4] = sto_fac
#     println(scen)
#     run_model(scen...)
# end

scen = [false, "Base", 1, 0, Time_steps]
run_model(scen...)
scen = [true, "Base", 1, 0, Time_steps]
run_model(scen...)
scen = [true, "Minus_50", 1 ,0, Time_steps]
run_model(scen...)
scen = [true, "Atomkraft_Nein_Danke", 1, 0, Time_steps]
run_model(scen...)
scen = [true, "Mixed_FR_PL_exept", 1, 0, Time_steps]
run_model(scen...)


scenario_list = [[true, "Base", 1, 0, Time_steps], #2
                [true, "Mixed_FR_PL_exept", 1, 0, Time_steps], #5
                ]
for ntc_fac in [1,0.8,1.2], sto_fac in [0,0.5,1,3], scen in copy(scenario_list)
    scen[3] = ntc_fac
    scen[4] = sto_fac
    println(scen)
    run_model(scen...)
end

