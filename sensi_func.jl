function run_model(enable_2030::Bool, remove_fuels_list::Vector, ntc_factor,extra_storage_factor)
    results_path  = joinpath(@__DIR__, "Results")
    if !isdir(results_path)
        mkdir(results_path)
    end
    remove_string = "WO_"
    for x in remove_fuels_list
        remove_string = remove_string*"_"*x
    end

    #  Preprocessing
    ### data load ###
    datapath = joinpath(@__DIR__, "data")
    # datapath = "data"

    timeseries_path = joinpath(datapath, "time_series")
    static_path = joinpath(datapath, "static")

    plants =  CSV.read(joinpath(static_path,"plants.csv"), DataFrame)
    ntc_data =  CSV.read(joinpath(static_path,"ntc.csv"), DataFrame)
    ntc_data[!,"ntc"] .= ntc_data[!,"ntc"] .* ntc_factor

    ############### switch this on for 2030 renewable capacity values #######################
    if enable_2030 == true
        renewables_2030 = CSV.read(joinpath(static_path,"renewables_2030.csv"), DataFrame)
        for x in eachrow(plants)
            for i in eachrow(renewables_2030)
                if (x.country == i.country) & (x.tech == i.tech) 
                    x["g_max"] =  i[:g_max_2030]
                end
            end
        end
    end

    ##################### Removing certain technologies (add to List "turn_of_fuel" to remove) ######################
    turn_of_fuel = remove_fuels_list
    for x in eachrow(plants)
        if x.fuel in turn_of_fuel
            x["g_max"] =  0
        end
        
    end

    ############### Add new storage tech #######################
    if extra_storage_factor != 0
        mc_el_sto = 20 # sets the marginal costs of the new storage in EUR per MWh
        g_max_fac = 1 # sets the factor to multiply the max generation by (times cumulated generation capacity in country  divided by 100) 
        eta_sto = 0.9 # sets the efficiency of storing in and out of new storage
        d_max_fac = 1 # sets the factor to multiply the max demand by (times cumulated generation capacity in country divided by 100) 
        storage_capacity_fac = 1 #  # sets the factor to multiply the max demand by (times cumulated generation capacity in country) 
        for zone in unique(ntc_data[:,:from_country])
            sum_gen_cap = sum(plants[plants[:,:country] .== zone, :g_max])
            println(sum_gen_cap)
            push!(plants,[(zone*"_storage_NA") 
                zone 
                "undetermined Storage" 
                "NA" 
                mc_el_sto 
                g_max_fac*(sum_gen_cap/100) * extra_storage_factor
                eta_sto 
                1 
                d_max_fac*(sum_gen_cap/100) * extra_storage_factor
                storage_capacity_fac*sum_gen_cap]) * extra_storage_factor
        end
    end
    #########################################################################################


    timeseries = Dict(splitext(files)[1] => CSV.read(joinpath(timeseries_path, files), DataFrame)
        for files in readdir(timeseries_path))



    ######################################################################
    ### create sets based on input data ###

    ### sets ###
    T = 1:size(timeseries["demand_2015"], 1) |> collect
    # T = 1:4000 |> collect
    TECH = unique(plants[:,:id] |> Vector)
    DISP = unique(plants[plants[:,:disp] .== 1, :id])
    NONDISP =  unique(plants[plants[:,:disp] .== 0, :id])
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
table2_name = string(enable_2030)*"_"*remove_string*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"FromTo_ntc_utilisation.csv"

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

table3_name = string(enable_2030)*"_"*remove_string*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"Country_ntc_import_export.csv"
CSV.write(joinpath(results_path, table3_name), table3)

ldc_plot = plot_ldc("DE")
ldc_name = string(enable_2030)*"_"*remove_string*"_"*string(ntc_factor)*"_"*string(extra_storage_factor)*"_"*"ldc_DE.png"
savefig(ldc_plot,ldc_name)
end
# set first as true for 2030 renewable values, second: list of fuels to remove, third: ntc faktor to be multiplied, fourth: extra storage to be added (factor)
# setting the basic scenarios for fuel removal
scenario_list = [[false, [], 1, 0], #1
                [true, [], 1, 0], #2
                [true, ["uran"], 1 ,0], #3
                [true, ["uran", "lignite"], 1, 0], #4
                [true, ["uran", "lignite", "hard coal"], 1, 0], #5
                ]
# itterating through ntc and storage changes 
# !!!!!!!!!! takes approx 2 to 3h !!!!!!!!!!!!!!!!!!
for ntc_fac in [0.8,1,1.2], sto_fac in [0,0.5,1,3], scen in copy(scenario_list[2:end])
    scen[3] = ntc_fac
    scen[4] = sto_fac
    println(scen)
    # run_model(scen...)``
end