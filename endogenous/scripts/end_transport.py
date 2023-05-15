#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:01:01 2023

@author: au731993
"""

import pypsa 
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import logging

# Add multilink functionality to PyPSA
override_component_attrs = pypsa.descriptors.Dict(
    {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
)
override_component_attrs["Link"].loc["bus2"] = [
    "string",
    np.nan,
    np.nan,
    "2nd bus",
    "Input (optional)",
]
override_component_attrs["Link"].loc["bus3"] = [
    "string",
    np.nan,
    np.nan,
    "3rd bus",
    "Input (optional)",
]
override_component_attrs["Link"].loc["efficiency2"] = [
    "static or series",
    "per unit",
    1.0,
    "2nd bus efficiency",
    "Input (optional)",
]
override_component_attrs["Link"].loc["efficiency3"] = [
    "static or series",
    "per unit",
    1.0,
    "3rd bus efficiency",
    "Input (optional)",
]
override_component_attrs["Link"].loc["p2"] = [
    "series",
    "MW",
    0.0,
    "2nd bus output",
    "Output",
]
override_component_attrs["Link"].loc["p3"] = [
    "series",
    "MW",
    0.0,
    "3rd bus output",
    "Output",
]

def annuity(n,r):
    """ Calculate the annuity factor 
    for an asset with lifetime n years and 
    discount rate of f, e.g. annuity(20,0.05)*20 = 1.6 """
    
    if r > 0:
        return r/(1. - 1./(1.+r)**n)
    else:
        return 1/n

def calculate_EV_timeseries(n,alpha):
    cars_driving = 0.5 # Assume a certain share of the fleet that can be driving at the same time
    L_t = n.loads_t.p_set['land_transport'] # load time series
    L_norm_t = L_t/max(L_t) # normalize time series
    EV_c = alpha*(1-L_norm_t) # we can charge all alpha [%] of parked cars  
    EV_d = cars_driving*L_norm_t # 70 % of EV fleet is available for driving 
    return EV_c, EV_d

def transport(string, costs):
    """
    Define endogenous transport PyPSA model
    """
    
    network = pypsa.Network(override_component_attrs=override_component_attrs)
    
    hours_in_2013 = pd.date_range('2013-01-01T00:00Z','2013-12-31T23:00Z', freq='H')
    network.set_snapshots(hours_in_2013)
    
    # Add oil bus and generator
    #discount rate 0.07
    oil_co2 = costs.at['oil', 'CO2 intensity']
    
    network.add("Carrier", "oil")
    
    network.add("Bus", 
                "oil bus", 
                carrier="oil")
    
    network.add("Generator",
                'oil',
                bus="oil bus",
                #p_nom=100,
                marginal_cost = costs.at['oil', 'fuel'], # EUR/MWh
                capital_cost = 0.,
                carrier = "oil",
                p_nom_extendable=True)
                #committable=True,
                #p_min_pu=0,
                #p_max_pu=1) 
    
    
    
    network.add("Carrier", "electricity")
    
    network.add("Bus", 
                "electricity bus", 
                carrier="electricity")
    network.add("Carrier", "solar")
    #add solar PV generator with data from one country
    df_solar = pd.read_csv(string + 'resources/pv_optimal.csv', sep=';', index_col=0)
    df_solar.index = pd.to_datetime(df_solar.index)
    CF_solar = df_solar[country][[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
    
    
    
    network.add("Generator",
                "solar",
                bus="electricity bus",
                p_nom_extendable=True,
                carrier="solar",
                #p_nom_max=1000, #maximum capacity can be limited due to environmental constraints
                capital_cost = 0, #1000*costs.at['solar', 'investment']*(annuity(costs.at['solar', 'lifetime'], 0.07)+costs.at['solar', 'FOM']/100),
                marginal_cost = costs.at['solar', 'VOM'],
                p_max_pu = CF_solar)
    
    # Add EV links
    network.add("Bus", "EV battery bus")

    if options["bev_dsm"]:
        network.add("Store", "EV battery storage",
                    bus="EV battery bus",
                    e_nom_extendable=True,
                    e_cyclic = True,
                    #capital_cost=costs.at['battery storage', 'investment'],
                    lifetime = costs.at['battery storage', 'lifetime'])
    
    EV_c,EV_d = calculate_EV_timeseries(network,alpha=0.5)

    network.add("Carrier",'EV battery')
    bev_charge_efficiency = options['bev_charge_efficiency']
    network.add("Link",
                "EV battery charger",
                bus0 = "electricity bus",
                bus1 = 'EV battery bus',
                carrier = 'EV battery charger',
                p_nom_extendable = 'True',
                p_max_pu = EV_c,
                efficiency = bev_charge_efficiency, #costs.at['battery inverter', 'efficiency']**0.5,
                #capital_cost = costs.at['battery inverter', 'investment']*(annuity(costs.at['battery inverter', 'lifetime'], 0.07)+costs.at['battery inverter', 'FOM']/100),
                lifetime = costs.at['battery inverter', 'lifetime'])

    
    
    # Only when v2g option is configured, add the v2g link
    if options["v2g"]:
        network.add("Link",
                    "EV battery discharger",
                    bus0 = "EV battery bus",
                    bus1 = 'electricity bus',
                    carrier = 'EV battery discharger',
                    p_max_pu = EV_c,
                    p_nom_extendable = 'True',
                    efficiency = costs.at['battery inverter', 'efficiency']**0.5,
                    capital_cost = costs.at['battery inverter', 'investment']*(annuity(costs.at['battery inverter', 'lifetime'], 0.07)+costs.at['battery inverter', 'FOM']/100),
                    lifetime = costs.at['battery inverter', 'lifetime'])
    
    network.add("Carrier", "land transport demand")
    
    network.add("Bus", 
                "land transport bus",
                carrier="land transport demand")
    
    network.add("Carrier", "co2", co2_emissions=1.)

    network.add("Bus", 
                "co2 atmosphere", 
                carrier="co2")
    
    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[country_short, 'road non-elec']
    co2_limit = co2_limit * limit
    
    network.add("Store",
                "co2 atmosphere",
                e_initial= 0, 
                e_nom = co2_limit,
                e_min_pu = -1,
                bus="co2 atmosphere",
                carrier="co2")

    # Network add ICEV
    ice_efficiency = options['transport_internal_combustion_efficiency']
    network.add(
        "Link",
        "ICE Vehicle",
        bus0="oil bus",  
        bus1="land transport bus",
        bus2 = "co2 atmosphere",                           
        carrier = "land transport demand",
        capital_cost = snakemake.config['costs']['capital_cost']['ICE_vehicle']/options['energy_to_cars'] * ice_efficiency,
        efficiency = ice_efficiency,  
        efficiency2 = oil_co2*ice_efficiency,
        p_nom_extendable=True,
    )
    
    # Add EV link
    ev_efficiency = options['bev_bat_to_wheel_efficiency']
    network.add(
        "Link",
        "EV",
        bus0="EV battery bus",                               
        bus1 = "land transport bus",
        carrier = "land transport demand",
        capital_cost = snakemake.config['costs']['capital_cost']['E_vehicle']/options['energy_to_cars']*(ev_efficiency*bev_charge_efficiency),
        efficiency = ev_efficiency,#bev_charge_efficiency,        
        p_nom_extendable=True,
    )

    # Add H2 cars 
    network.add(
        'Carrier',
        'H2',
    )

    network.add(
        'Bus',
        'H2 bus',
    )
    H2_vehicle_efficiency = options['H2_car_efficiency']
    network.add("Generator",
            'H2',
            bus="H2 bus",
            #p_nom=100,
            marginal_cost = snakemake.config['costs']['marginal_cost']['H2']/(options['energy_to_cars'])*H2_vehicle_efficiency, # EUR/MWh
            capital_cost = 0.,
            carrier = "H2",
            p_nom_extendable=True,
            efficiency=H2_vehicle_efficiency)
            #committable=True,
            #p_min_pu=0,
            #p_max_pu=1) 

    network.add(
        "Link",
        "H2 Vehicle",
        bus0="H2 bus",                               
        bus1 = "land transport bus",
        carrier = "land transport demand",
        efficiency= snakemake.config['sector']['H2_car_efficiency'],
        capital_cost = snakemake.config['costs']['capital_cost']['H2_vehicle'],#bev_charge_efficiency,        
        p_nom_extendable=True,
    )

    # Add loads

    #df_load  = pd.read_csv(snakemake.input.transport_demand, index_col=0, parse_dates=True)
    df_load = pd.read_csv(string + 'resources/transport_demand_s_45.csv', sep=',', index_col=0) 
    df_load.index = pd.to_datetime(df_load.index, utc=True)
    df_load.index.name = 'utc_time'
    load_p = df_load['DK1 0']
    # data for load: resources/transport_demand_s{simpl}_{clusters}.csv 
    #print(network.loads_t.p_set)
    print(network.carriers)
    network.add("Load",
                "load2",
                bus = "land transport bus",
                carrier = "land transport demand",
                p_set = load_p)
        
    
    network.add(
        "GlobalConstraint",
        "co2_limit",
        sense = "<=",
        carrier_attribute = "co2_emissions",
        constant = co2_limit,
    )
    return network


def add_v2g_constraint():
    lhs1 = network.model.variables['Link-p_nom'].sel({'Link-ext':'EV battery discharger'})
    lhs2 = network.model.variables['Link-p_nom'].sel({'Link-ext':'EV battery charger'})
    lhs = lhs1-lhs2
    rhs = 0
    network.model.add_constraints(lhs==rhs, name="constraint_v2g")

def add_EV_storage_constraint():
    lhs1 = network.model.variables['Link-p_nom'].sel({'Link-ext':'EV battery charger'})/options['EV_charge_rate']
    lhs2 = network.model.variables['Store-e_nom'].sel({'Store-ext':'EV battery storage'})/options['bev_energy']
    lhs = lhs1-lhs2
    rhs = 0
    network.model.add_constraints(lhs==rhs, name="constraint_EV_storage")

def add_EV_number_constraint():
    lhs1 = network.model.variables['Link-p_nom'].sel({'Link-ext':'EV battery charger'})/options['EV_charge_rate']
    lhs2 = network.model.variables['Link-p_nom'].sel({'Link-ext':'EV'})/options['EV_consumption_1car']
    lhs = lhs1-lhs2
    rhs = 0
    network.model.add_constraints(lhs==rhs, name="constraint_EV_number")

def extra_functionality(network, snapshots):
    
    #m = network.optimize.create_model() #for debugging
    print(network.model.variables['Link-p_nom'])
    print(network.model.variables['Link-p_nom']['EV battery charger'])
    network.model.variables['Link-p_nom'].sel({'Link-ext':'EV battery charger'})
    #network.model.variables['Link-p_nom'].sel('battery_charger') #use dictionary for network.model.variables().sel(key:value)
    if options["bev_dsm"]:
        add_EV_storage_constraint()
    add_EV_number_constraint()
    if options["v2g"]:
        add_v2g_constraint()
    
    # add_EV_timeconstraint(n)
    
    
def solve_network(network):

    status, condition = network.optimize(
            solver_name='gurobi',
            extra_functionality=extra_functionality)
    # status, status2 = network.lopf(network.snapshots, 
    #              pyomo=False,
    #              solver_name='gurobi')
    logger = logging.getLogger(__name__)
    pypsa.pf.logger.setLevel(logging.WARNING)

    if status != "ok":
        logger.warning(
            f"Solving status '{status}' with termination condition '{condition}'"
        )
    if "infeasible" in condition:
        raise RuntimeError("Solving status 'infeasible'")
        
    if(status == 'ok'):
        print(network)
        print(network.links)
        print(network.generators_t.p.sum())
        # print objective value
        print("Objective value: %f" % network.objective)
    
        # print oil generator value
        print("Oil generator value: %f" % network.generators.p_nom_opt["oil"])
    
        # print solar generator value
        print("Solar generator value: %f" % network.generators.p_nom_opt["solar"])
    
        network.generators_t.p.plot()

        plt.savefig(snakemake.output.results, dpi=1200, bbox_inches='tight')
        #plt.show()

        #print('year: ', Nyear)
        #print('co2 limit: ', co2_limit)

        network.export_to_netcdf(snakemake.output.network)

        return network
    else:
        print('model infeasable')
        return 1

def prepare_costs(cost_file, config, nyears):
    # set all asset costs and other parameters
    costs = pd.read_csv(cost_file, index_col=[0, 1]).sort_index()

    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = (
        costs.loc[:, "value"].unstack(level=1).groupby("technology").sum(min_count=1)
    )

    costs = costs.fillna(config["fill_values"])

    def annuity_factor(v):
        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100

    costs["fixed"] = [
        annuity_factor(v) * v["investment"] * nyears for i, v in costs.iterrows()
    ]

    return costs
    #get costs for technology
    df_costs = pd.read_csv(string + 'resources/costs_' + str(Nyear) + '.csv', sep=',', index_col=0)
    #set first column as second index and keep original index
    df_costs.set_index(df_costs.columns[0],append=True,inplace=True)
    #set parameter/value names as column name, min_count=1 fills not existing values with NaN
    costs = df_costs.loc[:,"value"].unstack(level=1).groupby("technology").sum(min_count=1)    

if __name__ == "__main__":
    string = '' 
    if 'snakemake' not in globals():
        print('run with snakemake')
        sys.exit(1)
        string = '../'
        
        class snakemake:
            output= pd.Series()
            output['results'] = string +'results/generators_profile.png'
            input = pd.Series()
            input['co2_totals_name'] = string + 'resources/co2_totals.csv'
            config = pd.Series()
            config['transport_internal_combustion_efficiency'] = 0.2
            config['bev_charge_efficiency'] = 0.8
    country = 'DNK'
    country_short = 'DK'

    Nyears = snakemake.config['scenario']['planning_horizons']
    print(Nyears)
    #for x in range(len(Nyears)):
    Nyear = snakemake.wildcards.planning_horizons #Nyears[x]

    Nyear = int(Nyear)
    limit = snakemake.config['co2_budget'].get(Nyear)
    cost_file = string + 'resources/costs_' + str(Nyear) + '.csv'
    nyears = 1
    costs = prepare_costs(
        cost_file,
        snakemake.config["costs"],
        nyears,
    )
    options = snakemake.config["sector"]
    network = transport(string, costs)
    solve_network(network)

