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

#makes it possible to add more busses to Links
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

def transport(string):
    network = pypsa.Network(override_component_attrs=override_component_attrs)
    
    hours_in_2013 = pd.date_range('2013-01-01T00:00Z','2013-12-31T23:00Z', freq='H')
    network.set_snapshots(hours_in_2013)
    
    df_tech = pd.read_csv(string + 'resources/costs_' + str(Nyear) + '.csv', sep=',', index_col=0)
    df_oil = df_tech.loc['oil']
    oil = df_oil.set_index(['parameter'])
    #discount rate 0.07
    oil_capital_cost = oil.at['investment','value']*1000*(annuity(oil.at['lifetime','value'], 0.07)+oil.at['FOM','value']/100) # in €/MW
    oil_marginal_cost = oil.at['VOM','value']
    oil_co2 = oil.at['CO2 intensity','value']
    
    network.add("Carrier", "co2_oil", co2_emissions=oil_co2)
    
    network.add("Bus", 
                "oil bus", 
                carrier="co2_oil")
    
    network.add("Generator",
                'oil',
                bus="oil bus",
                p_nom=100,
                marginal_cost = oil_marginal_cost, #million EUR/MWh
                capital_cost = 0,
                carrier = "co2_oil",
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
    solar = df_tech.loc['solar']
    solar = df_tech.loc['solar'].set_index(['parameter'])
    capital_cost_solar = solar.at['investment','value']*1000*(annuity(solar.at['lifetime','value'], 0.07)+solar.at['FOM','value']/100) # in €/MW
    # capital_cost_solar = annuity(25, 0.07)*425000*(1+0.03) # in €/MW
    network.add("Generator",
                "solar",
                bus="electricity bus",
                p_nom_extendable=True,
                carrier="solar",
                #p_nom_max=1000, #maximum capacity can be limited due to environmental constraints
                capital_cost = capital_cost_solar,
                marginal_cost = solar.at['VOM','value'],
                p_max_pu = CF_solar)
    
    network.add("Bus", "battery bus")

    df_costs = pd.read_csv(string + 'resources/costs_' + str(Nyear) + '.csv', sep=',', index_col=0)
    #set first column as second index and keep original index
    df_costs.set_index(df_costs.columns[0],append=True,inplace=True)
    #set parameter/value names as column name, min_count=1 fills not existing values with NaN
    costs = df_costs.loc[:,"value"].unstack(level=1).groupby("technology").sum(min_count=1)
    network.add("Store", "battery storage",
                bus="battery bus",
                e_nom_extendable=True,
                e_cyclic = True,
                capital_cost=costs.at['battery storage', 'investment'],
                lifetime = costs.at['battery storage', 'lifetime'])
    network.add("Carrier",'battery')
    network.add("Link",
                "battery charger",
                bus0 = "electricity bus",
                bus1 = 'battery bus',
                carrier = 'battery charger',
                p_nom_extendable = 'True',
                efficiency = costs.at['battery inverter', 'efficiency']**0.5,
                capital_cost = costs.at['battery inverter', 'investment']*(annuity(costs.at['battery inverter', 'lifetime'], 0.07)+costs.at['battery inverter', 'FOM']/100),
                lifetime = costs.at['battery inverter', 'lifetime'])
    network.add("Link",
                "battery discharger",
                bus0 = "battery bus",
                bus1 = 'electricity bus',
                carrier = 'battery discharger',
                p_nom_extendable = 'True',
                efficiency = costs.at['battery inverter', 'efficiency']**0.5,
                capital_cost = costs.at['battery inverter', 'investment']*(annuity(costs.at['battery inverter', 'lifetime'], 0.07)+costs.at['battery inverter', 'FOM']/100),
                lifetime = costs.at['battery inverter', 'lifetime'])
    network.add("Carrier", "vehicle")
    
    network.add("Bus", 
                "vehicle bus",
                carrier="vehicle")
    
    network.add("Carrier", "co2", co2_emissions=1.)

    network.add("Bus", 
                "co2 atmosphere", 
                carrier="co2")
    
    # convert Mt to tCO2
    co2_totals = 1e6 * pd.read_csv(snakemake.input.co2_totals_name, index_col=0)

    co2_limit = co2_totals.loc[country_short, 'road non-elec']
    co2_limit = co2_limit * limit
    #co2_limit = co2_limit * Nyear
    target = 0#3.392052e+03/8760
    
    network.add("Store",
                "co2 atmosphere",
                e_initial= 0, 
                e_nom = co2_limit,
                e_min_pu = -1,
                bus="co2 atmosphere", 
                carrier="co2")

    ice_efficiency = snakemake.config['transport_internal_combustion_efficiency']
    bev_charge_efficiency = snakemake.config['bev_charge_efficiency']
    network.add(
        "Link",
        "ICE Vehicle",
        bus0="oil bus",  
        bus1="vehicle bus",
        bus2 = "co2 atmosphere",                               
        carrier = "vehicle",
        efficiency = ice_efficiency,  
        efficiency2 = oil_co2*ice_efficiency,
        p_nom_extendable=True,
        #p_nom_max= ,
        #lifetime = ,
        #capital_cost =  
    )
    
    network.add(
        "Link",
        "EV",
        bus0="electricity bus",                               
        bus1 = "vehicle bus",
        carrier = "vehicle",
        efficiency = bev_charge_efficiency,        
        p_nom_extendable=True,
        #p_nom_max= ,
        #lifetime = ,
        #capital_cost =  
    )
    #pop_layout = pd.read_csv('../resources/pop_layout_elec_s_45.csv', index_col=0)
    #df_load  = pd.read_csv(snakemake.input.transport_demand, index_col=0, parse_dates=True)
    df_load = pd.read_csv(string + 'resources/transport_demand_s_45.csv', sep=',', index_col=0) 
    df_load.index = pd.to_datetime(df_load.index, utc=True)
    df_load.index.name = 'utc_time'
    load_p = df_load['DK1 0']
    # data for load: resources/transport_demand_s{simpl}_{clusters}.csv 

    network.add("Load",
                "load2",
                bus = "vehicle bus",
                carrier = "vehicle",
                p_set = load_p)
        
    
    network.add(
        "GlobalConstraint",
        "co2_limit",
        sense = "<=",
        carrier_attribute = "co2_emissions",
        constant = co2_limit,
    )
    
    #print(network.loads_t.p_set)
    print(network.carriers)
    
    status, status2 = network.lopf(network.snapshots, 
                 pyomo=False,
                 solver_name='gurobi')
    
    if(status == 'ok'):
        # print objective value
        print("Objective value: %f" % network.objective)
    
        # print oil generator value
        print("Oil generator value: %f" % network.generators.p_nom_opt["oil"])
    
        # print solar generator value
        print("Solar generator value: %f" % network.generators.p_nom_opt["solar"])
    
        network.generators_t.p.plot()

        plt.savefig(snakemake.output.results, dpi=1200, bbox_inches='tight')
        #plt.show()

        print('year: ', Nyear)
        print('co2 limit: ', co2_limit)
        return network
    else:
        print('model infeasable')
        return 1
    

if __name__ == "__main__":
    string = '' 
    if 'snakemake' not in globals():
        print('run with snakemake')
        #sys.exit(1)
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
    print(Nyear)
    Nyear = int(Nyear)
    limit = snakemake.config['co2_budget'].get(Nyear)
    
    network = transport(string)

