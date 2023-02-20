#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 16:01:01 2023

@author: au731993
"""

import pypsa 
import pandas as pd
import numpy as np

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
    
    df_tech = pd.read_csv(string + 'resources/costs_2025.csv', sep=',', index_col=0)
    df_oil = df_tech.loc['oil']
    oil = df_oil.set_index(['parameter'])
    #discount rate 0.07
    oil_capital_cost = annuity(oil.at['lifetime','value'], 0.07)*oil.at['investment','value']*1000*(1+oil.at['FOM','value']/100) # in €/MW
    oil_marginal_cost = oil.at['VOM','value']
    
    network.add("Carrier", "oil", co2_emissions=1.)
    
    network.add("Bus", 
                "oil bus", 
                carrier="oil")
    
    network.add("Generator",
                'oil',
                bus="oil bus",
                p_nom=100,
                marginal_cost = oil_marginal_cost, #million EUR/MWh
                capital_cost = oil_capital_cost,
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
    CF_solar = df_solar['ESP'][[hour.strftime("%Y-%m-%dT%H:%M:%SZ") for hour in network.snapshots]]
    solar = df_tech.loc['solar']
    solar = df_tech.loc['solar'].set_index(['parameter'])
    capital_cost_solar = annuity(solar.at['lifetime','value'], 0.07)*solar.at['investment','value']*1000*(1+solar.at['FOM','value']/100) # in €/MW
    # capital_cost_solar = annuity(25, 0.07)*425000*(1+0.03) # in €/MW
    network.add("Generator",
                "solar",
                bus="electricity bus",
                p_nom_extendable=True,
                carrier="solar",
                #p_nom_max=1000, #maximum capacity can be limited due to environmetal constraints
                capital_cost = capital_cost_solar,
                marginal_cost = 0.01,
                p_max_pu = CF_solar)
    
    
    network.add("Carrier", "vehicle")
    
    network.add("Bus", 
                "vehicle bus",
                carrier="vehicle")
    
    
    network.add(
        "Link",
        "ICE Vehicle",
        bus0="oil bus",  
        bus1="vehicle bus",                               
        carrier = "vehicle",
        efficiency = 1.,        
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
        efficiency = 1.,        
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
    load_p = df_load['AT1 0']
    # data for load: resources/transport_demand_s{simpl}_{clusters}.csv 
    
    network.add("Load",
                "load2",
                bus = "vehicle bus",
                carrier = "vehicle",
                p_set = load_p)
    
    #print(network.loads_t.p_set)
    
    network.lopf(network.snapshots, 
                 pyomo=False,
                 solver_name='gurobi')
    
    network.generators_t.p.plot()

if __name__ == "__main__":
    string = ''
    if 'snakemake' not in globals():
        string = '../'
    transport(string)

