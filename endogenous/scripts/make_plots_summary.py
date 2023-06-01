#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 10:06:17 2023

@author: au731993
"""

import pypsa 
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
import logging

def plot_EV_timeseries(n):
     # read EV time-dependent variables
     EV_charging = -n.links_t.p1['EV battery charger'] # from EV battery point-of-view
     EV_discharging = n.links_t.p0['EV'] # from EV battery point-of-view
     EV_soc = n.stores_t.e["EV battery storage"]/n.stores.loc["EV battery storage"].e_nom_opt

     # read general transportation time dependent variables
     ICE_discharging_p1 = -n.links_t.p1['ICE Vehicle'] # from land transport bus point-of-view
     EV_discharging_p1 = -n.links_t.p1['EV'] # from land transport bus point-of-view
     H2_p1 = -n.links_t.p1['H2 Vehicle'] # from land transport bus point-of-view
     load_t = n.loads_t.p_set['land transport'] # land transport load
     t_df = pd.DataFrame(index=n.snapshots) # collecting all cars in a pandas dataframe
     t_df['ICE'] = ICE_discharging_p1
     t_df['EV'] = EV_discharging_p1
     t_df['H2 Vehicle'] = H2_p1
     t_df[t_df < 0] = 0

     # read capacities
     EV_c_p_nom_opt = n.links.query('carrier == "EV battery charger"').p_nom_opt.sum() # EV charging power capacity
     EV_d_p_nom_opt = n.links.loc['EV'].p_nom_opt.sum() # EV "driving" power capacity

     # charge efficiency
     EV_c_eta = n.links.query('carrier == "EV battery charger"').efficiency # efficiency of EV battery chargers
          
     # normalize power with opt capacities to test constraints
     EV_charging_norm = EV_charging/(EV_c_p_nom_opt*EV_c_eta).item() # normalized charging (i.e., the fraction of cars being charged in each hour)
     EV_discharging_norm = EV_discharging/EV_d_p_nom_opt # normalized "driving" (i.e., the fraction of cars being used in each hour)

     # make plot of balancing of EV battery bus
     fig,ax = plt.subplots(figsize=(10,5))
     EV_charging_norm.plot(ax=ax,label='charging') # plotting the fraction of cars being charged
     EV_discharging_norm.plot(ax=ax,label='driving') # plotting the fraction of cars being used 
     ax.set_xlim([pd.to_datetime('5/5/2013'),pd.to_datetime('14/5/2013')])
     ax.set_ylabel('Fraction of EVs [-]')
     ax.legend()

     # make plot of state of charge of EV battery (full year)
     fig1,ax1 = plt.subplots(figsize=(10,5))
     EV_soc.plot(ax=ax1,label='SOC')
     ax1.legend()

     # make plot of state of charge of EV battery (1 day)
     fig1_1,ax1_1 = plt.subplots(figsize=(10,5))
     EV_soc.plot(ax=ax1_1,label='SOC')
     ax1_1.set_xlim([pd.to_datetime('5/5/2013'),pd.to_datetime('5/6/2013')])
     ax1_1.legend()

     # make plot of how land transport demand is met (full year)
     fig2,ax2 = plt.subplots(figsize=(10,5))
     t_df.resample('w').sum().plot.area(ax=ax2,stacked=True,alpha=0.5,lw=0)
     load_t.resample('w').sum().plot(ax=ax2,ls='--',lw=1, color='k',label='Load',zorder=10,alpha=0.5)
     ax2.legend()

     # make plot of how land transport demand is met (specified period)
     fig3,ax3 = plt.subplots(figsize=(10,5))
     t_df.plot.area(ax=ax3,stacked=True,alpha=0.5,lw=0)
     load_t.plot(ax=ax3,ls='--',lw=1, color='k',label='Load',zorder=10,alpha=0.5)
     ax3.set_xlim([pd.to_datetime('5/5/2013'),pd.to_datetime('14/5/2013')])
     ax3.legend()

     return fig, fig1, fig1_1, fig2, fig3

def plot_networks_together(networks, year):
    list_gen_p = list()
    list_gen_name = list()
    list_year = list()
    list_vehicle_name = list()
    nb_vehicle = list()
    count_year = -1
    list_vehicles = []
    # read data from all networks
    for network in networks:
        count_year+=1
        # get all links directly connected to land transport
        links = network.links.query('bus1 == "land transport bus"')
        vehicles = list()
        for i in links.index:            
             if i.find('Vehicle'):
                  ind = i.partition(' ')[0]
             #calculate number of vehicles
             vehicles.append(network.links.p_nom_opt[i]/snakemake.config["sector"][ind + '_consumption_1car'])
             if i == 'EV':
                print(vehicles)
                vehicles[-1] = vehicles[-1]*snakemake.config["sector"]['increase_nb_cars']
                print(vehicles)
                print(network.links_t.p1.sum())
             list_year.append(year[count_year])
             list_vehicle_name.append(i)
        veh = np.array(vehicles)
        veh = veh/sum(veh)
        vehicles = veh.tolist()
        nb_vehicle.append(vehicles)
    
        for index in network.generators.index: 
            #if list_gen_name.count(index) == 0:
            list_gen_name.append(index)
            #list_year.append(str(label[1]))
            print(network.generators.p_nom_opt[index])
            list_gen_p.append(network.generators_t.p[index].sum()/network.generators.p_nom_opt[index])
            print(network.generators_t.p[index].sum())
        # nb_EV = network.links.p_nom_opt['EV']/snakemake.config["sector"]['EV_consumption_1car']
        # nb_ICE = network.links.p_nom_opt['ICE Vehicle']/snakemake.config["sector"]['ICE_consumption_1car']
        # sum_cars = nb_EV + nb_ICE
        # nb_EV = nb_EV*snakemake.config["sector"]['increase_nb_cars']/sum_cars
        # nb_ICE = nb_ICE/sum_cars
        # print(network.loads_t.p.sum())
        # print(network.links_t.p1['EV'].sum())
        # print(network.links_t.p1['ICE Vehicle'].sum())
        # print('Ev, ICe', nb_EV, nb_ICE, nb_EV+nb_ICE)
        #nb_vehicle[str(label[1])]=[nb_ICE, nb_EV]
            
    index = pd.MultiIndex.from_tuples(tuple(zip(list_year,list_gen_name)))   
    generators = pd.DataFrame(list_gen_p,index=index)    
    print(generators) 
    generators = generators.unstack(level=-1)
    #plt.bar(snakemake.config["scenario"]["planning_horizons"], generators)
    for element in nb_vehicle:
        for nb in element:
            list_vehicles.append(nb)
    nb_vehicle = list_vehicles
    #nb_vehicle['index']=['EV', 'ICE']
    #nb_vehicle = nb_vehicle.set_index('index')
    index = pd.MultiIndex.from_tuples(tuple(zip(list_year,list_vehicle_name))) 
    print(index, nb_vehicle)
    nb_vehicles = pd.DataFrame(nb_vehicle, index=index)
    nb_vehicles = nb_vehicles.unstack(level=-1)
    print(nb_vehicles)
    plt.figure()
    nb_vehicles.plot.bar(stacked=True)
    plt.savefig(snakemake.output.vehiclenb, dpi=1200, bbox_inches='tight')

    plt.figure()
    generators.plot.bar(stacked=True)
    plt.savefig(snakemake.output.result, dpi=1200, bbox_inches='tight')



def overrides():
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

logger = logging.getLogger(__name__)


networks_dict = {
    (simpl, planning_horizon): "results/" +
    f"elec_s{simpl}_{planning_horizon}.nc"
    for simpl in snakemake.config["scenario"]["simpl"]
    for planning_horizon in snakemake.config["scenario"]["planning_horizons"]
}

networks = list()
year = list()
for label, filename in networks_dict.items():
        logger.info(f"Make summary for scenario {label}, using {filename}")

        network = pypsa.Network(filename, override_component_attrs=overrides())
        print('network', filename, label)
        networks.append(network)
        year.append(label[1])

plot_networks_together(networks,year)







