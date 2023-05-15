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
     # read time-dependent variables
     EV_charging = -n.links_t.p1['EV battery charger'] # from EV battery point-of-view
     EV_discharging = n.links_t.p0['EV'] # from EV battery point-of-view

     ICE_discharging_p1 = -n.links_t.p1['EV'] # from land transport bus point-of-view
     EV_discharging_p1 = -n.links_t.p1['EV'] # from EV battery point-of-view
     load_t = n.loads_t.p_set['land transport']

     t_df = pd.DataFrame(index=n.snapshots)
     t_df['ICE'] = ICE_discharging_p1
     t_df['EV'] = EV_discharging_p1

     # read capacities
     EV_c_p_nom_opt = n.links.query('carrier == "EV battery charger"').p_nom_opt.sum()
     EV_d_p_nom_opt = n.links.loc['EV'].p_nom_opt.sum()

     EV_c_eta = n.links.query('carrier == "EV battery charger"').efficiency
          
     # normalize power with opt capacities
     EV_charging_norm = EV_charging/(EV_c_p_nom_opt*EV_c_eta).item()
     EV_discharging_norm = EV_discharging/EV_d_p_nom_opt

     # make plot of balancing of EV battery bus
     fig,ax = plt.subplots(figsize=(10,5))
     EV_charging_norm.plot(ax=ax,label='charging')
     EV_discharging_norm.plot(ax=ax,label='driving')
     ax.set_xlim([pd.to_datetime('5/5/2013'),pd.to_datetime('14/5/2013')])
     ax.legend()

     # make plot of how land transport demand is met
     fig1,ax1 = plt.subplots(figsize=(10,5))
     t_df.plot.area(ax=ax1,stacked=True,alpha=0.5,lw=0)
     load_t.plot(ax=ax1,ls='--',lw=1, color='k',label='Load',zorder=10,alpha=0.5)
     ax1.set_xlim([pd.to_datetime('5/5/2013'),pd.to_datetime('14/5/2013')])
     ax1.legend()
     return fig, fig1



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

list_gen_p = list()
list_gen_name = list()
list_year = list()
nb_vehicle = pd.DataFrame()
for label, filename in networks_dict.items():
        logger.info(f"Make summary for scenario {label}, using {filename}")

        network = pypsa.Network(filename, override_component_attrs=overrides())
        print('network', filename, label)
        for index in network.generators.index: 
            #if list_gen_name.count(index) == 0:
            list_gen_name.append(index)
            list_year.append(str(label[1]))
            print(network.generators.p_nom_opt[index])
            list_gen_p.append(network.generators_t.p[index].sum()/network.generators.p_nom_opt[index])
            print(network.generators_t.p[index].sum())
        nb_EV = network.links.p_nom_opt['EV']/snakemake.config["sector"]['EV_consumption_1car']
        nb_ICE = network.links.p_nom_opt['ICE Vehicle']/snakemake.config["sector"]['ICE_consumption_1car']
        sum_cars = nb_EV + nb_ICE
        nb_EV = nb_EV/sum_cars
        nb_ICE = nb_ICE/sum_cars
        print(network.loads_t.p.sum())
        print(network.links_t.p1['EV'].sum())
        print(network.links_t.p1['ICE Vehicle'].sum())
        print('Ev, ICe', nb_EV, nb_ICE, nb_EV+nb_ICE)
        nb_vehicle[str(label[1])]=[nb_ICE, nb_EV]
        
index = pd.MultiIndex.from_tuples(tuple(zip(list_year,list_gen_name)))   
generators = pd.DataFrame(list_gen_p,index=index)    
print(generators) 
generators = generators.unstack(level=-1)
#plt.bar(snakemake.config["scenario"]["planning_horizons"], generators)

plt.figure()
generators.plot.bar(stacked=True)
plt.savefig(snakemake.output.result, dpi=1200, bbox_inches='tight')

fig,fig1 = plot_EV_timeseries(network)
fig.savefig(snakemake.output.EV_times1, dpi=1200, bbox_inches='tight')
fig1.savefig(snakemake.output.EV_times2, dpi=1200, bbox_inches='tight')


nb_vehicle['index']=['EV', 'ICE']
nb_vehicle = nb_vehicle.set_index('index')
print(nb_vehicle)
plt.figure()
nb_vehicle.T.plot.bar(stacked=True)
plt.savefig(snakemake.output.vehiclenb, dpi=1200, bbox_inches='tight')

