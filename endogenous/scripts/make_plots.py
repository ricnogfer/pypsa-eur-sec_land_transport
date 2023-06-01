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


nb_vehicle['index']=['EV', 'ICE']
nb_vehicle = nb_vehicle.set_index('index')
print(nb_vehicle)
plt.figure()
nb_vehicle.T.plot.bar(stacked=True)
plt.savefig(snakemake.output.vehiclenb, dpi=1200, bbox_inches='tight')

