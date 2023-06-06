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
     ax.set_xlim([pd.to_datetime('5/5/2013', dayfirst=True),pd.to_datetime('14/5/2013', dayfirst=True)])
     ax.set_ylabel('Fraction of EVs [-]')
     ax.legend()

     # make plot of state of charge of EV battery (full year)
     fig1,ax1 = plt.subplots(2, figsize=(10,5))
     EV_soc.plot(ax=ax1[0],label='SOC')
     EV_charging_norm.plot(ax=ax1[1],label='charging') # plotting the fraction of cars being charged
     EV_discharging_norm.plot(ax=ax1[1],label='driving') # plotting the fraction of cars being used 
     for ax in ax1:
          ax.legend()

     # make plot of state of charge of EV battery (1 day)
     fig1_1,ax1_1 = plt.subplots(figsize=(10,5))
     EV_soc.plot(ax=ax1_1,label='SOC')
     #ax1_1.plot(n.generators_t.marginal_cost['solar'], label='marginal cost solar')
     ax1_1.set_xlim([pd.to_datetime('5/5/2013', dayfirst=True),pd.to_datetime('5/6/2013', dayfirst=True)])
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
     ax3.set_xlim([pd.to_datetime('5/5/2013', dayfirst=True),pd.to_datetime('14/5/2013', dayfirst=True)])
     ax3.legend()
     
     #compare charging/discharging with constraints 
     if n.links_t.p_min_pu.empty:
          EV_L_norm = pd.Series(np.full(len(network.snapshots), n.links.p_min_pu['EV']), index=network.snapshots)
     else:
          EV_L_norm = n.links_t.p_min_pu['EV']
     max_charge = n.links_t.p_max_pu['EV battery charger']
     fig4,ax4 = plt.subplots(2, sharex=True, sharey=True)
     EV_charging_norm.plot(ax=ax4[0],label='charging') # plotting the fraction of cars being charged
     EV_discharging_norm.plot(ax=ax4[1],label='driving') # plotting the fraction of cars being used 
     EV_L_norm.plot(ax=ax4[1],ls='-.',alpha=0.7,label='minimum driving required') #plotting the minimum amount of normalized load that needs to be fulfiled 
     max_charge.plot(ax=ax4[0], ls=':', alpha=0.7, label='maximum charging') #plot maximum of battery charge allowed
     ax4[0].set_ylabel('Charging of EVs [-], p_max_pu of EV battery')
     ax4[1].set_ylabel('Discharging of EVs [-], p_min_pu of EV')
     for ax in ax4:
        ax.set_xlim([pd.to_datetime('5/5/2013', dayfirst=True),pd.to_datetime('14/5/2013', dayfirst=True)])   
        ax.legend()

     return fig, fig1, fig1_1, fig2, fig3, fig4


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


# networks_dict = {
#     (simpl, planning_horizon): "results/" +
#     f"elec_s{simpl}_{planning_horizon}.nc"
#     for simpl in snakemake.config["scenario"]["simpl"]
#     for planning_horizon in snakemake.config["scenario"]["planning_horizons"]
# }
# dict_test = {key:value}
networks_dict = {snakemake.input.network: snakemake.input.network,}
networks = list()
for label, filename in networks_dict.items():
        logger.info(f"Make summary for scenario {label}, using {filename}")

        network = pypsa.Network(filename, override_component_attrs=overrides())
        fig,fig1,fig2,fig3,fig4, fig5 = plot_EV_timeseries(network)
        fig.savefig(snakemake.output.EV_balance, dpi=1200, bbox_inches='tight')
        fig1.savefig(snakemake.output.EV_soc_1year, dpi=1200, bbox_inches='tight')
        fig2.savefig(snakemake.output.EV_soc_1day, dpi=1200, bbox_inches='tight')
        fig3.savefig(snakemake.output.transport_balance_1year, dpi=1200, bbox_inches='tight')
        fig4.savefig(snakemake.output.transport_balance_period, dpi=1200, bbox_inches='tight')
        fig5.savefig(snakemake.output.transport_balance_and_load, dpi=1200, bbox_inches='tight')

        print('network', filename, label)
        networks.append(network)






