#!/usr/bin/env python3


# This Python module models an exogenous-driven land transport scenario involving ICE/BE-based vehicles (graphically represented in
# https://github.com/ricnogfer/pypsa-eur-sec_land_transport/blob/master/resources/docs/exogenous_endogenous_land_transport.pdf)



# import necessary packages
import sys
import yaml
import pypsa



# function to create PyPSA-based network from scratch
def create_network(oil_capital_cost, oil_marginal_cost, solar_capital_cost, solar_marginal_cost, ice_share, bev_share):
    """
    Parameters
    ----------
    oil_capital_cost : float
        Capital cost of oil generator
    oil_marginal_cost : float
        Marginal cost of oil generator
    solar_capital_cost : float
        Capital cost of solar generator
    solar_marginal_cost : float
        Marginal cost of solar generator
    ice_share : Panda dataframe
        ICE-vehicle share values for horizon planning
    bev_share : Panda dataframe
        BE-vehicle share values forhorizon planning

    Returns
    -------
    A Network object representing a PyPSA-based network.
    """

    # create (empty) network
    network = pypsa.Network()


    # add bus "oil" to network
    network.add("Bus",
                "oil")


    # add generator "oil" to bus "oil" with given capital and marginal costs
    network.add("Generator",
                "oil",
                bus = "oil",
                capital_cost = oil_capital_cost,
                marginal_cost = oil_marginal_cost)


    # add load "ICE" to bus "oil" with given ICE-vehicle share values
    network.add("Load",
                "ICE",
                bus = "oil",
                p_set = ice_share)


    # add bus "electricity" to network
    network.add("Bus",
                "electricity")


    # add generator "solar" to bus "electricity" with given capital and marginal costs
    network.add("Generator",
                "solar",
                bus = "electricity",
                capital_cost = solar_capital_cost,
                marginal_cost = solar_marginal_cost)


    # add load "BEV" to bus "electricity" with given BE-vehicle share values
    network.add("Load",
                "BEV",
                bus = "electricity",
                p_set = bev_share)


    # return network
    return network



# run code (if present module is not used by another one)
if __name__ == "__main__":


    # read "config.yaml" file
    with open("config.yaml", "r") as stream:
        try:
            config = yaml.safe_load(stream)
        except Exception as e:
            print(e)
            sys.exit(-1)   # exit unsuccessfully


    # get capital and marginal costs for oil and solar generators
    oil_capital_cost = config["costs"]["capital_cost"]["oil"]
    oil_marginal_cost = config["costs"]["marginal_cost"]["oil"]
    solar_capital_cost = config["costs"]["capital_cost"]["solar"]
    solar_marginal_cost = config["costs"]["marginal_cost"]["solar"]


    # get ICE/BE-vehicles shares for horizon planning
    ice_share = config["sector"]["land_transport_ice_share"]
    bev_share = config["sector"]["land_transport_electric_share"]


    # TODO: create Panda dataframe from ICE/BE-vehicles shares (so that PyPSA may consume it)


    # create network
    network = create_network(oil_capital_cost,
                             oil_marginal_cost,
                             solar_capital_cost,
                             solar_marginal_cost,
                             0,
                             0)

    # print network content summary
    print(network)


    # solve network
    network.lopf(pyomo = False, solver_name = "gurobi")


    # print objective value
    objective_value = network.objective
    print("Objective value=%d" % objective_value)


    # print optimal ICE-based vehicles
    optimal_ICE_vehicles = network.generators.p_nom_opt["oil"]
    print("Optimal ICE-based vehicles=%d" % optimal_ICE_vehicles)


    # print optimal BE-based vehicles
    optimal_BE_vehicles = network.generators.p_nom_opt["solar"]
    print("Optimal BE-based vehicles=%d" % optimal_BE_vehicles)


    # exit successfully
    sys.exit(0)

