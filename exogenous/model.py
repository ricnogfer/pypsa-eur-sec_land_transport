#!/usr/bin/env python3


# This Python module implements an exogenous-driven land transport (toy) model involving ICE/BE-based vehicles (represented in
# https://github.com/ricnogfer/pypsa-eur-sec_land_transport/blob/master/resources/docs/exogenous_endogenous_land_transport.pdf)



# TODO:
# 1) generate results for multiple years (as specified in the "sector:land_transport_electric_share" and "sector:land_transport_ice_share" options in the config.yaml file)
# 2) use CO2 emission values from CSV file (instead of being hard-coded to 1.0 as it is currently)



# import necessary modules
import os
import sys
import pypsa
import pandas



# declare global variables
_COUNTRY = "DNK"   # country to be analised (chosen arbitrarly)
_NODE = "DK1 0"   # node within country to be analised (chosen arbitrarly)



# function to create model (i.e. PyPSA network)
def create_model(parameters):
    """
    Parameters
    ----------
    parameters : a Python dictionary (dict())
        Parameters needed to create and configure the model (i.e. PyPSA network) properly.

    Returns
    -------
    network : a PyPSA network object (pypsa.Network())
        The created model (i.e. PyPSA network).
    """

    # create (empty) PyPSA network
    network = pypsa.Network()
    network.set_snapshots(parameters["snapshots"])


    # add carrier "oil" to network
    network.add("Carrier",
                "oil",
                co2_emissions = 1.0)


    # add bus "oil" to network
    network.add("Bus",
                "oil",
                carrier = "oil")


    # add generator "oil" to bus "oil" with associated costs
    network.add("Generator",
                "oil",
                bus = "oil",
                carrier = "oil",
                p_nom_extendable = True,
                capital_cost = parameters["oil_capital_cost"],
                marginal_cost = parameters["oil_marginal_cost"])


    if parameters["ICE_shares"][0] > 0:

        # add load "ICE" to bus "oil" with associated demand
        value = parameters["transport_demand"] * parameters["ICE_shares"][0] / parameters["ICE_efficiency"]
        network.add("Load",
                    "ICE",
                    bus = "oil",
                    p_set = value)


    # add bus "electricity" to network
    network.add("Bus",
                "electricity")


    # add generator "solar" to bus "electricity" with associated solar profile and costs
    solar_profile = parameters["solar_profile"][[snapshot.strftime("%Y-%m-%dT%H:%M:%SZ") for snapshot in network.snapshots]]   # select solar profile period based on the snapshots period
    network.add("Generator",
                "solar",
                bus = "electricity",
                p_nom_extendable = True,
                #p_max_pu = solar_profile,   # fix solar profile
                capital_cost = parameters["solar_capital_cost"],
                marginal_cost = parameters["solar_marginal_cost"])


    if parameters["BEV_shares"][0] > 0:

        # add bus "BEV" to network
        network.add("Bus",
                    "BEV")


        # add load "BEV" to bus "BEV" with associated demand
        value = parameters["transport_demand"] * parameters["BEV_shares"][0]
        network.add("Load",
                    "BEV",
                    bus = "BEV",
                    p_set = value)


        # add link "electricity_2_BEV" which connects bus "electricity" to bus "BEV" with associated availability and efficiency
        value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][0] * parameters["BEV_charge_rate"]
        network.add("Link",
                    "electricity_2_BEV",
                    bus0 = "electricity",
                    bus1 = "BEV",
                    p_nom = value,
                    p_max_pu = parameters["availability_profile"].values,
                    efficiency = parameters["BEV_charge_efficiency"])


        # enable BE-vehicle batteries to send power to the grid  (i.e. bus "electricity")
        if parameters["BEV_V2G"]:
            value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][0] * parameters["BEV_charge_rate"]
            network.add("Link",
                        "BEV_2_electricity",
                        bus0 = "BEV",
                        bus1 = "electricity",
                        p_nom = value,
                        p_max_pu = parameters["availability_profile"].values,
                        efficiency = parameters["BEV_charge_efficiency"])


        # enable BE-vehicle batteries based on demand-side management (DSM) profile
        if parameters["BEV_DSM"]:
            value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][0] * parameters["BEV_availability"] * parameters["BEV_energy"]
            network.add("Store",
                        "battery",
                        bus = "BEV",
                        e_cyclic = True,
                        e_nom = value,
                        e_max_pu = 1,
                        e_min_pu = parameters["DSM_profile"].values)


    # add global constraint "global_constraint_co2" for CO2 emissions
    network.add("GlobalConstraint",
                "global_constraint_co2",
                carrier_attribute = "co2_emissions",
                sense = "<=",
                constant = 10**8)   # high value so that solver can find a solution (this is just to illustrate the usage of the CO2 global constraint - no practical use for this (toy) model though)


    return network



# function to get Snakemake parameters
def get_snakemake_parameters():
    """
    Returns
    -------
    parameters : a Python dictionary (dict())
        Contains the Snakemake parameters needed to create the model (i.e. PyPSA network) parameters.
    """

    parameters = dict()

    # check how present module was launched
    if "snakemake" in globals():   # through Snakemake
        parameters["technology_costs_file"] = snakemake.input["technology_costs_file"]
        parameters["solar_profile_file"] = snakemake.input["solar_profile_file"]
        parameters["transport_data_file"] = snakemake.input["transport_data_file"][0]
        parameters["transport_demand_file"] = snakemake.input["transport_demand_file"][0]
        parameters["availability_profile_file"] = snakemake.input["availability_profile_file"][0]
        parameters["DSM_profile_file"] = snakemake.input["dsm_profile_file"][0]
        parameters["ICE_shares"] = tuple(snakemake.config["sector"]["land_transport_ice_share"].values())
        parameters["BEV_shares"] = tuple(snakemake.config["sector"]["land_transport_electric_share"].values())
        parameters["ICE_efficiency"] = snakemake.config["sector"]["transport_internal_combustion_efficiency"]
        parameters["BEV_DSM"] = snakemake.config["sector"]["bev_dsm"]
        parameters["BEV_availability"] = snakemake.config["sector"]["bev_availability"]
        parameters["BEV_energy"] = snakemake.config["sector"]["bev_energy"]
        parameters["BEV_charge_rate"] = snakemake.config["sector"]["bev_charge_rate"]
        parameters["BEV_charge_efficiency"] = snakemake.config["sector"]["bev_charge_efficiency"]
        parameters["BEV_V2G"] = snakemake.config["sector"]["v2g"]
        parameters["result_txt_file"] = snakemake.output["result_txt_file"]
        parameters["result_png_file"] = snakemake.output["result_png_file"]
        parameters["snapshots"] = pandas.date_range("%sT00:00Z" % snakemake.config["snapshots"]["start"], "%sT23:00Z" % snakemake.config["snapshots"]["end"], freq = "H")
    else:   # through terminal
        resource_path = "../resources/data"
        result_path = "results"
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        parameters["technology_costs_file"] = "%s/costs_2025.csv" % resource_path
        parameters["solar_profile_file"] = "%s/solar_profile_1979_2017.csv" % resource_path
        parameters["transport_data_file"] = "%s/transport_data_s_37.csv" % resource_path
        parameters["transport_demand_file"] = "%s/transport_demand_s_37.csv" % resource_path
        parameters["availability_profile_file"] = "%s/avail_profile_s_37.csv" % resource_path
        parameters["DSM_profile_file"] = "%s/dsm_profile_s_37.csv" % resource_path
        parameters["ICE_shares"] = (0.9, 0.75, 0.4, 0.0)
        parameters["BEV_shares"] = (0.1, 0.25, 0.6, 1.0)
        parameters["ICE_efficiency"] = 0.3
        parameters["BEV_DSM"] = True
        parameters["BEV_availability"] = 0.5
        parameters["BEV_energy"] = 0.05
        parameters["BEV_charge_efficiency"] = 0.9
        parameters["BEV_charge_rate"] = 0.011
        parameters["BEV_V2G"] = True
        parameters["result_txt_file"] = "%s/%s" % (result_path, "results.txt")
        parameters["result_png_file"] = "%s/%s" % (result_path, "results.png")
        parameters["snapshots"] = pandas.date_range("2013-01-01T00:00Z", "2013-12-31T23:00Z", freq = "H")


    return parameters



# function to get parameters that the model (i.e. PyPSA network) will be based upon
def get_model_parameters(snakemake_parameters):
    """
    Parameters
    ----------
    snakemake_parameters : a Python dictionary (dict())
        Contains the Snakemake parameters needed to create the model (i.e. PyPSA network) parameters.

    Returns
    -------
    parameters : a Python dictionary (dict())
        Contains the parameters that the model (i.e. PyPSA network) will be based upon.
    """



    def calculate_annuity(period, discount_rate):
        """
        Parameters
        ----------
        period : a Python integer (int)
            Contains the fixed length of time to receive payments from an investment.
        discount_rate : a Python float (float)
            Contains the assumed rate of return that is used to determine the present value of future payments.

        Returns
        -------
        None.
        """

        if discount_rate > 0:
            return discount_rate / (1.0 - 1.0 / (1.0 + discount_rate) ** period)

        return 1 / period



    parameters = dict()


    # read technology costs (CSV) file
    technology_costs = pandas.read_csv(snakemake_parameters["technology_costs_file"], index_col = "technology")


    # get oil generator costs
    oil = technology_costs.loc["oil"].set_index("parameter")
    parameters["oil_capital_cost"] = 0   # set to 0 since oil is imported to Europe (i.e. no relevant capital cost associated with oil generators)
    parameters["oil_marginal_cost"] = oil.at["fuel", "value"]


    # get solar generator costs
    solar = technology_costs.loc["solar"].set_index("parameter")
    solar_investment = solar.at["investment", "value"] * 1000.0   # in EUR/MW
    solar_annuity = calculate_annuity(solar.at["lifetime", "value"], 0.05) + solar.at["FOM", "value"] / 100.0   # in percentage
    parameters["solar_capital_cost"] = solar_investment * solar_annuity
    parameters["solar_marginal_cost"] = solar.at["VOM", "value"]


    # read solar profile CSV file and get profile for chosen country
    solar_profile = pandas.read_csv(snakemake_parameters["solar_profile_file"], index_col = 0, sep = ";")
    solar_profile.index = pandas.to_datetime(solar_profile.index)
    parameters["solar_profile"] = solar_profile[_COUNTRY]


    # read transport data (CSV) file and get data for chosen node
    transport_data = pandas.read_csv(snakemake_parameters["transport_data_file"], index_col = "name")
    parameters["transport_data"] = transport_data.loc[_NODE]


    # read transport demand (CSV) file and get demand for chosen node
    transport_demand = pandas.read_csv(snakemake_parameters["transport_demand_file"], index_col = 0)
    transport_demand.index = pandas.to_datetime(transport_demand.index, utc = True)
    parameters["transport_demand"] = transport_demand[_NODE]


    # read availability CSV file and get availability for chosen node
    availability_profile = pandas.read_csv(snakemake_parameters["availability_profile_file"], index_col = 0, parse_dates = True)
    parameters["availability_profile"] = availability_profile[_NODE]


    # read DSM profile CSV file and get profile for chosen node
    dsm_profile = pandas.read_csv(snakemake_parameters["DSM_profile_file"], index_col = 0, parse_dates = True)
    parameters["DSM_profile"] = dsm_profile[_NODE]


    # get snapshots information
    parameters["snapshots"] = snakemake_parameters["snapshots"]


    # copy certain Snakemake parameters that are needed by the model (i.e. PyPSA network)
    parameters["ICE_shares"] = snakemake_parameters["ICE_shares"]
    parameters["BEV_shares"] = snakemake_parameters["BEV_shares"]
    parameters["ICE_efficiency"] = snakemake_parameters["ICE_efficiency"]
    parameters["BEV_DSM"] = snakemake_parameters["BEV_DSM"]
    parameters["BEV_availability"] = snakemake_parameters["BEV_availability"]
    parameters["BEV_energy"] = snakemake_parameters["BEV_energy"]
    parameters["BEV_charge_rate"] = snakemake_parameters["BEV_charge_rate"]
    parameters["BEV_charge_efficiency"] = snakemake_parameters["BEV_charge_efficiency"]
    parameters["BEV_V2G"] = snakemake_parameters["BEV_V2G"]


    return parameters



# run code (if present module is not imported by another one)
if __name__ == "__main__":


    # get Snakemake parameters
    snakemake_parameters = get_snakemake_parameters()


    # check that ICE/BE-vehicle shares are properly specified
    if not snakemake_parameters["ICE_shares"] or not snakemake_parameters["BEV_shares"]:
        print("The ICE vehicle shares and/or BE vehicle shares is/are not specified!")
        sys.exit(-1)   # exit unsuccessfully
    if len(snakemake_parameters["ICE_shares"]) != len(snakemake_parameters["BEV_shares"]):
        print("The number of ICE vehicle shares does not match the number of BE vehicle shares!")
        sys.exit(-1)   # exit unsuccessfully
    for i in range(len(snakemake_parameters["ICE_shares"])):
        if snakemake_parameters["ICE_shares"][i] + snakemake_parameters["BEV_shares"][i] != 1:
            print("The sum of ICE vehicle share with BE vehicle share is not equal to 1!")
            sys.exit(-1)   # exit unsuccessfully


    # get model parameters
    model_parameters = get_model_parameters(snakemake_parameters)


    # create model (i.e. PyPSA network) based on the parameters
    network = create_model(model_parameters)


    # solve network using Gurobi
    status = network.lopf(pyomo = False, solver_name = "gurobi")
    if status[0] != "ok":
        print("The solver did not reach a solution (infeasible or unbounded)!")
        sys.exit(-1)   # exit unsuccessfully


    # print results
    print("Objective value=%.2f M€" % (network.objective / 10**6))
    print("Optimal nominal power oil generator=%.2f MW" % network.generators.p_nom_opt["oil"])
    print("Optimal nominal power solar generator=%.2f MW" % network.generators.p_nom_opt["solar"])


    # write results to file
    with open(snakemake_parameters["result_txt_file"], "w") as handle:
        try:
            handle.write("Objective value=%.2f M€\n" % (network.objective / 10**6))
            handle.write("Optimal nominal power oil generator=%.2f MW\n" % network.generators.p_nom_opt["oil"])
            handle.write("Optimal nominal power solar generator=%.2f MW\n" % network.generators.p_nom_opt["solar"])
        except:
            sys.exit(-1)   # exit unsuccessfully


    # get generator plot
    plot = network.generators_t.p.plot()


    # plot store results (if environment allows it - e.g. JupyterLab, Spyder IDE)
    #network.stores_t.p.plot()


    # save generator plot to file (in png format)
    figure = plot.get_figure()
    figure.savefig(snakemake_parameters["result_png_file"])


    sys.exit(0)   # exit successfully

