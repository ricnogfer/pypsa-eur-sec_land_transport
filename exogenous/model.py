#!/usr/bin/env python3



# This Python module implements an exogenous-driven land transport (toy) model involving ICE/BE-based vehicles (represented in
# https://raw.githubusercontent.com/ricnogfer/pypsa-eur-sec_land_transport/resources/docs/exogenous_land_transport.svg)



# import necessary modules
import os
import sys
import pypsa
import pandas
import matplotlib.pyplot as plt



# declare global variables
_COUNTRY = "DNK"   # country to be analised (chosen arbitrarly)
_NODE = "DK1 0"   # node within country to be analised (chosen arbitrarly)



# function to create model (i.e. PyPSA network)
def create_model(parameters, index):
    """
    Parameters
    ----------
    parameters : a Python dictionary (dict())
        Parameters needed to create and configure the model (i.e. PyPSA network) properly.
    index : a Python integer (int)
        Index representing the parameters of a certain planning horizon.

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


    # add load "ICEV" to bus "oil" with associated demand
    value = parameters["transport_demand"] * parameters["ICEV_shares"][index] / parameters["ICEV_efficiency"]
    network.add("Load",
                "ICEV",
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


    # add bus "BEV" to network
    network.add("Bus",
                "BEV")


    # add load "BEV" to bus "BEV" with associated demand
    value = parameters["transport_demand"] * parameters["BEV_shares"][index]
    network.add("Load",
                "BEV",
                bus = "BEV",
                p_set = value)


    # add link "electricity_2_BEV" which connects bus "electricity" to bus "BEV" with associated availability and efficiency
    value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][index] * parameters["BEV_charge_rate"]
    network.add("Link",
                "charge",
                bus0 = "electricity",
                bus1 = "BEV",
                p_nom = value,
                p_max_pu = parameters["availability_profile"].values,
                efficiency = parameters["BEV_charge_efficiency"])


    # enable BE-vehicle batteries to send power to the grid  (i.e. bus "electricity")
    if parameters["BEV_V2G"]:
        value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][index] * parameters["BEV_charge_rate"]
        network.add("Link",
                    "discharge",
                    bus0 = "BEV",
                    bus1 = "electricity",
                    p_nom = value,
                    p_max_pu = parameters["availability_profile"].values,
                    efficiency = parameters["BEV_charge_efficiency"])


    # enable BE-vehicle batteries based on demand-side management (DSM) profile
    if parameters["BEV_DSM"]:
        value = parameters["transport_data"]["number cars"] * parameters["BEV_shares"][index] * parameters["BEV_availability"] * parameters["BEV_energy"]
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
        parameters["technology_costs_file"] = tuple(snakemake.input["technology_costs_file"])
        parameters["solar_profile_file"] = snakemake.input["solar_profile_file"]
        parameters["transport_data_file"] = snakemake.input["transport_data_file"][0]
        parameters["transport_demand_file"] = snakemake.input["transport_demand_file"][0]
        parameters["availability_profile_file"] = snakemake.input["availability_profile_file"][0]
        parameters["DSM_profile_file"] = snakemake.input["dsm_profile_file"][0]
        parameters["planning_horizons"] = tuple(snakemake.config["scenario"]["planning_horizons"])
        parameters["ICEV_shares"] = tuple(snakemake.config["sector"]["land_transport_ice_share"].values())
        parameters["BEV_shares"] = tuple(snakemake.config["sector"]["land_transport_electric_share"].values())
        parameters["ICEV_efficiency"] = snakemake.config["sector"]["transport_internal_combustion_efficiency"]
        parameters["BEV_DSM"] = snakemake.config["sector"]["bev_dsm"]
        parameters["BEV_availability"] = snakemake.config["sector"]["bev_availability"]
        parameters["BEV_energy"] = snakemake.config["sector"]["bev_energy"]
        parameters["BEV_charge_rate"] = snakemake.config["sector"]["bev_charge_rate"]
        parameters["BEV_charge_efficiency"] = snakemake.config["sector"]["bev_charge_efficiency"]
        parameters["BEV_V2G"] = snakemake.config["sector"]["v2g"]
        parameters["ICEV_consumption_per_unit"] = snakemake.config["sector"]["icev_consumption_per_unit"]
        parameters["BEV_consumption_per_unit"] = snakemake.config["sector"]["bev_consumption_per_unit"]
        parameters["technology_costs"] = snakemake.config["costs"]
        parameters["summary_file"] = snakemake.output["summary_file"]
        parameters["generators_file"] = tuple(snakemake.output["generators_file"])
        parameters["stores_file"] = tuple(snakemake.output["stores_file"])
        parameters["network_file"] = tuple(snakemake.output["network_file"])
        parameters["vehicle_file"] = snakemake.output["vehicle_file"]
        parameters["snapshots"] = pandas.date_range("%sT00:00Z" % snakemake.config["snapshots"]["start"], "%sT23:00Z" % snakemake.config["snapshots"]["end"], freq = "H")
    else:   # through terminal
        resource_path = "../resources/data"
        result_path = "results"
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        parameters["technology_costs_file"] = ("%s/costs_2025.csv" % resource_path, )
        parameters["solar_profile_file"] = "%s/solar_profile_1979_2017.csv" % resource_path
        parameters["transport_data_file"] = "%s/transport_data_s_37.csv" % resource_path
        parameters["transport_demand_file"] = "%s/transport_demand_s_37.csv" % resource_path
        parameters["availability_profile_file"] = "%s/avail_profile_s_37.csv" % resource_path
        parameters["DSM_profile_file"] = "%s/dsm_profile_s_37.csv" % resource_path
        parameters["planning_horizons"] = (2025, )
        parameters["ICEV_shares"] = (0.9, )
        parameters["BEV_shares"] = (0.1, )
        parameters["ICEV_efficiency"] = 0.3
        parameters["BEV_DSM"] = True
        parameters["BEV_availability"] = 0.5
        parameters["BEV_energy"] = 0.05
        parameters["BEV_charge_efficiency"] = 0.9
        parameters["BEV_charge_rate"] = 0.011
        parameters["BEV_V2G"] = True
        parameters["ICEV_consumption_per_unit"] = 0.05   # arbitrary value
        parameters["BEV_consumption_per_unit"] = 0.01   # arbitrary value
        parameters["technology_costs"] = {"year": 2030, "version": "v0.5.0", "rooftop_share": 0.14, "fill_values": {"FOM": 0, "VOM": 0, "efficiency": 1, "fuel": 0, "investment": 0, "lifetime": 25, "CO2 intensity": 0, "discount rate": 0.07}, "marginal_cost": {"solar": 0.01, "onwind": 0.015, "offwind": 0.015, "hydro": 0.0, "H2": 0.0, "electrolysis": 0.0, "fuel cell": 0.0, "battery": 0.0, "battery inverter": 0.0}, "emission_prices": {"co2": 0.0}}
        parameters["summary_file"] = "%s/summary.txt" % result_path
        parameters["generators_file"] = ("%s/generators_2025.png" % result_path, )
        parameters["stores_file"] = ("%s/stores_2025.png" % result_path, )
        parameters["network_file"] = ("%s/network_2025.nc" % result_path, )
        parameters["vehicle_file"] = "%s/vehicle.png" % result_path
        parameters["snapshots"] = pandas.date_range("2013-01-01T00:00Z", "2013-12-31T23:00Z", freq = "H")


    return parameters



# function to get parameters that the model (i.e. PyPSA network) will be based upon
def get_model_parameters(snakemake_parameters, index):
    """
    Parameters
    ----------
    snakemake_parameters : a Python dictionary (dict())
        Contains the Snakemake parameters needed to create the model (i.e. PyPSA network) parameters.
    index : a Python integer (int)
        Index representing the parameters of a certain planning horizon.

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
        a Python integer (int)
            Contains the annuity value.
        """

        if discount_rate > 0:
            return discount_rate / (1.0 - 1.0 / (1.0 + discount_rate) ** period)

        return 1 / period



    parameters = dict()


    # read technology costs (CSV) file and prepare (i.e. process) these for proper consumption afterwards
    technology_costs = prepare_costs(snakemake_parameters["technology_costs_file"][index], snakemake_parameters["technology_costs"], 1)


    # get oil generator costs
    oil = technology_costs.loc["oil"]
    parameters["oil_capital_cost"] = 0   # set to 0 since oil is imported to Europe (i.e. no relevant capital cost associated with oil generators)
    parameters["oil_marginal_cost"] = oil.at["fuel"]


    # get solar generator costs
    solar = technology_costs.loc["solar"]
    solar_investment = solar.at["investment"]
    solar_FOM = calculate_annuity(solar.at["lifetime"], 0.07) + solar.at["FOM"] / 100.0   # in percentage
    parameters["solar_capital_cost"] = solar_investment + solar_investment * solar_FOM
    parameters["solar_marginal_cost"] = solar.at["VOM"]


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
    parameters["ICEV_shares"] = snakemake_parameters["ICEV_shares"]
    parameters["BEV_shares"] = snakemake_parameters["BEV_shares"]
    parameters["ICEV_efficiency"] = snakemake_parameters["ICEV_efficiency"]
    parameters["BEV_DSM"] = snakemake_parameters["BEV_DSM"]
    parameters["BEV_availability"] = snakemake_parameters["BEV_availability"]
    parameters["BEV_energy"] = snakemake_parameters["BEV_energy"]
    parameters["BEV_charge_rate"] = snakemake_parameters["BEV_charge_rate"]
    parameters["BEV_charge_efficiency"] = snakemake_parameters["BEV_charge_efficiency"]
    parameters["BEV_V2G"] = snakemake_parameters["BEV_V2G"]


    return parameters



def solve_model(snakemake_parameters):
    """
    Parameters
    ----------
    snakemake_parameters : a Python dictionary (dict())
        Contains the Snakemake parameters needed to create the model (i.e. PyPSA network) parameters.

    Returns
    -------
    a Python integer (int)
        Contains the status of solving the model (0 for successful; not 0 for unsuccessful).
    """

    # declare lists to hold ICE/BE-based vehicle units throughout planning horizons
    ICEV_units = list()
    BEV_units = list()


    # open file to write results
    try:
        handle = open(snakemake_parameters["summary_file"], "w")
    except:
        print("Error when creating/opening file '%s'" % snakemake_parameters["summary_file"])
        return -1   # return unsuccessfully


    # loop through planning horizons
    for i in range(len(snakemake_parameters["planning_horizons"])):

        # get model parameters
        model_parameters = get_model_parameters(snakemake_parameters, i)


        # create model (i.e. PyPSA network) based on the parameters of the planning horizon
        network = create_model(model_parameters, i)


        # solve network using Gurobi
        status = network.lopf(pyomo = False, solver_name = "gurobi")
        if status[0] != "ok":
            print("The solver did not reach a solution (unfeasible or unbounded)!")
            handle.close()
            return -1   # return unsuccessfully


        # save network (in netcdf format)
        network.export_to_netcdf(snakemake_parameters["network_file"][i])


        # calculate ICE/BE-based vehicle units of the planning horizon
        ICEV_units.append(network.loads_t.p_set["ICEV"].max() / snakemake_parameters["ICEV_consumption_per_unit"])   # TODO: check if logic is correct
        BEV_units.append(network.loads_t.p_set["BEV"].max() / snakemake_parameters["BEV_consumption_per_unit"])   # TODO: check if logic is correct


        # print results
        print("Year of the almighty Lord %d AC:" % snakemake_parameters["planning_horizons"][i])
        print("   Objective value=%.2f M€" % (network.objective / 10**6))
        print("   Optimal nominal power of oil generator=%.2f MW" % network.generators.p_nom_opt["oil"])
        print("   Optimal nominal power of solar generator=%.2f MW" % network.generators.p_nom_opt["solar"])
        print("   Optimal nominal energy of battery=%.2f MWh" % network.stores.e_nom_opt["battery"])
        print("   Number of ICEV=%d units" % ICEV_units[-1])
        print("   Number of BEV=%d units" % BEV_units[-1])


        # write results to file
        try:
            handle.write("Year of the almighty Lord %d AC:\n" % snakemake_parameters["planning_horizons"][i])
            handle.write("   Objective value=%.2f M€\n" % (network.objective / 10**6))
            handle.write("   Optimal nominal power of oil generator=%.2f MW\n" % network.generators.p_nom_opt["oil"])
            handle.write("   Optimal nominal power of solar generator=%.2f MW\n" % network.generators.p_nom_opt["solar"])
            handle.write("   Optimal nominal energy of battery=%.2f MWh\n" % network.stores.e_nom_opt["battery"])
            handle.write("   Number of ICEV=%d units\n" % ICEV_units[-1])
            handle.write("   Number of BEV=%d units\n" % BEV_units[-1])
            handle.write("\n")
        except:
            print("Error when writing to file '%s'" % snakemake_parameters["summary_file"])
            handle.close()
            return -1   # return unsuccessfully


        # save generators plot to file (in png format)
        plot = network.generators_t.p.plot()
        figure = plot.get_figure()
        figure.savefig(snakemake_parameters["generators_file"][i])


        # save stores plot to file (in png format)
        plot = network.stores_t.p.plot()
        figure = plot.get_figure()
        figure.savefig(snakemake_parameters["stores_file"][i])


    # close file handle
    handle.close()


    # save vehicle plot to file (in png format)
    figure, axes = plt.subplots()
    axes.set_title("Vehicle Usage")
    axes.set_xlabel("Year")
    axes.set_ylabel("Unit")
    axes.set_xlim(2015, 2055)
    axes.plot(snakemake_parameters["planning_horizons"], ICEV_units, "o-", label = "ICEV")
    axes.plot(snakemake_parameters["planning_horizons"], BEV_units, "x-", label = "BEV")
    axes.legend()
    figure.savefig(snakemake_parameters["vehicle_file"])


    return 0   # return successfully



def prepare_costs(cost_file, config, nyears):
    """
    Parameters
    ----------
    cost_file : a Python string (str)
        Lorem Ipsum.
    config : a Python dictionary (dict())
        Lorem Ipsum.
    nyears : a Python integer (int)
        Lorem Ipsum.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """


    def annuity(n,r):
        """ Calculate the annuity factor
        for an asset with lifetime n years and
        discount rate of f, e.g. annuity(20,0.05)*20 = 1.6 """

        if r > 0:
            return r/(1. - 1./(1.+r)**n)
        else:
            return 1/n


    def annuity_factor(v):

        return annuity(v["lifetime"], v["discount rate"]) + v["FOM"] / 100.0


    # set all asset costs and other parameters
    costs = pandas.read_csv(cost_file, index_col = [0, 1]).sort_index()


    # correct units to MW and EUR
    costs.loc[costs.unit.str.contains("/kW"), "value"] *= 1e3

    # min_count=1 is important to generate NaNs which are then filled by fillna
    costs = (costs.loc[:, "value"].unstack(level = 1).groupby("technology").sum(min_count = 1))


    costs = costs.fillna(config["fill_values"])


    costs["fixed"] = [annuity_factor(v) * v["investment"] * nyears for i, v in costs.iterrows()]


    return costs



# run code (if present module is not imported by another one)
if __name__ == "__main__":


    # get Snakemake parameters
    snakemake_parameters = get_snakemake_parameters()


    # check that Snakemake parameters "planning_horizons", "ICEV_shares" and "BEV_shares" are properly specified
    if not snakemake_parameters["planning_horizons"] or not snakemake_parameters["ICEV_shares"] or not snakemake_parameters["BEV_shares"]:
        print("The Snakemake parameters 'planning_horizons', 'ICEV_shares' and/or 'BEV_shares' is/are not specified!")
        sys.exit(-1)   # exit unsuccessfully
    if len(snakemake_parameters["planning_horizons"]) != len(snakemake_parameters["ICEV_shares"]) != len(snakemake_parameters["BEV_shares"]):
        print("The number of Snakemake parameters in 'planning_horizons', 'ICEV_shares' and 'BEV_shares' do not match!")
        sys.exit(-1)   # exit unsuccessfully


    # check that ICE/BE-based vehicle shares are properly specified
    for i in range(len(snakemake_parameters["ICEV_shares"])):
        if snakemake_parameters["ICEV_shares"][i] + snakemake_parameters["BEV_shares"][i] != 1:
            print("The sum of ICE vehicle share with BE vehicle share is not equal to 1!")
            sys.exit(-1)   # exit unsuccessfully


    # solve model
    status = solve_model(snakemake_parameters)


    # exit with returned status value
    sys.exit(status)

