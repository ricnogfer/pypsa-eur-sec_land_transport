#!/usr/bin/env python3


# This Python module implements an exogenous-driven land transport (toy) model involving ICE/BE-based vehicles (represented in
# https://github.com/ricnogfer/pypsa-eur-sec_land_transport/blob/master/resources/docs/exogenous_endogenous_land_transport.pdf)



# TODO:
# 1) fix solar profile
# 2) implement CO2 constraint
# 3) enable BE-vehicles to behave like stores (i.e. to send back electricity to the bus "electricity")



# import necessary packages
import os
import sys
import pypsa
import pandas



# function to create PyPSA network
def create_network(parameters, snapshots):
    """
    Parameters
    ----------
    parameters : a Python dictionary (dict())
        Parameters needed to create the PyPSA network.
    snapshots : a Pandas date/time series (pandas.core.indexes.datetimes.DatetimeIndex())
        Number of snapshots that the PyPSA network has/will solve.

    Returns
    -------
    network : a PyPSA network object (pypsa.Network())
        The created PyPSA network.
    """

    # create (empty) PyPSA network
    network = pypsa.Network()
    network.set_snapshots(snapshots)


    # add bus "oil" to network
    network.add("Bus",
                "oil")


    # add generator "oil" to bus "oil" with associated costs
    network.add("Generator",
                "oil",
                bus = "oil",
                p_nom_extendable = True,
                capital_cost = parameters["oil_capital_cost"],
                marginal_cost = parameters["oil_marginal_cost"])


    # add load "ICE" to bus "oil" with associated ICE-vehicle demand
    network.add("Load",
                "ICE",
                bus = "oil",
                p_set = parameters["ICE_vehicle_demand"])


    # add bus "electricity" to network
    network.add("Bus",
                "electricity")


    # add generator "solar" to bus "electricity" with associated costs
    network.add("Generator",
                "solar",
                bus = "electricity",
                p_nom_extendable = True,
                #p_max_pu = parameters["solar_profile"],
                capital_cost = parameters["solar_capital_cost"],
                marginal_cost = parameters["solar_marginal_cost"])


    # add load "BEV" to bus "electricity" with associated BE-vehicle demand
    network.add("Load",
                "BEV",
                bus = "electricity",
                p_set = parameters["BE_vehicle_demand"])


    return network



# function to get parameters that the PyPSA network will be based upon
def get_parameters(technology_costs_file, solar_profile_file, transport_demand_file, ICE_vehicle_shares, BE_vehicle_shares, snapshots):
    """   
    Parameters
    ----------
    technology_costs_file : a Python string (str)
        Contains the name of the (CSV) file with technology costs.
    solar_profile_file : a Python string (str)
        Contains the name of the (CSV) file with solar profiles.
    transport_demand_file : a Python string (str)
        Contains the name of the (CSV) file with transport demand data.
    ICE_vehicle_shares : a Python tuple (tuple())
        Contains the shares of ICE-vehicles (e.g. (1.0, 0.7, 0.25, 0.0)).
    BE_vehicle_shares : a Python tuple (tuple)
        Contains the shares of BE-vehicles (e.g. (0.0, 0.3, 0.75, 1.0)).
    snapshots : a Pandas date/time series (pandas.core.indexes.datetimes.DatetimeIndex())
        Number of snapshots that the PyPSA network has/will solve.

    Returns
    -------
    parameters : a Python dictionary (dict())
        Contains the parameters that the PyPSA network will be based upon.
    """

    parameters = dict()


    # read technology costs (CSV) file
    technology_costs = pandas.read_csv(technology_costs_file, index_col = "technology")


    # get oil generator costs
    oil = technology_costs.loc["oil"].set_index("parameter")
    parameters["oil_capital_cost"] = 0   # set to 0 since oil is imported to Europe
    parameters["oil_marginal_cost"] = oil.at["VOM", "value"]


    # get solar generator costs
    solar = technology_costs.loc["solar"].set_index("parameter")
    parameters["solar_capital_cost"] = ((solar.at["investment", "value"] * 1000.0) / solar.at["lifetime", "value"]) * (1 + solar.at["FOM", "value"] / 100.0)   # TODO: check if the formula is correct
    parameters["solar_marginal_cost"] = solar.at["VOM", "value"]


    # read solar profile from CSV file
    solar_profile = pandas.read_csv(solar_profile_file, sep = ';', index_col = 0)
    solar_profile.index = pandas.to_datetime(solar_profile.index)
    parameters["solar_profile"] = solar_profile["DNK"]   # solar profile for Denmark (chosen arbitrarly)


    # read transport demand (CSV) file
    transport_demand = pandas.read_csv(transport_demand_file, index_col = 0)


    # get transport demand for "DK1 0" (chosen arbitrarly)
    transport_demand.index = pandas.to_datetime(transport_demand.index, utc = True)
    parameters["transport_demand"] = transport_demand["DK1 0"]


    # get ICE/BE-vehicle demand based on their demand
    # TODO: only the first specified share/year is considered for now; in a subsequent iteration, all specified shares/years will be considered
    parameters["ICE_vehicle_demand"] = parameters["transport_demand"] * ICE_vehicle_shares[0]
    parameters["BE_vehicle_demand"] = parameters["transport_demand"] * BE_vehicle_shares[0]


    return parameters



# run code (if present module is not imported by another one)
if __name__ == "__main__":

    # check how present module was launched
    if "snakemake" in globals():   # through Snakemake
        technology_costs_file = snakemake.input["technology_costs_file"]
        solar_profile_file = snakemake.input["solar_profile_file"]
        transport_demand_file = snakemake.input["transport_demand_file"][0]
        ICE_vehicle_shares = tuple(snakemake.config["sector"]["land_transport_ice_share"].values())
        BE_vehicle_shares = tuple(snakemake.config["sector"]["land_transport_electric_share"].values())
        result_txt_file = snakemake.output["result_txt_file"]
        result_png_file = snakemake.output["result_png_file"]
        snapshots = pandas.date_range("%sT00:00Z" % snakemake.config["snapshots"]["start"], "%sT23:00Z" % snakemake.config["snapshots"]["end"], freq = "H")
    else:   # through terminal
        resource_path = "../resources/technology_data"
        result_path = "results"
        if not os.path.exists(result_path):
            os.makedirs(result_path)
        technology_costs_file = "%s/costs_2025.csv" % resource_path
        solar_profile_file = "%s/solar_profile_1979_2017.csv" % resource_path
        transport_demand_file = "%s/transport_demand_s_37.csv" % resource_path
        ICE_vehicle_shares = (0.9, 0.75, 0.4, 0.0)
        BE_vehicle_shares = (0.1, 0.25, 0.6, 1.0)
        result_txt_file = "%s/%s" % (result_path, "results.txt")
        result_png_file = "%s/%s" % (result_path, "results.png")
        snapshots = pandas.date_range("2013-01-01T00:00Z", "2013-12-31T23:00Z", freq = "H")


    # check that ICE/BE-vehicle shares are properly specified
    if not ICE_vehicle_shares or not BE_vehicle_shares:
        print("The ICE vehicle shares and/or BE vehicle shares is/are not specified!")
        sys.exit(-1)   # exit unsuccessfully
    if len(ICE_vehicle_shares) != len(BE_vehicle_shares):
        print("The number of ICE vehicle shares does not match the number of BE vehicle shares!")
        sys.exit(-1)   # exit unsuccessfully
    for i in range(len(ICE_vehicle_shares)):
        if ICE_vehicle_shares[i] + BE_vehicle_shares[i] != 1:
            print("The sum of ICE vehicle share with BE vehicle share is not equal to 1!")
            sys.exit(-1)   # exit unsuccessfully


    # get (model) parameters
    parameters = get_parameters(technology_costs_file, solar_profile_file, transport_demand_file, ICE_vehicle_shares, BE_vehicle_shares, snapshots)


    # create network based on (model) parameters with snapshots information
    network = create_network(parameters, snapshots)


    # solve network using Gurobi
    status = network.lopf(pyomo = False, solver_name = "gurobi")
    if status[0] != "ok":
        print("The solver did not reach a solution (infeasible or unbounded)!")
        sys.exit(-1)   # exit unsuccessfully


    # print results
    # TODO: print the units too
    print("Objective value=%f" % network.objective)
    print("Optimal nominal power oil generator=%f" % network.generators.p_nom_opt["oil"])
    print("Optimal nominal power solar generator=%f" % network.generators.p_nom_opt["solar"])


    # write results to file
    with open(result_txt_file, "w") as handle:
        try:
            # TODO: write the units too
            handle.write("Objective value=%f\n" % network.objective)
            handle.write("Optimal nominal power oil generator=%f\n" % network.generators.p_nom_opt["oil"])
            handle.write("Optimal nominal power solar generator=%f\n" % network.generators.p_nom_opt["solar"])
        except Exception as e:
            sys.exit(-1)   # exit unsuccessfully


    # plot results (only if environment allows it - e.g. in JupyterLab, Spyder IDE)
    plot = network.generators_t.p.plot()


    # save plot to file (in png format)
    figure = plot.get_figure()
    figure.savefig(result_png_file)


    sys.exit(0)   # exit successfully

