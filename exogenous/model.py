#!/usr/bin/env python3

# This Python module models a exogenous-driven land transport scenario (as graphically described in https://github.com/ricnogfer/pypsa-eur-sec_land_transport/blob/master/resources/exogenous_endogenous_land_transport.pdf)


# import PyPSA package
import pypsa


# create (empty) network
network = pypsa.Network()


# add bus "oil" to network
network.add("Bus", "oil")


# add generator "oil" to bus "oil"
network.add("Generator", "oil", bus = "oil")


# add load "ICE" to bus "oil"
network.add("Load", "ICE", bus = "oil")


# add bus "electricity" to network
network.add("Bus", "electricity")


# add generator "solar" to bus "electricity"
network.add("Generator", "solar", bus = "electricity")


# add load "BEV" to bus "electricity"
network.add("Load", "BEV", bus = "electricity")


# run code (if the present module is not used by another one)
if __name__ == "__main__":

	# print network content summary
	print(network)


