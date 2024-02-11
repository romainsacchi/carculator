"""
inventory.py contains Inventory which provides all methods to solve inventories.
"""

import warnings

import numpy as np
from carculator_utils.inventory import Inventory, format_array

from . import DATA_DIR

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

IAM_FILES_DIR = DATA_DIR / "IAM"


class InventoryCar(Inventory):
    """
    Build and solve the inventory for results
    characterization and inventory export

    """

    def fill_in_A_matrix(self):
        """
        Fill-in the A matrix. Does not return anything. Modifies in place.
        Shape of the A matrix (values, products, activities).

        :param array: :attr:`array` from :class:`CarModel` class
        """

        # Glider
        self.A[
            :,
            self.find_input_indices(("market for glider, passenger car",)),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="glider base mass") * -1
        )

        self.A[
            :,
            self.find_input_indices(("glider lightweighting",)),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="lightweighting")
            * self.array.sel(parameter="glider base mass")
            * -1
        )

        self.A[
            :,
            self.find_input_indices(("maintenance, passenger car",)),
            self.find_input_indices(("transport, car, ",)),
        ] = (
            self.array.sel(parameter="curb mass") / 1240 / 150000 * -1
        )

        # Fuel tank EoL
        self.A[
            :,
            self.find_input_indices(
                ("market for waste plastic, industrial electronics",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = self.array.sel(parameter="fuel tank mass")

        # EoL Glider
        # the EoL treatment of the glider is already
        # taken into account in the glider production dataset
        # we do not consider the reduced EoL due to lightweighting

        # Powertrain components
        self.A[
            :,
            self.find_input_indices(("market for charger, electric passenger car",)),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="charger mass") * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for converter, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="converter mass") * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for electric motor, electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="electric engine mass") * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for inverter, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="inverter mass") * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for power distribution unit, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter="power distribution unit mass") * -1
        )

        l_elec_pt = [
            "charger mass",
            "converter mass",
            "inverter mass",
            "power distribution unit mass",
            # "powertrain mass",
            "electric engine mass",
            "fuel cell stack mass",
            "fuel cell ancillary BoP mass",
            "fuel cell essential BoP mass",
        ]

        self.A[
            :,
            self.find_input_indices(
                (
                    "market for used powertrain from electric passenger car, manual dismantling",
                )
            ),
            self.find_input_indices(("Car, ",)),
        ] = self.array.sel(parameter=l_elec_pt).sum(dim="parameter")

        self.A[
            :,
            self.find_input_indices(
                ("market for internal combustion engine, passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array.sel(parameter=["combustion engine mass", "powertrain mass"]).sum(
                dim="parameter"
            )
            * -1
        )

        # EoL internal combustion engine
        # is already taken into account
        # in the ICE production dataset

        # Energy storage
        self.add_fuel_cell_stack()
        self.add_hydrogen_tank()
        self.add_battery()

        self.A[
            :,
            self.find_input_indices(
                contains=("polyethylene production, high density, granulate",)
            ),
            self.find_input_indices(
                contains=("Car, ",),
            ),
        ] = (
            self.array.sel(parameter="fuel tank mass")
            * (self.array.sel(parameter="combustion power") > 0)
            * (self.array.sel(parameter="CNG tank mass intercept") == 0)
            * -1
        )

        self.add_cng_tank()

        # END of vehicle building

        # Add vehicle dataset to transport dataset
        self.add_vehicle_to_transport_dataset()

        self.display_renewable_rate_in_mix()

        self.add_electricity_to_electric_vehicles()

        self.add_hydrogen_to_fuel_cell_vehicles()

        self.add_fuel_to_vehicles("methane", ["ICEV-g"], "EV-g")

        # Gas leakage to air
        self.A[
            :,
            self.find_input_indices(("fuel supply for methane vehicles",)),
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] *= 1 + self.array.sel(parameter="CNG pump-to-tank leakage")

        # Gas leakage to air
        self.A[
            :,
            self.inputs[("Methane, fossil", ("air",), "kilogram")],
            self.find_input_indices((f"transport, {self.vm.vehicle_type}, ",)),
        ] *= self.array.sel(parameter="CNG pump-to-tank leakage")

        self.add_fuel_to_vehicles("diesel", ["ICEV-d", "PHEV-d", "HEV-d"], "EV-d")

        self.add_fuel_to_vehicles("petrol", ["ICEV-p", "PHEV-p", "HEV-p"], "EV-p")

        self.add_abrasion_emissions()

        self.add_road_construction()

        self.add_road_maintenance()

        self.add_exhaust_emissions()

        self.add_noise_emissions()

        self.add_refrigerant_emissions()

        print("*********************************************************************")
