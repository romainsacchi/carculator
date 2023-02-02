"""
inventory.py contains Inventory which provides all methods to solve inventories.
"""

import numpy as np
from carculator_utils.inventory import Inventory

from . import DATA_DIR

np.warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)

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
        ] = (self.array[self.array_inputs["glider base mass"], :]) * -1

        self.A[
            :,
            self.find_input_indices(("Glider lightweighting",)),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["lightweighting"]]
            * self.array[self.array_inputs["glider base mass"]]
        ) * -1

        self.A[
            :,
            self.find_input_indices(("maintenance, passenger car",)),
            self.find_input_indices(("transport, car, ",)),
        ] = (
            self.array[self.array_inputs["curb mass"]] / 1240 / 150000 * -1
        )

        # Glider EoL + fuel tank
        self.A[
            :,
            self.find_input_indices(
                ("treatment of used glider, passenger car, shredding",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["glider base mass"]]
            * (1 - self.array[self.array_inputs["lightweighting"]])
        ) + self.array[
            self.array_inputs["fuel tank mass"]
        ]

        # Combustion engine EoL
        self.A[
            :,
            self.find_input_indices(
                (
                    "treatment of used internal combustion engine, passenger car, shredding",
                )
            ),
            self.find_input_indices(("Car, ",)),
        ] = self.array[
            [
                self.array_inputs[l]
                for l in ["combustion engine mass", "powertrain mass"]
            ],
            :,
        ].sum(
            axis=0
        )

        # Powertrain components
        self.A[
            :,
            self.find_input_indices(("market for charger, electric passenger car",)),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["charger mass"], :] * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for converter, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["converter mass"], :] * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for electric motor, electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["electric engine mass"], :] * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for inverter, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["inverter mass"], :] * -1
        )

        self.A[
            :,
            self.find_input_indices(
                ("market for power distribution unit, for electric passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[self.array_inputs["power distribution unit mass"], :] * -1
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
        ] = self.array[[self.array_inputs[l] for l in l_elec_pt], :].sum(axis=0)

        self.A[
            :,
            self.find_input_indices(
                ("market for internal combustion engine, passenger car",)
            ),
            self.find_input_indices(("Car, ",)),
        ] = (
            self.array[
                [
                    self.array_inputs[l]
                    for l in ["combustion engine mass", "powertrain mass"]
                ],
                :,
            ].sum(axis=0)
        ) * -1

        # Energy storage
        self.add_fuel_cell_stack()
        self.add_hydrogen_tank()
        self.add_battery()

        index = self.get_index_vehicle_from_array(
            ["ICEV-p", "ICEV-d", "HEV-p", "HEV-d", "PHEV-p", "PHEV-d"]
        )

        self.A[
            :,
            self.find_input_indices(
                contains=("polyethylene production, high density, granulate",)
            ),
            self.find_input_indices(
                contains=("Car, ",), excludes=("BEV", "ICEV-g", "FCEV")
            ),
        ] = (
            self.array[self.array_inputs["fuel tank mass"], :, index] * -1
        )

        self.add_cng_tank()

        # END of vehicle building

        # Add vehicle dataset to transport dataset
        self.add_vehicle_to_transport_dataset()

        self.display_renewable_rate_in_mix()

        self.add_electricity_to_electric_vehicles()

        self.add_hydrogen_to_fuel_cell_vehicles()

        self.add_fuel_to_vehicles("cng", ["ICEV-g"], "EV-g")

        for year in self.scope["year"]:
            cng_idx = self.get_index_vehicle_from_array(
                [
                    "ICEV-g",
                ],
                [
                    year,
                ],
                method="and",
            )

            self.A[
                :,
                self.find_input_indices(("fuel supply for cng vehicles", str(year))),
                self.find_input_indices(
                    (f"transport, {self.vm.vehicle_type}, ", "ICEV-g", str(year))
                ),
            ] *= (
                1
                + self.array[self.array_inputs["CNG pump-to-tank leakage"], :, cng_idx]
            )

        # Gas leakage to air
        cng_idx = self.get_index_vehicle_from_array(
            [
                "ICEV-g",
            ]
        )
        self.A[
            :,
            self.inputs[("Methane, fossil", ("air",), "kilogram")],
            self.find_input_indices(
                (
                    f"transport, {self.vm.vehicle_type}, ",
                    "ICEV-g",
                )
            ),
        ] *= self.array[self.array_inputs["CNG pump-to-tank leakage"], :, cng_idx]

        self.add_fuel_to_vehicles("diesel", ["ICEV-d", "PHEV-d", "HEV-d"], "EV-d")

        self.add_fuel_to_vehicles("petrol", ["ICEV-p", "PHEV-p", "HEV-p"], "EV-p")

        self.add_abrasion_emissions()

        self.add_road_construction()

        self.add_road_maintenance()

        self.add_exhaust_emissions()

        self.add_noise_emissions()

        self.add_refrigerant_emissions()

        print("*********************************************************************")
