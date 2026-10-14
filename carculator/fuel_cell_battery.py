"""Matched passenger-car FCEV pack reference; charging remains uncalibrated."""

from dataclasses import dataclass

import numpy as np

PROFILE_SOURCE = (
    "Toyota second-generation Mirai, manufacturer announcement of 30 November "
    "2020: 310.8 V, 4 Ah (1.2432 kWh), complete pack 44.6 kg, 31.5 kW "
    "discharge for 10 seconds. "
    "https://newsroom.toyota.eu/introducing-the-all-new-toyota-mirai/"
)


@dataclass(frozen=True)
class FuelCellBatteryProfile:
    """Pack-level proxy: nominal kWh, complete kg and 10-second terminal kW.

    Applied from 2020 onward without an assumed technology improvement. The
    existing cell-mass fraction allocates the pack internally; it is not a
    measured Mirai cell fraction. Chemistry remains a manufacturing/cost proxy.
    No charging, continuous-output, SOC, thermal or lifetime rating is inferred.
    """

    nominal_kwh: float = 1.2432
    pack_kg: float = 44.6
    discharge_10s_kw: float = 31.5

    def validate(self):
        for name, value in vars(self).items():
            if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
                raise ValueError(f"FCEV battery {name} must be finite and positive.")
