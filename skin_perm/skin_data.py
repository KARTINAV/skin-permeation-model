"""
SkinProp.java, EnvOpt.java, DosageDat.java, VehicleDat.java,
OutputParameters.java, DosageScenario.java → Python dataclasses.
"""
from dataclasses import dataclass, field
from typing import List
from .constants import *


@dataclass
class SkinProp:
    hydration: str = PARTIALLY_HYDRATED
    sc_thickness: float = DEF_SC_THICKNESS_PARTIALLY_HYDRATED
    ve_thickness: float = DEF_VE_THICKNESS
    dermis_thickness: float = DEF_DERMIS_IN_VIVO
    sc_ph: float = DEF_SC_PH
    ve_ph: float = DEF_VE_PH
    dermis_ph: float = DEF_DERMIS_PH
    in_vivo_vitro: int = IN_VIVO
    is_sc: bool = True
    is_ve: bool = True
    is_dermis: bool = True
    species: str = "human"

    def set_hydration(self, s: str):
        if s == FULLY_HYDRATED:
            self.hydration = FULLY_HYDRATED
            self.sc_thickness = DEF_SC_THICKNESS_FULLY_HYDRATED
        elif s == PARTIALLY_HYDRATED:
            self.hydration = PARTIALLY_HYDRATED
            self.sc_thickness = DEF_SC_THICKNESS_PARTIALLY_HYDRATED

    def set_in_vitro_vivo(self, v: int):
        self.in_vivo_vitro = v
        self.dermis_thickness = DEF_DERMIS_IN_VIVO if v == IN_VIVO else DEF_DERMIS_IN_VITRO


@dataclass
class EnvOpt:
    temperature: float = DEF_SURFACE_TEMP
    wind_speed: float = DEF_WIND_INDOOR


@dataclass
class DosageScenario:
    permeant_dose: float = -1.0
    nvv_dose: float = 0.0
    time_v: float = 0.0
    removal_scenario: int = 0  # 0=none


@dataclass
class DosageDat:
    permeant_amount_applied: float = 0.0
    nvv_amount_applied: float = 0.0
    vv_amount_applied: float = 0.0
    area: float = 0.0
    is_permeant_aa: bool = False
    is_nvv_aa: bool = False
    is_vv_aa: bool = False
    is_area: bool = False
    is_terminal_wash: bool = False
    is_multiple_doses: bool = False
    scenarios: List[DosageScenario] = field(default_factory=lambda: [DosageScenario()])

    def set_permeant_amount(self, v: float):
        self.permeant_amount_applied = v
        self.is_permeant_aa = True
        if self.scenarios:
            self.scenarios[0].permeant_dose = v

    def set_area(self, v: float):
        self.area = v; self.is_area = True


@dataclass
class VehicleDat:
    volatile_type: int = NON_VOLATILE_VEHICLE
    name: str = ""
    is_present: bool = False
    vehicle_ph: float = DEF_VEHICLE_PH
    density: float = 0.0
    molecular_weight: float = 0.0
    _is_name: bool = False
    _is_ph: bool = True
    _is_mw: bool = False
    _is_density: bool = False
    is_water: bool = False

    def is_setup(self) -> bool:
        if not self.is_present:
            return False
        if not self._is_name or not self._is_ph:
            return False
        if self.volatile_type == NON_VOLATILE_VEHICLE and (not self._is_mw or not self._is_density):
            return False
        return True


@dataclass
class OutputParameters:
    max_duration: float = 1000.0
    set_concentration: float = 0.0
    min_output_step_size: float = 10.0
    max_step_size: float = 5.0
    output_flag: int = 1
    is_flux_plot: bool = False
    is_cumulative_plot: bool = False
    is_conc_profiles: bool = False
    is_systemic: bool = False
