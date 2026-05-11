"""
ChemFormula.java + Element.java + ChemDat.java → Python data classes.
Faithful line-by-line translation from decompiled Java.
"""
import re
import math
from dataclasses import dataclass, field
from .constants import (ATOMIC_NAMES, ATOMIC_MASSES, VA_COEFFICIENTS,
                         MINUS_BIG, MIN_MASS, MAX_MASS, MIN_LOGKOW, MAX_LOGKOW,
                         MIN_DENSITY, MAX_DENSITY, MAX_VAPOUR_PRESSURE,
                         MAX_WATER_SOLUBILITY, MIN_PKA, MAX_PKA,
                         MAX_PERMEABILITY_COEFFICIENT, MIN_FRACTIONS, MAX_FRACTIONS,
                         MIN_TEMP, MAX_TEMP, MIN_WATER_TEMP, MAX_WATER_TEMP,
                         GRAIN_OTHER)


class Element:
    """Single atomic element with mass and Schroeder Va coefficient."""
    def __init__(self, symbol: str, count: int):
        self.atom = symbol
        self.natom = count
        self.e_mass = 0.0
        self.va_coeff = 0.0
        self.is_valid = False
        self.is_other = False
        for i, name in enumerate(ATOMIC_NAMES):
            if self.atom == name:
                self.e_mass = ATOMIC_MASSES[i]
                self.va_coeff = VA_COEFFICIENTS[i]
                if self.va_coeff == 0.0:
                    self.is_other = True
                self.is_valid = True
                break
        self.t_mass = self.natom * self.e_mass if self.is_valid else 0.0
        self.total_va = self.natom * self.va_coeff if self.is_valid else 0.0


class ChemFormula:
    """Parse a molecular formula string like 'C19H28O2' into elements."""
    def __init__(self, formula_str: str = ""):
        self.chemform = ""
        self.mol_mass = 0.0
        self.va_factor = 0.0
        self.elements = []
        self.other_atoms = False
        self.is_valid = False
        if formula_str:
            self.set_formula(formula_str)

    def set_formula(self, s: str) -> bool:
        self.mol_mass = 0.0
        self.va_factor = 0.0
        self.elements = []
        self.other_atoms = False
        try:
            sf = s.strip()
            tokens = re.findall(r'[A-Z][a-z]*\d*', sf)
            self.is_valid = True
            for token in tokens:
                m_el = re.match(r'[A-Z][a-z]*', token)
                m_num = re.search(r'\d+', token)
                el_name = m_el.group() if m_el else ""
                el_count = int(m_num.group()) if m_num else 1
                ae = Element(el_name, el_count)
                self.mol_mass += ae.t_mass
                self.va_factor += ae.total_va
                if ae.is_other:
                    self.other_atoms = True
                if not ae.is_valid:
                    self.is_valid = False
                    return False
                self.elements.append(ae)
            self.chemform = sf
            return True
        except Exception:
            self.is_valid = False
            return False

    def get_total_mass(self) -> float:
        return self.mol_mass

    def get_total_atomic_va(self) -> float:
        return self.va_factor

    def has_other(self) -> bool:
        return self.other_atoms


@dataclass
class ChemDat:
    """Chemical property data — the INPUT interface of the model.
    Mirrors every field in ChemDat.java."""
    # Identification
    chem_name: str = ""
    cas: str = ""
    chem_id: str = ""
    smiles: str = ""
    formula: ChemFormula = field(default_factory=ChemFormula)

    # Required physical properties
    mw: float = 0.0
    logP: float = MINUS_BIG
    melting_point: float = MINUS_BIG
    boiling_point: float = MINUS_BIG
    vapour_pressure: float = 0.0  # stored in Pa internally

    # Optional properties
    density: float = 0.0
    density_temperature: float = MINUS_BIG
    water_solubility: float = 0.0  # stored in g/L (matching Java internally)
    water_solubility_temperature: float = MINUS_BIG
    permeability_coefficient: float = 0.0
    pka_HA: float = 0.0
    pka_BHplus: float = 0.0

    # Non-ionic / unbound fractions (user overrides)
    f_non_vt: float = 0.0
    f_u_vt: float = 0.0
    f_non_veh: float = 0.0
    f_u_sc: float = 0.0
    f_non_sc: float = 0.0
    non_volatile_vehicle_solubility: float = 0.0

    # Structural descriptors
    num_double_bonds: int = 0
    num_triple_bonds: int = 0
    num_rings: int = 0
    grain_class: int = GRAIN_OTHER
    has_pharmacophore: bool = False
    K_v_w: float = MINUS_BIG

    # ── Flags tracking which fields have been set ──
    _is_logkow: bool = field(default=False, repr=False)
    _is_mw: bool = field(default=False, repr=False)
    _is_mp: bool = field(default=False, repr=False)
    _is_bp: bool = field(default=False, repr=False)
    _is_vp: bool = field(default=False, repr=False)
    _is_density: bool = field(default=False, repr=False)
    _is_density_temp: bool = field(default=False, repr=False)
    _is_ws: bool = field(default=False, repr=False)
    _is_ws_temp: bool = field(default=False, repr=False)
    _is_kp: bool = field(default=False, repr=False)
    _is_pka_HA: bool = field(default=False, repr=False)
    _is_pka_BHplus: bool = field(default=False, repr=False)
    _is_formula: bool = field(default=False, repr=False)
    _is_f_non_vt: bool = field(default=False, repr=False)
    _is_f_u_vt: bool = field(default=False, repr=False)
    _is_f_non_veh: bool = field(default=False, repr=False)
    _is_f_non_sc: bool = field(default=False, repr=False)
    _is_f_u_sc: bool = field(default=False, repr=False)
    _is_nvv_sol: bool = field(default=False, repr=False)
    _is_kvw: bool = field(default=False, repr=False)

    # ── Convenience setters that also validate ──
    def set_formula(self, s: str) -> bool:
        self.formula = ChemFormula(s)
        if self.formula.is_valid:
            ok = self.set_mw(self.formula.get_total_mass())
            self._is_formula = ok
            return ok
        self._is_formula = False
        return False

    def set_mw(self, m: float) -> bool:
        if MIN_MASS < m < MAX_MASS:
            self.mw = m; self._is_mw = True; return True
        self._is_mw = False; return False

    def set_logkow(self, b: float) -> bool:
        if MIN_LOGKOW <= b <= MAX_LOGKOW:
            self.logP = b; self._is_logkow = True; return True
        self._is_logkow = False; return False

    def set_melting_point(self, t: float) -> bool:
        if MIN_TEMP < t < MAX_TEMP:
            self.melting_point = t; self._is_mp = True; return True
        return False

    def set_boiling_point(self, t: float) -> bool:
        if MIN_TEMP < t < MAX_TEMP:
            self.boiling_point = t; self._is_bp = True; return True
        return False

    def set_vapour_pressure(self, val: float, scale: float = 1.0) -> bool:
        """val in user units, scale converts to Pa. Java stores in Pa."""
        press = val * scale
        if 0 <= press < MAX_VAPOUR_PRESSURE:
            self.vapour_pressure = press; self._is_vp = True; return True
        self._is_vp = False; return False

    def set_density(self, d: float) -> bool:
        if MIN_DENSITY < d < MAX_DENSITY:
            self.density = d; self._is_density = True; return True
        return False

    def set_density_temperature(self, t: float) -> bool:
        if MIN_TEMP < t < MAX_TEMP:
            self.density_temperature = t; self._is_density_temp = True; return True
        return False

    def set_water_solubility(self, v: float, scale: float = 1.0) -> bool:
        """Set water solubility. Internal unit is g/L (same as Java).
        If passing mg/L, use scale=0.001. If passing g/100mL, use scale=10."""
        val = v * scale
        if 0 <= val <= MAX_WATER_SOLUBILITY:
            self.water_solubility = val; self._is_ws = True; return True
        return False

    def set_water_solubility_temperature(self, t: float) -> bool:
        if 0 < t < 100:
            self.water_solubility_temperature = t; self._is_ws_temp = True; return True
        return False

    def set_permeability_coefficient(self, t: float) -> bool:
        if 0 < t < MAX_PERMEABILITY_COEFFICIENT:
            self.permeability_coefficient = t; self._is_kp = True; return True
        return False

    def set_pka_HA(self, v: float) -> bool:
        if MIN_PKA < v < MAX_PKA:
            self.pka_HA = v; self._is_pka_HA = True; return True
        return False

    def set_pka_BHplus(self, v: float) -> bool:
        if MIN_PKA < v < MAX_PKA:
            self.pka_BHplus = v; self._is_pka_BHplus = True; return True
        return False

    def set_nvv_solubility(self, v: float, scale: float = 1.0) -> bool:
        val = v * scale
        if 0 <= val <= MAX_WATER_SOLUBILITY:
            self.non_volatile_vehicle_solubility = val; self._is_nvv_sol = True; return True
        return False

    def set_f_non_vt(self, v): self.f_non_vt = v; self._is_f_non_vt = True
    def set_f_u_vt(self, v): self.f_u_vt = v; self._is_f_u_vt = True
    def set_f_non_veh(self, v): self.f_non_veh = v; self._is_f_non_veh = True
    def set_f_non_sc(self, v): self.f_non_sc = v; self._is_f_non_sc = True
    def set_f_u_sc(self, v): self.f_u_sc = v; self._is_f_u_sc = True

    # ── Va (Schroeder molar volume) ──
    def get_struct_va(self) -> float:
        return self.num_rings * (-7) + self.num_double_bonds * 7.0 + self.num_triple_bonds * 14.0

    def get_va(self) -> float:
        return self.formula.get_total_atomic_va() + self.get_struct_va()
