"""
SystemConstants.java + ChemCons.java → Python constants.
Every value copied verbatim from the decompiled Java source.
"""
import math

# ─── Physical constants ─────────────────────────────────────────────
PI = math.pi
MINUS_BIG = -9.99999999999999e14
R_GAS = 0.082056  # ideal-gas constant used in Kv, ProgramOpt

# ─── Unit conversions ───────────────────────────────────────────────
TORR_TO_PA = 133.32237
PA_TO_TORR = 1.0 / TORR_TO_PA  # 0.007500616738211299

# ─── Validation bounds ──────────────────────────────────────────────
MIN_PH, MAX_PH = 0.0, 14.0
MIN_DENSITY, MAX_DENSITY = 1e-4, 13.0
MIN_MASS, MAX_MASS = 16.0, 800.0
IMP_MASS_DIFF = 0.5
MIN_LOGKOW, MAX_LOGKOW = -2.0, 5.5
MIN_WATER_TEMP, MAX_WATER_TEMP = 0.0, 100.0
MIN_TEMP, MAX_TEMP = -273.0, 5000.0
MAX_VAPOUR_PRESSURE = 1e9
MAX_AMOUNT_APPLIED = 1e6
MAX_PERMEABILITY_COEFFICIENT = 1e7
MAX_PKA, MIN_PKA = 100.0, -50.0
MIN_FRACTIONS, MAX_FRACTIONS = 0.0, 1.0
MAX_WATER_SOLUBILITY = 20000.0
MAX_DURATION_UPPER_LIMIT = 1e7
MAX_DOUBLE_BONDS = 100
MAX_TRIPLE_BONDS = 100
MAX_RINGS = 50

# ─── Default layer thicknesses (µm) ─────────────────────────────────
DEF_VE_THICKNESS = 100.0
DEF_SC_THICKNESS_PARTIALLY_HYDRATED = 13.365
DEF_SC_THICKNESS_FULLY_HYDRATED = 43.365
DEF_DERMIS_IN_VIVO = 2000.0
DEF_DERMIS_IN_VITRO = 400.0
VE_MAX_THICKNESS = 5000.0
DERMIS_MAX_THICKNESS = 20000.0

# ─── Default pH values ──────────────────────────────────────────────
DEF_SC_PH = 5.0
DEF_VE_PH = 7.4
DEF_DERMIS_PH = 7.4
DEF_VEHICLE_PH = 7.4

# ─── Default environmental parameters ───────────────────────────────
DEF_SURFACE_TEMP = 32.0
MIN_SURF_TEMP, MAX_SURF_TEMP = 0.0, 50.0
DEF_WIND_INDOOR = 0.165
DEF_WIND_OUTDOOR = 0.72
MIN_WIND, MAX_WIND = 0.1, 10.0

# ─── Grain-class identifiers ────────────────────────────────────────
GRAIN_OTHER = 1
GRAIN_ALCOHOL = 2
GRAIN_HYDROCARBON = 3

# ─── Hydration strings ──────────────────────────────────────────────
PARTIALLY_HYDRATED = "Partially Hydrated"
FULLY_HYDRATED = "Fully Hydrated"
ENV_INDOORS = "Indoors"
ENV_OUTDOORS = "Outdoors"
ENV_CUSTOM = "Custom defined"

# ─── Volatility / in-vivo flags ──────────────────────────────────────
VOLATILE_VEHICLE = 1
NON_VOLATILE_VEHICLE = 0
IN_VIVO = 0
IN_VITRO = 1

# ─── Vehicle names ──────────────────────────────────────────────────
VV_NAMES = ["Other", "Water", "Ethanol", "Acetone"]

# ─── Atomic data tables (from ChemCons.java) ────────────────────────
ATOMIC_NAMES = [
    "H", "C", "N", "P", "O", "F", "Cl", "Br", "I", "S",
    "Hg", "Si", "Pb", "Zn", "Cu", "As", "Na", "K", "Li", "Ca",
    "Mg", "Be", "Fe", "Ni", "Cr", "Cd", "Al", "Se", "Te", "Sn", "B"
]
ATOMIC_MASSES = [
    1.00794, 12.0107, 14.0067, 30.973762, 15.9994, 18.9984032,
    35.453, 79.904, 126.90447, 32.065,
    200.59, 28.0855, 207.2, 65.38, 63.546, 74.9216,
    22.98976928, 39.0983, 6.941, 40.078, 24.305, 9.012182,
    55.845, 58.6934, 51.9961, 112.411, 26.9815386,
    78.96, 127.6, 118.71, 10.811
]
VA_COEFFICIENTS = [
    7.0, 7.0, 7.0, 20.5, 7.0, 10.5, 24.5, 31.5, 38.5, 21.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
]
