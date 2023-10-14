"""It is a file to create component data to insert in database"""

NAME = "name"
CRITICAL_TEMPERATURE = "critical_temperature"
CRITICAL_PRESSURE = "critical_pressure"
ACENTRIC_FACTOR = "acentric_factor"
AAT1 = "aat1"
AAT2 = "aat2"
AAT3 = "aat3"
R = "r"
Q = "q"

COMPONENTS = {
    "H2S": {
        NAME: "H2S",
        CRITICAL_TEMPERATURE: 373.53,
        CRITICAL_PRESSURE: 89.629,
        ACENTRIC_FACTOR: 0.09417,
        AAT1: 0.44736,
        AAT2: 1.01502,
        AAT3: 0.07205,
        R: 1.00000,
        Q: 1.00000
    },
    "N2": {
        NAME: "N2",
        CRITICAL_TEMPERATURE: 126.60,
        CRITICAL_PRESSURE: 34.000,
        ACENTRIC_FACTOR: 0.03772,
        AAT1: 0.37729,
        AAT2: 0.98587,
        AAT3: 0.05244,
        R: 1.08800,
        Q: 1.04153
    },
    "CO2": {
        NAME: "CO2",
        CRITICAL_TEMPERATURE: 304.21,
        CRITICAL_PRESSURE: 73.830,
        ACENTRIC_FACTOR: 0.22362,
        AAT1: 0.44521,
        AAT2: 0.94117,
        AAT3: 0.17451,
        R: 1.28000,
        Q: 1.32000
    },
    "C1": {
        NAME: "C1",
        CRITICAL_TEMPERATURE: 190.56,
        CRITICAL_PRESSURE: 45.990,
        ACENTRIC_FACTOR: 0.01155,
        AAT1: 0.34674,
        AAT2: 0.97715,
        AAT3: 0.04375,
        R: 1.15200,
        Q: 1.12393
    },
    "C2": {
        NAME: "C2",
        CRITICAL_TEMPERATURE: 305.32,
        CRITICAL_PRESSURE: 48.720,
        ACENTRIC_FACTOR: 0.09949,
        AAT1: 0.40750,
        AAT2: 0.94471,
        AAT3: 0.07863,
        R: 1.70000,
        Q: 1.80000
    },
    "C3": {
        NAME: "C3",
        CRITICAL_TEMPERATURE: 369.83,
        CRITICAL_PRESSURE: 42.480,
        ACENTRIC_FACTOR: 0.15229,
        AAT1: 0.53120,
        AAT2: 0.99911,
        AAT3: 0.06650,
        R: 2.24000,
        Q: 2.48000
    },
    "iC4": {
        NAME: "iC4",
        CRITICAL_TEMPERATURE: 408.14,
        CRITICAL_PRESSURE: 36.480,
        ACENTRIC_FACTOR: 0.18077,
        AAT1: 0.60219,
        AAT2: 0.99617,
        AAT3: 0.06256,
        R: 2.77000,
        Q: 3.15000
    },
    "nC4": {
        NAME: "nC4",
        CRITICAL_TEMPERATURE: 425.12,
        CRITICAL_PRESSURE: 37.960,
        ACENTRIC_FACTOR: 0.20016,
        AAT1: 0.60219,
        AAT2: 0.99617,
        AAT3: 0.06256,
        R: 2.77600,
        Q: 3.15100
    },
    "iC5": {
        NAME: "iC5",
        CRITICAL_TEMPERATURE: 460.43,
        CRITICAL_PRESSURE: 33.810,
        ACENTRIC_FACTOR: 0.22746,
        AAT1: 0.61593,
        AAT2: 0.96228,
        AAT3: 0.08320,
        R: 3.31200,
        Q: 3.82460
    },
    "nC5": {
        NAME: "nC5",
        CRITICAL_TEMPERATURE: 469.70,
        CRITICAL_PRESSURE: 33.700,
        ACENTRIC_FACTOR: 0.25151,
        AAT1: 0.61593,
        AAT2: 0.96228,
        AAT3: 0.08320,
        R: 3.31600,
        Q: 3.82540
    },
    "C6": {
        NAME: "C6",
        CRITICAL_TEMPERATURE: 507.60,
        CRITICAL_PRESSURE: 30.250,
        ACENTRIC_FACTOR: 0.30126,
        AAT1: 0.71923,
        AAT2: 0.98457,
        AAT3: 0.07047,
        R: 3.85600,
        Q: 4.49980
    },
    "C7+": {
        NAME: "C7+",
        CRITICAL_TEMPERATURE: 617.70,
        CRITICAL_PRESSURE: 21.100,
        ACENTRIC_FACTOR: 0.49233,
        AAT1: 1.03502,
        AAT2: 1.00688,
        AAT3: 0.03594,
        R: 6.01600,
        Q: 7.19740
    },
}
