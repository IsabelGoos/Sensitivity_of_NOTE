import numpy as np

# Basic settings
DPI = 300

# Atomic numbers Z (first value) and atomic weights A (second value).
# The A values are obtained by averaging over the possible isotopes 
# and using the abundances given in the "NUCLEAR WALLET CARDS".
elements = {
    "H":  (1,   1.01), "He": (2,   4.00),
    "Li": (3,   6.92), "Be": (4,   9.00),
    "B":  (5,  10.80), "C":  (6,  12.01),
    "N":  (7,  14.01), "O":  (8,  16.00),
    "F":  (9,  19.00), "Ne": (10, 20.19),
    "Na": (11, 22.99), "Mg": (12, 24.31),
    "Al": (13, 26.98), "Si": (14, 28.09),
    "P":  (15, 30.97), "S":  (16, 32.07),
    "Cl": (17, 35.45), "Ar": (18, 39.95),
    "K":  (19, 39.10), "Ca": (20, 40.08),
    "Sc": (21, 45.00), "Ti": (22, 47.92),
    "V":  (23, 51.00), "Cr": (24, 52.00),
    "Mn": (25, 54.94), "Fe": (26, 55.85),
    "Co": (27, 59.00), "Ni": (28, 58.69),
    "Cu": (29, 63.62), "Zn": (30, 65.47),
    "Ga": (31, 69.80), "Ge": (32, 72.72),
    "As": (33, 75.00), "Se": (34, 79.04),
    "Br": (35, 79.90), "Kr": (36, 83.89),
    "Rb": (37, 85.56), "Sr": (38, 87.71), 
}

# Elements present in the Earth (with abundances above 0.01%).
labels_Earth = [
    "H", "C", "N", "O", "Na", "Mg", "Al", "Si", "P", "S", 
    "Cl", "Ar", "K", "Ca", "Cr", "Mn", "Fe", "Ni", "Br"
]
Z_Earth = np.array([elements[e][0] for e in labels_Earth])
A_Earth = np.array([elements[e][1] for e in labels_Earth])

# Element fractions for different composites.
# BE   = bulk Earth
# BSE  = bulk silicate Earth
# MORB = mid-ocean ridge basalt
# BE, core, BSE from:
# McDonough, W. F. (2014). 3.16–Compositional model for the Earth’s core.
# Treatise on geochemistry, 2, 559-577.
# MORB from:
# Gale, A., et al. (2013). The chemical compositions of global MORB and their implications for mantle properties. In AGU Fall Meeting Abstracts (Vol. 2013, pp. V32A-01).
wt_rocks = {
    "H":  {"BE":  0.2600, "core":  0.0600, "BSE":  0.1000, "MORB":  0.2000},
    "O":  {"BE": 30.0000, "core":  1.0000, "BSE": 44.0000, "MORB": 44.3800},
    "Na": {"BE":  0.1800, "core":  0.0000, "BSE":  0.0030, "MORB":  2.7900},
    "Mg": {"BE": 15.4000, "core":  0.0000, "BSE": 22.8000, "MORB":  4.5700},
    "Al": {"BE":  1.5900, "core":  0.0000, "BSE":  2.5300, "MORB":  7.7800},
    "Si": {"BE": 15.3000, "core":  3.5000, "BSE": 21.0000, "MORB": 23.5900},
    "P":  {"BE":  0.1500, "core":  0.4300, "BSE":  0.0100, "MORB":  0.0000},
    "S":  {"BE":  1.1000, "core":  3.3400, "BSE":  0.0300, "MORB":  0.0000},
    "K":  {"BE":  0.0016, "core":  0.0000, "BSE":  0.0000, "MORB":  0.1300},
    "Ca": {"BE":  1.7100, "core":  0.0000, "BSE":  2.5300, "MORB":  8.1400},
    "Ti": {"BE":  0.0810, "core":  0.0000, "BSE":  0.0010, "MORB":  1.0100},
    "Cr": {"BE":  0.4200, "core":  0.7500, "BSE":  0.2600, "MORB":  0.0000},
    "Fe": {"BE": 31.9000, "core": 85.3000, "BSE":  6.3000, "MORB":  8.1100}, 
    "Co": {"BE":  0.0900, "core":  0.2600, "BSE":  0.0100, "MORB":  0.0000},
    "Ni": {"BE":  1.8300, "core":  5.2400, "BSE":  0.2000, "MORB":  0.0000},
}



