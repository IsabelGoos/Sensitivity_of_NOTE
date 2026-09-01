import numpy as np

# Basic settings
DPI = 300

# Atomic numbers (Z) and atomic weights (A)
elements = {
    "H":  (1,  1.01),   "He": (2,  4.0),
    "Li": (3,  6.92),   "Be": (4,  9.0),
    "B":  (5,  10.8),   "C":  (6,  12.01),
    "N":  (7,  14.01),  "O":  (8,  16.0),
    "F":  (9,  19.0),   "Ne": (10, 20.19),
    "Na": (11, 22.99),  "Mg": (12, 24.31),
    "Al": (13, 26.98),  "Si": (14, 28.09),
    "P":  (15, 30.97),  "S":  (16, 32.07),
    "Cl": (17, 35.45),  "Ar": (18, 39.95),
    "K":  (19, 39.1),   "Ca": (20, 40.08),
    "Sc": (21, 45.0),   "Ti": (22, 47.92),
    "V":  (23, 51.0),   "Cr": (24, 52.0),
    "Mn": (25, 54.94),  "Fe": (26, 55.85),
    "Co": (27, 59.0),   "Ni": (28, 58.69),
    "Cu": (29, 63.62),  "Zn": (30, 65.47),
    "Ga": (31, 69.8),   "Ge": (32, 72.72),
    "As": (33, 75.0),   "Se": (34, 79.04),
    "Br": (35, 79.904), "Kr": (36, 83.89),
    "Rb": (37, 85.56),  "Sr": (38, 87.71), 
}

# Elements present in the Earth
labels_Earth = [
    "H", "C", "N", "O", "Na", "Mg", "Al", "Si", "P", "S", 
    "Cl", "Ar", "K", "Ca", "Cr", "Mn", "Fe", "Ni", "Br"
]
Z_Earth = np.array([elements[e][0] for e in labels_Earth])
A_Earth = np.array([elements[e][1] for e in labels_Earth])

# Rock elements
wt_rocks = {
    "H":  {"BE": 0.003, "core": 0.001, "BSE": 0.001, "MORB": 0.002},
    "O":  {"BE": 0.150, "core": 0.005, "BSE": 0.220, "MORB": 0.222},
    "Na": {"BE": 0.000, "core": 0.000, "BSE": 0.000, "MORB": 0.005},
    "Mg": {"BE": 0.076, "core": 0.000, "BSE": 0.113, "MORB": 0.023},
    "Al": {"BE": 0.008, "core": 0.000, "BSE": 0.012, "MORB": 0.037},
    "Si": {"BE": 0.076, "core": 0.017, "BSE": 0.105, "MORB": 0.118},
    "P":  {"BE": 0.001, "core": 0.002, "BSE": 0.000, "MORB": 0.000},
    "S":  {"BE": 0.005, "core": 0.017, "BSE": 0.000, "MORB": 0.000},
    "Ca": {"BE": 0.009, "core": 0.000, "BSE": 0.013, "MORB": 0.041},
    "Cr": {"BE": 0.002, "core": 0.003, "BSE": 0.001, "MORB": 0.000},
    "Fe": {"BE": 0.149, "core": 0.397, "BSE": 0.014, "MORB": 0.038}, 
    "Co": {"BE": 0.000, "core": 0.001, "BSE": 0.000, "MORB": 0.000},
    "Ni": {"BE": 0.009, "core": 0.025, "BSE": 0.001, "MORB": 0.000},
}

wt_BE   = np.array([v["BE"]   for v in wt_rocks.values()])
wt_core = np.array([v["core"] for v in wt_rocks.values()])
wt_BSE  = np.array([v["BSE"]  for v in wt_rocks.values()])
wt_MORB = np.array([v["MORB"] for v in wt_rocks.values()])

Z_rocks = np.array([elements[element][0] for element in wt_rocks])
A_rocks = np.array([elements[element][1] for element in wt_rocks])

Ye_BE   = np.average(Z_rocks/A_rocks, weights=wt_BE)
Ye_core = np.average(Z_rocks/A_rocks, weights=wt_core)
Ye_BSE  = np.average(Z_rocks/A_rocks, weights=wt_BSE)
Ye_MORB = np.average(Z_rocks/A_rocks, weights=wt_MORB)
