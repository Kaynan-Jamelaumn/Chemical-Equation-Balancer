from collections import defaultdict
import re
from sympy import Matrix
from sympy import Rational

class EquationBalancer:
    """
    A class to balance chemical equations by parsing compounds, building a matrix of element counts,
    solving for stoichiometric coefficients, and constructing the balanced equation.
    """
    def __init__(self, elements_list=None):
        elements = [
            {"symbol": "H",  "name": "Hydrogen",    "atomic_number": 1,  "molar_mass": 1.008,      "nox": [+1, -1], "electrons_to_neutral": 1},
            {"symbol": "He", "name": "Helium",      "atomic_number": 2,  "molar_mass": 4.0026,     "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Li", "name": "Lithium",     "atomic_number": 3,  "molar_mass": 6.94,       "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Be", "name": "Beryllium",   "atomic_number": 4,  "molar_mass": 9.0122,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "B",  "name": "Boron",       "atomic_number": 5,  "molar_mass": 10.81,      "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "C",  "name": "Carbon",      "atomic_number": 6,  "molar_mass": 12.011,     "nox": [+4, -4], "electrons_to_neutral": 4},
            {"symbol": "N",  "name": "Nitrogen",    "atomic_number": 7,  "molar_mass": 14.007,     "nox": [+5, -3], "electrons_to_neutral": 3},
            {"symbol": "O",  "name": "Oxygen",      "atomic_number": 8,  "molar_mass": 15.999,     "nox": [-2],     "electrons_to_neutral": 2},
            {"symbol": "F",  "name": "Fluorine",    "atomic_number": 9,  "molar_mass": 18.998,     "nox": [-1],     "electrons_to_neutral": 1},
            {"symbol": "Ne", "name": "Neon",        "atomic_number": 10, "molar_mass": 20.180,     "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Na", "name": "Sodium",      "atomic_number": 11, "molar_mass": 22.990,     "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Mg", "name": "Magnesium",   "atomic_number": 12, "molar_mass": 24.305,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Al", "name": "Aluminum",    "atomic_number": 13, "molar_mass": 26.982,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Si", "name": "Silicon",     "atomic_number": 14, "molar_mass": 28.085,     "nox": [+4, -4], "electrons_to_neutral": 4},
            {"symbol": "P",  "name": "Phosphorus",  "atomic_number": 15, "molar_mass": 30.974,     "nox": [+5, -3], "electrons_to_neutral": 3},
            {"symbol": "S",  "name": "Sulfur",      "atomic_number": 16, "molar_mass": 32.06,      "nox": [+6, -2], "electrons_to_neutral": 2},
            {"symbol": "Cl", "name": "Chlorine",    "atomic_number": 17, "molar_mass": 35.45,      "nox": [+7, -1], "electrons_to_neutral": 1},
            {"symbol": "Ar", "name": "Argon",       "atomic_number": 18, "molar_mass": 39.948,     "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "K",  "name": "Potassium",   "atomic_number": 19, "molar_mass": 39.098,     "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Ca", "name": "Calcium",     "atomic_number": 20, "molar_mass": 40.078,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Sc", "name": "Scandium",    "atomic_number": 21, "molar_mass": 44.956,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Ti", "name": "Titanium",    "atomic_number": 22, "molar_mass": 47.867,     "nox": [+4, +3], "electrons_to_neutral": 4},
            {"symbol": "V",  "name": "Vanadium",    "atomic_number": 23, "molar_mass": 50.942,     "nox": [+5, +4], "electrons_to_neutral": 5},
            {"symbol": "Cr", "name": "Chromium",    "atomic_number": 24, "molar_mass": 51.996,     "nox": [+6, +3], "electrons_to_neutral": 6},
            {"symbol": "Mn", "name": "Manganese",   "atomic_number": 25, "molar_mass": 54.938,     "nox": [+7, +2], "electrons_to_neutral": 7},
            {"symbol": "Fe", "name": "Iron",        "atomic_number": 26, "molar_mass": 55.845,     "nox": [+3, +2], "electrons_to_neutral": 3},
            {"symbol": "Co", "name": "Cobalt",      "atomic_number": 27, "molar_mass": 58.933,     "nox": [+3, +2], "electrons_to_neutral": 3},
            {"symbol": "Ni", "name": "Nickel",      "atomic_number": 28, "molar_mass": 58.693,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Cu", "name": "Copper",      "atomic_number": 29, "molar_mass": 63.546,     "nox": [+2, +1], "electrons_to_neutral": 2},
            {"symbol": "Zn", "name": "Zinc",        "atomic_number": 30, "molar_mass": 65.38,      "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Ga", "name": "Gallium",     "atomic_number": 31, "molar_mass": 69.723,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Ge", "name": "Germanium",   "atomic_number": 32, "molar_mass": 72.630,     "nox": [+4, -4], "electrons_to_neutral": 4},
            {"symbol": "As", "name": "Arsenic",     "atomic_number": 33, "molar_mass": 74.922,     "nox": [+5, -3], "electrons_to_neutral": 3},
            {"symbol": "Se", "name": "Selenium",    "atomic_number": 34, "molar_mass": 78.971,     "nox": [+6, -2], "electrons_to_neutral": 2},
            {"symbol": "Br", "name": "Bromine",     "atomic_number": 35, "molar_mass": 79.904,     "nox": [+7, -1], "electrons_to_neutral": 1},
            {"symbol": "Kr", "name": "Krypton",     "atomic_number": 36, "molar_mass": 83.798,     "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Rb", "name": "Rubidium",    "atomic_number": 37, "molar_mass": 85.468,     "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Sr", "name": "Strontium",   "atomic_number": 38, "molar_mass": 87.62,      "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Y",  "name": "Yttrium",     "atomic_number": 39, "molar_mass": 88.906,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Zr", "name": "Zirconium",   "atomic_number": 40, "molar_mass": 91.224,     "nox": [+4],     "electrons_to_neutral": 4},
            {"symbol": "Nb", "name": "Niobium",     "atomic_number": 41, "molar_mass": 92.906,     "nox": [+5],     "electrons_to_neutral": 5},
            {"symbol": "Mo", "name": "Molybdenum",  "atomic_number": 42, "molar_mass": 95.95,      "nox": [+6],     "electrons_to_neutral": 6},
            {"symbol": "Tc", "name": "Technetium",  "atomic_number": 43, "molar_mass": 98,         "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Ru", "name": "Ruthenium",   "atomic_number": 44, "molar_mass": 101.07,     "nox": [+8, +3], "electrons_to_neutral": 8},
            {"symbol": "Rh", "name": "Rhodium",     "atomic_number": 45, "molar_mass": 102.91,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Pd", "name": "Palladium",   "atomic_number": 46, "molar_mass": 106.42,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Ag", "name": "Silver",      "atomic_number": 47, "molar_mass": 107.87,     "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Cd", "name": "Cadmium",     "atomic_number": 48, "molar_mass": 112.41,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "In", "name": "Indium",      "atomic_number": 49, "molar_mass": 114.82,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Sn", "name": "Tin",         "atomic_number": 50, "molar_mass": 118.71,     "nox": [+4, +2], "electrons_to_neutral": 4},
            {"symbol": "Sb", "name": "Antimony",    "atomic_number": 51, "molar_mass": 121.76,     "nox": [+5, -3], "electrons_to_neutral": 3},
            {"symbol": "Te", "name": "Tellurium",   "atomic_number": 52, "molar_mass": 127.60,     "nox": [+6, -2], "electrons_to_neutral": 2},
            {"symbol": "I",  "name": "Iodine",      "atomic_number": 53, "molar_mass": 126.90,     "nox": [+7, -1], "electrons_to_neutral": 1},
            {"symbol": "Xe", "name": "Xenon",       "atomic_number": 54, "molar_mass": 131.29,     "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Cs", "name": "Cesium",      "atomic_number": 55, "molar_mass": 132.91,     "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Ba", "name": "Barium",      "atomic_number": 56, "molar_mass": 137.33,     "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "La", "name": "Lanthanum",   "atomic_number": 57, "molar_mass": 138.91,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Ce", "name": "Cerium",      "atomic_number": 58, "molar_mass": 140.12,     "nox": [+4, +3], "electrons_to_neutral": 4},
            {"symbol": "Pr", "name": "Praseodymium","atomic_number": 59, "molar_mass": 140.91,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Nd", "name": "Neodymium",   "atomic_number": 60, "molar_mass": 144.24,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Pm", "name": "Promethium",  "atomic_number": 61, "molar_mass": 145,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Sm", "name": "Samarium",    "atomic_number": 62, "molar_mass": 150.36,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Eu", "name": "Europium",    "atomic_number": 63, "molar_mass": 151.96,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Gd", "name": "Gadolinium",  "atomic_number": 64, "molar_mass": 157.25,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Tb", "name": "Terbium",     "atomic_number": 65, "molar_mass": 158.93,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Dy", "name": "Dysprosium",  "atomic_number": 66, "molar_mass": 162.50,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Ho", "name": "Holmium",     "atomic_number": 67, "molar_mass": 164.93,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Er", "name": "Erbium",      "atomic_number": 68, "molar_mass": 167.26,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Tm", "name": "Thulium",     "atomic_number": 69, "molar_mass": 168.93,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Yb", "name": "Ytterbium",   "atomic_number": 70, "molar_mass": 173.05,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Lu", "name": "Lutetium",    "atomic_number": 71, "molar_mass": 174.97,     "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Hf", "name": "Hafnium",     "atomic_number": 72, "molar_mass": 178.49,     "nox": [+4],     "electrons_to_neutral": 4},
            {"symbol": "Ta", "name": "Tantalum",    "atomic_number": 73, "molar_mass": 180.95,     "nox": [+5],     "electrons_to_neutral": 5},
            {"symbol": "W",  "name": "Tungsten",    "atomic_number": 74, "molar_mass": 183.84,     "nox": [+6],     "electrons_to_neutral": 6},
            {"symbol": "Re", "name": "Rhenium",     "atomic_number": 75, "molar_mass": 186.21,     "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Os", "name": "Osmium",      "atomic_number": 76, "molar_mass": 190.23,     "nox": [+8],     "electrons_to_neutral": 8},
            {"symbol": "Ir", "name": "Iridium",     "atomic_number": 77, "molar_mass": 192.22,     "nox": [+4],     "electrons_to_neutral": 4},
            {"symbol": "Pt", "name": "Platinum",    "atomic_number": 78, "molar_mass": 195.08,     "nox": [+4, +2], "electrons_to_neutral": 4},
            {"symbol": "Au", "name": "Gold",        "atomic_number": 79, "molar_mass": 196.97,     "nox": [+3, +1], "electrons_to_neutral": 3},
            {"symbol": "Hg", "name": "Mercury",     "atomic_number": 80, "molar_mass": 200.59,     "nox": [+2, +1], "electrons_to_neutral": 2},
            {"symbol": "Tl", "name": "Thallium",    "atomic_number": 81, "molar_mass": 204.38,     "nox": [+3, +1], "electrons_to_neutral": 3},
            {"symbol": "Pb", "name": "Lead",        "atomic_number": 82, "molar_mass": 207.2,      "nox": [+4, +2], "electrons_to_neutral": 4},
            {"symbol": "Bi", "name": "Bismuth",     "atomic_number": 83, "molar_mass": 208.98,     "nox": [+5, +3], "electrons_to_neutral": 5},
            {"symbol": "Po", "name": "Polonium",    "atomic_number": 84, "molar_mass": 209,        "nox": [+4, +2], "electrons_to_neutral": 4},
            {"symbol": "At", "name": "Astatine",    "atomic_number": 85, "molar_mass": 210,        "nox": [-1],     "electrons_to_neutral": 1},
            {"symbol": "Rn", "name": "Radon",       "atomic_number": 86, "molar_mass": 222,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Fr", "name": "Francium",    "atomic_number": 87, "molar_mass": 223,        "nox": [+1],     "electrons_to_neutral": 1},
            {"symbol": "Ra", "name": "Radium",      "atomic_number": 88, "molar_mass": 226,        "nox": [+2],     "electrons_to_neutral": 2},
            {"symbol": "Ac", "name": "Actinium",    "atomic_number": 89, "molar_mass": 227,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Th", "name": "Thorium",     "atomic_number": 90, "molar_mass": 232.04,     "nox": [+4],     "electrons_to_neutral": 4},
            {"symbol": "Pa", "name": "Protactinium","atomic_number": 91, "molar_mass": 231.04,     "nox": [+5],     "electrons_to_neutral": 5},
            {"symbol": "U",  "name": "Uranium",     "atomic_number": 92, "molar_mass": 238.03,     "nox": [+6],     "electrons_to_neutral": 6},
            {"symbol": "Np", "name": "Neptunium",   "atomic_number": 93, "molar_mass": 237,        "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Pu", "name": "Plutonium",   "atomic_number": 94, "molar_mass": 244,        "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Am", "name": "Americium",   "atomic_number": 95, "molar_mass": 243,        "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Cm", "name": "Curium",      "atomic_number": 96, "molar_mass": 247,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Bk", "name": "Berkelium",   "atomic_number": 97, "molar_mass": 247,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Cf", "name": "Californium", "atomic_number": 98, "molar_mass": 251,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Es", "name": "Einsteinium", "atomic_number": 99, "molar_mass": 252,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Fm", "name": "Fermium",     "atomic_number": 100,"molar_mass": 257,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Md", "name": "Mendelevium", "atomic_number": 101,"molar_mass": 258,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "No", "name": "Nobelium",    "atomic_number": 102,"molar_mass": 259,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Lr", "name": "Lawrencium",  "atomic_number": 103,"molar_mass": 262,        "nox": [+3],     "electrons_to_neutral": 3},
            {"symbol": "Rf", "name": "Rutherfordium","atomic_number": 104,"molar_mass": 267,       "nox": [+4],     "electrons_to_neutral": 4},
            {"symbol": "Db", "name": "Dubnium",     "atomic_number": 105,"molar_mass": 268,        "nox": [+5],     "electrons_to_neutral": 5},
            {"symbol": "Sg", "name": "Seaborgium",  "atomic_number": 106,"molar_mass": 269,        "nox": [+6],     "electrons_to_neutral": 6},
            {"symbol": "Bh", "name": "Bohrium",     "atomic_number": 107,"molar_mass": 270,        "nox": [+7],     "electrons_to_neutral": 7},
            {"symbol": "Hs", "name": "Hassium",     "atomic_number": 108,"molar_mass": 269,        "nox": [+8],     "electrons_to_neutral": 8},
            {"symbol": "Mt", "name": "Meitnerium",  "atomic_number": 109,"molar_mass": 278,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Ds", "name": "Darmstadtium","atomic_number": 110,"molar_mass": 281,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Rg", "name": "Roentgenium", "atomic_number": 111,"molar_mass": 282,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Cn", "name": "Copernicium", "atomic_number": 112,"molar_mass": 285,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Nh", "name": "Nihonium",    "atomic_number": 113,"molar_mass": 286,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Fl", "name": "Flerovium",   "atomic_number": 114,"molar_mass": 289,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Mc", "name": "Moscovium",   "atomic_number": 115,"molar_mass": 290,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Lv", "name": "Livermorium", "atomic_number": 116,"molar_mass": 293,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Ts", "name": "Tennessine",  "atomic_number": 117,"molar_mass": 294,        "nox": [0],      "electrons_to_neutral": 0},
            {"symbol": "Og", "name": "Oganesson",   "atomic_number": 118,"molar_mass": 294,        "nox": [0],      "electrons_to_neutral": 0},
        ]

        self.elements_list = elements
        self.element_nox_map = {elem['symbol']: elem['nox'] for elem in self.elements_list}
        self.polyatomic_ions = {
            "NH4": {"N": -3, "H": +1},
            "C2H3O2": {"C": +3, "H": +1, "O": -2},  # Acetate
            "CN": {"C": +2, "N": -3},               # Cyanide
            "OH": {"O": -2, "H": +1},               # Hydroxide
            "NO2": {"N": +3, "O": -2},              # Nitrite
            "NO3": {"N": +5, "O": -2},              # Nitrate
            "MnO4": {"Mn": +7, "O": -2},            # Permanganate
            "SCN": {"S": -2, "C": +4, "N": -3},     # Thiocyanate
            "ClO": {"Cl": +1, "O": -2},             # Hypochlorite
            "ClO2": {"Cl": +3, "O": -2},            # Chlorite
            "ClO3": {"Cl": +5, "O": -2},            # Chlorate
            "ClO4": {"Cl": +7, "O": -2},            # Perchlorate
            "HCO3": {"H": +1, "C": +4, "O": -2},    # Hydrogen Carbonate
            "H2PO4": {"H": +1, "P": +5, "O": -2},   # Dihydrogen Phosphate
            "HSO3": {"H": +1, "S": +4, "O": -2},    # Hydrogen Sulfite
            "HSO4": {"H": +1, "S": +6, "O": -2},    # Hydrogen Sulfate
            "CO3": {"C": +4, "O": -2},              # Carbonate
            "SO3": {"S": +4, "O": -2},              # Sulfite
            "SO4": {"S": +6, "O": -2},              # Sulfate
            "CrO4": {"Cr": +6, "O": -2},            # Chromate
            "Cr2O7": {"Cr": +6, "O": -2},           # Dichromate
            "C2O4": {"C": +3, "O": -2},             # Oxalate
            "O2": {"O": -1},                        # Peroxide
            "S2O3": {"S": +2, "O": -2},             # Thiosulfate
            "PO4": {"P": +5, "O": -2},              # Phosphate
            "PO3": {"P": +3, "O": -2},              # Phosphite
            "AsO4": {"As": +5, "O": -2},            # Arsenate
            "BO3": {"B": +3, "O": -2},              # Borate
            "SiO3": {"Si": +4, "O": -2},            # Silicate
            "AlO2": {"Al": +3, "O": -2},            # Aluminate
            "Fe(CN)6": {"Fe": +3, "C": +2, "N": -3},# Ferricyanide
            "Fe(CN)6": {"Fe": +2, "C": +2, "N": -3},# Ferrocyanide
            "MnO4": {"Mn": +6, "O": -2},            # Manganate
            "TcO4": {"Tc": +7, "O": -2},            # Pertechnetate
            "ReO4": {"Re": +7, "O": -2},            # Perrhenate
            "P2O7": {"P": +5, "O": -2},             # Pyrophosphate
            "S4O6": {"S": +2.5, "O": -2},           # Tetrathionate
        }


    def gcd(self, a, b):
        """
        Compute the greatest common divisor (GCD) of two numbers using the Euclidean algorithm.

        Args:
            a (int): The first number.
            b (int): The second number.

        Returns:
            int: The GCD of `a` and `b`.
        """
        while b:
            a, b = b, a % b
        return a
    def compute_oxidation_states(self, elements_dict, charge):
        """
        Compute the oxidation states of elements in a compound based on its composition and charge.

        Args:
            elements_dict (dict): A dictionary where keys are element symbols and values are their counts in the compound.
            charge (int): The total charge of the compound.

        Returns:
            dict: A dictionary where keys are element symbols and values are their computed oxidation states.

        Raises:
            ValueError: If the oxidation states cannot be determined due to insufficient information.

        Description:
        This function calculates the oxidation states of elements in a compound using the following rules:
        1. Polyatomic ions are handled first, using predefined oxidation states for their constituent elements.
        2. Oxygen is assigned an oxidation state of -2 (except in special cases like peroxides).
        3. Hydrogen is assigned an oxidation state of +1 (except in metal hydrides).
        4. For other elements, the first possible oxidation state from the element's oxidation state map is used.
        5. If the total oxidation does not match the compound's charge, the function adjusts the oxidation state
           of an element with variable oxidation states to balance the charge.

        Example:
            Consider the compound SO4^2- (sulfate ion).

            Step 1: Handle polyatomic ions.
            - SO4 is a polyatomic ion with predefined oxidation states: S = +6, O = -2.

            Step 2: Assign oxidation states for oxygen and hydrogen.
            - Oxygen is assigned -2 (default rule).
            - Hydrogen is not present in this compound.

            Step 3: Assign oxidation states for other elements.
            - Sulfur (S) is assigned +6 (from the polyatomic ion rules).

            Step 4: Check total oxidation.
            - Total oxidation = (1 * +6) + (4 * -2) = +6 - 8 = -2.
            - This matches the charge of the compound (-2), so no adjustment is needed.

            Final output:
            - oxidation_states: {'S': 6, 'O': -2}
        """
        # Initialize a dictionary to store the oxidation states of each element in the compound.
        oxidation_states = {}
        total_oxidation = 0  # Tracks the sum of (oxidation state * count) for all elements.

        # Handle polyatomic ions.
        # Polyatomic ions have fixed oxidation states for their constituent elements.
        for ion, ion_nox in self.polyatomic_ions.items():
            if ion in elements_dict:
                # Assign the predefined oxidation states for elements in the polyatomic ion.
                for element, nox in ion_nox.items():
                    oxidation_states[element] = nox
                    # Update the total oxidation by adding (oxidation state * count).
                    total_oxidation += nox * elements_dict[element]

        # Assign oxidation states for oxygen and hydrogen.
        # Oxygen typically has an oxidation state of -2 (except in peroxides, superoxides, etc.).
        if "O" in elements_dict and "O" not in oxidation_states:
            oxidation_states["O"] = -2
            total_oxidation += -2 * elements_dict["O"]

        # Hydrogen typically has an oxidation state of +1 (except in metal hydrides).
        if "H" in elements_dict and "H" not in oxidation_states:
            oxidation_states["H"] = +1
            total_oxidation += +1 * elements_dict["H"]

        # Assign oxidation states for other elements.
        # For elements not in polyatomic ions and not oxygen or hydrogen,
        # use the first possible oxidation state from the element's oxidation state map.
        for element, count in elements_dict.items():
            if element not in oxidation_states:
                # Get the possible oxidation states for the element (default to [0] if not found).
                possible_nox = self.element_nox_map.get(element, [0])
                # Assign the first possible oxidation state.
                oxidation_states[element] = possible_nox[0]
                # Update the total oxidation by adding (oxidation state * count).
                total_oxidation += possible_nox[0] * count

        # Adjust for variable oxidation states if the total oxidation does not match the charge.
        if total_oxidation != charge:
            # Find an element with variable oxidation states (i.e., more than one possible oxidation state).
            variable_element = None
            for element in elements_dict:
                if len(self.element_nox_map.get(element, [0])) > 1 and element not in self.polyatomic_ions:
                    variable_element = element
                    break

            if variable_element:
                # Calculate the difference between the total oxidation and the compound's charge.
                delta = charge - total_oxidation
                # Adjust the oxidation state of the variable element to balance the charge.
                # Use Rational to handle fractional adjustments precisely.
                oxidation_states[variable_element] += Rational(delta, elements_dict[variable_element])
            else:
                # If no variable element is found, raise an error.
                raise ValueError(
                    f"Could not determine oxidation states for compound with elements {elements_dict} and charge {charge}"
                )

        return oxidation_states

    def parse_compound(self, compound):
        """
        Parse a chemical compound into its constituent elements and their counts, including handling charges.

        Args:
            compound (str): The chemical formula of the compound (e.g., "H2O", "Fe^3+", "SO4^2-").

        Returns:
            tuple: A tuple containing:
                - elements_dict (defaultdict): A dictionary mapping each element to its count in the compound.
                - charge (int): The net charge of the compound.
                
        Example:
            Consider the compound "Fe2(SO4)3^2-".

            Step 1: Handle charge notation.
            - The compound has a charge notation: "^2-".
            - Extract the charge: -2.

            Step 2: Parse the formula part "Fe2(SO4)3".
            - Use a regex to match elements and polyatomic ions:
            - "Fe2" -> Fe (count = 2)
            - "(SO4)3" -> SO4 (count = 3)

            Step 3: Recursively parse polyatomic ions.
            - Parse "SO4":
            - S (count = 1)
            - O (count = 4)

            Step 4: Combine counts.
            - Fe: 2
            - S: 3 (from (SO4)3)
            - O: 12 (from (SO4)3, 4 O per SO4 * 3)

            Final output:
            - elements_dict: {'Fe': 2, 'S': 3, 'O': 12}
            - charge: -2
        """
        formula_part = compound
        charge = 0

       # Handle charge notation (e.g., Fe^3+ or SO4^2-)
        if '^' in compound:
            # Split the compound into the formula part and the charge string using '^' as the delimiter.
            # Example: "Fe^3+" -> formula_part = "Fe", charge_str = "3+"
            formula_part, charge_str = compound.split('^', 1)

            # Initialize the charge sign to +1 (default for positive charges).
            sign = 1  # Default to positive charge
            digits = ''  # Initialize an empty string to store the numeric part of the charge.

            # Check if the charge string ends with a '-' (negative charge).
            if charge_str.endswith('-'):
                sign = -1  # Set the sign to -1 for negative charges.
                digits = charge_str[:-1]  # Extract the numeric part by removing the '-' sign.
            # Check if the charge string ends with a '+' (positive charge).
            elif charge_str.endswith('+'):
                digits = charge_str[:-1]  # Extract the numeric part by removing the '+' sign.
            else:
                # If the charge string does not end with '+' or '-', assume it's a positive charge.
                digits = charge_str

            # Calculate the charge value:
            # - If digits is not empty, convert it to an integer and multiply by the sign.
            # - If digits is empty (e.g., "^+"), assume the charge is +1 or -1 based on the sign.
            charge = sign * int(digits) if digits else sign * 1
        else:
            # Handle charge notation without '^' (e.g., Fe+3 or SO4-2).
            # Use a regex to search for a charge at the end of the formula part.
            # The regex looks for:
            # - A '+' or '-' sign ([+-]).
            # - Zero or more digits (\d*) after the sign.
            charge_match = re.search(r'([+-])(\d*)$', formula_part)
            if charge_match:
                # Extract the sign and digits from the regex match.
                sign_str, digits = charge_match.groups()

                # Calculate the charge value:
                # - Combine the sign and digits (e.g., "+3" or "-2").
                # - If digits is empty (e.g., "+"), assume the charge is +1 or -1.
                charge = int(sign_str + (digits if digits else '1'))

                # Remove the charge notation from the formula part to avoid parsing it as part of the element.
                formula_part = formula_part[:charge_match.start()]

        # Parse the formula part into elements and their counts.
        # Use a defaultdict to store the count of each element, defaulting to 0 for missing keys.
        elements_dict = defaultdict(int)  # Default value is 0 for missing keys

        # Use a regex to find all matches of elements or polyatomic ions in the formula part.
        # The regex matches:
        # - Polyatomic ions enclosed in parentheses (e.g., (SO4)) and their optional subscripts.
        # - Single elements (e.g., H, O, Fe) and their optional subscripts.
        matches = re.finditer(r'(\(([^)]+)\)|([A-Z][a-z]*))(\d*)', formula_part)

        # Iterate over all matches found by the regex.
        for match in matches:
            # Extract the polyatomic ion content (if present) and the element symbol.
            poly_content, elem = match.group(2, 3)  # poly_content is for polyatomic ions, elem is for single elements

            # Extract the subscript (count) of the element or polyatomic ion.
            subscript = match.group(4)  # The subscript (count) of the element or polyatomic ion

            # Convert the subscript to an integer. If no subscript is present, default to 1.
            count = int(subscript) if subscript else 1  # Default count is 1 if no subscript

            if poly_content:
                # If the match is a polyatomic ion (e.g., (SO4)), recursively parse its content.
                sub_elements, _ = self.parse_compound(poly_content)

                # Add the counts of the elements in the polyatomic ion to the main elements_dict.
                # Multiply by the subscript of the polyatomic ion (e.g., (SO4)2 -> multiply counts by 2).
                for el, el_count in sub_elements.items():
                    elements_dict[el] += el_count * count
            else:
                # If the match is a single element (e.g., H2, O2), add its count to the elements_dict.
                elements_dict[elem] += count

        # Return the dictionary of element counts and the calculated charge.
        return elements_dict, charge


    def parse_equation(self, equation):
        """
        Parse a chemical equation into reactants and products, and further parse each compound.

        Args:
            equation (str): The chemical equation (e.g., "H2 + O2 -> H2O").

        Returns:
            tuple: A tuple containing:
                - A tuple of parsed data for reactants and products.
                - A tuple of reactant and product formulas.
        Example:
            Consider the equation: "H2 + O2 -> H2O"

            Step 1: Split the equation into reactants and products.
            - Reactants: "H2 + O2"
            - Products: "H2O"

            Step 2: Parse the reactants side.
            - Split into compounds: ["H2", "O2"]
            - Parse each compound:
            - H2: {'elements': {'H': 2}, 'charge': 0, 'oxidation_states': {'H': 1}}
            - O2: {'elements': {'O': 2}, 'charge': 0, 'oxidation_states': {'O': 0}}

            Step 3: Parse the products side.
            - Split into compounds: ["H2O"]
            - Parse the compound:
            - H2O: {'elements': {'H': 2, 'O': 1}, 'charge': 0, 'oxidation_states': {'H': 1, 'O': -2}}

            Step 4: Return the parsed data and formulas.
            - Parsed data: ([{'H2': ...}, {'O2': ...}], [{'H2O': ...}])
            - Formulas: (["H2", "O2"], ["H2O"])
        """
        reactants_str, products_str = equation.split('->')

        def parse_side(side_str):
            """
            Parse one side of the equation (reactants or products) into individual compounds.

            Args:
                side_str (str): The string representation of one side of the equation.

            Returns:
                tuple: A tuple containing:
                    - parsed_data (list): A list of dictionaries with element counts and charges for each compound.
                    - formulas (list): A list of chemical formulas for each compound.
            """
            # Split only on plus signs that are surrounded by spaces
            # This ensures that we correctly split the chemical equation into individual compounds.
            # Example: "H2 + O2 -> H2O" will split into ["H2", "O2"] for the reactants.
            formulas = [c.strip() for c in re.split(r'\s+\+\s+', side_str)]

            # Initialize an empty list to store parsed data for each compound.
            # Each entry in this list will be a dictionary containing:
            # - 'elements': A dictionary of element counts for the compound.
            # - 'charge': The net charge of the compound.
            parsed_data = []

            # Iterate over each compound formula in the list.
            for formula in formulas:
                # Parse the compound formula into its constituent elements and charge.
                # The `parse_compound` method returns:
                # - elements_dict: A dictionary mapping each element to its count in the compound.
                # - charge: The net charge of the compound.
                elements_dict, charge = self.parse_compound(formula)

                oxidation_states = self.compute_oxidation_states(elements_dict, charge)

                # Append the parsed data for this compound to the `parsed_data` list.
                parsed_data.append({'elements': elements_dict, 'charge': charge, 'oxidation_states': oxidation_states})

            # Return the parsed data and the original list of formulas.
            # - `parsed_data`: A list of dictionaries containing element counts and charges for each compound.
            # - `formulas`: A list of the original compound formulas (e.g., ["H2", "O2"]).
            return parsed_data, formulas

        # Parse the reactants side of the equation.
        # Example: If `reactants_str` is "H2 + O2", this will:
        # - Split it into ["H2", "O2"].
        # - Parse each compound to get element counts and charges.
        reactants_data, reactant_formulas = parse_side(reactants_str)

        # Parse the products side of the equation.
        # Example: If `products_str` is "H2O", this will:
        # - Split it into ["H2O"].
        # - Parse the compound to get element counts and charges.
        products_data, product_formulas = parse_side(products_str)

        # Return the parsed data for both reactants and products, along with their original formulas.
        # The returned value is a tuple containing:
        # - A tuple of parsed data for reactants and products.
        # - A tuple of reactant and product formulas.
        return (reactants_data, products_data), (reactant_formulas, product_formulas)

    def build_matrix(self, reactants_data, products_data):
        """
        Build a matrix representing the element counts and charges for the reactants and products.

        Args:
            reactants_data (list): Parsed data for reactants.
            products_data (list): Parsed data for products.
            oxidation_states: A dictionary of oxidation states for each element in the compound.

        Returns:
            list: A matrix where each row represents an element (or charge), and each column represents a compound.


        Example:
            Consider the unbalanced equation: H2 + O2 -> H2O

            Step 1: Collect all unique elements in the equation.
            - Elements: H, O

            Step 2: Sort the elements alphabetically.
            - Sorted elements: H, O

            Step 3: Build the matrix rows for each element.
            - For H:
            - Reactants: H2 has 2 H, O2 has 0 H.
            - Products: H2O has 2 H.
            - Row: [2, 0, -2]  # 2 for H2, 0 for O2, -2 for H2O
            - For O:
            - Reactants: H2 has 0 O, O2 has 2 O.
            - Products: H2O has 1 O.
            - Row: [0, 2, -1]  # 0 for H2, 2 for O2, -1 for H2O

            Step 4: Add a charge balance row.
            - No charges in this equation, so the row is [0, 0, 0].

            Step 5: Add an electron transfer row (for redox reactions).
            - No redox in this equation, so the row is [0, 0, 0].

            Final matrix:
            [
                [2, 0, -2],  # Row for Hydrogen (H)
                [0, 2, -1],  # Row for Oxygen (O)
                [0, 0, 0],   # Row for charge balance
                [0, 0, 0]    # Row for electron transfer balance
            ]

            This matrix represents the system of equations:
            - 2a + 0b - 2c = 0  (Hydrogen balance)
            - 0a + 2b - 1c = 0  (Oxygen balance)
            - 0a + 0b + 0c = 0  (Charge balance)
            - 0a + 0b + 0c = 0  (Electron transfer balance)

            Solving this system gives the coefficients for the balanced equation:
            - a = 2 (H2)
            - b = 1 (O2)
            - c = 2 (H2O)

            Balanced equation: 2H2 + O2 -> 2H2O
        """
        # Initialize a set to store all unique elements present in the reactants and products.
        # Using a set ensures that each element is included only once, even if it appears in multiple compounds.
        all_elements = set()

        # Iterate over all compounds in both reactants and products.
        # `reactants_data + products_data` combines the two lists into a single list for iteration.
        for comp in reactants_data + products_data:
            # Update the set with the keys from the 'elements' dictionary of the current compound.
            # `comp['elements'].keys()` gives the list of elements in the compound.
            all_elements.update(comp['elements'].keys())

        # Sort the elements alphabetically to ensure consistent ordering in the matrix.
        # This is important because the order of elements affects the structure of the matrix.
        all_elements = sorted(all_elements)  # Sort elements for consistent ordering

        # Initialize an empty list to store the matrix rows.
        # Each row in the matrix will represent an element and its counts in the reactants and products.
        matrix = []

        # Iterate over each element in the sorted list of elements.
        for el in all_elements:
            # Initialize an empty row for the current element.
            row = []

            # Add counts for the current element in all reactants.
            for comp in reactants_data:
                # Get the count of the current element in the reactant compound.
                # If the element is not present in the compound, use 0 as the default value.
                row.append(comp['elements'].get(el, 0))  # Default to 0 if element is not present

            # Subtract counts for the current element in all products.
            # This is done because products are on the right side of the equation, and their counts are subtracted
            # to balance the equation (reactants = products).
            for comp in products_data:
                # Get the count of the current element in the product compound.
                # If the element is not present in the compound, use 0 as the default value.
                # Multiply by -1 to represent the subtraction.
                row.append(-comp['elements'].get(el, 0))

            # Append the completed row to the matrix.
            matrix.append(row)

        # Add a charge balance row to the matrix.
        # This row ensures that the total charge on the reactants' side equals the total charge on the products' side.
        # Charge balance is a critical part of balancing chemical equations, especially for ionic reactions.
        charge_row = []
        # Iterate over all reactant compounds and add their charges to the charge row.
        for comp in reactants_data:
            # Append the charge of the current reactant compound to the charge row.
            # The charge is stored in the 'charge' key of the compound's dictionary.
            charge_row.append(comp['charge'])

        
        # Iterate over all product compounds and subtract their charges from the charge row.
        # This is done because products are on the right side of the equation, and their charges are subtracted
        # to balance the equation (reactants' total charge = products' total charge).
        for comp in products_data:
            # Append the negative charge of the current product compound to the charge row.
            # Multiplying by -1 ensures that the charge is subtracted.
            charge_row.append(-comp['charge'])

        # Append the charge row to the matrix.
        # This row represents the charge balance equation:
        # (Sum of charges of reactants) - (Sum of charges of products) = 0
        matrix.append(charge_row)

        # Add an electron transfer row to the matrix.
        # This row ensures that the total number of electrons lost in oxidation equals the total number of electrons gained in reduction.
        electron_row = []
        # Iterate over all compounds (both reactants and products).
        for comp in reactants_data + products_data:
            # Initialize a variable to store the contribution of the current compound to electron transfer.
            contrib = 0

            # Iterate over each element in the current compound.
            for el, count in comp['elements'].items():
                # If the compound is a reactant, calculate the change in oxidation state (delta) for the element.
                if comp in reactants_data:
                    # Get the oxidation state of the element in the reactant.
                    reactant_nox = comp['oxidation_states'][el]
                    # Find the oxidation state of the same element in the products.
                    # If the element is not present in any product, use the reactant's oxidation state as a default.
                    product_noxs = [p['oxidation_states'].get(el, reactant_nox) for p in products_data if el in p['elements']]
                    product_nox = product_noxs[0] if product_noxs else reactant_nox
                    # Calculate the change in oxidation state (delta).
                    delta = product_nox - reactant_nox
                # If the compound is a product, calculate the change in oxidation state (delta) for the element.
                else:
                    # Get the oxidation state of the element in the product.
                    product_nox = comp['oxidation_states'][el]
                    # Find the oxidation state of the same element in the reactants.
                    # If the element is not present in any reactant, use the product's oxidation state as a default.
                    reactant_noxs = [r['oxidation_states'].get(el, product_nox) for r in reactants_data if el in r['elements']]
                    reactant_nox = reactant_noxs[0] if reactant_noxs else product_nox
                    # Calculate the change in oxidation state (delta).
                    delta = product_nox - reactant_nox

                # Multiply the delta by the count of the element in the compound to get the total contribution.
                contrib += delta * count

            # If the compound is a reactant, append its contribution to the electron row.
            if comp in reactants_data:
                electron_row.append(contrib)
            # If the compound is a product, append the negative of its contribution to the electron row.
            else:
                electron_row.append(-contrib)

        # Append the electron row to the matrix.
        # This row represents the electron transfer balance equation:
        # (Sum of electrons lost by reactants) - (Sum of electrons gained by products) = 0
        matrix.append(electron_row)

        return matrix


    def balance_equation(self, reactants_data, products_data):
        """
        Balance the chemical equation by solving the matrix for the nullspace and finding the coefficients.

        Args:
            reactants_data (list): Parsed data for reactants.
            products_data (list): Parsed data for products.

        Returns:
            tuple: A tuple containing:
                - react_coeffs (list): Coefficients for the reactants.
                - prod_coeffs (list): Coefficients for the products.
        Example:
            Consider the equation: "H2 + O2 -> H2O"

            Step 1: Build the matrix.
            - Matrix represents element counts and charge balance.

            Step 2: Find the nullspace of the matrix.
            - Nullspace gives the stoichiometric coefficients.

            Step 3: Scale the solution to integer coefficients.
            - Example solution: [2, 1, 2] (for H2, O2, H2O).

            Step 4: Split the solution into reactant and product coefficients.
            - react_coeffs: [2, 1] (for H2, O2)
            - prod_coeffs: [2] (for H2O)

            Step 5: Return the coefficients.
        """
        # Build the matrix representing the system of linear equations for the chemical equation.
        # The matrix is constructed using the `build_matrix` method, which takes the parsed data
        # for reactants and products and returns a matrix where each row corresponds to an element
        # and each column corresponds to a compound.
        matrix = self.build_matrix(reactants_data, products_data)

        # Convert the matrix to a SymPy Matrix object.
        # This allows us to use SymPy's linear algebra capabilities, such as finding the nullspace.
        sympy_matrix = Matrix(matrix)

        # Find the nullspace of the matrix.
        # The nullspace represents the set of solutions to the homogeneous system of equations (Ax = 0).
        # For balancing chemical equations, the nullspace gives the stoichiometric coefficients.
        nullspace = sympy_matrix.nullspace()  # Find the nullspace of the matrix

        # Check if the nullspace is empty.
        # If no nullspace is found, it means there is no valid solution to balance the equation.
        if not nullspace:
            return None  # No solution found

        # Extract the first vector from the nullspace.
        # The nullspace may contain multiple vectors, but for balancing chemical equations,
        # we typically use the first vector (smallest non-trivial solution).
        solution = nullspace[0]

        # Extract the denominators of the rational numbers in the solution.
        # This is necessary to convert the solution to integer coefficients.
        denominators = [val.q for val in solution]  # Extract denominators from the solution

        # Compute the least common multiple (LCM) of the denominators.
        # The LCM is used to scale the solution so that all coefficients are integers.
        lcm = 1
        for d in denominators:
            if d != 0:
                # Update the LCM using the formula: LCM(a, b) = (a * b) // GCD(a, b)
                lcm = lcm * d // self.gcd(lcm, d)  # Compute the least common multiple (LCM)

        # Scale the solution by the LCM to convert it to integer coefficients.
        # The `abs` function ensures that all coefficients are positive.
        solution = [abs(val * lcm) for val in solution]  # Scale the solution by the LCM

        # Split the solution into coefficients for reactants and products.
        # The first `len(reactants_data)` values correspond to the reactants.
        # The remaining values correspond to the products.
        react_coeffs = solution[:len(reactants_data)]  # Coefficients for reactants
        prod_coeffs = solution[len(reactants_data):]  # Coefficients for products

        # Return the coefficients for reactants and products.
        # These coefficients can be used to construct the balanced chemical equation.
        return react_coeffs, prod_coeffs

    def construct_balanced_equation(self, reactant_formulas, product_formulas, react_coeffs, prod_coeffs):
        """
        Construct the balanced chemical equation from the coefficients and formulas.

        Args:
            reactant_formulas (list): List of reactant formulas.
            product_formulas (list): List of product formulas.
            react_coeffs (list): Coefficients for the reactants.
            prod_coeffs (list): Coefficients for the products.

        Returns:
            str: The balanced chemical equation as a string.

        Example:
            Consider the coefficients and formulas:
            - react_coeffs: [2, 1] (for H2, O2)
            - prod_coeffs: [2] (for H2O)
            - reactant_formulas: ["H2", "O2"]
            - product_formulas: ["H2O"]

            Step 1: Format reactant terms.
            - "2H2", "O2"

            Step 2: Format product terms.
            - "2H2O"

            Step 3: Combine into a balanced equation.
            - "2H2 + O2 -> 2H2O"
        """
        # Initialize an empty list to store the formatted terms for the reactants.
        # Each term will represent a reactant compound with its stoichiometric coefficient.
        reactant_terms = []

        # Iterate over the reactant coefficients and their corresponding formulas.
        # `zip(react_coeffs, reactant_formulas)` pairs each coefficient with its formula.
        for coeff, formula in zip(react_coeffs, reactant_formulas):
            # Check if the coefficient is zero.
            # A zero coefficient indicates an error in balancing, as it implies the compound is not part of the reaction.
            if coeff == 0:
                raise ValueError("Zero coefficient detected; balancing failed.")

            # Format the term for the reactant:
            # - If the coefficient is 1, omit it for simplicity (e.g., "H2" instead of "1H2").
            # - If the coefficient is greater than 1, include it (e.g., "2H2").
            term = f"{coeff if coeff != 1 else ''}{formula}"  # Omit coefficient if it's 1

            # Append the formatted term to the list of reactant terms.
            reactant_terms.append(term)

        # Initialize an empty list to store the formatted terms for the products.
        # Each term will represent a product compound with its stoichiometric coefficient.
        product_terms = []

        # Iterate over the product coefficients and their corresponding formulas.
        # `zip(prod_coeffs, product_formulas)` pairs each coefficient with its formula.
        for coeff, formula in zip(prod_coeffs, product_formulas):
            # Check if the coefficient is zero.
            # A zero coefficient indicates an error in balancing, as it implies the compound is not part of the reaction.
            if coeff == 0:
                raise ValueError("Zero coefficient detected; balancing failed.")

            # Format the term for the product:
            # - If the coefficient is 1, omit it for simplicity (e.g., "H2O" instead of "1H2O").
            # - If the coefficient is greater than 1, include it (e.g., "2H2O").
            term = f"{coeff if coeff != 1 else ''}{formula}"  # Omit coefficient if it's 1

            # Append the formatted term to the list of product terms.
            product_terms.append(term)

        # Combine the reactant and product terms into a balanced chemical equation.
        # The reactants and products are joined with " + ", and the two sides are separated by " -> ".
        # Example: "2H2 + O2 -> 2H2O"
        return " + ".join(reactant_terms) + " -> " + " + ".join(product_terms)
    def get_balanced_equation(self, equation):
        """
        Get the balanced chemical equation from the input equation string.

        Args:
            equation (str): The unbalanced chemical equation (e.g., "H2 + O2 -> H2O").

        Returns:
            str: The balanced chemical equation.

        Raises:
            ValueError: If the equation cannot be balanced.

        Example:
            Input: "H2 + O2 -> H2O"
            Output: "2H2 + O2 -> 2H2O"

        """
        # Parse the input chemical equation into structured data.
        # The `parse_equation` method splits the equation into reactants and products,
        # and further parses each compound into its constituent elements and charges.
        # It returns two tuples:
        # - The first tuple contains parsed data for reactants and products.
        # - The second tuple contains the original formulas for reactants and products.
        (reactants_data, products_data), (reactant_formulas, product_formulas) = self.parse_equation(equation)

        # Balance the chemical equation using the parsed data.
        # The `balance_equation` method calculates the stoichiometric coefficients
        # for the reactants and products to satisfy the law of conservation of mass.
        result = self.balance_equation(reactants_data, products_data)

        # Check if the balancing process was successful.
        # If `result` is `None`, it means the equation could not be balanced.
        if not result:
            raise ValueError("Failed to balance the equation.")

        # Unpack the result into coefficients for reactants and products.
        # The `balance_equation` method returns a tuple containing:
        # - `react_coeffs`: A list of coefficients for the reactants.
        # - `prod_coeffs`: A list of coefficients for the products.
        react_coeffs, prod_coeffs = result

        # Construct the balanced chemical equation using the coefficients and formulas.
        # The `construct_balanced_equation` method formats the reactants and products
        # into a readable string representation of the balanced equation.
        balanced = self.construct_balanced_equation(reactant_formulas, product_formulas, react_coeffs, prod_coeffs)

        # Return the balanced chemical equation as a string.
        # Example: "2H2 + O2 -> 2H2O"
        return balanced




    def is_acid_or_base(self, compound):
        """
        Determines if a compound is an acid, base, or neutral.
        Handles common acids, bases, and edge cases.

        Args:
            compound (str): The chemical formula of the compound (e.g., "HCl", "NaOH", "H2O").

        Returns:
            str: "Acid", "Base", or "Neutral".

        Example:
            Consider the compound "HCl":
            - It starts with "H" and has no charge, so it is identified as an acid.
            - Output: "Acid"

            Consider the compound "NaOH":
            - It contains "OH", so it is identified as a base.
            - Output: "Base"

            Consider the compound "H2O":
            - It does not match any acid or base patterns, so it is identified as neutral.
            - Output: "Neutral"
        """
        # Parse the compound into its constituent elements and charge.
        elements_dict, charge = self.parse_compound(compound)

        # Check for common acid patterns.
        if self._is_acid(compound, elements_dict, charge):
            return "Acid"

        # Check for common base patterns.
        if self._is_base(compound, elements_dict, charge):
            return "Base"

        # If neither acid nor base, return neutral.
        return "Neutral"
    def _is_acid(self, compound, elements_dict, charge):
        """
        Checks if the compound is an acid.

        Args:
            compound (str): The chemical formula of the compound.
            elements_dict (dict): A dictionary mapping element symbols to their counts in the compound.
            charge (int): The net charge of the compound.

        Returns:
            bool: True if the compound is an acid, False otherwise.

        Example:
            Consider the compound "HCl":
            - It starts with "H" and has no charge, so it is identified as an acid.
            - Output: True

            Consider the compound "CH3COOH":
            - It contains "COOH", so it is identified as an acid.
            - Output: True

            Consider the compound "H2SO4":
            - It is in the list of common acid polyatomics, so it is identified as an acid.
            - Output: True
        """
        # Acids typically start with H (e.g., HCl, H2SO4).
        if 'H' in elements_dict and elements_dict['H'] > 0 and charge == 0:
            return True

        # Check for organic acids (e.g., CH3COOH).
        if 'COOH' in compound or 'CO2H' in compound:
            return True

        # Check for polyatomic ions that are acids (e.g., H3PO4, HNO3).
        common_acid_polyatomics = ['H3PO4', 'HNO3', 'H2SO4', 'HClO4']
        if compound in common_acid_polyatomics:
            return True

        # If none of the above conditions are met, the compound is not an acid.
        return False
    def _is_base(self, compound, elements_dict, charge):
        """
        Checks if the compound is a base.

        Args:
            compound (str): The chemical formula of the compound.
            elements_dict (dict): A dictionary mapping element symbols to their counts in the compound.
            charge (int): The net charge of the compound.

        Returns:
            bool: True if the compound is a base, False otherwise.

        Example:
            Consider the compound "NaOH":
            - It contains "OH", so it is identified as a base.
            - Output: True

            Consider the compound "NH3":
            - It is ammonia, so it is identified as a base.
            - Output: True

            Consider the compound "CO3^2-":
            - It contains "CO3", so it is identified as a base.
            - Output: True
        """
        # Bases typically contain OH (e.g., NaOH, KOH).
        if 'OH' in compound:
            return True

        # Check for ammonia and amines (e.g., NH3, CH3NH2).
        if compound == 'NH3' or 'NH2' in compound:
            return True

        # Check for polyatomic ions that are bases (e.g., CO3^2-, OH-).
        common_base_polyatomics = ['CO3', 'OH', 'HCO3']
        for poly in common_base_polyatomics:
            if poly in compound:
                return True

        # If none of the above conditions are met, the compound is not a base.
        return False
    

    def calculate_molar_mass(self, compound):
        """
        Calculate the molar mass of a chemical compound.

        Args:
            compound (str): The chemical formula of the compound (e.g., "H2O", "C6H12O6").

        Returns:
            float: The molar mass of the compound in grams per mole (g/mol).

        Example:
            Consider the compound "H2O":
            - H: 2 atoms * 1.008 g/mol = 2.016 g/mol
            - O: 1 atom * 15.999 g/mol = 15.999 g/mol
            - Total molar mass = 2.016 + 15.999 = 18.015 g/mol

            Output: 18.015
        """
        # Parse the compound into its constituent elements and their counts.
        elements_dict, _ = self.parse_compound(compound)

        # Initialize the total molar mass to 0.
        total_molar_mass = 0.0

        # Iterate over each element in the compound.
        for element, count in elements_dict.items():
            # Find the element in the elements_list to get its molar mass.
            element_data = next((elem for elem in self.elements_list if elem['symbol'] == element), None)
            if element_data:
                # Multiply the element's molar mass by its count in the compound.
                total_molar_mass += element_data['molar_mass'] * count
            else:
                # If the element is not found, raise an error.
                raise ValueError(f"Element '{element}' not found in the elements list.")

        # Return the total molar mass.
        return total_molar_mass