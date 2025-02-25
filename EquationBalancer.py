from collections import defaultdict
import re
from sympy import Matrix

class EquationBalancer:
    """
    A class to balance chemical equations by parsing compounds, building a matrix of element counts,
    solving for stoichiometric coefficients, and constructing the balanced equation.
    """

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

    def parse_compound(self, compound):
        """
        Parse a chemical compound into its constituent elements and their counts, including handling charges.

        Args:
            compound (str): The chemical formula of the compound (e.g., "H2O", "Fe^3+", "SO4^2-").

        Returns:
            tuple: A tuple containing:
                - elements_dict (defaultdict): A dictionary mapping each element to its count in the compound.
                - charge (int): The net charge of the compound.
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

                # Append the parsed data for this compound to the `parsed_data` list.
                parsed_data.append({'elements': elements_dict, 'charge': charge})

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

        Returns:
            list: A matrix where each row represents an element (or charge), and each column represents a compound.
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




elements = [
    {"symbol": "H",  "name": "Hydrogen",    "nox": [+1, -1], "electrons_to_neutral": 1},
    {"symbol": "He", "name": "Helium",      "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Li", "name": "Lithium",     "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Be", "name": "Beryllium",   "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "B",  "name": "Boron",       "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "C",  "name": "Carbon",      "nox": [+4, -4], "electrons_to_neutral": 4},
    {"symbol": "N",  "name": "Nitrogen",    "nox": [+5, -3], "electrons_to_neutral": 3},
    {"symbol": "O",  "name": "Oxygen",      "nox": [-2],     "electrons_to_neutral": 2},
    {"symbol": "F",  "name": "Fluorine",    "nox": [-1],     "electrons_to_neutral": 1},
    {"symbol": "Ne", "name": "Neon",        "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Na", "name": "Sodium",      "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Mg", "name": "Magnesium",   "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Al", "name": "Aluminum",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Si", "name": "Silicon",     "nox": [+4, -4], "electrons_to_neutral": 4},
    {"symbol": "P",  "name": "Phosphorus",  "nox": [+5, -3], "electrons_to_neutral": 3},
    {"symbol": "S",  "name": "Sulfur",      "nox": [+6, -2], "electrons_to_neutral": 2},
    {"symbol": "Cl", "name": "Chlorine",    "nox": [+7, -1], "electrons_to_neutral": 1},
    {"symbol": "Ar", "name": "Argon",       "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "K",  "name": "Potassium",   "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Ca", "name": "Calcium",     "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Sc", "name": "Scandium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Ti", "name": "Titanium",    "nox": [+4, +3], "electrons_to_neutral": 4},
    {"symbol": "V",  "name": "Vanadium",    "nox": [+5, +4], "electrons_to_neutral": 5},
    {"symbol": "Cr", "name": "Chromium",    "nox": [+6, +3], "electrons_to_neutral": 6},
    {"symbol": "Mn", "name": "Manganese",   "nox": [+7, +2], "electrons_to_neutral": 7},
    {"symbol": "Fe", "name": "Iron",        "nox": [+3, +2], "electrons_to_neutral": 3},
    {"symbol": "Co", "name": "Cobalt",      "nox": [+3, +2], "electrons_to_neutral": 3},
    {"symbol": "Ni", "name": "Nickel",      "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Cu", "name": "Copper",      "nox": [+2, +1], "electrons_to_neutral": 2},
    {"symbol": "Zn", "name": "Zinc",        "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Ga", "name": "Gallium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Ge", "name": "Germanium",   "nox": [+4, -4], "electrons_to_neutral": 4},
    {"symbol": "As", "name": "Arsenic",     "nox": [+5, -3], "electrons_to_neutral": 3},
    {"symbol": "Se", "name": "Selenium",    "nox": [+6, -2], "electrons_to_neutral": 2},
    {"symbol": "Br", "name": "Bromine",     "nox": [+7, -1], "electrons_to_neutral": 1},
    {"symbol": "Kr", "name": "Krypton",     "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Rb", "name": "Rubidium",    "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Sr", "name": "Strontium",   "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Y",  "name": "Yttrium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Zr", "name": "Zirconium",   "nox": [+4],     "electrons_to_neutral": 4},
    {"symbol": "Nb", "name": "Niobium",     "nox": [+5],     "electrons_to_neutral": 5},
    {"symbol": "Mo", "name": "Molybdenum",  "nox": [+6],     "electrons_to_neutral": 6},
    {"symbol": "Tc", "name": "Technetium",  "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Ru", "name": "Ruthenium",   "nox": [+8, +3], "electrons_to_neutral": 8},
    {"symbol": "Rh", "name": "Rhodium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Pd", "name": "Palladium",   "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Ag", "name": "Silver",      "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Cd", "name": "Cadmium",     "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "In", "name": "Indium",      "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Sn", "name": "Tin",         "nox": [+4, +2], "electrons_to_neutral": 4},
    {"symbol": "Sb", "name": "Antimony",    "nox": [+5, -3], "electrons_to_neutral": 3},
    {"symbol": "Te", "name": "Tellurium",   "nox": [+6, -2], "electrons_to_neutral": 2},
    {"symbol": "I",  "name": "Iodine",      "nox": [+7, -1], "electrons_to_neutral": 1},
    {"symbol": "Xe", "name": "Xenon",       "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Cs", "name": "Cesium",      "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Ba", "name": "Barium",      "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "La", "name": "Lanthanum",   "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Ce", "name": "Cerium",      "nox": [+4, +3], "electrons_to_neutral": 4},
    {"symbol": "Pr", "name": "Praseodymium","nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Nd", "name": "Neodymium",   "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Pm", "name": "Promethium",  "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Sm", "name": "Samarium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Eu", "name": "Europium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Gd", "name": "Gadolinium",  "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Tb", "name": "Terbium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Dy", "name": "Dysprosium",  "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Ho", "name": "Holmium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Er", "name": "Erbium",      "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Tm", "name": "Thulium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Yb", "name": "Ytterbium",   "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Lu", "name": "Lutetium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Hf", "name": "Hafnium",     "nox": [+4],     "electrons_to_neutral": 4},
    {"symbol": "Ta", "name": "Tantalum",    "nox": [+5],     "electrons_to_neutral": 5},
    {"symbol": "W",  "name": "Tungsten",    "nox": [+6],     "electrons_to_neutral": 6},
    {"symbol": "Re", "name": "Rhenium",     "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Os", "name": "Osmium",      "nox": [+8],     "electrons_to_neutral": 8},
    {"symbol": "Ir", "name": "Iridium",     "nox": [+4],     "electrons_to_neutral": 4},
    {"symbol": "Pt", "name": "Platinum",    "nox": [+4, +2], "electrons_to_neutral": 4},
    {"symbol": "Au", "name": "Gold",        "nox": [+3, +1], "electrons_to_neutral": 3},
    {"symbol": "Hg", "name": "Mercury",     "nox": [+2, +1], "electrons_to_neutral": 2},
    {"symbol": "Tl", "name": "Thallium",    "nox": [+3, +1], "electrons_to_neutral": 3},
    {"symbol": "Pb", "name": "Lead",        "nox": [+4, +2], "electrons_to_neutral": 4},
    {"symbol": "Bi", "name": "Bismuth",     "nox": [+5, +3], "electrons_to_neutral": 5},
    {"symbol": "Po", "name": "Polonium",    "nox": [+4, +2], "electrons_to_neutral": 4},
    {"symbol": "At", "name": "Astatine",    "nox": [-1],     "electrons_to_neutral": 1},
    {"symbol": "Rn", "name": "Radon",       "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Fr", "name": "Francium",    "nox": [+1],     "electrons_to_neutral": 1},
    {"symbol": "Ra", "name": "Radium",      "nox": [+2],     "electrons_to_neutral": 2},
    {"symbol": "Ac", "name": "Actinium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Th", "name": "Thorium",     "nox": [+4],     "electrons_to_neutral": 4},
    {"symbol": "Pa", "name": "Protactinium","nox": [+5],     "electrons_to_neutral": 5},
    {"symbol": "U",  "name": "Uranium",     "nox": [+6],     "electrons_to_neutral": 6},
    {"symbol": "Np", "name": "Neptunium",   "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Pu", "name": "Plutonium",   "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Am", "name": "Americium",   "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Cm", "name": "Curium",      "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Bk", "name": "Berkelium",   "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Cf", "name": "Californium", "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Es", "name": "Einsteinium", "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Fm", "name": "Fermium",     "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Md", "name": "Mendelevium", "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "No", "name": "Nobelium",    "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Lr", "name": "Lawrencium",  "nox": [+3],     "electrons_to_neutral": 3},
    {"symbol": "Rf", "name": "Rutherfordium","nox": [+4],    "electrons_to_neutral": 4},
    {"symbol": "Db", "name": "Dubnium",     "nox": [+5],     "electrons_to_neutral": 5},
    {"symbol": "Sg", "name": "Seaborgium",  "nox": [+6],     "electrons_to_neutral": 6},
    {"symbol": "Bh", "name": "Bohrium",     "nox": [+7],     "electrons_to_neutral": 7},
    {"symbol": "Hs", "name": "Hassium",     "nox": [+8],     "electrons_to_neutral": 8},
    {"symbol": "Mt", "name": "Meitnerium",  "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Ds", "name": "Darmstadtium","nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Rg", "name": "Roentgenium", "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Cn", "name": "Copernicium", "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Nh", "name": "Nihonium",   "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Fl", "name": "Flerovium",   "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Mc", "name": "Moscovium",   "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Lv", "name": "Livermorium", "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Ts", "name": "Tennessine",  "nox": [0],      "electrons_to_neutral": 0},
    {"symbol": "Og", "name": "Oganesson",   "nox": [0],      "electrons_to_neutral": 0},
]
