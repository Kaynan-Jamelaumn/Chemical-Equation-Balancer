# EquationBalancer

`EquationBalancer` is a Python class designed to balance chemical equations or detect if it's acid or base.

---

## Features

- **Balances Chemical Equations**: Automatically balances chemical equations, including those with polyatomic ions and charges.
- **Handles Charges**: Supports balancing of ionic equations by accounting for the charge of each species.
- **Error Handling**: Detects and reports errors in balancing, such as invalid input or unsolvable equations.
- **Balances NOX**: Automatically balances nox.

---

## Installation

To use `EquationBalancer`, ensure you have Python 3.x installed. The class requires the `sympy` library for matrix operations. You can install it using pip:

```bash
pip install sympy
```

import the EquationBalancer  class

```bash
from EquationBalancer import EquationBalancer
```

Create an instance of the EquationBalancer class:

```bash
balancer = EquationBalancer()
```

Use the get_balanced_equation method to balance a chemical equation. The method takes a string representing the unbalanced equation and returns the balanced equation as a string.
```bash
equation = "Fe^3+ + SO4^2- -> Fe2(SO4)3" 
#equation = "Fe^2+ + Cr2O7^2- + H^+ -> Fe^3+ + Cr^3+ + H2O" #equation to test that balances nox
balanced_equation = balancer.get_balanced_equation(equation)
print("Balanced Equation:")
print(balanced_equation)
```

Here’s a complete example script that demonstrates how to use the EquationBalancer class:
```bash
from EquationBalancer import EquationBalancer

def print_balanced_equation(equation: str):
    try:
        # Create an instance of EquationBalancer
        balancer = EquationBalancer()
        
        # Get the balanced equation
        balanced_equation = balancer.get_balanced_equation(equation)
        
        # Print the balanced equation
        print("Balanced Equation:")
        print(balanced_equation)
    except Exception as e:
        # Handle any errors that occur during balancing
        print(f"Error: {e}")

if __name__ == "__main__":
    # Example chemical equation
    equation = "Fe^3+ + SO4^2- -> Fe2(SO4)3"

    # Print the balanced equation
    print_balanced_equation(equation)


    print(balancer.is_acid_or_base("HCl"))          # Acid
    print(balancer.is_acid_or_base("H2SO4"))       # Acid
    print(balancer.is_acid_or_base("CH3COOH"))     # Acid (organic acid)
    print(balancer.is_acid_or_base("NaOH"))        # Base
    print(balancer.is_acid_or_base("NH3"))         # Base
    print(balancer.is_acid_or_base("Na2CO3"))      # Base (contains CO3^2-)
    print(balancer.is_acid_or_base("NaCl"))        # Neutral
    print(balancer.is_acid_or_base("H2O"))         # Neutral
```