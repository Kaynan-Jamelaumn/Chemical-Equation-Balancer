from EquationBalancer import EquationBalancer

def print_balanced_equation(equation: str):
    try:
        balanced_equation = balancer.get_balanced_equation(equation)
        print("Balanced Equation:")
        print(balanced_equation)
    except Exception as e:
        print(f"Error: {e}")



if __name__ == "__main__":
    balancer = EquationBalancer()
    #equation = "Fe^3+ + SO4^2- -> Fe2(SO4)3"
    equation = "Fe^2+ + Cr2O7^2- + H^+ -> Fe^3+ + Cr^3+ + H2O"
    print_balanced_equation(equation)
    equation = 'H2O'
    print(balancer.calculate_molar_mass(equation))