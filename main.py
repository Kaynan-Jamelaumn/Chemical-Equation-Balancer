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
    equation = "Fe^3+ + SO4^2- -> Fe2(SO4)3"
    print_balanced_equation(equation)