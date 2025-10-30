import sympy as sym

def linearize_constraint(constraint, nx, nu, np, var_type, var_ind):
    #LINEARIZE_CONSTRAINT Summary of this function goes here
    #   Detailed explanation goes here

    t_sym = sym.symbols("t");
    x_sym = sym.MatrixSymbol("x", nx, 1);
    u_sym = sym.MatrixSymbol("u", nu, 1);
    p_sym = sym.MatrixSymbol("p", np, 1);

    if var_type is "x":
        x_selected = [x_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(x_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x_ref[:, k], u[:, k], p) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (x[var_ind, k] - x_ref[var_ind, k]);
    elif var_type is "u":
        u_selected = [u_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(u_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x[:, k], u_ref[:, k], p) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (u[var_ind, k] - u_ref[var_ind, k]);
    elif var_type is "p":
        p_selected = [p_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(p_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x[:, k], u[:, k], p_ref) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (p[var_ind, k] - p_ref[var_ind, k]);    
    
    return linearized_constraint

