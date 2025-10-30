import numpy as np
import cvxpy as cp
from cvxpygen import cpg
import time
import pickle

#from linearize_constraint import linearize_constraint

import sympy as sym

def linearize_constraint(constraint, nx, nu, np, var_type, var_ind):
    #LINEARIZE_CONSTRAINT Summary of this function goes here
    #   Detailed explanation goes here

    t_sym = sym.symbols("t");
    x_sym = sym.MatrixSymbol("x", nx, 1);
    u_sym = sym.MatrixSymbol("u", nu, 1);
    p_sym = sym.MatrixSymbol("p", np, 1);

    if var_type == "x":
        x_selected = [x_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(x_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x_ref[:, k], u[:, k], p) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (x[var_ind, k] - x_ref[var_ind, k]);
    elif var_type == "u":
        u_selected = [u_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(u_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x[:, k], u_ref[:, k], p) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (u[var_ind, k] - u_ref[var_ind, k]);
    elif var_type == "p":
        p_selected = [p_sym[i, 0] for i in var_ind]
        constraint_jacobian = sym.lambdify([t_sym, x_sym, u_sym, p_sym], sym.Matrix([constraint(t_sym, x_sym, u_sym, p_sym)]).jacobian(p_selected));
        linearized_constraint = lambda t, x, u, p, x_ref, u_ref, p_ref, k : constraint(t, x[:, k], u[:, k], p_ref) + constraint_jacobian(t, x_ref[:, k], u_ref[:, k], p_ref) * (p[var_ind, k] - p_ref[var_ind, k]);    
    
    return linearized_constraint




#from Code.Discretization.discretize import discretize_dynamics_ZOH

save = False;
load = False;

if not load:
    # Input from matlab
    x_ref = np.array(x_ref)
    u_ref = np.array(u_ref)

    x_0 = np.array(x_0).reshape(-1, 1)
    x_f = np.array(x_f).reshape(-1, 1)

    # Problem parameters
    # T_min, T_max, tau_max, alpha, glideslope_max_angle, gimbal_max_angle, pitch_max, angvel_max, time_min_max_thrust, max_gimbal_rate
    T_min = np.array(params[0])
    T_max = np.array(params[1])
    tau_max = np.array(params[2])
    alpha = np.array(params[3])
    glideslope_max_angle = np.array(params[4])
    gimbal_max_angle = np.array(params[5])
    pitch_max = np.array(params[6])
    angvel_max = np.array(params[7])
    time_min_max_thrust = np.array(params[8])
    max_gimbal_rate = np.array(params[9])
    S_x = np.array(params[10 : 24]).T
    S_u = np.array(params[24 : 28]).T
    c_x = np.array(params[28 : 42]).T
    c_u = np.array(params[42 : 46]).T
    N = int(N);
    delta_t = np.array(delta_t)

    # PTR parameters
    w_vc = np.array(w_vc)
    w_tr = np.array(w_tr)

    ### Discretization
    A_k = np.array(A_k)
    B_k_minus = np.array(B_k_minus)
    B_k_plus = np.array(B_k_plus)
    c_k = np.array(c_k)
else:
    loaded_data = np.load(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_input_FOH.npz', allow_pickle=True)
    x_ref = loaded_data["x_ref"]
    u_ref = loaded_data["u_ref"]

    x_0 = loaded_data["x_0"]
    x_f = loaded_data["x_f"]

    # Problem parameters
    # T_min, T_max, alpha, glideslope_max_angle, gimbal_max_angle
    T_min = loaded_data["T_min"]
    T_max = loaded_data["T_max"]
    tau_max = loaded_data["tau_max"]
    alpha = loaded_data["alpha"]
    glideslope_max_angle = loaded_data["glideslope_max_angle"]
    gimbal_max_angle = loaded_data["gimbal_max_angle"]
    pitch_max = loaded_data["pitch_max"]
    angvel_max = loaded_data["angvel_max"]
    time_min_max_thrust = loaded_data["time_min_max_thrust"]
    max_gimbal_rate = loaded_data["max_gimbal_rate"]
    N = loaded_data["N"]
    delta_t = loaded_data["delta_t"]
    S_x = loaded_data["S_x"]
    S_u = loaded_data["S_u"]
    c_x = loaded_data["c_x"]
    c_u = loaded_data["c_u"]

    # PTR parameters
    w_vc = loaded_data["w_vc"]# / 1000
    w_tr = loaded_data["w_tr"]

    # Discretization
    A_k = loaded_data["A_k"]
    B_k_minus = loaded_data["B_k_minus"]
    B_k_plus = loaded_data["B_k_plus"]
    c_k = loaded_data["c_k"]

    # Saved Problem
    dbfile = open(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_prob_FOH', 'rb')    
    problem = pickle.load(dbfile)

    dbfile.close()

t_k = np.linspace(0, N * delta_t, N + 1)

nx = 14
nu = 4
n_r = 3

def virtual_control_cost_func(w_vc, V, v_0, v_N):
    #J_vc = w_vc * (cp.sum(cp.sum_squares(V)) + cp.sum(cp.sum_squares(v_0)) + cp.sum(cp.sum_squares(v_N)))
    J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.norm(v_0, 1) + cp.norm(v_N, 1))
    return J_vc

def trust_region_cost_func(w_tr, eta):
    J_tr = w_tr @ eta.T
    return J_tr

def trust_region_cost_func_explicit(w_tr, x, x_ref, u, u_ref):
    J_tr = (cp.sum(cp.square(x - x_ref), 0) + cp.sum(cp.square(u - u_ref), 0)) @ w_tr.T
    return J_tr


def virtual_control_cost_func_np(w_vc, V, v_0, v_N):
    J_vc = w_vc * (np.sum(np.linalg.norm(V, 1, 0)) + np.linalg.norm(v_0, 1) + np.linalg.norm(v_N, 1))
    return J_vc

def scale_x(x, S_x, c_x):
    x_scaled = (x - c_x) / S_x
    return x_scaled 

def unscale_x(x, S_x, c_x):
    x_unscaled = cp.multiply(S_x, x) + c_x
    return x_unscaled

def unscale_u(u, S_u, c_u):
    u_unscaled = cp.multiply(S_u, u) + c_u
    return u_unscaled

def fine_gimbal(U, u_ref_param):
    fine_gimbal_constraint = cp.square(u_ref_param[1, :]) + cp.square(u_ref_param[2, :]) + cp.multiply(2 * u_ref_param[1, :], (U[1, :] - u_ref_param[1, :])) + cp.multiply(2 * u_ref_param[2, :], (U[2, :] - u_ref_param[2, :]));
    return fine_gimbal_constraint

def thrust_rate(U, u_ref_param):
    thrust_rate_constraint = cp.abs((cp.norm2(u_ref_param[0:2, k + 1]) + u_ref_param[0:2, k + 1].T / cp.norm2(u_ref_param[0:2, k + 1]) * (U[0:2, k + 1] - u_ref_param[0:2, k + 1])) - (cp.norm2(u_ref_param[0:2, k]) + u_ref_param[0:2, k].T / cp.norm2(u_ref_param[0:2, k]) * (U[0:2, k] - u_ref_param[0:2, k]))) / delta_t
    return thrust_rate_constraint

def gimbal_rate(U, u_ref_param):
    gimbal_rate_constraint = -(u_ref_param[0:2, k].T * u_ref_param[0:2, k + 1] + u_ref_param[0:2, k].T * (U[0:2, k + 1] - u_ref_param[0:2, k + 1]) + u_ref_param[0:2, k + 1].T * (U[0:2, k] - u_ref_param[0:2, k])) + (cp.norm2(u_ref_param[0:2, k + 1]) * cp.norm2(u_ref_param[0:2, k]) + cp.norm2(u_ref_param[0:2, k]) * U[0:2, k + 1].T / cp.norm2(u_ref_param[0:2, k + 1]) * (U[0:2, k + 1] - u_ref_param[0:2, k + 1]) + cp.norm2(u_ref_param[0:2, k + 1]) * u_ref_param[0:2, k].T / cp.norm2(u_ref_param[0:2, k]) * (U[0:2, k] - u_ref_param[0:2, k])) * np.cos(max_gimbal_rate * delta_t)
    return gimbal_rate_constraint

## Save input data so it can be debugged in python
#if save:
#    np.savez(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_input_FOH.npz', A_k = A_k, B_k_minus = B_k_minus, B_k_plus = B_k_plus, c_k = c_k, w_vc = w_vc, w_tr = w_tr, T_min = T_min, T_max = T_max, tau_max = tau_max, alpha = alpha, glideslope_max_angle = glideslope_max_angle, gimbal_max_angle = gimbal_max_angle, pitch_max = pitch_max, angvel_max = angvel_max, time_min_max_thrust = time_min_max_thrust, max_gimbal_rate = max_gimbal_rate, delta_t = delta_t, N = N, x_0 = x_0, x_f = x_f, x_ref = x_ref, u_ref = u_ref, S_x = S_x, S_u = S_u, c_x = c_x, c_u = c_u)
    
#    #dbfile = open(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_prob_FOH', 'ab')
    
#    ## source, destination
#    #pickle.dump(problem, dbfile)                    
#    #dbfile.close()

## Create Convex Problem
# Define variables
if str(type(problem)) == "<class 'array.array'>" or str(type(problem)) == "<class 'numpy.ndarray'>":
    #print("Creating problem")
    X = cp.Variable((nx, N + 1), name='X')
    U = cp.Variable((nu, N + 1), name='U')
    #eta = cp.Variable((1, N + 1), name = 'eta')
    V = cp.Variable((nx, N + 0), name = 'V')
    v_0 = cp.Variable((nx, 1), name = 'v_0')
    v_N = cp.Variable((nx - 1, 1), name = 'v_N')
    #v_prime = cp.Variable((, 1), name = "v_prime")

    # Define parameters
    #Ak_param = cp.Parameter((nx, nx, N + 0), name = 'Ak')
    #Bk_param = cp.Parameter((nx, nu, N + 0), name = 'Bk')
    ck_param = cp.Parameter((nx, N + 0), name = 'ck')
    x_0_param = cp.Parameter((nx, 1), name = 'x_0')
    x_f_param = cp.Parameter((nx - 1, 1), name = 'x_f')

    x_ref_param = cp.Parameter((nx, N + 1), name = 'x_ref')
    u_ref_param = cp.Parameter((nu, N + 1), name = 'u_ref')
    p_ref_param = 0 # Haven't implemented parameters for this

    # Define objective
    objective = cp.Minimize(-(S_x[13] * X[13, -1] + c_x[13]) / x_0[13, -1] + cp.sum(cp.abs(S_u[3] * U[3, :] + c_u[3])))
    virtual_control_cost = cp.Minimize(virtual_control_cost_func(w_vc, V, v_0, v_N))
    #trust_region_cost = cp.Minimize(trust_region_cost_func(w_tr, eta)) # 
    trust_region_cost = cp.Minimize(trust_region_cost_func_explicit(w_tr, X, x_ref_param, U, u_ref_param))
    augmented_objective = objective + virtual_control_cost + trust_region_cost

    # Define constraints
    m_i = 6

    Ak_params = []
    Bk_minus_params = []
    Bk_plus_params = []
    dynamics_constraints = []
    for k in range(N):
        Ak_params.append(cp.Parameter((nx, nx), name = "Ak_" + str(k)))
        Bk_plus_params.append(cp.Parameter((nx, nu), name = "Bk_plus_" + str(k)))
        Bk_minus_params.append(cp.Parameter((nx, nu), name = "Bk_minus_" + str(k)))

        dynamics_constraints.append(X[:, k + 1] == scale_x(Ak_params[k] @ unscale_x(X[:, k], S_x, c_x) + Bk_minus_params[k] @ unscale_u(U[:, k], S_u, c_u) + Bk_plus_params[k] @ unscale_u(U[:, k + 1], S_u, c_u) + ck_param[:, k] + V[:, k], S_x, c_x))


    #dynamics_constraints = []
    #for k in range(N):
    #    dynamics_constraints.append(X[:, k + 1] == Ak_param[:, :, k] @ X[:, k] + Bk_param[:, :, k] @ U[:, k] + ck_param[:, 0, k] + V[:, k])

    # Linearize constraints
    #min_thrust_constraint = lambda t, x, u, p : T_min ^ 2 - (u[0] ** 2 + u[1] ** 2 + u[2] ** 2)
    #min_thrust_constraint_linearized = linearize_constraint(min_thrust_constraint, nx, nu, 0, "u", np.arange(3))
    
    #pitch_constraint = lambda t, x, u, p : 1 - (2 * (x[6] * x[8] - X[7] * X[9])) ** 2 - np.sin(pitch_max) ** 2
    #pitch_constraint_linearized = linearize_constraint(pitch_constraint, nx, nu, 0, "x", np.arange(6, 10))

    #gimbal_constraint = lambda t, x, u, p : (u[1] ** 2 + u[2] ** 2)
    #gimbal_constraint_linearized = linearize_constraint(gimbal_constraint, nx, nu, 0, "u", np.arange(1, 3))


    # Min thrust, max thrust, glideslope, max gimbal, max roll torque, max angular velocity, thrust rate
    constraints = dynamics_constraints + [ # Dynamics constraint
                   T_min <= (S_u[0] * U[0, :] + c_u[0]), # Thrust min constraint
                   (S_u[0] * U[0, :] + c_u[0]) <= T_max, # Thrust max constraint
                   cp.norm2(cp.multiply(np.reshape(S_x[0:3], [3, 1]), X[0:3, :]) + np.reshape(c_x[0:3], [3, 1]), 0) - (S_x[2] * X[2, :] + c_x[2]) / np.cos(glideslope_max_angle) <= 0, # Glideslope constraint
                   #cp.abs((cp.multiply(np.reshape(S_u[1:3], [2, 1]), U[1:3, :]) + np.reshape(c_u[1:3], [2, 1]))) <= gimbal_max_angle, # Rough gimbal constraint (needed?)
                   cp.norm2(cp.multiply(np.reshape(S_u[1:3], [2, 1]), U[1:3, :]) + np.reshape(c_u[1:3], [2, 1]), 0) <= gimbal_max_angle,
                   #fine_gimbal(cp.multiply(np.reshape(S_u, [4, 1]), U) + np.reshape(c_u, [4, 1]), cp.multiply(np.reshape(S_u, [4, 1]), u_ref_param) + np.reshape(c_u, [4, 1])) <= gimbal_max_angle ** 2, #gimbal_constraint_linearized <= gimbal_max_angle ** 2, # Fine gimbal constraint
                   cp.abs((S_u[3] * U[3, :] + c_u[3])) <= tau_max, # Max roll torque constraint
                   cp.norm_inf((cp.multiply(np.reshape(S_x[10:13], [3, 1]), X[10:13, :]) + np.reshape(c_x[10:13], [3, 1]))) <= angvel_max, # Angular velocity constraint
                   #pitch_constraint_linearized(t, x, u, p, x_ref_param, u_ref_param, p_ref_param, k) <= 0, # Pitch constraint
                   cp.abs(S_u[0] * (U[0, 1 : N] - U[0, N - 1])) / delta_t <= (T_max - T_min) / time_min_max_thrust, # Thrust rate
                   #thrust_rate(cp.multiply(np.reshape(S_u, [4, 1]), U) + np.reshape(c_u, [4, 1]), cp.multiply(np.reshape(S_u, [4, 1]), u_ref_param) + np.reshape(c_u, [4, 1])) <= (T_max - T_min) / time_min_max_thrust, # Thrust rate constraint
                   #gimbal_rate(cp.multiply(np.reshape(S_u, [4, 1]), U) + np.reshape(c_u, [4, 1]), cp.multiply(np.reshape(S_u, [4, 1]), u_ref_param) + np.reshape(c_u, [4, 1])) <= 0, # Gimbal rate constraint
                   (cp.multiply(S_x, X[:, 0]) + c_x) + v_0 == x_0_param.flatten(), # Initial condition constraint
                   (cp.multiply(S_x[0:-1], X[0:-1, N]) + c_x[0:-1]) + v_N == x_f_param.flatten()] # Terminal condition constraint
                   #cp.sum(cp.sum_squares(X - x_ref_param)) + cp.sum(cp.sum_squares(U - u_ref_param)) <= eta]

    # Define problem
    problem = cp.Problem(augmented_objective, constraints)



# Save input data so it can be debugged in python
if save:
    np.savez(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_input_FOH.npz', A_k = A_k, B_k_minus = B_k_minus, B_k_plus = B_k_plus, c_k = c_k, w_vc = w_vc, w_tr = w_tr, T_min = T_min, T_max = T_max, tau_max = tau_max, alpha = alpha, glideslope_max_angle = glideslope_max_angle, gimbal_max_angle = gimbal_max_angle, pitch_max = pitch_max, angvel_max = angvel_max, time_min_max_thrust = time_min_max_thrust, max_gimbal_rate = max_gimbal_rate, delta_t = delta_t, N = N, x_0 = x_0, x_f = x_f, x_ref = x_ref, u_ref = u_ref, S_x = S_x, S_u = S_u, c_x = c_x, c_u = c_u)
    
    dbfile = open(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Stochastic-SCP-Rocket-Landing\Deterministic 6DoF CVXPY\Deterministic 6DoF CVXPY\Deterministic_6DoF_fixed_prob_FOH', 'ab')
    
    # source, destination
    pickle.dump(problem, dbfile)                    
    dbfile.close()


# Set parameters
#Ak_param.value = A_k[:, :, 0:N]
#Bk_param.value = B_k[:, :, 0:N]
#ck_param.value = c_k[:, :, 0:N]

# Set parameters
for k in range(N):
    problem.param_dict["Ak_" + str(k)].value = A_k[:, :, k]
    problem.param_dict["Bk_minus_" + str(k)].value = B_k_minus[:, :, k]
    problem.param_dict["Bk_plus_" + str(k)].value = B_k_plus[:, :, k]

problem.param_dict["ck"].value = c_k[:, 0, 0:N]

problem.param_dict["x_0"].value = x_0
problem.param_dict["x_f"].value = x_f

problem.param_dict["x_ref"].value = x_ref
problem.param_dict["u_ref"].value = u_ref

# Solve
#t0 = time.time()
#val = problem.solve(solver = "QOCO", verbose = False)#solver = "ECOS")
#t1 = time.time()
##print('\nCVXPY Clarabel\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))
# val = problem.solve(solver = "QOCO", verbose = False)#solver = "ECOS")
# 
# t0 = time.time()
# val = problem.solve(solver = "QOCO", verbose = False)#solver = "ECOS")
# t1 = time.time()
# print('\nCVXPY Clarabel\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))

# t0 = time.time()
# val = problem.solve(solver = "osqp", verbose = False)#solver = "ecos")
# t1 = time.time()
#print('\ncvxpy osqp\nsolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))


#t0 = time.time()
#val = problem.solve(solver = "OSQP", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY OSQP\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))

#t0 = time.time()
#val = problem.solve(solver = "OSQP", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY OSQP\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))


#t0 = time.time()
#val = problem.solve(solver = "ECOS", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY ECOS\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))

#t0 = time.time() 
#val = problem.solve(solver = "QOCO", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY QOCO\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val)) # don't know why QOCO is so slow for this... :(


#t0 = time.time() 
#val = problem.solve(solver = "QOCO", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY QOCO\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val)) # don't know why QOCO is so slow for this... :(


#cpg.generate_code(problem, code_dir='Deterministic_6DoF_fixed_FOH_ECOS', solver = "ECOS")
#cpg.generate_code(problem, code_dir='Deterministic_6DoF_fixed_FOH_QOCO_N5', solver = "QOCO")
##cpg.generate_code(problem, code_dir='Ast2Ast_fixed_FOH_PTR_SCS', solver = "SCS")


#from Deterministic_6DoF_fixed_FOH_ECOS.cpg_solver import cpg_solve
from Deterministic_6DoF_fixed_FOH_QOCO.cpg_solver import cpg_solve
# 
problem.register_solve('cpg', cpg_solve)
# 
# 
t0 = time.time()
val = problem.solve(method='cpg')
t1 = time.time()
print('\ncvxpy ecos_gen \nsolve time: %.3f ms with %.3f and %.5f ms solve' % (1000 * (t1 - t0), val, 1000 * problem.solution.attr["solve_time"]))

# Extract Solution
X_sol = problem.var_dict['X'].value
U_sol = problem.var_dict['U'].value

obj_val = -(S_x[13] * X_sol[13, -1] + c_x[13]) / x_0[13, -1] + np.sum(np.abs(S_u[3] * U_sol[3, :] + c_u[3]))

eta = np.sum(np.square(X_sol - x_ref), 0) + np.sum(np.square(U_sol - u_ref), 0) #problem.var_dict['eta'].value
V = problem.var_dict['V'].value
v_0 = problem.var_dict['v_0'].value
v_N = problem.var_dict['v_N'].value

solve_status = problem.status

#dyn_err = np.zeros([nx, N])
#dyn_err_ref = np.zeros([nx, N])
#for k in range(N):
#    dyn_err[:, k] = X_sol[:, k + 1] - (A_k[:, :, k] @ X_sol[:, k] + B_k[:, :, k] @ U_sol[:, k] + c_k[:, 0, k])
#    dyn_err_ref[:, k] = x_ref[:, k + 1] - (A_k[:, :, k] @ x_ref[:, k] + B_k[:, :, k] @ u_ref_param[:, k] + c_k[:, 0, k])


# Debugs
#print(problem.var_dict['v_N'].value)
#print(virtual_control_cost_func_np(w_vc, problem.var_dict['V'].value, problem.var_dict['v_0'].value, problem.var_dict['v_N'].value))
#print(trust_region_cost_func(w_tr, problem.var_dict['eta'].value))
a = 1;