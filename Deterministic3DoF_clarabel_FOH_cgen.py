import numpy as np
import cvxpy as cp
from cvxpygen import cpg
import time
import pickle


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
    # T_min, T_max, alpha, glideslope_max_angle, gimbal_max_angle
    T_min = np.array(params[0])
    T_max = np.array(params[1])
    alpha = np.array(params[2])
    glideslope_max_angle = np.array(params[3])
    gimbal_max_angle = np.array(params[4])
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
    loaded_data = np.load(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Applied-Control-in-Astronautics\HW6\Parts\Q1\saved_3DoF_input_FOH.npz', allow_pickle=True)
    x_ref = loaded_data["x_ref"]
    u_ref = loaded_data["u_ref"]

    x_0 = loaded_data["x_0"]
    x_f = loaded_data["x_f"]

    # Problem parameters
    # T_min, T_max, alpha, glideslope_max_angle, gimbal_max_angle
    T_min = loaded_data["T_min"]
    T_max = loaded_data["T_max"]
    alpha = loaded_data["alpha"]
    glideslope_max_angle = loaded_data["glideslope_max_angle"]
    gimbal_max_angle = loaded_data["gimbal_max_angle"]
    N = loaded_data["N"]
    delta_t = loaded_data["delta_t"]

    # PTR parameters
    w_vc = loaded_data["w_vc"] / 1000
    w_tr = loaded_data["w_tr"]

    # Discretization
    A_k = loaded_data["A_k"]
    B_k_minus = loaded_data["B_k_minus"]
    B_k_plus = loaded_data["B_k_plus"]
    c_k = loaded_data["c_k"]

    # Saved Problem
    dbfile = open(r'C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 590ACA\AAE590ACA-Applied-Control-in-Astronautics\HW6\Parts\Q1\det_3DoF_prob_FOH', 'rb')    
    problem = pickle.load(dbfile)

    dbfile.close()

t_k = np.linspace(0, N * delta_t, N + 1)

nx = 7#np.size(x_ref, 1)
nu = 3#np.size(u_ref, 1)
n_r = 2;


# Save input data so it can be debugged in python
if save:
    np.savez("saved_3DoF_input_FOH.npz", A_k = A_k, B_k_minus = B_k_minus, B_k_plus = B_k_plus, c_k = c_k, w_vc = w_vc, w_tr = w_tr, delta_t = delta_t, N = N, gimbal_max_angle = gimbal_max_angle, glideslope_max_angle = glideslope_max_angle, alpha = alpha, T_min = T_min, T_max = T_max, x_0 = x_0, x_f = x_f, x_ref = x_ref, u_ref = u_ref)
    
    dbfile = open('det_3DoF_prob_FOH', 'ab')
    
    # source, destination
    pickle.dump(problem, dbfile)                    
    dbfile.close()

    
def virtual_control_cost_func(w_vc, V, v_0, v_N):
    J_vc = w_vc * (cp.sum(cp.norm(V, 1, 0)) + cp.norm(v_0, 1) + cp.norm(v_N, 1))
    return J_vc

def trust_region_cost_func(w_tr, eta):
    J_tr = w_tr @ eta.T
    return J_tr


def virtual_control_cost_func_np(w_vc, V, v_0, v_N):
    J_vc = w_vc * (np.sum(np.linalg.norm(V, 1, 0)) + np.linalg.norm(v_0, 1) + np.linalg.norm(v_N, 1))
    return J_vc


## Create Convex Problem
# Define variables
if str(type(problem)) == "<class 'array.array'>":
    #print("Creating problem")
    X = cp.Variable((nx, N + 1), name='X')
    U = cp.Variable((nu, N + 1), name='U')
    eta = cp.Variable((1, N + 1), name = 'eta')
    V = cp.Variable((nx, N + 0), name = 'V')
    v_0 = cp.Variable((nx, 1), name = 'v_0')
    v_N = cp.Variable((nx - 1, 1), name = 'v_N')
    v_extra = cp.Variable((2))

    # Define parameters
    #Ak_param = cp.Parameter((nx, nx, N + 0), name = 'Ak')
    #Bk_param = cp.Parameter((nx, nu, N + 0), name = 'Bk')
    ck_param = cp.Parameter((nx, N + 0), name = 'ck')
    x_0_param = cp.Parameter((nx, 1), name = 'x_0')

    x_ref_param = cp.Parameter((nx, N + 1), name = 'x_ref')
    u_ref_param = cp.Parameter((nu, N + 1), name = 'u_ref')

    # Define objective
    objective = cp.Minimize(cp.sum(U[2,:]) * delta_t + cp.sum(v_extra))
    virtual_control_cost = cp.Minimize(virtual_control_cost_func(w_vc, V, v_0, v_N))
    trust_region_cost = cp.Minimize(trust_region_cost_func(w_tr, eta))
    augmented_objective = objective + virtual_control_cost + trust_region_cost

    # Define constraints
    z_lb = lambda t : np.log(np.exp(x_0[6]) - alpha * T_max * t)
    z_lb_k = z_lb(t_k)

    #plt.plot(z_lb_k)
    #plt.show()

    z_i = 6
    sigma_i = 2

    Ak_params = []
    Bk_minus_params = []
    Bk_plus_params = []
    dynamics_constraints = []
    for k in range(N):
        Ak_params.append(cp.Parameter((nx, nx), name = "Ak_" + str(k)))
        Bk_plus_params.append(cp.Parameter((nx, nu), name = "Bk_plus_" + str(k)))
        Bk_minus_params.append(cp.Parameter((nx, nu), name = "Bk_minus_" + str(k)))

        dynamics_constraints.append(X[:, k + 1] == Ak_params[k] @ X[:, k] + Bk_minus_params[k] @ U[:, k] + Bk_plus_params[k] @ U[:, k + 1] + ck_param[:, k])


    #dynamics_constraints = []
    #for k in range(N):
    #    dynamics_constraints.append(X[:, k + 1] == Ak_param[:, :, k] @ X[:, k] + Bk_param[:, :, k] @ U[:, k] + ck_param[:, 0, k] + V[:, k])

    constraints = dynamics_constraints + [ # Dynamics constraint
                   T_min * cp.multiply(np.exp(-z_lb_k), (1 - (X[z_i, :] - z_lb_k) + 0.5 * (X[z_i, :] - z_lb_k) ** 2)) <= cp.minimum(U[sigma_i, :], U[0, :]), # Thrust min constraint
                   U[sigma_i, :] <= T_max * cp.multiply(np.exp(-z_lb_k), (1 - (X[z_i, :] - z_lb_k))), # Thrust min constraint
                   cp.norm2(U[0:2, :], 0) <= U[sigma_i, :], # Lcvx constraint
                   cp.norm2(X[0:2, :], 0) - X[1, :] / np.cos(glideslope_max_angle) <= 0, # Glideslope constraint
                   cp.norm2(U[0:2, :], 0) - U[0, :] / np.cos(gimbal_max_angle) <= 0, # Gimbal constraint
                   np.deg2rad(45) <= X[4, :] + v_extra[0],
                   X[4, :] <= np.deg2rad(135) + v_extra[1],
                   0 <= v_extra,
                   X[:, 0] + v_0 == x_0_param.flatten(), # Initial condition constraint
                   X[0:2, N] + v_N[0:2] == x_f[0:2], X[2:4, N] + v_N[2:4] == x_f[2:4], X[4, N] + v_N[4] == x_f[4], X[5, N] + v_N[5] == x_f[5], # Terminal condition constraint
                   cp.norm(X - x_ref_param, 2, 0) + cp.norm(U - u_ref_param, 2, 0) <= eta]

    # Define problem
    problem = cp.Problem(augmented_objective, constraints)

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

problem.param_dict["x_ref"].value = x_ref
problem.param_dict["u_ref"].value = u_ref

# Solve
#t0 = time.time()
#val = problem.solve(solver = "MOSEK", verbose = False)#solver = "ECOS")
#t1 = time.time()
#print('\nCVXPY MOSEK\nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))

#cpg.generate_code(problem, code_dir='Deterministic_3DoF_FOH_PTR_QOCO2', solver = "QOCO")

#from Deterministic_3DoF_FOH_PTR_QOCO2.cpg_solver import cpg_solve

#problem.register_solve('CPG', cpg_solve)


t0 = time.time()
val = problem.solve(method='CPG')
t1 = time.time()
print('\nCVXPY QOCO_gen \nSolve time: %.3f ms with %.3f' % (1000 * (t1 - t0), val))

# Extract Solution
X_sol = problem.var_dict['X'].value
U_sol = problem.var_dict['U'].value

eta = problem.var_dict['eta'].value
V = problem.var_dict['V'].value
v_0 = problem.var_dict['v_0'].value
v_N = problem.var_dict['v_N'].value

solve_status = problem.status

#dyn_err = np.zeros([nx, N])
#dyn_err_ref = np.zeros([nx, N])
#for k in range(N):
#    dyn_err[:, k] = X_sol[:, k + 1] - (A_k[:, :, k] @ X_sol[:, k] + B_k[:, :, k] @ U_sol[:, k] + c_k[:, 0, k])
#    dyn_err_ref[:, k] = x_ref[:, k + 1] - (A_k[:, :, k] @ x_ref[:, k] + B_k[:, :, k] @ u_ref[:, k] + c_k[:, 0, k])


# Debugs
#print(problem.var_dict['v_N'].value)
#print(virtual_control_cost_func_np(w_vc, problem.var_dict['V'].value, problem.var_dict['v_0'].value, problem.var_dict['v_N'].value))
#print(trust_region_cost_func(w_tr, problem.var_dict['eta'].value))
