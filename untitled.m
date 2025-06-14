%% Call Convex Optimizer
[X_py, U_py] = pyrunfile("MarsLanding3D_3DoF.py", ["X_sol", "U_sol"], x_0 = x_0, A_k = A_k, B_k = B_k, c_k = c_k);
X = double(X_py);
U = double(U_py);