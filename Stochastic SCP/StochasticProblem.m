classdef StochasticProblem
    %DETERMINISTICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        xf 
        P0
        Pf
        N
        Nu
        n % Has .x, .u, .p, .cvx, .ncvx, .w, .y
        tf
        u_hold string {mustBeMember(u_hold, ["ZOH"])} = "ZOH"
        guess % Has to have values .x, .u, .p
        cont % Originally .f and .G then after prob.linearize() it has .A, .B, .E, .c
        disc % Has .A_k, .B_k, .E_k, .c_k, .G_k, .C_k, .D_k, .Ptilde_minus_k, .Ptilde_k, .L_k
        filter % Has f_0, g_0, .C, .D
        stoch % Has .w, .delta_t
        convex_constraints % Cell array of convex constraint functions @(x, u, p)
        nonconvex_constraints % Cell array of nonconvex constraint functions @(x, u, p, x_ref, u_ref, p_ref)
        initial_bc % Has to be @(x, p)
        terminal_bc % Has to be @(x, p)
        objective % Has to be @(x, u, p)
        scale
        scaling
        sol
        tolerances
    end
    
    methods
        function obj = DeterministicProblem(x0, xf, P0, Pf, N, u_hold, tf, f, guess, convex_constraints, objective, options)
            arguments
                x0
                xf
                P0
                Pf
                N
                u_hold
                tf
                f
                guess % Has to have values .x, .u, .p
                convex_constraints % Cell array of constraint functions @(x, u, p)
                objective % Has to be @(x, p)s
                options.initial_bc = @(x, p) x - x0 % Has to be @(x, p)
                options.terminal_bc = @(x, p) x - xf % Has to be @(x, p)
                options.integration_tolerance = 1e-12
                options.scale = true
                options.nonconvex_constraints = [] % Cell array of constraint functions @(x, u, p, x_ref, u_ref, p_ref)
                options.w
                options.v
            end
            %DETERMINISTICPROBLEM Construct an instance of this class
            %   Detailed explanation goes here

            obj.x0 = x0;
            obj.xf = xf;
            obj.P0 = P0;
            obj.Pf = Pf;
            obj.N = N;
            obj.Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;
            obj.n.x = numel(x0);
            obj.n.u = size(guess.u, 1);
            obj.n.p = size(guess.p, 1);
            obj.n.cvx = numel(convex_constraints);
            obj.n.ncvx = numel(options.nonconvex_constraints);
            obj.tf = tf;
            obj.u_hold = u_hold;
            obj.guess = guess;
            obj.cont.f = f;
            obj = linearize(obj);
            obj.convex_constraints = convex_constraints;
            obj.nonconvex_constraints = options.nonconvex_constraints;
            obj.initial_bc = options.initial_bc;
            obj.terminal_bc = options.terminal_bc;
            obj.objective = objective;
            obj.scale = options.scale;
            obj.scaling = obj.compute_scaling();
            obj.tolerances = odeset(RelTol=options.integration_tolerance, AbsTol=options.integration_tolerance);
        end

        function prob = linearize(prob)
            %LINEARIZE 

            t_sym = sym("t");
            x_sym = sym("x", [prob.n.x, 1]);
            u_sym = sym("u", [prob.n.u, 1]);
            p_sym = sym("p", [prob.n.p, 1]);
            
            % Linearize Dynamics
            prob.cont.A = matlabFunction(jacobian(prob.cont.f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
            prob.cont.B = matlabFunction(jacobian(prob.cont.f(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
            prob.cont.E = matlabFunction(jacobian(prob.cont.f(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

            prob.cont.c = @(t, x, u, p) prob.cont.f(t, x, u, p) - prob.cont.A(t, x, u, p) * x - prob.cont.B(t, x, u, p) * u - zero_if_empty(prob.cont.E(t, x, u, p) * p);

            % Linearize Measurement
            prob.filter.C = matlabFunction(jacobian(prob.filter.f_0(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
            prob.filter.D = matlabFunction(jacobian(prob.filter.g_0(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

            % Linearize Nonconvex Constraints (Needed??)

            % Linearize Boundary Conditions (6DoF Quaternion)
            
        end
        
        function [prob, Delta] = discretize(prob, x_ref, u_ref, p_ref)
            %DISCRETIZE Summary of this method goes here
            %   Detailed explanation goes here
            
            % Discretize Dynamics
            [prob.disc.A_k, prob.disc.B_k, prob.disc.E_k, prob.disc.c_k, prob.disc.G_k, Delta] = discretize_stochastic_dynamics_ZOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.cont.G, prob.N, [0, prob.tf], x_ref, u_ref, p_ref, prob.tolerances);

            % Discretize Measurement
            t_k = linspace(0, prob.tf, prob.N);
            prob.disc.C_k = prob.filter.C(t_k, x_ref, u_ref, p_ref); 
            prob.disc.D_k = prob.filter.D(t_k, x_ref, u_ref, p_ref);

            % Precompute Kalman Filter Gain and 
            [prob.disc.Ptilde_minus_k, prob.disc.Ptilde_k, prob.disc.L_k] = compute_Kalman_matrices_apriori(prob.disc.A_k, prob.disc.G_k, prob.disc.C_k, prob.disc.D_k);
        end

        function [scaling] = compute_scaling(obj)
            z_ub = 1;
            z_lb = 0;

            x_max = max(obj.guess.x, [], 2);
            u_max = max(obj.guess.u, [], 2);
            p_max = max(obj.guess.p);

            x_min = min(obj.guess.x, [], 2);
            u_min = min(obj.guess.u, [], 2);
            p_min = min(obj.guess.p);

            if ~obj.scale
                scaling.S_x = eye(obj.n.x);%diag(make_not_zero(x_max - x_min) / (z_ub - z_lb));
                scaling.S_u = eye(obj.n.u);%diag(make_not_zero(u_max - u_min) / (z_ub - z_lb));
                scaling.S_p = eye(obj.n.p);%diag(make_not_zero(p_max - p_min) / (z_ub - z_lb));
    
                scaling.c_x = zeros([obj.n.x, 1]);%x_min - scaling.S_x * ones([obj.n.x, 1]) * z_lb;
                scaling.c_u = zeros([obj.n.u, 1]);%u_min - scaling.S_u * ones([obj.n.u, 1]) * z_lb;
                scaling.c_p = zeros([obj.n.p, 1]);%p_min - scaling.S_p * ones([obj.n.p, 1]) * z_lb;
            else
                scaling.S_x = diag(make_not_zero(x_max - x_min) / (z_ub - z_lb));
                scaling.S_u = diag(make_not_zero(u_max - u_min) / (z_ub - z_lb));
                scaling.S_p = diag(make_not_zero(p_max - p_min) / (z_ub - z_lb));
    
                scaling.c_x = x_min - scaling.S_x * ones([obj.n.x, 1]) * z_lb;
                scaling.c_u = u_min - scaling.S_u * ones([obj.n.u, 1]) * z_lb;
                scaling.c_p = p_min - scaling.S_p * ones([obj.n.p, 1]) * z_lb;
            end
            
            function [not_zero] = make_not_zero(maybe_zero)
                % If number is too close to zero, make it 1 so that the
                % scaling matrix stays invertible

                tol = 1e-2;
                not_zero = (maybe_zero < tol) + (maybe_zero >= tol) .* maybe_zero;
            end
        end

        function [x] = unscale_x(prob, xhat)
            if numel(size(xhat)) == 3
                x = pagemtimes(prob.scaling.S_x, xhat) + prob.scaling.c_x;
            else
                x = prob.scaling.S_x * xhat + repmat(prob.scaling.c_x, 1, size(xhat, 2));
            end
        end

        function [u] = unscale_u(prob, uhat)
            if numel(size(uhat)) == 3
                u = pagemtimes(prob.scaling.S_u, uhat) + prob.scaling.c_u;
            else
                u = prob.scaling.S_u * uhat + repmat(prob.scaling.c_u, 1, size(uhat, 2));
            end
        end

        function [p] = unscale_p(prob, phat)
            if numel(size(phat)) == 2
                if isempty(phat)
                    p = phat;
                else
                    p = prob.scaling.S_p * phat + prob.scaling.c_p;
                end
            else
                p = prob.scaling.S_p * phat + prob.scaling.c_p;
            end
        end

        function [xhat] = scale_x(prob, x)
            xhat = pagemldivide(prob.scaling.S_x, x - prob.scaling.c_x);
        end

        function [uhat] = scale_u(prob, u)
            uhat = pagemldivide(prob.scaling.S_u, u - prob.scaling.c_u);
        end

        function [phat] = scale_p(prob, p)
            phat = prob.scaling.S_p \ (p - prob.scaling.c_p); % shape????
        end

        function [t_cont, x_cont, xhat_cont, Phat_cont, u_cont] = cont_prop(prob, x_ref, u_ref, p, K)
            %DISC_PROP Summary of this function goes here
            %   Detailed explanation goes here
            t_k = linspace(0, prob.tf, prob.N);
            x_ref_func = @(t) linspace(t_k, x_ref, t);

            w_k = w(numel(t_k));
            v_k = v(numel(t_k));

            u_ref_func = @(t) interp1(t_k(1:prob.Nu), u_ref', t, "previous", "extrap")';
            K_func = @(t) interp1(t_k(1:prob.Nu), K', t, "previous", "extrap")';
            u_func = @(t, xhat) u_ref_func(t) + K_func(t) * (xhat - x_ref_func(t));

            N_sub = 15;

            [t_cont, x_cont, xhat_cont, Phat_cont] = propagate_cont_kalman_filter(prob.x0, prob.P0, u_func, prob.cont.f, prob.cont.G, A_ref, B_ref, c_ref, G_ref, prob.disc.L_k, prob.disc.C_k, prob.disc.D_k, prob.stoch.f_0, prob.stoch.g_0, t_k, [0, prob.tf], N_sub, w_k, v_k, prob.tolerances);

            u_cont = u_func(t_cont(1:(numel(t_cont) - 1)));
        end

        function [t_k, x_disc, xhat_disc, Phat_disc, u_disc] = disc_prop(prob, x_ref, u_ref, p, K)
            %DISC_PROP Summary of this function goes here
            %   Detailed explanation goes here
            x_disc = zeros([prob.n.x, prob.N]);
            x_disc(:, 1) = prob.x0;
            xhat_disc = x_disc;

            u_disc = zeros([prob.n.u, prob.N - 1]);

            Phat_disc = zeros([prob.n.x, prob.n.x, prob.N - 1]);
            Phat_disc(:, :, 1) = prob.P0;

            t_k = linspace(0, prob.tf, prob.N);
            w_k = w(numel(t_k));
            v_k = v(numel(t_k));

            %u_ref_func = @(t) interp1(t_k(1:prob.Nu), u_ref', t, "previous", "extrap")';
            %K_func = @(t) interp1(t_k(1:prob.Nu), K', t, "previous", "extrap")';

            for k = 1:(prob.N - 1)
                % Compute control
                u_disc(:, k) = u_ref(:, k) + K * (xhat_disc(:, k) - x_ref(:, k));

                % Time update
                % True state
                u_func = @(t, x) u_disc(:, k);%u_ref_func(t) + K_func(t) * (xhat_disc(:, k) - x_ref(t));
                [~, x_cont_k] = sode45(@(t, x) prob.cont.f, @(t, x) prob.cont.G, u_func, p, prob.stoch.w, [t_k(k), t_k(k + 1)], prob.stoch.delta_t, x_disc(:, k), prob.tolerances);
                x_disc(:, k + 1) = x_cont_k(end, :)';
                
                % Estimated state
                xhat_disc(:, k + 1) = prob.disc.A_k(:, :, k) * xhat_disc(:, k) ...
                             + prob.disc.B_k(:, :, k) * u_disc(:, k) ...
                             + zero_if_empty(prob.disc.E_k(:, :, k) * p) ...
                             + prob.disc.c_k(:, :, k);
                Phat_disc(:, :, k + 1) = covariance_time_update(prob.disc.A_k(:, :, k), Phat_disc(:, :, k), prob.disc.G_k(:, :, k));

                % Measurement update
                y_k = f_0(t, x_disc(:, k), u_disc(:, k)) + g_0(x_m, u_disc(:, k)) * v_k;
                ytilde_minus_k = innovation_process(y_k, C_k(k), xhat_disc(:, :, k + 1));
                xhat_disc(:, k + 1) = estimate_measurement_update(xhat_disc(:, k + 1), L_k(k), ytilde_minus_k);
                Phat_disc(:, :, k + 1) = covariance_measurement_update(L_k(k), C_k(k), Phat_disc(:, :, k + 1), D_k(k));
            end
        end


    end
end

