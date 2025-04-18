classdef StochasticProblem
    %STOCHASTICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        xf 
        P0
        Pf
        N
        Nu
        n % Has .x, .u, .p, .cvx, .ncvx, .y, .w
        tf
        u_hold string {mustBeMember(u_hold, ["ZOH"])} = "ZOH"
        guess % Has to have values .x, .u, .p
        cont % Originally .f and .G then after prob.linearize() it has .A, .B, .E, .c
        disc % Has .A_k, .B_k, .E_k, .c_k, .G_k, .C_k, .D_k, .Ptilde_minus_k, .Ptilde_k, .L_k
        filter % Has f_0, g_0, .C
        stoch % Has .w, .v, .delta_t
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
        function obj = StochasticProblem(x0, xf, P0, Pf, N, u_hold, tf, f, G, f_0, g_0, guess, convex_constraints, objective, options)
            arguments
                x0
                xf
                P0
                Pf
                N
                u_hold
                tf
                f
                G
                f_0
                g_0
                guess % Has to have values .x, .u, .p, .X_k, .S_k
                convex_constraints % Cell array of constraint functions @(x, u, p)
                objective % Has to be @(x, p)
                options.initial_bc = @(x, p) x - x0 % Has to be @(x, p)
                options.terminal_bc = @(x, p) x - xf % Has to be @(x, p)
                options.integration_tolerance = 1e-12
                options.scale = true
                options.nonconvex_constraints = [] % Cell array of constraint functions @(x, u, p, x_ref, u_ref, p_ref)
                options.w = @(n) randn([size(G(0, x0, guess.u(:, 1), guess.p), 2), n])
                options.v = @(n) randn([numel(f_0(0, x0, guess.u(:, 1), guess.p)), n])
                options.delta_t = 1e0
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
            obj.n.y = numel(f_0(0, x0, guess.u(:, 1), guess.p));
            obj.tf = tf;
            obj.u_hold = u_hold;
            obj.guess = guess;
            obj.cont.f = f;
            obj.cont.G = G;
            obj.filter.f_0 = f_0;
            obj.filter.g_0 = g_0;
            obj = linearize(obj);
            obj.convex_constraints = convex_constraints;
            obj.nonconvex_constraints = options.nonconvex_constraints;
            obj.initial_bc = options.initial_bc;
            obj.terminal_bc = options.terminal_bc;
            obj.objective = objective;
            obj.scale = options.scale;
            obj.scaling = obj.compute_scaling();
            obj.tolerances = odeset(RelTol=options.integration_tolerance, AbsTol=options.integration_tolerance);
            obj.stoch.w = options.w;
            obj.stoch.v = options.v;
            obj.stoch.delta_t = options.delta_t;
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

            prob.disc.C_k = zeros(prob.n.y, prob.n.x);
            prob.disc.D_k = zeros(prob.n.y, prob.n.y);
            for k = 1:prob.N
                prob.disc.C_k(:, :, k) = prob.filter.C(t_k(k), x_ref(:, k), u_ref(:, min(k, prob.Nu)), p_ref); 
                prob.disc.D_k(:, :, k) = prob.filter.g_0(t_k(k), x_ref(:, k), u_ref(:, min(k, prob.Nu)), p_ref);
            end

            % Precompute Kalman Filter Gain and 
            [prob.disc.Ptilde_minus_k, prob.disc.Ptilde_k, prob.disc.L_k] = compute_Kalman_matrices_apriori(prob.P0, prob.disc.A_k, prob.disc.G_k, prob.disc.C_k, prob.disc.D_k);
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

        function [t_cont, x_cont, u_cont] = cont_prop_without_feedback_control(prob, u_ref, p, tspan)
            %CONT_PROP_WITHOUT_FEEDBACK_CONTROL Summary of this function goes here
            %   Detailed explanation goes here
            t_k = linspace(0, prob.tf, prob.N);

            u_func = @(t, x) interp1(t_k(1:prob.Nu), u_ref', t, "previous", "extrap")';

            [t_cont, x_cont] = sode45(prob.cont.f, prob.cont.G, u_func, p, w, tspan, prob.stoch.delta_t, prob.x0, prob.tolerances);

            u_cont = u_func(t_cont(1:(numel(t_cont) - 1)));
        end

        function [t_cont, x_cont, u_cont] = cont_prop_feedback_no_kalman_filter(prob, x_ref, u_ref, p, K_k, N_sub)
            %CONT_PROP Summary of this function goes here
            %   Detailed explanation goes here
            t_k = linspace(0, prob.tf, prob.N);

            [t_cont, x_cont, u_cont] = propagate_cont_feedback_no_kalman_filter(prob.x0, x_ref, u_ref, K_k, prob.cont.f, prob.cont.G, t_k, N_sub, prob.stoch.w, prob.stoch.delta_t, prob.tolerances);
        end

        function [t_cont, x_cont, xhat_cont, Phat_cont, u_cont] = cont_prop(prob, x_ref, u_ref, p, K, options)
            arguments
                prob
                x_ref
                u_ref
                p
                K
                options.x_0 = [prob.sample_initial_condition(), prob.sample_initial_condition()]
                options.w_k_func = prob.create_w_func()
                options.v_k = prob.stoch.v(prob.N)
                options.N_sub = 15
            end
            %CONT_PROP Summary of this function goes here
            %   Detailed explanation goes here
            t_k = linspace(0, prob.tf, prob.N);

            [t_cont, x_cont, xhat_cont, Phat_cont, u_cont] = propagate_cont_feedback_kalman_filter(options.x_0, prob.P0, p, prob.cont.f, prob.cont.G, prob.cont.A, prob.cont.B, prob.cont.c, x_ref, u_ref, K, prob.disc.L_k, prob.disc.C_k, prob.disc.D_k, prob.filter.f_0, prob.filter.g_0, t_k, [0, prob.tf], options.N_sub, options.w_k_func, options.v_k, prob.tolerances);
        end

        function [t_k, x_disc, xhat_disc, Phat_disc, u_disc] = disc_prop(prob, x_ref, u_ref, p, K, options)
            arguments
                prob
                x_ref
                u_ref
                p
                K
                options.x_0 = [prob.sample_initial_condition(), prob.sample_initial_condition()]
                options.w_k_func = prob.create_w_func()
                options.v_k = prob.stoch.v(prob.N)
            end
            %DISC_PROP Summary of this function goes here
            %   Detailed explanation goes here
            x_disc = zeros([prob.n.x, prob.N]);
            x_disc(:, 1) = options.x_0(:, 1);
            xhat_disc = x_disc;
            xhat_disc(:, 1) = options.x_0(:, 2);

            u_disc = zeros([prob.n.u, prob.N - 1]);

            Phat_disc = zeros([prob.n.x, prob.n.x, prob.N]);
            Phat_disc(:, :, 1) = prob.P0;

            t_k = linspace(0, prob.tf, prob.N);
            w_k_func = options.w_k_func;
            v_k = options.v_k;

            %u_ref_func = @(t) interp1(t_k(1:prob.Nu), u_ref', t, "previous", "extrap")';
            %K_func = @(t) interp1(t_k(1:prob.Nu), K', t, "previous", "extrap")';

            for k = 1:(prob.N - 1)
                % Compute control
                u_disc(:, k) = u_ref(:, k) + K(:, :, k) * (xhat_disc(:, k) - x_ref(:, k));

                % Time update
                % True state
                u_func = @(t, x) u_disc(:, k);%u_ref_func(t) + K_func(t) * (xhat_disc(:, k) - x_ref(t));
                [~, x_cont_k] = sode45(prob.cont.f,prob.cont.G, u_func, p, prob.stoch.w, [t_k(k), t_k(k + 1)], prob.stoch.delta_t, x_disc(:, k), prob.tolerances, w_k_func = w_k_func);
                x_disc(:, k + 1) = x_cont_k(:, end);
                
                % Estimated state
                xhat_disc(:, k + 1) = prob.disc.A_k(:, :, k) * xhat_disc(:, k) ...
                             + prob.disc.B_k(:, :, k) * u_disc(:, k) ...
                             + zero_if_empty(prob.disc.E_k(:, :, k) * p) ...
                             + prob.disc.c_k(:, :, k);
                Phat_disc(:, :, k + 1) = covariance_time_update(prob.disc.A_k(:, :, k), Phat_disc(:, :, k), prob.disc.G_k(:, :, k));

                % Measurement update
                y_k = prob.filter.f_0(t_k(k + 1), x_disc(:, k + 1), u_disc(:, k)) + prob.filter.g_0(xhat_disc(:, k + 1), u_disc(:, k)) * v_k(:, k + 1);
                ytilde_minus_k = innovation_process(y_k, prob.disc.C_k(:, :, k + 1), xhat_disc(:, k + 1));
                xhat_disc(:, k + 1) = estimate_measurement_update(xhat_disc(:, k + 1), prob.disc.L_k(:, :, k + 1), ytilde_minus_k);
                Phat_disc(:, :, k + 1) = covariance_measurement_update(prob.disc.L_k(:, :, k + 1), prob.disc.C_k(:, :, k + 1), Phat_disc(:, :, k + 1), prob.disc.D_k(:, :, k + 1));
            end
        end

        function [x_0] = sample_initial_condition(prob)
            x_0 = prob.x0 + chol(prob.P0, "lower") * randn([prob.n.x, 1]);
        end

        function [w_func] = create_w_func(prob)
            t_k = 0:prob.stoch.delta_t:prob.tf;
            w_k = prob.stoch.w(numel(t_k));

            w_func = @(t) interp1(t_k, w_k', t, "previous", "extrap")';
        end
    end
    methods(Static)
        function [stoch_prob] = stochastify_discrete_problem(disc_prob, G, f_0, g_0, P0, Pf, options)
            arguments
                disc_prob
                G
                f_0
                g_0
                P0
                Pf
                options.w = @(n) randn([size(G(0, disc_prob.x0, disc_prob.guess.u(:, 1), disc_prob.guess.p), 2), n])
                options.v = @(n) randn([numel(f_0(0, disc_prob.x0, disc_prob.guess.u(:, 1), disc_prob.guess.p)), n])
                options.delta_t = 1e0
                options.sol = []
            end
            
            if ~isempty(options.sol)
                if numel(size(options.sol.x)) == 3 % Check if ptr_sol which stores solutions from all iterations
                    stoch_guess.x = options.sol.x(:, :, options.sol.converged_i);
                    stoch_guess.u = options.sol.u(:, :, options.sol.converged_i);
                    stoch_guess.p = options.sol.p(:, options.sol.converged_i);
                else
                    stoch_guess = options.sol;
                end
            else
                stoch_guess = disc_prob.guess;
            end

    
            stoch_prob = StochasticProblem(disc_prob.x0, disc_prob.xf, P0, Pf, disc_prob.N, disc_prob.u_hold, disc_prob.tf, ...
                disc_prob.cont.f, G, f_0, g_0, stoch_guess, disc_prob.convex_constraints, ...
                disc_prob.objective, initial_bc = disc_prob.initial_bc, terminal_bc = disc_prob.terminal_bc, ...
                integration_tolerance = disc_prob.tolerances.AbsTol, scale = disc_prob.scale, ...
                nonconvex_constraints = disc_prob.nonconvex_constraints, w = options.w, v = options.v, delta_t = options.delta_t);
        end
    end
end

