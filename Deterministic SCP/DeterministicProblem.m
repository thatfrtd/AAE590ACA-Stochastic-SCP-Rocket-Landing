classdef DeterministicProblem
    %DETERMINISTICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        N
        Nu
        n
        tf
        u_hold string {mustBeMember(u_hold, ["ZOH", "FOH"])} = "ZOH"
        guess
        cont
        disc
        convex_constraints
        nonconvex_constraints
        initial_bc
        terminal_bc
        objective
        scaling
        sol
        tolerances
    end
    
    methods
        function obj = DeterministicProblem(x0, N, u_hold, tf, f, guess, convex_constraints, objective, options)
            arguments
                x0
                N
                u_hold
                tf
                f
                guess % Has to have values .x, .u, .p
                convex_constraints % Cell array of constraint functions @(x, u, p)
                objective
                options.initial_bc = @(x, p) x - x0 % Has to be @(x, p)
                options.terminal_bc = @(x, p) x - guess.x(:, end) % Has to be @(x, p)
                options.integration_tolerance = 1e-12
                options.nonconvex_constraints = [] % Cell array of constraint functions @(x, u, p, x_ref, u_ref, p_ref)
            end
            %DETERMINISTICPROBLEM Construct an instance of this class
            %   Detailed explanation goes here

            obj.x0 = x0;
            obj.N = N;
            obj.Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;
            obj.n.x = numel(x0);
            obj.n.u = size(guess.u, 1);
            obj.n.p = size(guess.p, 1);
            obj.n.cvx = numel(convex_constraints);
            obj.n.ncvx = numel(nonconvex_constraints);
            obj.tf = tf;
            obj.u_hold = u_hold;
            obj.guess = guess;
            obj.cont.f = f;
            obj = linearize(obj);
            obj.convex_constraints = convex_constraints;
            obj.nonconvex_constraints = nonconvex_constraints;
            obj.initial_bc = initial_bc;
            obj.terminal_bc = terminal_bc;
            obj.objective = objective;
            obj.scaling = obj.compute_scaling();
            obj.tolerances = odeset(RelTol=options.integration_tolerance, AbsTol=options.integration_tolerance);
        end

        function prob = linearize(prob)
            %LINEARIZE 

            x_sym = sym("x");
            u_sym = sym("u");
            p_sym = sym("p");
            
            % Linearize Dynamics
            prob.cont.A = matlabFunction(jacobian(prob.cont.f(x_sym, u_sym, p_sym), x_sym));
            prob.cont.B = matlabFunction(jacobian(prob.cont.f(x_sym, u_sym, p_sym), u_sym));
            prob.cont.E = matlabFunction(jacobian(prob.cont.f(x_sym, u_sym, p_sym), p_sym));

            prob.cont.c = @(x, u, p) prob.cont.f(x, u, p) - A(x, u, p) * x - B(x, u, p) * u - E(x, u, p) * p;

            % Linearize Nonconvex Constraints (Needed??)

            % Linearize Boundary Conditions (6DoF Quaternion)
            
        end
        
        function [prob, Delta] = discretize(prob, x_ref, u_ref, p_ref)
            %DISCRETIZE Summary of this method goes here
            %   Detailed explanation goes here
            
            % Discretize Dynamics
            if prob.u_hold == "ZOH"
                [prob.disc.A_k, prob.disc.B_k, prob.disc.E_k, prob.disc.c_k, Delta] = discretize_dynamics_ZOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, [0, prob.tf], x_ref, u_ref, p_ref, prob.tolerances);
            elseif prob.u_hold == "FOH"
                [prob.disc.A_k, prob.disc.B_plus_k, prob.disc.B_minus_k, prob.disc.E_k, prob.disc.c_k, Delta] = discretize_dynamics_FOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, [0, prob.tf], x_ref, u_ref, p_ref, prob.tolerances);
            end
        end

        function [obj] = compute_scaling(obj)
            z_ub = 1;
            z_lb = 0;

            x_max = max(obj.guess.x, [], 1);
            u_max = max(obj.guess.u, [], 1);
            p_max = max(obj.guess.p);

            x_min = min(obj.guess.x, [], 1);
            u_min = min(obj.guess.u, [], 1);
            p_min = min(obj.guess.p);

            obj.S_x = diag((x_max - x_min) / (z_ub - z_lb));
            obj.S_u = diag((u_max - u_min) / (z_ub - z_lb));
            obj.S_p = diag((p_max - p_min) / (z_ub - z_lb));

            obj.c_x = x_min - S_x * ones([obj.n.x, 1]) * z_lb;
            obj.c_u = u_min - S_u * ones([obj.n.u, 1]) * z_lb;
            obj.c_p = p_min - S_p * ones([obj.n.p, 1]) * z_lb;
        end

        function [xhat] = scale_x(prob, x)
            xhat = pagemtimes(prob.scaling.S_x, x) + prob.scaling.c_x;
        end

        function [uhat] = scale_u(prob, u)
            uhat = pagemtimes(prob.scaling.S_u, u) + prob.scaling.c_u;
        end

        function [phat] = scale_p(prob, p)
            phat = prob.scaling.S_p * p + prob.scaling.c_p; % shape????
        end

        function [x] = unscale_x(prob, xhat)
            x = pagemldivide(prob.scaling.S_x, xhat - prob.scaling.c_x);
        end

        function [u] = unscale_u(prob, uhat)
            u = pagemldivide(prob.scaling.S_u, uhat - prob.scaling.c_u);
        end

        function [p] = unscale_p(prob, phat)
            p = prob.scaling.S_p \ (phat - prob.scaling.c_p); % shape????
        end
    end
end

