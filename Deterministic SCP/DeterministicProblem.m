classdef DeterministicProblem
    %DETERMINISTICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        xf
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
        scale
        scaling
        sol
        tolerances
        phase_transition
        phase_transition_k
        t_k
    end
    
    methods
        function obj = DeterministicProblem(x0, xf, N, u_hold, tf, f, guess, convex_constraints, objective, options)
            arguments
                x0
                xf
                N
                u_hold
                tf
                f
                guess % Has to have values .x, .u, .p
                convex_constraints % Cell array of constraint functions @(t, x, u, p)
                objective
                options.initial_bc = @(x, p, x_ref, p_ref) x - x0 % Has to be @(x, p, x_ref, p_ref)
                options.terminal_bc = @(x, p, x_ref, p_ref) x - xf % Has to be @(x, p, x_ref, p_ref)
                options.integration_tolerance = 1e-12
                options.scale = true
                options.nonconvex_constraints = [] % Cell array of constraint functions @(t, x, u, p, x_ref, u_ref, p_ref)
                options.phase_transition = []
                options.phase_transition_k = []
                options.t_k = []
            end
            %DETERMINISTICPROBLEM Construct an instance of this class
            %   Detailed explanation goes here

            obj.x0 = x0;
            obj.xf = xf;
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
            obj.n.pncvx = 0;
            obj.n.phase = numel(options.phase_transition) + 1;
            if ~isempty(options.phase_transition)
                for e = 1:numel(options.phase_transition)
                    options.phase_transition{e}.n.cvx = numel(options.phase_transition{e}.convex_constraints);
                    options.phase_transition{e}.n.ncvx = numel(options.phase_transition{e}.nonconvex_constraints);

                    obj.n.pncvx = obj.n.pncvx + options.phase_transition{e}.n.ncvx;
                end
            end
            obj.phase_transition = options.phase_transition;
            obj.phase_transition_k = options.phase_transition_k;
            if ~isempty(options.t_k)
                obj.t_k = options.t_k;
            else
                obj.t_k = linspace(0, tf, N);
            end
        end     

        function prob = linearize(prob)
            %LINEARIZE 
            
            t_sym = sym("t");
            x_sym = sym("x", [prob.n.x, 1]);
            u_sym = sym("u", [prob.n.u, 1]);
            p_sym = sym("p", [prob.n.p, 1]);
            
            % Linearize Dynamics
            for k = 1:(prob.N - 1)
                prob.cont.A{k} = matlabFunction(jacobian(prob.cont.f{k}(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
                prob.cont.B{k} = matlabFunction(jacobian(prob.cont.f{k}(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
                prob.cont.E{k} = matlabFunction(jacobian(prob.cont.f{k}(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

                prob.cont.c{k} = @(t, x, u, p) prob.cont.f{k}(t, x, u, p) - prob.cont.A{k}(t, x, u, p) * x - prob.cont.B{k}(t, x, u, p) * u - zero_if_empty(prob.cont.E{k}(t, x, u, p) * p);
            end
        end
        
        function [prob, Delta] = discretize(prob, x_ref, u_ref, p_ref)
            %DISCRETIZE Summary of this method goes here
            %   Detailed explanation goes here
            
            % Discretize Dynamics
            if prob.u_hold == "ZOH"
                [prob.disc.A_k, prob.disc.B_k, prob.disc.E_k, prob.disc.c_k, Delta] = discretize_dynamics_ZOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, prob.t_k, x_ref, u_ref, p_ref, prob.tolerances);
            elseif prob.u_hold == "FOH"
                [prob.disc.A_k, prob.disc.B_plus_k, prob.disc.B_minus_k, prob.disc.E_k, prob.disc.c_k, Delta] = discretize_dynamics_FOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, prob.t_k, x_ref, u_ref, p_ref, prob.tolerances);
            end
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
            if numel(size(x)) == 3
                xhat = pagemldivide(prob.scaling.S_x, x - prob.scaling.c_x);
            else
                xhat = prob.scaling.S_x \ (x - prob.scaling.c_x);
            end
        end

        function [uhat] = scale_u(prob, u)
            uhat = pagemldivide(prob.scaling.S_u, u - prob.scaling.c_u);
        end

        function [phat] = scale_p(prob, p)
            phat = prob.scaling.S_p \ (p - prob.scaling.c_p); % shape????
        end

        function [t_cont, x_cont, u_cont] = cont_prop(prob, u, p, options)
            arguments
                prob
                u
                p
                options.tspan = [0, prob.tf]
                options.x_k = []
            end
            %DISC_PROP Summary of this function goes here
            %   Detailed explanation goes here
            if prob.u_hold == "ZOH"
                u_func = @(t) interp1(prob.t_k(1:prob.Nu), u', t, "previous", "extrap")';
            elseif prob.u_hold == "FOH"
                u_func = @(t) interp1(prob.t_k, u', t)';
            end

            x_cont = prob.x0;
            t_cont = 0;
            for e = 1:prob.n.phase
                if e == 1
                    x_0 = x_cont(:, end);
                    if prob.n.phase == 1
                        k_0 = 1;
                        k_f = prob.N;
                    else
                        k_0 = 1;
                        k_f = prob.phase_transition_k(e);
                    end
                else
                    t_0 = prob.t_k(prob.phase_transition_k(e - 1) + 1);
                    x_0 = prob.phase_transition{e - 1}.transition_func(t_0, x_cont(:, end), options.x_k(:, prob.phase_transition_k(e - 1)), options.x_k(:, prob.phase_transition_k(e - 1) + 1));

                    if e == prob.n.phase
                        k_0 = prob.phase_transition_k(e - 1) + 1;
                        k_f = prob.N;
                    else
                        k_0 = prob.phase_transition_k(e - 1) + 1;
                        k_f = prob.phase_transition_k(e);
                    end
                end

                for k = k_0 : (k_f - 1)
                    tspan_k = [prob.t_k(k) + 1e-12, prob.t_k(k + 1)];
                    if k ~= k_0 
                        x_0 = x_cont(:, end);
                    end
                        [t_cont_k, x_cont_k] = ode45(@(t, x) prob.cont.f{k}(t, x, u_func(t), p), tspan_k, x_0, prob.tolerances);
    
                    if k == 1
                        t_cont = t_cont_k;
                        x_cont = x_cont_k';
                    else
                        t_cont = [t_cont; t_cont_k];
                        x_cont = [x_cont, x_cont_k'];
                    end
                end
            end

            if prob.u_hold == "ZOH"
                u_cont = u_func(t_cont(1:(numel(t_cont) - 1)));
            elseif prob.u_hold == "FOH"
                u_cont = u_func(t_cont);
            end
        end

        function [x_disc] = disc_prop(prob, u, p)
            %DISC_PROP Summary of this function goes here
            %   Detailed explanation goes here
            x_disc = zeros([prob.n.x, prob.N - 1]);
            x_disc(:, 1) = prob.x0;

            if prob.u_hold == "ZOH"
                for k = 1:(prob.N - 1)
                    x_disc(:, k + 1) = prob.disc.A_k(:, :, k) * x_disc(:, k) ...
                                 + prob.disc.B_k(:, :, k) * u(:, k) ...
                                 + zero_if_empty(prob.disc.E_k(:, :, k) * p) ...
                                 + prob.disc.c_k(:, :, k);
                end
            elseif prob.u_hold == "FOH"
                for k = 1:(prob.N - 1)
                    x_disc(:, k + 1) = prob.disc.A_k(:, :, k) * x_disc(:, k) ...
                                 + prob.disc.B_minus_k(:, :, k) * u(:, k) ...
                                 + prob.disc.B_plus_k(:, :, k) * u(:, k + 1) ...
                                 + zero_if_empty(prob.disc.E_k(:, :, k) * p) ...
                                 + prob.disc.c_k(:, :, k);
                end
            end
        end


    end
end

