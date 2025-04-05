classdef DeterministicProblem
    %DETERMINISTICPROBLEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        x0
        N
        Nu
        tf
        u_hold
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
        function obj = DeterministicProblem(inputArg1,inputArg2)
            %DETERMINISTICPROBLEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;

            obj = linearize(obj);
        end

        function prob = linearize(prob)
            %LINEARIZE 

            % Linearize Dynamics

            % Linearize Nonconvex Constraints (?)

            % Linearize Boundary Conditions (?)

        end
        
        function prob = discretize(prob, x_ref, u_ref, p_ref)
            %DISCRETIZE Summary of this method goes here
            %   Detailed explanation goes here
            
            % Discretize Dynamics
            if prob.u_hold == "ZOH"
                [prob.disc.A_k, prob.disc.B_k, prob.disc.E_k, prob.disc.c_k] = discretize_dynamics_ZOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, [0, prob.tf], x_ref, u_ref, p_ref, prob.tolerances);
            elseif prob.u_hold == "FOH"
                [prob.disc.A_k, prob.disc.B_plus_k, prob.disc.B_minus_k, prob.disc.E_k, prob.disc.c_k] = discretize_dynamics_FOH(prob.cont.f, prob.cont.A, prob.cont.B, prob.cont.E, prob.cont.c, prob.N, [0, prob.tf], x_ref, u_ref, p_ref, prob.tolerances);
            end
        end
    end
end

