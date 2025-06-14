classdef SOC2_Ball < Cone
    %SOC2_BALL l2-norm ball intersected with type 2 SOC
    %   Detailed explanation goes here
    
    properties
        indices
        
        soc2
        ball
    end
    
    methods
        function obj = SOC2_Ball(indices, r, theta, u)
            %SOC2_BALL Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;

            obj.soc2 = SOC2(indices, theta, u);
            obj.ball = Ball(indices, r);
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here

            x = obj.soc2.project(k, z);
            x = obj.ball.project(k, x);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc_soc2 = obj.soc2.get_constraint_function;
            cfunc_ball = obj.ball.get_constraint_function;

            cfunc = @(z) [cfunc_soc2(z); cfunc_ball(z)];
        end

        function [obj] = shift_indices(obj, new_index_map)
            obj.soc2 = obj.soc2.shift_indices(new_index_map);
            obj.ball = obj.ball.shift_indices(new_index_map);
        end
    end
end

