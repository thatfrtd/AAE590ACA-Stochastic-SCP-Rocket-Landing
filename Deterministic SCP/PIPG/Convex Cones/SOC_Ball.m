classdef SOC_Ball < Cone
    %SOC_BALL l2-norm ball intersected with type 1 SOC
    %   Detailed explanation goes here
    
    properties
        indices

        soc
        ball
    end
    
    methods
        function obj = SOC_Ball(indices, r, alpha)
            %SOC_BALL Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;
            
            obj.soc = SOC(indices, alpha);
            obj.ball = Ball(indices, r);
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here

            x = obj.soc.project(k, z);
            x = obj.ball.project(k, x);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc_soc = obj.soc.get_constraint_function;
            cfunc_ball = obj.ball.get_constraint_function;

            cfunc = @(z) [cfunc_soc(z); cfunc_ball(z)];
        end

        function [obj] = shift_indices(obj, new_index_map)
            obj.soc = obj.soc.shift_indices(new_index_map);
            obj.ball = obj.ball.shift_indices(new_index_map);
        end
    end
end

