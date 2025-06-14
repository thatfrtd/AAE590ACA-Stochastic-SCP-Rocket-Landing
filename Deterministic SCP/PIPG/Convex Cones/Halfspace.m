classdef Halfspace < Cone
    %HALFSPACE Halfspace { x | dot(u, x) <= zeta }
    %   Detailed explanation goes here
    
    properties
        indices
        k_override = []

        u % Unit normal to supporting hyperplane
        zeta % scalar on RHS
        branching = false
    end
    
    methods
        function obj = Halfspace(indices, u, zeta)
            %HALFSPACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;

            obj.u = u;
            obj.zeta = zeta;
        end
        
        function x = project(obj, k, z)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(vecnorm(obj.u, 2, 1) ~= 1) == 0, "Normal vector u must have unit norm")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.u(:, k), 2));

            beta = dot(z1, obj.u(:, k)) - obj.zeta(k);
            
            if obj.branching
                num = numel(k);

                for i = 1:num
                    if beta(i) > 0 % Constraint violated
                        z1(:, i) = -beta(i) * obj.u(:, k(i)) + z1(:, i);
                    else % Constraint satisfied
                        z1(:, i) = z1(:, i);
                    end
                end
            else
                z1 = -beta .* obj.u(:, k) .* (beta > 0) + z1;
            end

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) dot(obj.u, z(obj.indices, :)) - obj.zeta;
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.u, 2));

            cval = dot(obj.u, z) - obj.zeta;
        end
        
        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

