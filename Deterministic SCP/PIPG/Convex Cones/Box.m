classdef Box < Cone
    %BOX Box linfinity-norm ball
    %   Detailed explanation goes here
    
    properties
        indices
        k_override = []

        l % lower bound
        u % upper bound
    end
    
    methods
        function obj = Box(indices, l, u)
            %BOX Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;

            obj.u = u;
            obj.l = l;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(obj.u <= obj.l, "all") == 0, "Upper bound must be greater than the lower bound")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.u(:, k), 2));

            z1 = min(obj.u(:, k), max(z1, obj.l(:, k)));

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.u, 2));

            cval = [z - obj.u ... % upper bound
                    obj.l - z]; % lower bound
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

