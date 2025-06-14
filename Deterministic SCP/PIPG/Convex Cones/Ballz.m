classdef Ballz
    %BALL Ball l2-norm ball
    %   Detailed explanation goes here
    
    properties
        indices
        k_override = []

        r % radius
    end
    
    methods
        function obj = Ballz(indices, r)
            %BALL Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;
            obj.r = r;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(obj.r <= 0) == 0, "Radius should be positive")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.r(k), 2));

            z1 = obj.r(k) ./ max(vecnorm(z1, 2, 1), obj.r(k)) .* z1;

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)

            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.r, 2));

            cval = norms(z, 2, 1) - obj.r;
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

