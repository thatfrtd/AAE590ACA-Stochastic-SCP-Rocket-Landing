classdef Singleton < Cone
    %SINGLETON Singleton {z | z = y} 
    %   Detailed explanation goes here
    
    properties
        indices
        k_override = []

        y % Value z is constrained to be
    end
    
    methods
        function obj = Singleton(indices, y)
            %SINGLETON Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;
            obj.y = y;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.y, 2));

            z1 = obj.y;

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.y, 2));

            cval = z - obj.y;
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);
            
            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

