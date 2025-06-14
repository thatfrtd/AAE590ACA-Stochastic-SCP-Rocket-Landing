classdef SOC < Cone
    %SOC Second-order cone
    %   Type 1 : cone axis along last coordinate and arctan(alpha) is cone
    %   angle
    
    properties
        indices
        k_override = []

        alpha % tan(cone angle)
        branching = false
    end
    
    methods
        function obj = SOC(indices, alpha)
            %SOC Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;
            
            obj.alpha = alpha;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(obj.alpha <= 0) == 0, "Cone angle actan(alpha) must be positive")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.alpha(k), 2));

            n = size(z1, 1);

            beta = vecnorm(z1(1 : (n - 1), :), 2, 1) + 1e-15;

            if obj.branching
                num = numel(k);
                
                for i = 1:num
                    if beta(i) <= obj.alpha(k(i)) * z1(n, i)
                        z1(:, i) = z1(:, i);
                    elseif obj.alpha(k(i)) * beta(i) + z1(n, i) <= 0
                        z1(:, i) = 0;
                    else
                        z1(:, i) = z1(:, i);
                        z1(n, i) = (beta(i) * obj.alpha(k(i)) + z1(n, i)) / (obj.alpha(k(i)) ^ 2 + 1);
                        z1(1 : (n - 1), i) = z1(1 : (n - 1), i) .* obj.alpha(k(i)) .* z1(n, i) ./ beta(i); 
                    end
                end
            else
                cond_1 = (beta <= obj.alpha(k) .* z1(n, :));
                cond_2 = (obj.alpha(k) .* beta + z1(n, :) <= 0) .* ~cond_1;
                cond_3 = ~cond_1 .* ~cond_2;

                x_n = (beta .* obj.alpha(k) + z1(n, :)) ./ (obj.alpha(k) .^ 2 + 1);

                z1 = cond_1 .* z1 ...
                  + cond_3 .*  [z1(1 : (n - 1), :) .* obj.alpha(k) .* x_n ./ beta; x_n];
            end

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.alpha, 2));

            n = size(z, 1);

            cval = norms(z(1 : (n - 1), :), 2, 1) - obj.alpha .* z(n, :);
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

