classdef SOC2 < Cone
    %SOC2 Second-order cone
    %   Type 2 : cone axis along unit vector u and arccos(theta) is cone
    %   angle
    
    properties
        indices
        k_override = []

        theta_c % cos(cone angle)
        theta_s % sin(cone angle)
        u % unit vector along cone axis
        branching = false
    end
    
    methods
        function obj = SOC2(indices, theta, u)
            %SOC2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;

            obj.theta_c = cos(theta);
            obj.theta_s = sin(theta);
            obj.u = u;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(vecnorm(obj.u, 2, 1) ~= 1) == 0, "Normal vector u must have unit norm")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.u(:, k), 2));

            beta = vecnorm(z1, 2, 1);
            gamma = dot(z1, obj.u(:, k));

            if obj.branching
                num = numel(k);

                for i = 1:num
                    if obj.theta_c(k(i)) * beta(i) <= gamma(i)
                        z1(:, i) = z1(:, i);
                    elseif obj.theta_s(k(i)) * beta(i) + gamma(i) <= 0
                        z1(:, i) = 0;
                    else
                        kappa = -gamma(i) * obj.u(:, k(i)) + z1(:, i);
                        x1 = obj.theta_s(k(i)) / norm(kappa) * kappa + obj.theta_c(k(i)) * obj.u(:, k(i));
                        z1(:, i) = dot(z1(:, i), x1) * x1;
                    end
                end
            else
                cond_1 = (obj.theta_c(k) .* beta <= gamma);
                cond_2 = (obj.theta_s(k) .* beta + gamma <= 0);
                cond_3 = ~cond_1 .* ~cond_2;

                kappa = -gamma .* obj.u(:, k) + z1(:, k);
                x1 = obj.theta_s(k) ./ vecnorm(kappa, 2, 1) .* kappa + obj.theta_c(k) .* obj.u(:, k);

                z1 = cond_1 .* z1 ...
                  + cond_3 .* dot(x(:, k), x1) .* x1;
            end

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.u, 2));

            cval = vecnorm(z, 2, 1) - 1 ./ obj.theta_c .* dot(obj.u, z);
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

