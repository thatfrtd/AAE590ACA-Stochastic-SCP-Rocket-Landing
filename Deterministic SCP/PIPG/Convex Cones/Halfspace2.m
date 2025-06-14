classdef Halfspace2 < Cone
    %HALFSPACE2 Halfspace2 { x | dot(u_1, x) <= zeta_1 } and { x | dot(u_2, x) <= zeta_2 }
    %   Intersection of two halfspaces
    
    properties
        indices
        k_override = []

        u_1 % Unit normal to supporting hyperplane 1
        u_2 % Unit normal to supporting hyperplane 2
        zeta_1 % scalar on RHS for 1
        zeta_2 % scalar on RHS for 1
        branching = false
    end
    
    methods
        function obj = Halfspace2(indices, u_1, u_2, zeta_1, zeta_2)
            %HALFSPACE2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.indices = indices;

            obj.u_1 = u_1;
            obj.u_2 = u_2;
            obj.zeta_1 = zeta_1;
            obj.zeta_2 = zeta_2;
        end
        
        function x = project(obj, k, z)
            %PROJECT Summary of this method goes here
            %   Detailed explanation goes here
            assert(sum(vecnorm(obj.u_1, 2, 1) ~= 1) == 0, "Normal vector u_1 must have unit norm")
            assert(sum(vecnorm(obj.u_2, 2, 1) ~= 1) == 0, "Normal vector u_2 must have unit norm")
            assert(sum(obj.zeta_1(k) + obj.zeta_2(k) <= 0) == 0, "The two halfspaces do not intersect")

            if ~isempty(obj.k_override)
                k = obj.k_override;
            end

            x = z;

            shape = size(z(obj.indices, :));
            z1 = reshape(z(obj.indices, :), [], size(obj.u_1(:, k), 2));

            beta = dot(obj.u_1(:, k), obj.u_2(:, k));
            alpha_1 = dot(z1(:, k), obj.u_1(:, k)) - obj.zeta_1(k);
            alpha_2 = dot(z1(:, k), obj.u_2(:, k)) - obj.zeta_2(k);
            
            if obj.branching
                num = numel(k);

                beta_2 = beta;

                for i = 1:num
                    if beta(i) == 1
                        beta_2(i) = min(obj.zeta_1(k(i)), obj.zeta_2(k(i)));
        
                        if alpha_1(i) > beta_2(i)
                            cond_1_1 = cond_1_1 + 1;

                            z1(obj.indices, i) = (beta_2(i) - alpha_1(i)) * obj.u_1(:, k(i)) + z1(obj.indices, i);
                        else
                            z1(obj.indices, i) = z1(obj.indices, i);
                        end
                    elseif beta(i) == -1
                        if alpha_1(i) < -obj.zeta_2(k(i))
                            z1(obj.indices, i) = (-alpha_1(i) - obj.zeta_2(k(i))) * obj.u_1(:, k(i)) + z1(obj.indices, i);
                        elseif alpha_1(i) > obj.zeta_1(k(i))
                            z1(obj.indices, i) = (-alpha_1(i) + obj.zeta_1(k(i))) * obj.u_1(:, k(i)) + z1(obj.indices, i);
                        else
                            z1(obj.indices, i) = z1(obj.indices, i);
                        end
                    else                        
                        if alpha_1(i) - obj.zeta_1(k(i)) + beta(i) * (obj.zeta_2(k(i)) - alpha_2(i)) > 0 && alpha_2(i) - obj.zeta_2(k(i)) + beta(i) * (obj.zeta_1(k(i)) - alpha_1(i)) > 0
                            z1(obj.indices, i) = (alpha_1(i) - obj.zeta_1(k(i)) + beta(i) * (obj.zeta_2(k(i)) - alpha_2(i))) / (beta(i) ^ 2 - 1) * obj.u_1(:, k(i)) + z1(obj.indices, i);
                            z1(obj.indices, i) = (alpha_2(i) - obj.zeta_2(k(i)) + beta(i) * (obj.zeta_1(k(i)) - alpha_1(i))) / (beta(i) ^ 2 - 1) * obj.u_2(:, k(i)) + z1(:, i);
                        elseif alpha_2(i) > obj.zeta_2(k(i)) && alpha_1(i) - obj.zeta_1(k(i)) + beta(i) * (obj.zeta_2(k(i)) - alpha_2(i)) <= 0
                            z1(obj.indices, i) = (obj.zeta_2(k(i)) - alpha_2(i)) * obj.u_2(:, k(i)) + z1(obj.indices, i);
                        elseif alpha_1(i) > obj.zeta_1(k(i)) || alpha_2(i) > obj.zeta_2(k(i))
                            z1(obj.indices, i) = (obj.zeta_1(k(i)) - alpha_1(i)) * obj.u_1(:, k(i)) + z1(obj.indices, i);
                        else 
                            z1(obj.indices, i) = z1(obj.indices, i);
                        end
                    end
                end
            else
                zeta_min = min(obj.zeta_1(k), obj.zeta_2(k));

                cond_1 = (beta == 1);
                cond_1_1 = (alpha_1 > zeta_min);

                cond_2 = (beta == -1);
                cond_2_1 = (alpha_1 < -obj.zeta_2(k));
                cond_2_2 = (alpha_1 > obj.zeta_1(k)) .* ~cond_2_1;

                cond_3 = ~(cond_1 .* cond_2);
                cond_3_1 = ((alpha_1 - obj.zeta_1(k) + beta .* (obj.zeta_2(k) - alpha_2) > 0) .* (alpha_2 - obj.zeta_2(k) + beta .* (obj.zeta_1(k) - alpha_1) > 0));
                cond_3_2 = ((alpha_2 > obj.zeta_2(k)) .* (alpha_1 - obj.zeta_1(k) + beta .* (obj.zeta_2(k) - alpha_2) <= 0)) .* ~cond_3_1;
                cond_3_3 = min(1, ((alpha_1 > obj.zeta_1(k)) + (alpha_2 > obj.zeta_2(k)))) .* ~cond_3_1 .* ~cond_3_2;

                z1 = cond_1 .* cond_1_1 .* (zeta_min - alpha_1) .* obj.u_1(:, k) ...
                  + cond_2 .* (cond_2_1 .* (-alpha_1 - obj.zeta_2(k)) .* obj.u_1(:, k) ...
                             + cond_2_2 .* (-alpha_1 + obj.zeta_1(k)) .* obj.u_1(:, k)) ...
                  + cond_3 .* (cond_3_1 .* (~cond_3_2 .* ~cond_3_3) ...
                                   .* ((alpha_1 - obj.zeta_1(k) + beta .* (obj.zeta_2(k) - alpha_2)) ./ (beta .^ 2 - 1) .* obj.u_1(:, k) ...
                                     + (alpha_2 - obj.zeta_2(k) + beta .* (obj.zeta_1(k) - alpha_1)) ./ (beta .^ 2 - 1) .* obj.u_2(:, k)) ...
                             + cond_3_2 .* (~cond_3_1 .* ~cond_3_3) ...
                                   .* (obj.zeta_2(k) - alpha_2) .* obj.u_2(:, k) ...
                             + cond_3_3 .* (~cond_3_1 .* ~cond_3_2) ...
                                   .* (obj.zeta_1(k) - alpha_1) .* obj.u_1(:, k)) ...
                  + z1;
            end

            x(obj.indices, :) = reshape(z1, shape);
        end

        function [cfunc] = get_constraint_function(obj)
            cfunc = @(z) obj.constraint_function(z);
        end

        function [cval] = constraint_function(obj, z)
            shape = size(z);
            z = reshape(z(obj.indices, :), [], size(obj.u_1, 2));

            cval = [dot(obj.u_1, z) - obj.zeta_1;
                    dot(obj.u_2, z) - obj.zeta_2];
        end

        function [obj] = shift_indices(obj, new_index_map)
            new_indices = new_index_map(obj.indices, :);

            obj.indices = new_indices(:);

            obj.k_override = 1:size(new_indices, 2);
        end
    end
end

