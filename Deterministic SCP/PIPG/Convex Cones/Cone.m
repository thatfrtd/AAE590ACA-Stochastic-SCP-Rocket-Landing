classdef Cone < matlab.mixin.Heterogeneous
    %CONE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

    end
    
    methods(Sealed)
        % This section is for methods that call the abstract methods of
        % arrays of subclasses.
        % Example:
        % function constraint_values = evaluate_constraint(aryIn, outputs)
        % 
        %     n = numel(aryIn);
        %     constraint_values = zeros(size(aryIn));
        % 
        %     for c = 1:n
        %         constraint_values(c) = evaluate(aryIn(c), outputs);
        %     end
        % end

        function [x] = project_onto_cones(aryIn, k, z)
            n = numel(aryIn);

            x = z;

            for c = 1:n
                x = project(aryIn(c), k, x);
            end
        end

        function [cfuncs] = get_constraint_functions(aryIn)
            n = numel(aryIn);

            cfuncs = {};

            for c = 1:n
                cfuncs{c} = get_constraint_function(aryIn(c));
            end
        end

        function [cvals] = eval_constraint_functions(aryIn, z)
            n = numel(aryIn);

            for c = 1:n
                cfunc = get_constraint_function(aryIn(c));
                
                %dot(aryIn.u_1, z(obj.indices, :)) - aryIn.zeta_1

                cvals{c} = cfunc(z);
            end
        end

        function [aryIn] = shift_all_indices(aryIn, new_index_map)
            n = numel(aryIn);

            for c = 1:n
                aryIn(c) = aryIn(c).shift_indices(new_index_map);
            end
        end
    end

    methods(Abstract)
        % Abstract methods that every subclass must implement

        x = project(obj, k, z)

        cfunc = get_constraint_function(obj);

        obj = shift_indices(obj, new_index_map);
    end
end

