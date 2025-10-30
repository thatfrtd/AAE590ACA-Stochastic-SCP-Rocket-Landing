function [multi_f] = multidynamic_f(f, k_set)
%MULTIDYNAMIC_F Summary of this function goes here
%   Detailed explanation goes here

multi_f = {};
for k = 1:numel(k_set)
    multi_f{k} = @(t, x, u, p) f(t, x, u, p, k_set(k));
end

end

