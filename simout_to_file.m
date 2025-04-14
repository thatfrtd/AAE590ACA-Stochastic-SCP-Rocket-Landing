function [] = simout_to_file(out, options)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% output_variables: output_variables = ["x_NED", "Quaternions"]
% simout_to_file(out, file_name = "testsimout.csv", rocket_length = rocket.NC_length + rocket.AF1_length)
arguments
    out
    options.output_variables = ["x_NED", "Quaternions", "COM"]
    options.file_name = "simout.csv"
    options.tstep = 0.1
    options.rocket_length = []
end

T = table;
T.Time = out.tout;

for s = options.output_variables
    d = squeeze(out.yout.get(s).Values.Data);
    if size(d, 1) ~= length(d)
        for i = 1:size(d, 1)
            T.(s + "_" + string(i)) = d(i, :)';
        end
    else
        T.(s) = d;
    end
end

if ~isempty(options.rocket_length)
    T.("COM") = options.rocket_length - T.("COM");
end

T_resamp = resample_data(T, 0:options.tstep:out.tout(end));

writetable(T_resamp, options.file_name)

end