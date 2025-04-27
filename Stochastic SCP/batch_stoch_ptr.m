function [ptr_sols] = batch_stoch_ptr(stoch_prob, ptr_ops, thread_number)
%BATCH_STOCH_PTR Summary of this function goes here
%   Detailed explanation goes here

clust = parcluster("Stoch2DoF");

prob_num = numel(stoch_prob);

if prob_num < thread_number
    thread_number = prob_num;
end

% Initialize jobs and job folders
for t = 1:thread_number
    % Create jobs so that tasks can be added to them
    jobs(t) = createJob(clust);

    w(t) = waitbar(0, sprintf("P#%d: Creating Job...", t));

    store{t} = jobs(t).ValueStore;
    store{t}.KeyUpdatedFcn = @(store,key) updateWaitbar(w(t),store,key);
end

for m = 1:prob_num    
    t = mod(m - 1, thread_number) + 1;
    tasks(m) = createTask(jobs(t), @Stochastic_ptr_Job, 1, {stoch_prob{m}, ptr_ops{m}, m},'CaptureDiary',false);
    waitbar(0, w(t), sprintf("P#%d: Solving PTR iter %d", t, 1));
end

% Submit Jobs
for t = 1:thread_number
    submit(jobs(t))
end

for t = 1:thread_number
    wait(jobs(t));
end

ptr_sols = {};

% Close loading bar and get results
for t = 1:thread_number
    job_output = fetchOutputs(jobs(t));
    ptr_sols{t} = job_output;

    close(w(t));
end

delete(jobs)
end

function updateWaitbar(w,store,key)
% Update a waitbar using the ValueStore property.
% Check if the waitbar is a reference to a deleted object.
if isvalid(w)
    % Update results from the job ValueStore object.
    result = store(key);

    if result.progress==1
        waitbar(result.progress,w,"Job Completed");
    else
        % Update the waitbar
        waitbar(result.progress, w, sprintf("P#%d: %s", result.prob_num, result.iter_info));
    end
end
end