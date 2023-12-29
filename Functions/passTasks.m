function [t2a, numPass] = passTasks(t2a, agents, tasks, simTh, numPass)
%This function performs reorganization of which agent works on which task
%
%Usage: [Task2Agent, NumberOfPasses] = pasTasks(Task2Agent, Agents, Tasks, SimilarityThreshold, NumberOfPasses)
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Perform pass calculation

% For each active task evaluate if the agent working on it should pass it or keep it
for tidx = find(not(isnan(t2a)))
    
    %Check if task was not passed too many times and allow further pass-evaluation only if passes are still allowed
    if numPass(tidx) > 0
        
        %Who is similar enough to be worthy for communication? -- asked the agent working to task tidx?
        simAgs = find(sqrt(sum((agents - agents(t2a(tidx), :)).^2, 2)) < simTh);
        %This could be made even faster by changing the calculation to a simple table-lookup

        %Which are the free similar agents?
        freeAgs = simAgs(not(ismember(simAgs, t2a)));
        
        %How good am I to solve my task?
        myFit = CalcSolvFitness(agents(t2a(tidx), :), tasks(tidx, :));
                
        %How good are free agents I talk to?
        theirFit = CalcSolvFitness(agents(freeAgs, :), tasks(tidx, :));
                
        %Who's best?
        [~, bestidx] = min([myFit; theirFit]);
        
        %If I'm not the best, re-assign my task to the best agent
        if bestidx ~= 1
            t2a(tidx) = freeAgs(bestidx-1);
            numPass(tidx) = numPass(tidx) - 1;
        end
    end
end
