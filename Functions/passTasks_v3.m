function [t2a, numPass] = passTasks_v3(t2a, agents, tasks, simTh, numPass, agdistmat)
%This function performs reorganization of which agent works on which task. This is version 3.
%
%In the current version all agents that communicate first evaluate which is better in solving a given task. This way for each task a list of agents
%are assigned. Second, tasks are assigned to agents in order of their fitness. In this version, if there are N agents for N tasks (a subset of all
%tasks), their fitness is evaluated within the group, and agents assigned accordingly.
%
%Usage: [Task2Agent, NumberOfPasses] = pasTasks(Task2Agent, Agents, Tasks, SimilarityThreshold, NumberOfPasses)
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Perform pass calculation

activeTasks = find(not(isnan(t2a)));
agFit  = cell(numel(t2a), 1);
simAgs = cell(numel(t2a), 1);

% First, for each active task evaluate how good agents are in solving them, taking into account communication
for tidx = activeTasks

    %Check if task was not passed too many times and allow further pass-evaluation only if passes are still allowed
    if numPass(tidx) > 0

        %Who is similar enough to be worthy for communication? -- asked the agent working on task tidx?
        simAgs{tidx} = find(agdistmat(t2a(tidx), :) <= simTh);

        %Calculate how good all these friends are in solving task tidx
        agFit{tidx} = CalcSolvFitness(agents(simAgs{tidx}, :), tasks(tidx, :));
    end
end

%Second, convert simAgs and agFit to matrices for faster computation
maxAg = max(cellfun(@length, simAgs));
tmpSA = NaN(numel(t2a), maxAg);
tmpAF = NaN(numel(t2a), maxAg);
for tidx = 1:numel(t2a)
    if ~isempty(simAgs{tidx})
        tmpSA(tidx, 1:length(simAgs{tidx})) = simAgs{tidx};
        tmpAF(tidx, 1:length(simAgs{tidx})) = agFit{tidx};
    end
end
simAgs = tmpSA;
agFit  = tmpAF;

% I had quite some trouble with how to implement this resulting in programming errors :) Hence, some debug variables
simAgs_orig = simAgs;
agFit_orig = agFit;
t2a_orig = t2a;
activeTasks_orig = activeTasks;

% Third, assign agents to tasks, based on how much real work they can do. First assign the one that can do the most, ... until all tasks have an agent
try
    while ~isempty(activeTasks)
        % Test if there are enough agents for a given list of tasks
        tmp  = simAgs;
        tmp(isnan(tmp)) = inf;
        [~, ~, c] = unique(sort(tmp, 2), 'rows');
        agNum = sum(not(isinf(tmp)), 2);
        t = tabulate(c);
        %*********************
        %This is not right yet. There can be multiple tasks that have the same agents assigned to but some agents might overlap...
        if any(t(c, 2)-agNum == 0) %There is no choice, only one selection is possible
            tmp = find(t(c, 2)-agNum == 0); %This is the subset of tasks that has to be assigned right now
            while ~isempty(tmp)
                [t2a, simAgs, agFit, numPass, activeTasks] = reassign(t2a, simAgs, agFit, numPass, tmp, activeTasks);
            end
        else %If there will be agents allowing a choice
            [t2a, simAgs, agFit, numPass, activeTasks] = reassign(t2a, simAgs, agFit, numPass, activeTasks, activeTasks);
        end
    end

    % If there are remaining tasks, assign them to available agents
    if exist('remTasks', 'var')
        avAgs = setdiff(1:size(agents, 1), t2a(setdiff(1:end, remTasks))); %These are the available agents (agents not in t2a)
        t2a(remTasks) = randsample(avAgs, length(remTasks)); %Randomly assign available agents to remaining tasks
    end
catch ME
    fprintf('\nDump:\n')
    fprintf('********************\n')
    disp(activeTasks_orig)
    disp(t2a_orig)
    disp(simAgs_orig)
    disp(agFit_orig)
    fprintf('********************\n')
    disp(activeTasks)
    disp(t2a)
    disp(simAgs)
    disp(agFit)
    rethrow ME
end

    function [nt2a, nsimAgs, nagFit, nnumPass, nactiveTasks] = reassign(tta, sA, aF, nP, cT, aT)
        nt2a = tta;
        nsimAgs = sA;
        nagFit = aF;
        nnumPass = nP;
        nactiveTasks = aT;

        [I, J] = find(min(min(aF(cT, :))) == aF(cT, :)); %Which agent(J) at which task(I) would do the most work
        I = I(1); J = J(1); %In case of a draw, pick the first one
        bestAg = sA(cT(I), J); %The #number of the best agent
        nt2a(cT(I)) = bestAg; %(Re-)assign the task
        nnumPass(cT(I)) = nP(cT(I)) - 1;
        ntmp = nsimAgs == bestAg; %Locate all tasks where this agent shows up
        nsimAgs(ntmp) = NaN; %Remove agents from the list of available agents -- this line is unnecessary
        nagFit(ntmp) = NaN; %Remove its fitnes
        nagFit(cT(I), :) = NaN; %Remove task from list of tasks-to-be-assigned
        nactiveTasks(nactiveTasks == cT(I)) = []; %Remove task from list of tasks-to-be-assigned
    end

end
