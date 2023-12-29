function [t2a, taskhist, stepNo] = PassingSolveTasks(agents, tasks, etc, agdistmat, simTh)
%This function runs an instance of the simulation. In this version the simulation runs until the maximal number of steps is reached or all tasks are
%solved. The following processes are performed:
%1. Randomly assign tasks to agents (also, keep track of the assignment across stimulation steps)
%2. Figure out if tasks are to be passed around using the current rule
%3. Agents not in refractory work on the task
%4. Free agents re-assigned to task
%5. Go to step 2 unless maximum number of steps reached or all tasks solved.
%
% Input arguments:
% agents: is the agents matrix
% tasks: is the task matrix
% etc: other parameters, a struct containing:
%   emStop: a number that specifies the maximum number of steps allowed (to stop simulation in case the system gets stuck)
%   maxPass: maximum number of allowable passes
%   passCost: the number of steps an agent has to wait to do any work after receiving a task
%   NoPassInRef: a switch to tell if passing during refractory is allowed or not
% agdistmat: distance matrix -- pre-calculated for effective computation
% simTh: a threshold of similarity (percent of max) above which agents don't pass
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Some initialization
stepNo = 0; %Counter to track number of simulation steps needed to solve all tasks
numPass = etc.maxPass*ones(size(tasks,1), 1); %To assign the same numpass to all tasks at the start of the simulation
agDelay = zeros(1, size(agents, 1)); %After receiving a task agents enter a refractory period of etc.passCost length

%% Solving the tasks

%Assign tasks to agents
t2a = -1*ones(etc.emStop+1, size(tasks, 1));
tloc = randperm(size(tasks,1), min(size(agents,1), size(tasks,1)));
t2a(1, tloc) = randperm(size(agents,1), min(size(tasks,1), size(agents,1)));

%This variable will store task history
taskhist = NaN(size(tasks, 1), size(tasks, 2), etc.emStop+1);
taskhist(:, :, 1) = tasks;

%To be able to handle passing/no passing during refractory, and being stuck, I have to initialize this here
DelayedAgentIdx = false(1, size(tasks, 1));
StuckAg = false(1, size(tasks, 1));

% Iterate until there are tasks or emergeny stop reached
while any(not(isnan(t2a(stepNo+1, :)))) && stepNo < etc.emStop

    switch etc.passVersion
        case 1
            % This is the version in which agents evaluate in every step if another agent is better in solving a task than themselves
            [t2a(stepNo+2, :), numPass] = ...
                passTasks_v2(t2a(stepNo + 1, :), agents, tasks, simTh, numPass, agdistmat, DelayedAgentIdx, etc.NoPassInRef);

        case 2
            % In this version only stuck agents evaluate if others are better in solving than themselves
            if sum(StuckAg) > 0
                [t2a(stepNo+2, :), numPass] = ...
                    passTasks_v4(t2a(stepNo + 1, :), agents, tasks, simTh, numPass, agdistmat, StuckAg);
            else
                t2a(stepNo+2, :) = t2a(stepNo+1, :);
            end

        otherwise
            error('Unknown passVersion')
    end

    %If agent received new task, start refractory, if kept task, decrease refractory time (negative values don't matter)
    %These are the agents that got a task or swapped tasks (appeared (in new position) in t2a relative to previous line)
    newAg = t2a(stepNo+2, (t2a(stepNo+2, :) - t2a(stepNo+1, :)) ~= 0 & not(isnan(t2a(stepNo+2, :))) & t2a(stepNo+2, :)  ~= -1);
    %These are the agents that kept working on their task (in the same positoon
    oldAg = t2a(stepNo+2, (t2a(stepNo+2, :) - t2a(stepNo+1, :)) == 0 & not(isnan(t2a(stepNo+2, :))) & t2a(stepNo+2, :)  ~= -1);
    agDelay(newAg) = etc.passCost; %Enter refractory
    agDelay(oldAg) = agDelay(oldAg) - 1; %Working on getting out of refractory

    %Need to find busy agents and active tasks because NaN t2a cannot be used for index
    BusyAgentIdx = not(isnan(t2a(stepNo+2, :))) & t2a(stepNo+2, :)  ~= -1; % Index into t2a | These are agents assigned to tasks in t2a
    tmp = t2a(stepNo+2, :); %Helper variable to circumvent issue caused by indexing using NaN
    tmp(isnan(tmp) |  tmp == -1) = 1; %Have to index into agDelay but cannot use NaN or -1 as an index
    DelayedAgentIdx = agDelay(tmp) > 0; %Also index into t2a | agents that just got new task and so are delayed
    DelayedAgentIdx(isnan(t2a(stepNo+2, :)) | t2a(stepNo+2, :) == -1) = false; %These are tasks that are completed, need to put back to show that these are not delayed agents
    t2aIdx = BusyAgentIdx & not(DelayedAgentIdx); %To select tasks and agents using the t2a assignment variable

    %Active and not delayed agents work on tasks, mark by NaN if a component is completed
    tasks(t2aIdx, :) = tasks(t2aIdx, :) - agents(t2a(stepNo+2, t2aIdx), :);
    tasks(tasks <= 0) = NaN;

    %Find finished tasks
    isFinished = all(isnan(tasks), 2);

    %Free agent
    t2a(stepNo+2, isFinished) = NaN;

    %Store current tasks
    taskhist(:, :, stepNo + 2) = tasks;

    %See if an agent is stuck
    %An agent is considered to be stuck if he cannot finish his task before the deadline
    tmp1 = taskhist(:, :, stepNo + 2);
    remTask = nansum(tmp1, 2); %#ok<NANSUM>
    tmp2 = taskhist(:, :, stepNo + 1);
    tmp1(isnan(tmp1)) = -100; %Have to track change of task from previous step but sum([NaN, anythin]) == NaN, hence I use the ad hoc number -100
    tmp2(isnan(tmp2)) = -100;
    tmp = tmp2 - tmp1;
    StuckAg = ((remTask ./ sum(tmp, 2)) > (etc.emStop - (stepNo    +2)))' & not(DelayedAgentIdx) & t2a(stepNo+2, :)  ~= -1;

    %Find idle agents
    IdleAg = setdiff(1:size(agents, 1), t2a(stepNo+2, :)); %Gives the number of agent(s)

    %Find pending tasks
    PendTasks = find(t2a(stepNo+2, :) == -1); %Gives a location in t2a

    %Assign pending task to idle agent(s)
    if not(isempty(IdleAg)) && not(isempty(PendTasks))
        tloc = PendTasks(randperm(length(PendTasks), min(length(PendTasks), length(IdleAg)))); %Which task(s) to assign
        t2a(stepNo+2, tloc) = IdleAg(randperm(length(IdleAg), min(length(PendTasks), length(IdleAg)))); %Choose idle agent and put to work
        DelayedAgentIdx(tloc) = true;
        agDelay(t2a(stepNo+2, tloc)) = etc.passCost;
    end

    %Increase step counter
    stepNo = stepNo + 1;
end

% No work was done using the initial task assignment, only passing happened
%t2a = t2a(2:end, :);
