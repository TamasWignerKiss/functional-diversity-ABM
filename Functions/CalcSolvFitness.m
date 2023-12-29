function fitness = CalcSolvFitness(agents, task)
%This function determines how good agents are in solving a given task.
%
%Usage: Fitness = CalsSolvFitness(AgentList, Task)
%
%Here AgentList is an NxF double array, where N is the number of agents being queried and Task is the task to test agents on. Fitness is the amount of
%work agents in AgentList can perform on the given task
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Calculate agents' task solving fitness
tmp = task - agents;
tmp(tmp < 0) = 0;
fitness = nansum(tmp - task, 2); %#ok<NANSUM> 

end
