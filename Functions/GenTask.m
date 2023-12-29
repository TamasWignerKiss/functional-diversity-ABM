function tasks = GenTask(nf, nt, no, es, type, agents)
%This function generates tasks
%
%Usage: tasks = GenTask(NoFunctions, NoTasks, NormalFact, EmergencyStop, Type, Agents)
%
%Input arguments:
%NoFunctions: integer, the number of sub-elements in a task. Should be the same as the skills of agents
%NoTasks: integer, the number of tasks to be generated
%NormalFact: double, the normalizing factor. Task sub-elements will be normalized such that their sum is equal to NormalFact
%EmergencyStop: integer, the number of maximum steps allowed in the simulation
%Type: 1/2: switches between the type of task to be generated as follows:
%  1: simply generate random numbers and normalize them
%  2: only allow task sub-functions that the system of agents as a whole can solve within the maximum allowed time steps
%
%Agents: N x NoFunctions double array, the agents (used for type 2 and 3 tasks)
%
%Output values
%tasks: NoTasks x NoFunctions doulbe array, the tasks
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Generate tasks

% The proto-tasks
tasks = rand(nt, nf);
tasks = tasks./sum(tasks, 2)*no;

switch type
    case 1
        %Proto tasks are good enough
    case 2
        mas = max(agents); %The fastet solver for a given function
        tmp = tasks - mas*es; %This would be the value of the task after es steps if always the fastet was solving. If positive, cannot be done
        tasks(tmp > 0) = NaN;
        tasks = tasks./nansum(tasks, 2)*no; %#ok<NANSUM> %Need to normalize again after removing some components
    otherwise
        error('Don know hot to generate these tasks.')
end

end
