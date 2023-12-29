function [t2a, numPass] = passTasks_v4(t2a, agents, tasks, simTh, numPass, agdistmat, stuckag)
%% Description
% This function performs reorganization of which agent works on which task.
%
% In this version stuck agents ask all their friends (not engaged or stuck, leaving out busy agent friends) if they can solve the task of the stuck
% agent. If not, bummer, nothing happens, agent does not pass and stays in the stuck state. If yes, they pass on their task. If the receiving agent is
% also a stuck agent, their tasks are swapped. If the receiving agent is an empty agent, the empty agent becomes the task solver and the original
% stuck agent becomes empty.
%
% Question to myself: what shall happen to a task when numPass (which is not used now) is reached? The agent is stuck forever? -- In current version
% this is what happens.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>
% Date: 2023-03-03

%% Pass calculation

%The stuck agents
sa = t2a(stuckag);

% Loop through all stuck agents
for sidx = 1:length(sa)

    % Identify stuck and empty friends. These can take the task of the current stuck agent.
    friends = find(agdistmat(sa(sidx), :) <= simTh); %These are the similar agents to the stuck agent (including the original stuck agent)
    sfidxv = ismember(friends, sa); %Which friends are stuck agents?
    efidxv = not(ismember(friends, t2a)); %Which friends are empty agents?
    friends = friends(sfidxv | efidxv); %Carry on with stuck or empty friends, leaving out busy, working friends

    % Figure out which friend is best in solving the task of the stuck agent
    tidxv = t2a == sa(sidx); %This is gona index the task of the stuck agent
    if sum(tidxv) ~= 1 %Just to be on the safe side
        error('Malformed t2a!')
    end
    if numPass(tidxv) <= 0 %Task was passed too many times, cannot be passed, go to next agent.
        continue
    end
    rt = tasks(tidxv,:) - agents(friends, :); %This is if friends reduce task of stuck agent
    wd = nansum(tasks(tidxv, :) - rt, 2); %#ok<NANSUM> %The amount of work done

    %Pass the task if possible
    if any(wd > 0) %If task was reduced, passing will happen. If not, nothing happens
        [~, tmp] = max(wd);
        bestfr = friends(tmp(1)); %The friend that can do the most work on the task; tmp(1) is taken in case there are multiple best agents
        if not(ismember(bestfr, t2a)) % If the best agent is an empty agent, it takes the task and sa(sidx) becomes empty
            t2a(tidxv) = bestfr;
            numPass(tidxv) = numPass(tidxv) - 1;

        else %If the best agent is a stuck agent, swap tasks, if both task has positive numPass values
            if bestfr ~= sa(sidx) %If the original stuck agent is the best, no change happens
                ridxv = t2a == bestfr; %Task of receiving agent
                if numPass(ridxv) > 0 %If task of receiving stuck agent is passed too many times, swapping cannot happen
                    t2a(tidxv) = bestfr;
                    t2a(ridxv) = sa(sidx);
                    numPass(tidxv) = numPass(tidxv) - 1;
                    numPass(ridxv) = numPass(ridxv) - 1; %If passing happen, both tasks are passed and numPass decreases
                end
            end
        end
    end
end
