%% Description
% This ABM simulates problem solving agents as they work on a set of tasks. Each agent has a trait vector describing their skills. Two agents are
% collaborators if their traits are similar enough. Tasks are defined by a component vector. Agents work on tasks by decreasing tasks' components
% vectors towards zero.
%
% The simulation generates agents with given intrapersonal functional diversity (IFD), and dominant function diversity (DFD), and a set of random
% tasks and evaluates how agents' traits influence their level of communication and the performance of a team of agents.
%
% Motivation for creating this model was to explain findings in papers by Bunderson & Sutcliffe (Academy of Management Journal, 45(5), 875-893, 2002)
% and Zhou et al (Review of Managerial Science, 17, 259-286, 2023)
%
% Author: Tamas Kiss <kiss.t@wigner.hu>
% Date: 2023-12-29

%% Parameters
% IO
par.functiondir = '/data/bognor/ATRAN/Functions';

% Model parameters
par.numfuncs = 9; %Number of domain functions
par.numagents = 10; %Number of agents in the system
par.numtasks = 7; %Number of tasks injected
par.anorm = 10; %Total skill of agents (agent normalizing factor)
par.tnorm = 10; %Total labour requirement of tasks (task normalizing factor)
par.TaskType = 1; %Sets what type of task to generate. 1 is random tasks; 2 is solvable tasks; For more details see help GenTask
par.MixSkillStrengths = false; %Whether to have all agents with skill strengths ordered, or have it mixed.
par.AgentGenMethod = 3; % See GenAgent.m for an explanation what methods are implemented
par.maxPass = Inf; %Each task can be passed this many times
par.similThresh = 90; %The threshold of similarity (as percentage of maximal distance) below which agents are willing to communicate w/ each other
par.passCost = 5; %Steps, the cost of communication: the receiving agent cannot work until this amount of time stamps
par.NoPassInRef = true; %If false, agent can pass even in refractory. If true no passing nor working are allowed in refractory
par.passVersion = 2; %Multiple ways of task passing were implemented. Choose here.

% External simulation control parameters
par.numrepeats = 10; %Number of independent realization of the process
par.EmergencyStop = 2.5e2; %Number of time steps allowed
par.IFDs = 0:0.05:0.9; %Number of agent diversity steps
par.DFDs = [0 0.20 0.36 0.47 0.54 0.61 0.65 0.70 0.74 0.79 0.83 0.86 0.88 0.90 0.92 0.95 0.97 0.99]; %Prescribed DFD values
par.IFDSLevel = 5; %Runs from 1 to par.numagents/2
par.IFDSFitRange = [0 0.99]; %The range of IFDS values to be covered by agents

% Analysis / Plotting
par.IFDbins = linspace(0, 1, 11); %For coarse-graining the values calculated
par.DFDbins = linspace(0, 1, 11);
par.DFDmeasd = [0.27 0.92]; %This was the range seen by B&S w/ mean = 0.66, sd = 0.17
par.IFDmeasd = [0.12 0.51]; %w/ mean = 0.28, sd = 0.09
par.PlotCG = false; %Whether to plot coarse-grained (true) or original (false) results
par.plotfigs = [true true true];

%% Set the path
warning('off', 'MATLAB:rmpath:DirNotFound')
restoredefaultpath;
rmpath('/home/umat/Documents/MATLAB')
clear RESTOREDEFAULTPATH_EXECUTED
addpath(par.functiondir)
rehash
warning('on', 'MATLAB:rmpath:DirNotFound')

%% Some initialization
DFD = NaN(length(par.IFDs), length(par.DFDs));
IFD = NaN(length(par.IFDs), length(par.DFDs));
meanstn = NaN(length(par.IFDs), length(par.DFDs));
meanca = NaN(length(par.IFDs), length(par.DFDs));
meansno = NaN(length(par.IFDs), length(par.DFDs));
simTh = NaN(length(par.IFDs), length(par.DFDs));
simAgs = NaN(length(par.IFDs), length(par.DFDs));
etc.maxPass = par.maxPass;
etc.emStop = par.EmergencyStop;
etc.passCost = par.passCost; %Steps
etc.passVersion = par.passVersion;
etc.NoPassInRef = par.NoPassInRef;

%% Loop through DFD and IFD values
tic
T = 0;
for ai = 1:length(par.IFDs)
    fprintf('%0.2f ', par.IFDs(ai))
    
    %The following piece of code runs in parallel on multiple cores. Remove "par" from parfor below to allow debugging of the inside.
    parfor gi = 1:length(par.DFDs) %This is for parallel computing
        
        %Generate agents -- this could be regenerated for every realization, but I don't think it makes a big difference. Better this way I think.
        agents = GenAgent(par.numfuncs, par.numagents, [], [], par.anorm, 'Method', par.AgentGenMethod, 'IFDSLevel', par.IFDSLevel, ...
            'IFDSFitRange', par.IFDSFitRange, 'DesiredIFD', par.IFDs(ai), 'DesiredDFD', par.DFDs(gi), 'SigEstMeth', 2, ...
            'MixSkillStrengths', par.MixSkillStrengths); %#ok<PFBNS> %This Code Analyzer Message is due to the parfor -> for change for gi
        agdistmat = squareform(pdist(agents));

        maxDist = sqrt(2*par.anorm^2);
        simTh(ai, gi) = maxDist*par.similThresh/100; %This is given as a percentage of the largest difference between any two agents
        simAgs(ai, gi) = sum(agdistmat <= simTh(ai, gi), "all") - par.numagents;

        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        
        %The following variables will only live for an [ai, gi] point
        t2a = NaN(par.EmergencyStop+1, par.numtasks, par.numrepeats); %Initialize variable keeping track of which agent works on which task
        taskhists = NaN(par.numtasks, par.numfuncs, par.EmergencyStop+1, par.numrepeats); %Initialize variable keeping track of how tasks progress
        stn = zeros(par.numrepeats, 1); %Initialize variable storing the amount of work done in repeats
        ca = NaN(par.numrepeats, 1); %Init variable for communication amount
        sno = NaN(par.numrepeats, 1); %Number of steps taken until completion

        for ridx = 1:par.numrepeats
            
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm, par.EmergencyStop, par.TaskType, agents);

            %PassingSolveTasks is the function that runs the while loop running the whole game of passing-working
            %It returns number of steps required to solve all tasks
            [t2a(:, :, ridx), taskhists(:, :, :, ridx), sno(ridx)] = PassingSolveTasks(agents, tasks, etc, agdistmat, simTh(ai, gi));

            % Calculate the amount of solved sub-tasks
            tmp1 = taskhists(:, :, 1, ridx);
            tmp2 = taskhists(:, :, par.EmergencyStop+1, ridx);
            tmp1(isnan(tmp1)) = 0; %There can be NaNs initially due to sub-tasks that the system cannot solve!
            tmp2(isnan(tmp2)) = 0;
            stn(ridx) = sum(sum(tmp1 - tmp2));

            % Calculate amount of communication
            tmp = diff(t2a(:, :, ridx));
            ca(ridx) = sum(sum(tmp ~= 0 & not(isnan(tmp)))) ./ sno(ridx);

        end
        
        % Average solved subtasks
        meanstn(ai, gi) = mean(stn);

        % Average communication amount
        meanca(ai, gi) = mean(ca);

        %Average number of steps needed to complete all tasks
        meansno(ai, gi) = mean(sno);

        fprintf('.')
    end

    % Step counter
    T_old = T;
    T = toc;
    fprintf(' (%0.2f sec)\n', T-T_old)
end
toc

% Interim cleanup
clear adiv agents ai gdiv gi stn T t2a T_old taskhists etc

%% Coarse-grain (average in bins) results and normalize max of amount of solved subtasks to 1
cgmstn = NaN(length(par.IFDbins)-1, length(par.DFDbins)-1);
cgca = NaN(size(cgmstn));
cgIFD = NaN(size(cgmstn));
cgDFD = NaN(size(cgmstn));
for iidx = 1:length(par.IFDbins)-1
    for didx = 1:length(par.DFDbins)-1
        cgmstn(iidx, didx) = mean(meanstn(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgca(iidx, didx)   = mean(meanca(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgIFD(iidx, didx)  = par.IFDbins(iidx);
        cgDFD(iidx, didx)  = par.DFDbins(didx);
    end
end

cgmstn = cgmstn/(par.numtasks*par.tnorm);

% Interim cleanup
clear iidx didx

%% Plot amount of work solved vs. IFD, DFD, also amount of communication and similarity-related measures
if par.plotfigs(1)
figure
p = panel();
p.pack(2,2)

% The work plot
p(1,1).select()
ax1 = gca;
if par.PlotCG
    surf(cgDFD, cgIFD, cgmstn)
else
    surf(DFD, IFD, meanstn/(par.numtasks*par.tnorm))
end
grid on
title(['TaskType=' num2str(par.TaskType) ' | MixSkills=' num2str(par.MixSkillStrengths) ' | PassCost=' num2str(par.passCost) ' | SimThresh=' ...
    num2str(par.similThresh) '% | NoAgs= ' num2str(par.numagents) ' | NoTasks= ' num2str(par.numtasks) ' | NoRepeats= ' num2str(par.numrepeats) ...
    ' | PassInRef=' num2str(not(par.NoPassInRef)) ' | PassVersion=' num2str(par.passVersion)])
xlabel('Dominant Function Diversity')
ylabel('Individual Functional Diversity')
zlabel('Percent Solved Subtasks')
xlim([0 1])
ylim([0 1])
zlim([0 1])
clim([0 1])

% The communication plot
p(1,2).select()
ax2 = gca;
if par.PlotCG
    surf(cgDFD, cgIFD, cgca)
else
    surf(DFD, IFD, meanca)
end
grid on
xlabel('Dominant Function Diversity')
ylabel('Individual Functional Diversity')
zlabel('Task Passing Per Step')
xlim([0 1])
ylim([0 1])
zlims = zlim;

% The similarity threshol plot
p(2,1).select()
ax3 = gca;
if par.PlotCG
    text(0.5, 0.5, 'Sorry, this is not implemented', 'HorizontalAlignment', 'center')
else
    surf(DFD, IFD, simTh)
end
grid on
xlabel('Dominant Function Diversity')
ylabel('Individual Functional Diversity')
zlabel('Similarity Threshold')
xlim([0 1])
ylim([0 1])

% The number of similar agents plot
p(2,2).select()
ax4 = gca;
if par.PlotCG
    text(0.5, 0.5, 'Sorry, this is not implemented', 'HorizontalAlignment', 'center')
else
    surf(DFD, IFD, simAgs/((par.numagents-1)*par.numagents))
end
grid on
xlabel('Dominant Function Diversity')
ylabel('Individual Functional Diversity')
zlabel('Percentage of Similar Agent Pairs')
xlim([0 1])
ylim([0 1])
zlim([0 1])

%Linking the four plots
Link = linkprop([ax1, ax2, ax3, ax4],{'View'});
setappdata(gcf, 'StoreTheLink', Link);
end

if par.plotfigs(2)
figure
surf(DFD, IFD, meansno)
xlabel('Dominant Function Diversity')
ylabel('Individual Functional Diversity')
zlabel('Number of Steps')
end

% Cleanup
clear p ax1 ax2 ax3 ax4

%% Some projections
if par.plotfigs(3)
figure
p = panel();
p.pack(2,2)

% The work plots
p(1,1).select()
if par.PlotCG
    plot(cgDFD(1,:), nanmean(cgmstn), 'b-*'); %#ok<*NANMEAN> 
else
    plot(DFD(1,:), nanmean(meanstn)/(par.numtasks*par.tnorm), 'b-*'); %#ok<*NANMEAN> 
end
xlim([0 1])
ylim([0 1])
vline(par.DFDmeasd, 'r')
grid on
xlabel('Dominant Function Diversity')
ylabel('Projected Percent Solved Subtasks')

p(1,2).select()
if par.PlotCG
    plot(cgIFD(:,1), nanmean(cgmstn, 2), 'b-*');
else
    plot(IFD(:,1), nanmean(meanstn, 2)/(par.numtasks*par.tnorm), 'b-*');
end
xlim([0 1])
ylim([0 1])
vline(par.IFDmeasd, 'r')
grid on
xlabel('Individual Functional Diversity')
ylabel('Projected Percent Solved Subtasks')

% The communication plots
p(2,1).select()
if par.PlotCG
    plot(cgDFD(1,:), nanmean(cgca), 'b-*');
else
    plot(DFD(1,:), nanmean(meanca), 'b-*');
end
xlim([0 1])
ylim(zlims)
vline(par.DFDmeasd, 'r')
grid on
xlabel('Dominant Function Diversity')
ylabel('Projected Amount of Communication')

p(2,2).select()
if par.PlotCG
    plot(cgIFD(:,1), nanmean(cgca, 2), 'b-*');
else
    plot(IFD(:,1), nanmean(meanca, 2), 'b-*');
end
xlim([0 1])
ylim(zlims)
vline(par.IFDmeasd, 'r')
grid on
xlabel('Individual Functional Diversity')
ylabel('Projected Amount of communication')

% Cleanup
clear p zlims
end

%% Some regressions here

% Create matrix of observations
obss = table('Size', [numel(meanstn), 4], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames', {'IFD', 'DFD', 'Comm', 'Perf'});
obssidx = 1;
for ai = 1:length(par.IFDs)
    for gi = 1:length(par.DFDs)
        obss(obssidx,:) = table(IFD(ai, gi), DFD(ai, gi), meanca(ai, gi), meanstn(ai, gi));
        obssidx = obssidx + 1;
    end
end

% Do a linear fit
disp(fitlm(obss, 'Comm~IFD+DFD'))

% Constrain to observed interval
fprintf('\n\n\nWith restricted IFD [%0.2f %0.2f] and DFD [%0.2f %0.2f] regions\n', par.IFDmeasd, par.DFDmeasd)
ConstRows = obss{:, 'IFD'} > par.IFDmeasd(1) & obss{:, 'IFD'} < par.IFDmeasd(2) & obss{:, 'DFD'} > par.DFDmeasd(1) & obss{:, 'DFD'} < par.DFDmeasd(2);
disp(fitlm(obss(ConstRows, :), 'Perf~IFD+DFD'))

clear ai gi obssidx
