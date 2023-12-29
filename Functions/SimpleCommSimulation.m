function [IFD, DFD, meansn, meanstn, meanca, cgIFD, cgDFD, cgmsn, cgmstn, cgmca] = SimpleCommSimulation(par)
%This function is a literal copy of the core content of the ATRAN_SimpleComm.m script to save on typing and the size of the paramter testing script.
%Its goal is to take paramters, run the simulation and return simulation results such as the amount of tasks solved, the time it took to get to the
%end of tasks (though this is of less usability) and the amount of communication (task passing) performed.
%
%Usage: [IFD, DFD, meansn, meanstn, meanca, cgIFD, cgDFD, cgmsn, cgmstn, cgmca] = SimpleCommSimulation(par) 
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Some initialization
adivvals = 10.^(linspace(0.1, 0.2, par.noavals).*(linspace(-6, 7, par.noavals)));
gdivvals = 10.^(0.2*linspace(-3, 5, par.nogvals));
DFD = NaN(par.noavals, par.nogvals);
IFD = NaN(par.noavals, par.nogvals);
meansn = NaN(par.noavals, par.nogvals);
meanstn = NaN(par.noavals, par.nogvals);
meanca = NaN(par.noavals, par.nogvals);
etc.maxPass = par.maxPass;
etc.similThresh = par.similThresh;
etc.emStop = par.EmergencyStop;
etc.passCost = par.passCost; %Steps
cll = matlab.desktop.commandwindow.size;
cll = cll(1)-4;

%% Loop through DFD and IFD values
ocl = 0;
for ai = 1:length(adivvals)
    adiv = adivvals(ai);
    oc = [num2str(ai) ' '];
    fprintf('%s', oc)
    ocl = ocl + length(oc);
    if ocl >= cll
        fprintf('\n')
    end
    
    %The following piece of code runs in parallel on multiple cores. Remove "par" to allow debugging of the inside.
    parfor gi = 1:length(gdivvals) %This is for parallel computing
        gdiv = gdivvals(gi);

        %Generate agents
        agents = GenAgent(par.numfuncs, par.numagents, [adiv, par.agspread], exp((-(0:par.numfuncs-1).^2)/gdiv), par.anorm); %#ok<PFBNS> 
        agdistmat = squareform(pdist(agents));

        %Calculate and store diversity values
        [DFD(ai, gi), IFD(ai, gi)] = CalcFD(agents);
        
        %The following variables will only live for an [ai, gi] point
        t2a = NaN(par.EmergencyStop+1, par.numtasks, par.numrepeats); %Initialize variable keeping track of which agent works on which task
        taskhists = NaN(par.numtasks, par.numfuncs, par.EmergencyStop+1, par.numrepeats); %Initialize variable keeping track of how tasks progress
        sn = NaN(par.numrepeats, 1); %Initialize number of steps nedded to solve all tasks
        stn = zeros(par.numrepeats, 1); %Initialize variable storing the amount of work done in repeats
        ca = NaN(par.numrepeats, 1); %Init variable for communication amount

        for ridx = 1:par.numrepeats
            
            %Generate tasks
            tasks = GenTask(par.numfuncs, par.numtasks, par.tnorm, par.EmergencyStop, par.TaskType, agents);

            %PassingSolveTasks is the function that runs the while loop running the whole game of passing-working
            [t2a(:, :, ridx), taskhists(:, :, :, ridx)] = PassingSolveTasks(agents, tasks, etc, agdistmat); %This returns number of steps required to solve all tasks

            % Calculate steps taken to solve all tasks
            tmp = find(all(isnan(t2a(:, :, ridx)), 2), 1, 'first')';
            if isempty(tmp)
                sn(ridx) = par.EmergencyStop;
            else
                sn(ridx) = tmp;
            end

            % Calculate the amount of solved sub-tasks
            tmp1 = taskhists(:, :, 1, ridx);
            tmp2 = taskhists(:, :, par.EmergencyStop+1, ridx);
            tmp1(isnan(tmp1)) = 0; %There can be NaNs initially due to sub-tasks that the system cannot solve!
            tmp2(isnan(tmp2)) = 0;
            stn(ridx) = sum(sum(tmp1 - tmp2));

            % Calculate amount of communication
            tmp = diff(t2a(:, :, ridx));
            ca(ridx) = sum(sum(tmp ~= 0 & not(isnan(tmp))));

        end
        
        % Average steps
        meansn(ai, gi) = mean(sn);

        % Average solved subtasks
        meanstn(ai, gi) = mean(stn);

        % Average communication amount
        meanca(ai, gi) = mean(ca);

    end
end

% Interim cleanup
clear adiv agents ai gdiv gi sn stn T t2a T_old taskhists

%% Coarse-grain (average in bins) results and normalize max of amount of solved subtasks to 1
cgmstn = NaN(length(par.IFDbins)-1, length(par.DFDbins)-1);
cgmsn = NaN(size(cgmstn));
cgmca = NaN(size(cgmstn));
cgIFD = NaN(size(cgmstn));
cgDFD = NaN(size(cgmstn));
for iidx = 1:length(par.IFDbins)-1
    for didx = 1:length(par.DFDbins)-1
        cgmstn(iidx, didx) = mean(meanstn(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgmsn(iidx, didx)  = mean(meansn(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgmca(iidx, didx)   = mean(meanca(par.IFDbins(iidx)<=IFD & IFD<par.IFDbins(iidx+1) & par.DFDbins(didx)<=DFD & DFD<par.DFDbins(didx+1)));
        cgIFD(iidx, didx)  = par.IFDbins(iidx);
        cgDFD(iidx, didx)  = par.DFDbins(didx);
    end
end

cgmstn = cgmstn/(par.numtasks*par.tnorm);

% Interim cleanup
clear iidx didx

