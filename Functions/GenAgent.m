function [agents, protoagents] = GenAgent(nf, na, diva, divg, no, varargin)
% Agents have par.nofunc functions. This is stored in an array and the agents' capability in any given function is the value in the array at a given
% position. Values can range from 0 to 10, but their sum is normalized to 10. Capabilities are generated by a Gaussian function. It's width specifies
% how universal the agent is: the wider the Gaussian the agent has high values (good capabilities) at multiple functions, if narrower, only at few
% functions. This will set the intra-personal functional diversity value.
%
% For the dominant functional diversity the maximum of the array is to be varied across agents: if multiple functions have the same probability, DFD
% will be high; if only a few functions have high probability, DFD will be low.
%
% Usage: [agents, protoagents] = GenAgent(nf, na, diva, divg, no, ...)
%
% Where nf is the number of functions (9 in the Bunderson and Sutcliffe (2002) paper), na is the number of agents, diva is a two element array
% determining intrapersonal functional diversity score (and ultimately IFD), divg determines dominant functional diversity, no is the normalizing
% factor to agents.
%
% There are optional arguments as well: 'Method' can be 1 or 2. In method #1 the agents generated (protoagents) are about the same, i.e. their IFDSs
% are similar. In method #2 the final IFD is generated by mixing agents with IFDS = 1 and IFDS = 0. Method #2 has a second optional argument,
% 'MixingRatio' that specifies the ratio of hiper-specialists and hiper-generalists in the final team.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Parse input arguments
p = inputParser;
addRequired(p, 'nf');
addRequired(p, 'na');
addRequired(p, 'diva');
addRequired(p, 'divg');
addRequired(p, 'no');
addParameter(p, 'Method', 1, @isnumeric);
addParameter(p, 'MixingRatio', 0.5, @(n)validateattributes(n,{'numeric'}, {'>=',0,'<=',1, 'size', [1,1]}));
addParameter(p, 'IFDSFitRange', [0.1 0.9], @(n)validateattributes(n,{'numeric'}, {'>=',0,'<=',1, 'size', [1,2]}))
addParameter(p, 'DesiredIFD', 0.5, @(n)validateattributes(n,{'numeric'}, {'>=',0,'<=',1, 'size',[1,1]}))
addParameter(p, 'DesiredDFD', NaN, @(n)validateattributes(n,{'numeric'}, {'>=',0,'<=',1, 'size',[1,1]}))
addParameter(p, 'DFDSigmas', 0.1:0.01:10, @isnumeric)
addParameter(p, 'IFDSLevel', 1, @isnumeric)
addParameter(p, 'SigEstMeth', 2, @isnumeric)
addParameter(p, 'SigEstLims', [0 20], @isnumeric)
addParameter(p, 'SigEstStop', 10000, @isnumeric)
addParameter(p, 'SigEstAccuracy', 1e-5, @isnumeric)
addParameter(p, 'SigEstStart', 2, @isnumeric)
addParameter(p, 'DominantFunctions', [], @isnumeric)
addParameter(p, 'MixSkillStrengths', true, @islogical)
addParameter(p, 'Debug', false)

parse(p, nf, na, diva, divg, no, varargin{:});

%% Initialization
% Create "functions" for Method 1
x = repmat(0:(nf-1), na, 1); 

% Assign dominant function to agents
if ~isempty(p.Results.DominantFunctions)
    domf = p.Results.DominantFunctions;
elseif ~isnan(p.Results.DesiredDFD)
    domf = sweepDFD(p.Results.DFDSigmas, p.Results.DesiredDFD);
else
    domf = randsample(x(1,:), na, true, divg); %Find a dominant function for agents
end

%% Generate agents
switch p.Results.Method
    case 1
        % In this version agents are generated such that their skill strengths follow a Gaussian distribution, and IFDS is set by the width of the
        % distribution. All agents have about the same IFDS. To make agents different the strength values are mixed across skills. The highest
        % strength (determining DFD) will be put into a skill position specified by domf
        agents = exp((-(x).^2)./(diva(1)+diva(1)*diva(2)/100*rand(na,1)-diva(2)/200)); %Initialize agents
        agents = agents./sum(agents,2)*no; %Normalize agents to have a total capability of 10

        protoagents = agents; %To be able to generate a number DFDs with given IFD

        %Now some final touches
        for aidx = 1:na

            %Randomly mix up functions so not all agents have the same strength for given functions
            if p.Results.MixSkillStrengths
                agents(aidx,:) = agents(aidx, randperm(length(agents(aidx,:))));
            end

            %And finally, put the strongest function to its place determined by domf
            swp = agents(aidx,domf(aidx)+1);
            [m, l] = max(agents(aidx,:));
            agents(aidx,domf(aidx)+1) = m;
            agents(aidx, l) = swp;
        end

    case 2
        % This version of agent generation mixes full generalists with full specialists to step through IFD values. Here only 0 and 1 IFDS exists, the
        % average of it (=IFD) depends on the ratio of agent types in the mix.
        sagent = [1 zeros(1,nf-1)]*no; %Hiper-specialist
        cagent = ones(1,nf)/nf*no; % Hiper-generalist
        cagent(1) = cagent(1) + eps; %To have a dominant function (even if slightly)

        amr = round(p.Results.MixingRatio*na);
        agents = [repmat(sagent, na-amr, 1); repmat(cagent, amr, 1)];
        protoagents = agents;

        for aidx = 1:na
            swp = agents(aidx,domf(aidx)+1);
            agents(aidx,domf(aidx)+1) = agents(aidx, 1);
            agents(aidx, 1) = swp;
        end

    case 3
        % This agent generation version creates agents such that IFDS of individual agents are different and follow a gaussian around a mean, which is
        % IFD by definition. By changing the width of the IFDS distribution (IFDSLevel), systems of increasingly more different agents can be created 
        % with the same IFD.

        % Define required functions
        pj = @(sigma) exp(-1/2*((1:nf)./sigma).^2); %Determining p_j (skill strenght for skill number j) for a given sigma
        switch p.Results.SigEstMeth
            case 1
                s = @s1;
            case 2
                s = @s2;
            otherwise
                fprintf('Unknown sigma calculating method #%i.\n', p.Results.SigEstMeth)
        end

        % Step 1 -- find sigma values for given IFDSs
        % Divide the smaller fillable region to N/2-1 locations where individual agents can be. # of locs used is determined by IFDSLevel
        reg = min(abs(p.Results.IFDSFitRange-p.Results.DesiredIFD));
        dl = reg/(na/2-1);
        if p.Results.DesiredIFD == 1
            agents = ones(na, nf);
            agents(:, 1) = agents(:, 1) + eps; %Otherwise DFD can only be one
        else
            locs = [p.Results.DesiredIFD - ((na/2-1):-1:1)*dl, p.Results.DesiredIFD, p.Results.DesiredIFD + (1:(na/2-1))*dl]; %These are the IFDS locations I want to put agents. These will not be acurately hit
            sigmas = s(locs); %I need to use these sigmas to target the above IFDS locations. The actual IFDSs generated by these will be populated.
            if p.Results.Debug
                fprintf('Desired IFDS locations:\n')
                disp(locs)
                fprintf('Required sigma values:\n')
                disp(sigmas)
            end

            % Step 2
            % Populate a given team of agents
            if p.Results.IFDSLevel > na/2 || p.Results.IFDSLevel < 1
                error('Cannot use this IFDS Level. It must be at least 1 but less than or equal to half of the number of agents')
            end
            agents = NaN(na, nf);
            agents(1:na-(p.Results.IFDSLevel-1)*2, :) = repmat(pj(sigmas(na/2)), na-(p.Results.IFDSLevel-1)*2, 1);
            for lidx = 1:p.Results.IFDSLevel-1
                agents(na-(p.Results.IFDSLevel-1)*2+(lidx-1)*2+1, :) = pj(sigmas(na/2+lidx));
                agents(na-(p.Results.IFDSLevel-1)*2+lidx*2, :)       = pj(sigmas(na/2-lidx));
            end
        end
        agents = agents./sum(agents,2)*no;
        protoagents = agents;

        %Now some final touches
        for aidx = 1:na

            %Randomly mix up functions so not all agents have the same strength for given functions
            if p.Results.MixSkillStrengths
                agents(aidx,:) = agents(aidx, randperm(length(agents(aidx,:))));
            end

            %And finally, put the strongest function to its place determined by domf
            swp = agents(aidx,domf(aidx)+1);
            [m, l] = max(agents(aidx,:));
            agents(aidx,domf(aidx)+1) = m;
            agents(aidx, l) = swp;
        end

        if p.Results.Debug
            fprintf('Resulting IFDS locations:\n')
            [~, ~, IFDS] = CalcFD(protoagents);
            disp(IFDS')
        end

    otherwise
        fprintf('Unknown agent generating method #%i.\n', p.Results.Method)
end

%% Functions to calculate sigma from desired IFDS
    function sigma = s1(dIFDS) %This version uses a fit I did with nf == 9
        %s = @(ifds) -1*log((ifds-0.9824)/-1.823)/1.07; %To calculate sigma for a desired IFDS -- actually, only tested for nf == 9
        sigma = -1*log((dIFDS-0.9824)/-1.823)/1.07;
    end

    function sigma = s2(dIFDS) %This is the iterative solution, which is much more accurate but slower (~5x).
        epsilon = p.Results.SigEstAccuracy;
        normfac = 1-1/nf;
        sigma = NaN(1,length(dIFDS));
        for iidx = 1:length(dIFDS)
            if dIFDS(iidx) > 1-eps
                sigma(iidx) = inf;
            else
                counter = 1;
                ie = Inf;
                ul = p.Results.SigEstLims(2);
                ll = p.Results.SigEstLims(1);
                sigma(iidx) = p.Results.SigEstStart;
                while ie > epsilon && counter < p.Results.SigEstStop
                    agent = pj(sigma(iidx));
                    cIFDS = (1 - sum((agent./sum(agent)).^2))/normfac;
                    ie = abs(dIFDS(iidx) - cIFDS);
                    if cIFDS < dIFDS(iidx)
                        ll = sigma(iidx);
                        sigma(iidx) = sigma(iidx) + (ul-sigma(iidx))/2;
                    else
                        ul = sigma(iidx);
                        sigma(iidx) = sigma(iidx) - (sigma(iidx)-ll)/2;
                    end
                    counter = counter + 1;
                end
                if counter == p.Results.SigEstStop
                    error('Sigma calculation did not converge.')
                end
            end
        end
    end

%% Function to sweep through DFD values to pick the one nearest to desired and return dominant functions per agent
    function bestdomfu = sweepDFD(sigmas, dDFD)

        % Some init
        DFDdist = Inf;
        bestdomfu = NaN;
        DFD = NaN(length(sigmas), 1);
        ctmp = NaN(length(sigmas), 1);
        maxDFD = 1-1/nf;

        % Swipe through DFD(sigma) function to find best sigma
        for sidx = 1:length(sigmas)
            tmp = exp(-1/2*((0:nf-1)./sigmas(sidx)).^2); %This is an exponential. Could be many other things, for example a parameter
            tmp = round(na*tmp/sum(tmp));
            ctmp(sidx) = sum(tmp);
            if ctmp(sidx) < na
                tmpidx = find(tmp>0, 1, 'last');
                tmp(tmpidx) = tmp(tmpidx) + 1;
                ctmp(sidx) = sum(tmp);
            end
            if ctmp(sidx) > na
                tmpidx = find(tmp>0, 1, 'last');
                tmp(tmpidx) = tmp(tmpidx) - 1;
                ctmp(sidx) = sum(tmp);
            end

            domfu = [];
            for k = 1:length(tmp)
                if tmp(k) > 0
                    domfu = [domfu, repmat(k-1, 1, tmp(k))]; %#ok<AGROW>
                end
            end
            T = tabulate(domfu);
            DFD(sidx) = (1 - sum((T(:,3)/100).^2))/maxDFD;

            tmp = abs(DFD(sidx)-dDFD);
            if tmp < DFDdist
                DFDdist = tmp;
                bestdomfu = domfu;
            end
        end

        if p.Results.Debug
            figure
            plot(DFD, '*-')
            hold on
            plot(ctmp/na, '*-')

            fprintf('Achievable DFD values are:\n')
            disp(unique(DFD))
        end

    end

end
