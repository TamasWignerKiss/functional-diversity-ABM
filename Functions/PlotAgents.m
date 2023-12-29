function fh = PlotAgents(agents)
%This is a helper function, to quickly visualize what agents of a simulation instance look like.
%
%Usage: PlotAgents(agents)
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Get XFD information
[DFD, IFD, IFDS] = CalcFD(agents);

%% Plot agents
fh = figure('WindowStyle','docked', 'Name','Agents');
bh = bar(agents, 'stacked');

set(bh, 'FaceColor', 'Flat')
colors =  mat2cell(distinguishable_colors(numel(bh)),ones(numel(bh),1), 3);
set(bh, {'CData'}, colors) % using jet colormap

xlabel('Agent''s IFDS')
ylabel('Skills')

grid on

set(gca, 'XTickLabel', num2str(IFDS, '%0.4f'))
legend(strcat("Skill #", num2str((1:numel(bh))')))
title(['DFD=' num2str(DFD, '%0.4f') ' | IFD=' num2str(IFD, '%0.4f')])
