function fh = PlotRepeatEvents(rn, tn, t2a, th, par)
% PlotRepeatEvents is a helper function to help see how task passing and task subfunctions changed during a given instance of the simulation.
%
% Author: Tamas Kiss <kiss.t@wigner.hu>

%% Plot evolution of task assigment and task functions
fh = figure;
ax1 = subplot(2,1,1);
plot(t2a(:,:, rn), '.-');
legend(strcat("task #", num2str((1:par.numtasks)')))
title(['Repeat #' num2str(rn)])
xlim([0, par.EmergencyStop])
grid on

% An example to plot how task #taskNum evolved in repeat repeatNum
ax2 = subplot(2,1,2);
plot(squeeze(th(tn, :, :, rn))', '.-')
legend(strcat('Func. #', num2str((1:par.numfuncs)')))
title(['Repeat #' num2str(rn) ', Task #' num2str(tn)])
xlim([0, par.EmergencyStop])
grid on

linkaxes([ax1, ax2], 'x')

end