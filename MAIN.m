%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Matlab version used to write the code: R2014b
%
%FOR ANY QUESTION RELATED TO THE MATLAB FILES, PLEASE CONTACT JOHANNES
%OVERVELDE AT J.T.B.OVERVELDE@GMAIL.COM (WWW.OVERVELDE.COM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc
format long
clearvars -global
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE PREDEFINED GEOMETRY, SIMULATION AND PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt=initOpt('inputType','individual', 'plot','result',... 
            'template','triangular prism','onlyUnitCell', 'off',...
            'createFig', 'on','saveFig','on','saveMovie', 'off',...
            'interval', 1,'relInterval', 1,'constAnglePerc',0.999,...
            'periodic','off','figDPI',200,'safeMovieAntiAlias', 0,...
            'folAlgor', 'active-set','relAlgor', 'active-set',...
            'gradDescStep', 1e-1, 'gradDescTol', 1e-9,...
            'readAngFile', 'on','gethistory', 'on',...
            'constrFace','on','constrEdge','off',...
            'Khinge',0.01,'Kedge',4,'Kface',100,'KtargetAngle',1,...
            'stepkHinge', 1, 'stepkTargetAngle', 3, 'stepkEdge', 1,...
            'stepMaxStrech', 1, 'maxStretch', nan,'maxHinges',inf,'minHinges',0);

%opt.saveFile = strcat('/',date,'_temp');
opt.saveFile = strcat('/','_Pres_Ventura')
%         opt.saveFile = '/13-Dec-2017_noAngleCnstr';

%-pi if its the extruded version, pi if its
%only the internal polyheron
%hingeSet = [3 8 13 17 21 26 22 30]; %Multistable cube
hingeSet = [3 8 12 16 17 21];
hingeSet = [1 2 3 13 41];
hingeSet = [1 7 2 16 21 17 8 12];
% 1 2 13 23 36
% 1 2 13 37 41
% 1 2 19 37
% 1 2 3 13 14
% 1 2 3 14
% 1 2 3 37
% 1 2 7 13 41
% 1 2 7 28 37
% 1 3 18 37
% 1 3 19
% 1 3 24 37 41
% 1 3 7 13 41
% 1 3 7 28 37
% 1 3 9 24 36
% 1 3 9 28 37
% 1 9 18 24
% 3 9 32];
%hingeSet = [1];
if strcmp(opt.onlyUnitCell, 'on')
    opt.angleConstrFinal(1).val=[ hingeSet(:) , (pi*opt.constAnglePerc) *ones(length(hingeSet), 1)];
else
    opt.angleConstrFinal(1).val=[ hingeSet(:) , (-pi*opt.constAnglePerc) *ones(length(hingeSet), 1)];
end

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-8,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 1000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...
                         'RelLineSrchBnd', 0.025, 'RelLineSrchBndDuration', 10e10);

%                          'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', eps^(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SELECT HINGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selectHinges(unitCell, extrudedUnitCell, opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findDeformation(unitCell,extrudedUnitCell,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(opt.plot, 'result') && strcmp(opt.createFig, 'off'))
    fprintf('Not ploting any results.\n');
else
    ReadAndPlot(unitCell, extrudedUnitCell, opt);
end

t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


