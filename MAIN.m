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
opt=initOpt('template','truncated tetrahedron','analysis','result','readHingeFile','off',...
            'createFig', 'off','saveFig','on','saveMovie', 'on',...
            'figDPI',200,'safeMovieAntiAlias', 0,...
            'folAlgor', 'active-set','relAlgor', 'active-set',...
            'gethistory', 'on',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',0.001,'Kedge',3,'Kdiag',3,'Kface',100,'KtargetAngle',100,...
            'maxStretch', 0.3,...
            'maxHinges',inf,'minHinges',0);    %Only work when readHingeFile is 'on'

opt.saveFile = strcat('/',date,'_EnergyAllAngles_3_24');
% opt.saveFile = strcat('/','13-Feb-2018_EnergyLandPath_24to3');

hingeSet = [3 24];
opt.angleConstrFinal(1).val=[ hingeSet(:) , (-pi*(opt.constAnglePerc-0.005)) *ones(length(hingeSet), 1)];

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.025, 'RelLineSrchBndDuration', 500);

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

if (strcmp(opt.analysis, 'result') && strcmp(opt.createFig, 'off'))
    fprintf('Not ploting any results.\n');
else
    ReadAndPlot(unitCell, extrudedUnitCell, opt);
end

t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


