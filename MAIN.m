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
SimCase=1;

switch SimCase
    case 1
        opt=initOpt('inputType','individual', 'plot','result',... 
                    'template','cube','onlyUnitCell', 'off',...
                    'createFig', 'on', 'showFig', 'on','saveFig','off','saveMovie', 'off',...
                    'interval', 1,'relInterval', 1,'constAnglePerc',nan,...
                    'periodic','off','figDPI',300,'safeMovieAntiAlias', 0,...
                    'folAlgor', 'sqp','relAlgor', 'sqp',...
                    'gradDescStep', 1e-1, 'gradDescTol', 1e-9,...
                    'readAngFile', 'off','gethistory', 'on',...
                    'constrFace','on','constrEdge','off',...
                    'Khinge',0.01,'Kedge',10^0.5,'Kface',1,'KtargetAngle',1,...
                    'stepkHinge', 1, 'stepkTargetAngle', 3, 'stepkEdge', 1,...
                    'stepMaxStrech', 1, 'maxStretch', nan);

%         opt.saveFile = strcat('/',date,'_temp');
        opt.saveFile = '/13-Dec-2017_noAngleCnstr';

        %-pi if its the extruded version, pi if its
        %only the internal polyheron
        hingeSet = [3 8 13 17 21 26 22 30];%[7 3 8 13 17 26 22 30];
        if strcmp(opt.onlyUnitCell, 'on')
            opt.angleConstrFinal(1).val=[ hingeSet(:) , (pi*0.985) *ones(length(hingeSet), 1)];
        else
            opt.angleConstrFinal(1).val=[ hingeSet(:) , (-pi*0.985) *ones(length(hingeSet), 1)];
        end

        
        
    case 2
        opt=initOpt('inputType','individual',...
                    'template','cuboctahedron',...
                    'plot','result',...
                    'interval',1,'saveFig','off','periodic','off',...
                    'constrFace','off','constrEdge','off',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.98);
        opt.angleConstrFinal(1).val=[16 -pi*0.975
                                     18 -pi*0.975
                                     30 -pi*0.975
                                     53 -pi*0.975];        
    case 3        
        opt=initOpt('inputType','preDefined',...
                    'template','#22',...
                    'plot','result',...
                    'interval',20,'saveFig','off','periodic','on',...
                    'constrFace','on','constrEdge','on',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.98); 
        opt.angleConstrFinal(1).val=[1  -pi*0.975];
        opt.angleConstrFinal(2).val=[1  -pi*0.975
                                     13 -pi*0.975];
    case 4
        opt=initOpt('inputType','preDefined',...
                    'template','#26',...
                    'plot','result',...
                    'interval',20,'saveFig','off','periodic','on',...
                    'constrFace','on','constrEdge','on',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.98);  
        opt.angleConstrFinal(1).val=[1  -pi*0.975];
    case 5
        opt=initOpt('inputType','individual',...
                    'template','rhombicuboctahedron',...
                    'plot','result',...
                    'interval',5,'saveFig','off','periodic','off',...
                    'constrFace','on','constrEdge','off',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.99);
        opt.angleConstrFinal(1).val=[111 -pi*0.985
                                     105 -pi*0.985
                                     116 -pi*0.985
                                     106 -pi*0.985];
        opt.angleConstrFinal(2).val=[131 -pi*0.985
                                     125 -pi*0.985
                                     126 -pi*0.985
                                     136 -pi*0.985];
        opt.angleConstrFinal(3).val=[];
    case 6
        opt=initOpt('inputType','preDefined',...
                    'template','#18',...
                    'plot','result',...
                    'interval',10,'saveFig','off','periodic','on',...
                    'constrFace','on','constrEdge','on',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.98);  
        opt.angleConstrFinal(1).val=[1  -pi*0.975];
    case 7
        opt=initOpt('inputType','user',...
                    'template','userInputFile1',...
                    'plot','result',...
                    'interval',1,'saveFig','off','periodic','on',...
                    'constrFace','on','constrEdge','off',...
                    'Khinge',0.001,'Kedge',1,'Kface',1,'KtargetAngle',1,...
                    'constAnglePerc',0.98);  
        opt.angleConstrFinal(1).val=[19  -pi*0.6
                                     47  -pi*0.6]; 
        opt.angleConstrFinal(2).val=[];
end


tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-5,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 1000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun);
%                          'RelLineSrchBnd', 0.1, 'RelLineSrchBndDuration', 10e10,...
%                          'TypicalX', 0.1*ones(length(extrudedUnitCell.node(:,1))*3,1),...    

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



