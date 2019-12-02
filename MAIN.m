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
opt=initOpt('inputType', 'origami','template','SingleVertex',...
            'numVert', 4, 'vertexType', '2NC',...
            'restang', pi/4, 'angDesign', [0 120*pi/180 240*pi/180],...
            'analysis','info','analysisType','randomPert',...
            'numIterations', 1000,'RandstDev', 0.2,...
            'createFig', 'on','saveFig','on','saveMovie', 'off',...
            'figDPI',200,'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',10^-2.875,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1,...
            'periodic', 'off');    

tic;
hingeSet = [1 2];
opt.angleConstrFinal(1).val=[ hingeSet(:) , [ones(1,size(hingeSet,2))*pi/2]'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[extrudedUnitCell,opt]=buildGeometry(opt);

opt.saveFile = strcat('/',date,'_',num2str(opt.angDesign*180/pi,'%.2f_'),'temp');
% opt.saveFile = strcat('/24-Jun-2019_vertexDesign2');
opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType,opt.saveFile);

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','off','GradObj','off',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
findDeformation(extrudedUnitCell,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kappas = logspace(-3,1,33);
% close all
% for i = kappas
%     opt.Khinge = i;
    if (strcmp(opt.analysis, 'result') && strcmp(opt.createFig, 'off'))
        fprintf('Not ploting any results.\n');
    else
        ReadAndPlot( extrudedUnitCell, opt);
    end
% end
t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


