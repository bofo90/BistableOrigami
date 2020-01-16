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
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE PREDEFINED GEOMETRY, SIMULATION AND PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt=initOpt('template','Tessellation','numVert', 4,'vertexType', "2CFF",...
            'tessellationType', '25','xrep', 2, 'yrep', 2,...
            'restang', 0.785, 'angDesign', [0.00 90 180.00 270]*pi/180,...
            'analysis','plot','analysisType','randomPert',...
            'numIterations', 1000,'RandstDev', 0.2,...
            'figDPI',200, 'saveMovie', 'off', 'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off', 'plotAll', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',0.01,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1,...
            'periodic', 'off');    

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','off','GradObj','off',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);
                     
opt.file = '/09-Jan-2020_2_2_';
switch opt.analysis
    case{'info'}
        [extrudedUnitCell,opt]=obtainOrigami(opt);
        outputResults(extrudedUnitCell,[],opt,'')
    case{'result'}
        %when using des = 'non' the opt.angDes need to be specified
        %possible des ["2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"];
        des = ["2C"];
        ang = linspace(0,pi,5);
        ang = ang(2:4);
        kap = logspace(-3,0,4);
        xrep = 1:5; %only used when having tessellations
        yrep = 1:5; %only used when having tessellations
        findDeformation(opt, des, xrep, yrep, ang, kap) 
    case{'savedata'}
        %!!!When tessellation, the design angles need to be given!!!!!
        des = ["2CFF"];
        ReadAndPlot(opt,'oneDes', des) %other option is 'allDes', 'oneDes'
    case{'plot'}
        opt.sim = 367;   %only used when selecting option 'oneRes'
        PlotResults(opt,'oneRes') %other options is 'allRes','ststRes','oneRes'
end


t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


