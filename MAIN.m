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
            'tessellationType', '25','xrep', 8, 'yrep', 8, 'periodic', 'off',...
            'restang', 3*pi/4, 'angDesign', [0.00 90 180.00 270]*pi/180,...%2.356 %1.571 %0.785
            'analysis','savedata','analysisType','randomPert',...
            'numIterations', 1000,'RandstDev', 0.2,...
            'figDPI',200, 'saveMovie', 'off', 'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off', 'plotAll', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',1,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...%10^-1.5
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1);    

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);
                     
opt.file = '/07-Jul-2020_';
switch opt.analysis
    case{'info'}
        [extrudedUnitCell,opt]=obtainOrigami(opt);
        outputResults(extrudedUnitCell,[],opt,'')
    case{'result'}
        %when using des = 'non' the opt.angDes need to be specified
        %possible des ["2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"];
        des = opt.vertexType;
%         ang = linspace(0,pi,17);
        ang = opt.restang; %ang(2:end-1);%ang(2:end-1);%
%         ang(1:4:end) = [];
%         kap = logspace(-3,0,7);191408
        kap = logspace(-3,0,16);% opt.Khinge; %[kap kap*10 kap*100]; %[10^-3, logspace(-1,0,2)];% 
        xrep = opt.xrep;%opt.xrep;%[1:15];% opt.xrep; %only used when having tessellations
        yrep = opt.yrep;%opt.yrep; %only used when having tessellations
        findDeformation(opt, des, xrep, yrep, ang, kap)
    case{'savedata'}
        %!!!When tessellation, the design angles need to be given!!!!!
%         for i = 1:15
%         opt.file = sprintf('/06-Jul-2020_1_%d_', [i]);
        des = ["2CFF"];
        ReadAndPlot(opt,'oneDes', des) %other option is 'allDes', 'oneDes'
%         end
    case{'plot'}
        opt.sim = 707;   %only used when selecting option 'oneRes'
        PlotResults(opt,'oneRes') %other options is 'allRes','ststRes','oneRes','everyRes'
end


t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


