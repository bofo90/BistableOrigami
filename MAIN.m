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
%%%%%% to use when camview variable is going to be used
% close all
% clc
% format long
% clearvars -global poop
% clearvars -except camview camLightPos
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE PREDEFINED GEOMETRY, SIMULATION AND PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt=initOpt('template','SingleVertex','numVert', 4,'vertexType', "2CFF",...
            'tessellationType', '25','xrep', 4, 'yrep', 4, 'periodic', 'off',...
            'restang', 1.933, 'angDesign', [0.00 90 180.00 270]*pi/180,...%2.356 %1.571 %0.785
            'analysis','info','analysisType','randomPert',...
            'numIterations', 1000,'RandstDev', 0.2,...
            'figDPI',200, 'saveMovie', 'off', 'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off', 'plotAll', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',0.03162 ,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...%10^-1.5
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1);    

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);
                     
opt.file = '/19-Mar-2020_0.00_ 90.00_180.00_270.00_';
switch opt.analysis
    case{'info'}
        [extrudedUnitCell,opt]=obtainOrigami(opt);
        outputResults(extrudedUnitCell,[],opt,'')
    case{'result'}
        %when using des = 'non' the opt.angDes need to be specified
        %possible des ["2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"];
        des = opt.vertexType;
        ang = linspace(0,pi,5);
        ang = ang(2:end-1);%opt.restang; %ang(2:end-1);%
%         ang(1:4:end) = [];
%         kap = logspace(-3,0,7);191408
        kap = opt.Khinge; %logspace(-3,0,16);% [kap kap*10 kap*100]; %[10^-3, logspace(-1,0,2)];% 
        xrep = [12:15];%opt.xrep;%opt.xrep;% opt.xrep; %only used when having tessellations
        yrep = opt.yrep;%opt.yrep; %only used when having tessellations
        findDeformation(opt, des, xrep, yrep, ang, kap)
    case{'savedata'}
        %!!!When tessellation, the design angles need to be given!!!!!
        for i = 12:15
        opt.file = sprintf('/08-May-2020_%d_%d_', [i,i]);
        des = ["2CFF"];
        ReadAndPlot(opt,'oneDes', des) %other option is 'allDes', 'oneDes'
        end
    case{'plot'}
        %%% camview is a variable that can be used to define the camera position for
        %%% a nice view of the structure. Needs to be defined by the user. You can
        %%% use the following code to get it from an already created image.
        %%% camview = get(gca,{'CameraPosition', 'CameraTarget', 'CameraUpVector', 'CameraViewAngle'})
        %%% camlight1 = findobj(gcf,'Type','Light'); camLightPos = camlight1.Position
%         opt.camview = camview;
%         opt.camlight = camLightPos;
        
        opt.sim = 134;   %only used when selecting option 'oneRes'
        PlotResults(opt,'oneRes') %other options is 'allRes','ststRes','oneRes','everyRes'
end


t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


