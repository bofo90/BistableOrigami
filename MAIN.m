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
opt=initOpt('inputType', 'origami','template','Tessellation',...
            'numVert', 4, 'vertexType', "2CFF", ...
            'tessellationType', '25','xrep', 2, 'yrep', 1,...
            'restang', pi/4, 'angDesign', [0 65 130 195]*pi/180,...
            'analysis','result','analysisType','randomPert',...
            'numIterations', 2000,'RandstDev', 0.2,...
            'figDPI',200, 'saveMovie', 'off', 'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off', 'plotAll', 'on',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',10^-2.875,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1,...
            'periodic', 'off');    

tic;
% hingeSet = [1 2];
% opt.angleConstrFinal(1).val=[ hingeSet(:) , [ones(1,size(hingeSet,2))*pi/2]'];

% "CZ", "X1", "3L", "2NC", "2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"
% for i = [ "CFF", "CY"]
%     
%     opt.vertexType = i;
%     
%     for j = 1:10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BUILD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [extrudedUnitCell,opt]=buildGeometry(opt);

%         opt.saveFile = strcat(date,'_temp');
        opt.saveFile = strcat('/10-Dec-2019_',num2str([opt.xrep,opt.yrep],'%d_'));
        % opt.saveFile = strcat('/02-Dec-2019_0.00_ 79.60_159.21_280.40_temp\RestAng_1.571\kappa_10.00000');
        opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType,opt.saveFile);

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
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kappas = logspace(-3,1,33);
% close all
% for i = kappas
%     opt.Khinge = i;
    if (strcmp(opt.analysis, 'result'))
        fprintf('Not ploting any results.\n');
    else
        ReadAndPlot(extrudedUnitCell, opt);
    end
% end
t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


