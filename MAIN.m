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
            'angDesign', [0 120*pi/180 240*pi/180],...
            'restang', pi/4, 'numVert', 4, 'numIterations', 1000,'RandstDev', 0.2,...
            'analysis','result','analysisType','randomPert',...
            'createFig', 'off','saveFig','on','saveMovie', 'off',...
            'figDPI',200,'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',10^-2.875,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1,...
            'periodic', 'off');    


% opt.saveFile = strcat('/',date,'_temp');
% opt.saveFile = strcat('/24-Jun-2019_vertexDesign2');
tic;


hingeSet = [757 306];
opt.angleConstrFinal(1).val=[ hingeSet(:) , [ones(1,size(hingeSet,2))*opt.restang]'];
for i = 90:5:90
    for j = i:5:90

    
%         if (i+j)<180 || ((i <=90 && i >=35) && (i+j) <195)
%             continue;
%         end
        opt.angDesign = [0, i, 180, 360-j]*pi/180;
        opt.saveFile = strcat('/',date,'_KandTheta0/Angles_',num2str(i),'_',num2str(j));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %BUILD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [extrudedUnitCell,opt]=buildGeometry(opt);

        %SOLVER OPTIONS
        opt.options=optimoptions('fmincon','GradConstr','off','GradObj','off',...
                                 'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                                 'Display','off','DerivativeCheck','off',...
                                 'maxfunevals',30000, 'MaxIterations', 2000,...
                                 'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                                 'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SELECT HINGES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % selectHinges(unitCell, extrudedUnitCell, opt);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %ANALYSIS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        findDeformation(extrudedUnitCell,opt);
    
    end
end

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


