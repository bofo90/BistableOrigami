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
opt=initOpt('inputType', 'origami','template','SquareTiling',...
            'numUnitCell', 3,'ratio', 2,'layers', 3,...
            'analysis','plot','readHingeFile','on',...
            'createFig', 'off','saveFig','off','saveMovie', 'off',...
            'figDPI',200,'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',0.001,'Kedge',10,'Kdiag',10,'Kface',100,'KtargetAngle',100,...
            'maxStretch', 0.3,'steps',20,...
            'maxHinges',3,'minHinges',0,... %Only work when readHingeFile is 'on'
            'periodic', 'on');    

% opt.saveFile = strcat('/',date,'_Mode1Gif_2');
opt.saveFile = strcat('/','09-Jul-2019_Mode2Gif_2');
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[extrudedUnitCell,opt]=buildGeometry(opt);

% beta = 14+12*(round(opt.numUnitCell/2)-1)+12*opt.numUnitCell*(opt.layers-1);
% alpha = beta + 2 + 8*(opt.layers-1);

hingeSet = [extrudedUnitCell.alpha extrudedUnitCell.beta];
alpha_value = 1.55;
beta_value = -1.55;%0.02;%pi/2;%-pi/2*0.985
opt.angleConstrFinal(1).val=[ hingeSet(:) , [ones(1,size(extrudedUnitCell.alpha,2))*alpha_value...
                                             ones(1,size(extrudedUnitCell.beta,2))*beta_value]'];

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(opt.analysis, 'result') && strcmp(opt.createFig, 'off'))
    fprintf('Not ploting any results.\n');
else
    ReadAndPlot( extrudedUnitCell, opt);
end

t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);


