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
            'numVert', 4, 'vertexType', "2CFF", ...
            'tessellationType', '25','xrep', 2, 'yrep', 1,...
            'restang', pi/4, 'angDesign', [0 65 130 195]*pi/180,...
            'analysis','savedata','analysisType','randomPert',...
            'numIterations', 1000,'RandstDev', 0.2,...
            'figDPI',200, 'saveMovie', 'off', 'safeMovieAntiAlias', 0,...
            'folAlgor', 'sqp','relAlgor', 'sqp',...
            'gethistory', 'off', 'plotAll', 'off',...
            'constrEdge','off','constrFace','on','constAnglePerc',0.99,... 
            'Khinge',10^-2.875,'Kedge',1,'Kdiag',1,'Kface',100,'KtargetAngle',1000,...
            'maxStretch', nan,'steps',3,...
            'maxArea', 0.1,...
            'periodic', 'off');    

tic;
% hingeSet = [1 2];
% opt.angleConstrFinal(1).val=[ hingeSet(:) , [ones(1,size(hingeSet,2))*pi/2]'];

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','off','GradObj','off',...
                         'tolfun',1e-6,'tolx',1e-9, 'tolcon',1e-9,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',30000, 'MaxIterations', 2000,...
                         'Algorithm', opt.folAlgor, 'OutputFcn',@outfun,...               
                         'RelLineSrchBnd', 0.01, 'RelLineSrchBndDuration', 5000);
                     
% "CZ", "X1", "3L", "2NC", "2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"

% opt.analysisType = 'randomPert1';
% for i = ["GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"]
%     
%     opt.vertexType = i;
%     
%     [extrudedUnitCell,opt]=buildGeometry(opt);
%     opt.saveFile = strcat('/02-Dec-2019_',num2str(opt.angDesign*180/pi,'%.2f_'));
%     opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType,opt.saveFile);
%     
%     findDeformation(extrudedUnitCell,opt);
%     
% end

% opt.analysisType = 'randomPert2';
% for i = ["CFF", "CY", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"]
%     
%     opt.vertexType = i;
%     for j = 1:10
%         [extrudedUnitCell,opt]=buildGeometry(opt);
%         opt.saveFile = strcat('/06-Dec-2019_',num2str(opt.angDesign*180/pi,'%.2f_'));
%         opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType,opt.saveFile);
% 
%         findDeformation(extrudedUnitCell,opt);
%     end
% end


% opt.analysisType = 'randomPert3';
% opt.template = 'Tessellation';
% opt.vertexType = '2CFF';
% for i = 2:5
%     
%     opt.xrep = i;
%     
%     for j = 1:5
%         
%         opt.yrep = j;
%         
%         if i*j >10
%             opt.numIterations = 10000;
%         else
%             opt.numIterations = i*j*1000;
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %BUILD
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         [extrudedUnitCell,opt]=buildGeometry(opt);
%         opt.saveFile = strcat('/10-Dec-2019_',num2str([opt.xrep,opt.yrep],'%d_'));
%         opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType,opt.saveFile);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %ANALYSIS
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         findDeformation(extrudedUnitCell,opt);
%     end
% end
% 
% opt.yrep = 1;
% for i = 6:10
%     opt.xrep = i;
%     opt.numIterations = i*j*1000;
% 
%     [extrudedUnitCell,opt]=buildGeometry(opt);
%     opt.saveFile = strcat('/10-Dec-2019_',num2str([opt.xrep,opt.yrep],'%d_'));
%     opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType,opt.saveFile);
% 
%     findDeformation(extrudedUnitCell,opt);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%When using savedata
for i = ["CZ", "X1", "3L", "2NC", "2C", "GFF", "2OFF","3S", "2OM1", "2OM2", "2NM1", "2NM2", "2OL", "2NL", "2OS", "2NS", "Z1", "Z2", "Y1", "Y2", "X2"]
% for i = ["2CFF"]    
    
    fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',i);
    allFiles = dir(fileContainer);
    for ct = 3:length(allFiles)
        opt.saveFile = strcat('/',allFiles(ct).name);
        opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',i,opt.saveFile);
        file = opt.file;
        
        opt.vertexType = 'non';
        opt.angDesign = getAngles(opt.saveFile);
        [extrudedUnitCell,opt]=buildGeometry(opt);

        angles = linspace(0,pi,5);
        for j = angles(2:4)
            opt.restang = j;
            extrudedUnitCell.theta = ones(size(extrudedUnitCell.theta,1),1)*opt.restang;
            
            allFiles2 = dir(strcat(file, sprintf('/RestAng_%.3f', opt.restang)));
            for k = 3:length(allFiles2)
                opt.Khinge = getKappa(allFiles2(k).name);
                opt.file = strcat(file,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));
                fprintf(strcat('Saving data: ',allFiles(ct).name, sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge), '\n'))
                ReadAndPlot(extrudedUnitCell, opt);
            end
        end
    end
end

%%%%When using plot !!carefull, it is still not so well estabilshed!!
% for i = ["2CFF"]
% 
%     opt.saveFile = '/02-Dec-2019_0.00_90.00_180.00_270.00_';
%     
%     opt.vertexType = 'non';
%     opt.angDesign = getAngles(opt.saveFile);
%     [extrudedUnitCell,opt]=buildGeometry(opt);
%     
%     opt.origin = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',i,opt.saveFile);
%     
%     allFiles = csvread(strcat(opt.origin,'/Images/InfoforStableStatesImages.csv'),1);
%     
%     for j = 1:size(allFiles,1)
%         close all
%         opt.Khinge = allFiles(j,1);
%         opt.restang = allFiles(j,4);
%         opt.sim = allFiles(j,2);
%         opt.file = strcat(opt.origin,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));
%         ReadAndPlot(extrudedUnitCell, opt);
%     end
% end

t = toc;
fprintf('The whole program lasted %.2f seconds\n', t);

function Angles = getAngles(fileName)
    parsedName = strsplit(fileName(1:end), '_');
    Angl1 = str2double(parsedName{3});
    Angl2 = str2double(parsedName{4});
    Angl3 = str2double(parsedName{5});
%     Angles = [0, Angl1, 180, 360-Angl2]*pi/180;
    Angles = [0, Angl1, Angl2, Angl3]*pi/180;
end

function kappa = getKappa(fileName)
    parsedName = strsplit(fileName(1:end), '_');
    kappa = str2double(parsedName{2});
end
