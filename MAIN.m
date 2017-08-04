%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Original Matlab version used to write the code: R2014b
%
%FOR ANY QUESTION RELATED TO THE MATLAB FILES, PLEASE CONTACT JOHANNES
%OVERVELDE AT J.T.B.OVERVELDE@GMAIL.COM (WWW.OVERVELDE.COM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all , clc, format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CHOOSE PREDEFINED GEOMETRY, SIMULATION AND PLOT OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SimCase=1

switch SimCase
    case 1
        opt=initOpt('inputType','individual',...
                    'template','truncated tetrahedron',...
                    'plot','result',...
                    'interval', 1,'saveFig','on','periodic','on','figDPI',300,...
                    'saveMovie', 'on', 'safeMovieAntiAlias', 0,...
                    'closeAlgor', 'interior-point','relAlgor', 'sqp',...
                    'gradDescStep', 1e-1, 'gradDescTol', 1e-9,...
                    'constrFace','off','constrEdge','off',...
                    'Khinge',0.0005,'Kedge',1,'Kface',1,'KtargetAngle',0.5,...
                    'relInterval', 1, 'constAnglePerc',0.99);
        opt.angleConstrFinal(1).val=[1  -pi*0.985
                                     2  -pi*0.985
                                     14  -pi*0.985
                                     ];
        opt.angleConstrFinal(2).val=[1  -pi*0.985
                                     2  -pi*0.985
                                     19  -pi*0.985
                                     ];
        opt.angleConstrFinal(3).val=[1  -pi*0.985
                                     2  -pi*0.985
                                     8  -pi*0.985
                                     ];
        opt.angleConstrFinal(4).val=[2  -pi*0.985
                                     13  -pi*0.985
                                     41  -pi*0.985
                                     ];
        opt.angleConstrFinal(5).val=[3  -pi*0.985
                                     8  -pi*0.985
                                     ];
        opt.angleConstrFinal(6).val=[3  -pi*0.985
                                     19  -pi*0.985
                                     37  -pi*0.985
                                     ];
        opt.angleConstrFinal(7).val=[3  -pi*0.985
                                     41  -pi*0.985
                                     ];
        opt.angleConstrFinal(8).val=[3  -pi*0.985
                                     ];
        opt.angleConstrFinal(9).val=[24  -pi*0.985
                                     48  -pi*0.985
                                     ];
        opt.angleConstrFinal(10).val=[36  -pi*0.985
                                     48  -pi*0.985
                                     ];

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

%SOLVER OPTIONS
opt.options=optimoptions('fmincon','GradConstr','on','GradObj','on',...
                         'tolfun',1e-5','tolx',1e-9,'tolcon',1e-5,...
                         'Display','off','DerivativeCheck','off',...
                         'maxfunevals',100000,'Algorithm', opt.closeAlgor);
%                          'FiniteDifferenceType', 'central', 'FiniteDifferenceStepSize', eps^(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BUILD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unitCell,extrudedUnitCell,opt]=buildGeometry(opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[extrudedUnitCell]=findDeformation(unitCell,extrudedUnitCell,opt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OUTPUT AND PLOT GEOMETRY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileFolder = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat');
allFiles = dir(fileFolder);
for ct = 1:length(allFiles)
    if allFiles(ct).isdir
%         disp('skip all directories...')
        continue;
    end
    
    % parse the file name to get back hinge set and exit flag
    fileName = allFiles(ct).name;
    parsedName = strsplit(fileName(1:end-4), '_');
    hingeSetStr = parsedName{2};
    hingeSetStr = strsplit(hingeSetStr(2:end-1), ' ');
    hingeSet = str2double(hingeSetStr)';
    extrudedUnitCell.angleConstr = [hingeSet(:), -pi*0.985 * ones(length(hingeSet), 1)];
    load(strcat(fileFolder,'/', fileName));
    outputResults(unitCell,extrudedUnitCell,result,opt);
    close all;
end

