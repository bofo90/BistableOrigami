function PlotResults(opt, flag)

if strcmp(opt.template,'Tessellation')
    fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType);
elseif strcmp(opt.template,'SingleVertex')
    fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType);
end
opt.file = strcat(fileContainer,opt.file);

opt.vertexType = 'non';
if strcmp(opt.template,'Tessellation')
%     [opt.xrep, opt.yrep] = getTessell(opt.file);
elseif strcmp(opt.template,'SingleVertex')
    opt.angDesign = getAngles(opt.file);
end
[extrudedUnitCell,opt]=obtainOrigami(opt);

switch flag
    case {'allRes', 'ststRes'}
        if strcmp(flag,'allRes')
            allFiles = csvread(strcat(opt.file,'/Images/InfoforAllImages.csv'),1);
        else
            allFiles = csvread(strcat(opt.file,'/Images/InfoforStableStatesImages.csv'),1);
        end

        for j = 1:size(allFiles,1)
            close all            
            opt.Khinge = allFiles(j,1);
            opt.restang = allFiles(j,4);
            opt.sim = allFiles(j,2);
            opt.StSt = allFiles(j,3);
            SelectResult(extrudedUnitCell, opt, true);
        end

    case {'oneRes'}
        opt.StSt = 0;       
        SelectResult(extrudedUnitCell, opt, true)
        
    case {'everyRes'}  
        opt.StSt = 0;
        opt.sim = 0;
        SelectResult(extrudedUnitCell, opt, false)
end


function SelectResult(extrudedUnitCell,opt, check)
%get file of results

fileResults = strcat(opt.file,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));

folderResults = strcat(fileResults,'/mat');
if ~exist(folderResults, 'dir')
    error('\n----------\nNo data directory found\n----------\n',[])
end 
allFiles = dir(folderResults);

for ct = 1:length(allFiles)
    if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:end-4), 'metadata')
        % skip all directories and metadata file
        continue;
    end
    
    resfilename = allFiles(ct).name;
    if check
        hingeSet = getSimulation(resfilename);
    else
        hingeSet = 0;
    end
    if (hingeSet ~= opt.sim)
        continue;
    end

    lofile = load(strcat(folderResults,'/', resfilename));
    outputResults(extrudedUnitCell,lofile.result,opt,strcat(sprintf('RestAng_%.3f_kappa_%2.5f_', opt.restang, opt.Khinge),resfilename(1:end-4)));
end
    
function Angles = getAngles(fileName)
parsedName = strsplit(fileName(1:end), '_');
Angl1 = str2double(parsedName{end-3});
Angl2 = str2double(parsedName{end-2});
Angl3 = str2double(parsedName{end-1});
Angles = [0, Angl1, Angl2, Angl3]*pi/180;

function [x,y] = getTessell(fileName)
parsedName = strsplit(fileName(1:end), '_');
x = str2double(parsedName{end-2});
y = str2double(parsedName{end-1});

function hingeSet = getSimulation(fileName)
parsedName = strsplit(fileName(1:end-4), '_');
Angl1 = parsedName{2};
hingeSet = str2double(Angl1);