function ReadAndPlot(extrudedUnitCell, opt)

switch opt.analysis
    case 'info'
        result = [];
        outputResults(extrudedUnitCell,result,opt);
    case {'result', 'savedata', 'plot'}
        %get file of results
        extraName = sprintf('/kh%2.3f_kta%2.3f_ke%2.3f_kf%2.3f', opt.Khinge,opt.KtargetAngle,opt.Kedge, opt.Kface);
        folderResults = strcat(pwd, '/Results/', opt.template,num2str(opt.numVert),'/',opt.relAlgor,'/mat', opt.saveFile, extraName);
        if ~exist(folderResults, 'dir')
            fprintf(['No folder with results:',folderResults,'\n']);
        else
            
            %create folder of data with files
            if strcmp(opt.analysis,'savedata')
                folderEnergy = strcat(pwd, '/Results/', opt.template,num2str(opt.numVert),'/',opt.relAlgor,'/energy', opt.saveFile, extraName);
                [fMassDist, fHinge, fEnergy, fAngles] = makeFileswHeaders(folderEnergy, folderResults);
                
            end
            
            allFiles = dir(folderResults);
            directories = 0;
            succesfullFiles = 0;
            
            for ct = 1:length(allFiles)
                if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:end-4), 'metadata')
                    % skip all directories and metadata file
                    directories = directories+1;
                    continue;
                end

                % parse the file name to get back hinge set
                resfilename = allFiles(ct).name;
                [hingeSet, ~] = getHingeSet(resfilename);
                if strcmp(opt.analysisType,'single')
                    if ~isequal(hingeSet(:,1), opt.angleConstrFinal(end).val(:,1))
                        continue;
%                     elseif ~strcmp(resfilename(1:end-4), '[8 3]_Ang1_18_Angl2_27')
%                         continue;
                    end
                end
                % load results from file
                load(strcat(folderResults,'/', allFiles(ct).name), 'result');
                succesfullFiles = succesfullFiles + 1;
                fprintf('Plot of Hinges number %d/%d\n', succesfullFiles, length(allFiles)-directories);
                
                if strcmp(opt.analysis, 'savedata')
%                     [lowerR, upperR] = getData(extrudedUnitCell, opt, result);
                    Energies = [ones(size(result.E,2),1)*(ct-directories), result.Eedge(2,:)',...
                        result.Ediag(2,:)', result.Eface(2,:)', result.Ehinge(2,:)',...
                        result.EtargetAngle(2,:)', result.exfl(2,:)'];
%                     PosStad = [(ct-directories), lowerR, upperR];
                    Hinges = [num2str(ct-directories),',',mat2str(result.kappa,5),',',...
                        mat2str(result.anglConstr(1,2),5),',', int2str(result.angNum(1)),',',...
                        int2str(result.angNum(2))];
                    AllAngles = zeros(size(extrudedUnitCell.theta));
                    for iter = 1:size(result.deform,2)
                        AllAngles = [AllAngles result.deform(iter).theta];
                    end
                    AllAngles = [ones(size(AllAngles,2),1)*(ct-directories) AllAngles'];
%                     dlmwrite(fMassDist, PosStad, 'delimiter', ',', '-append','precision',7);
                    dlmwrite(fHinge, Hinges, 'delimiter', '', '-append');
                    dlmwrite(fEnergy, Energies, 'delimiter', ',', '-append','precision',7);
                    dlmwrite(fAngles, AllAngles, 'delimiter', ',', '-append','precision',7);
                end
                
                if strcmp(opt.createFig, 'on') || strcmp(opt.analysis, 'plot')
                    
                    nameFolderPlot=[pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.relAlgor,'/images',...
                        opt.saveFile,extraName];
                    nameFilePlot = ['/',resfilename(1:end-4),'_AnglEv.png'];
                    
                    if ~exist(nameFolderPlot, 'dir')
                        mkdir(nameFolderPlot);
                    end
                    allangles = [];
                    for iter = 1:size(result.deform,2)
                        allangles = [allangles result.deform(iter).interV(:).theta];
                    end
                    p = plot(allangles', 'Color', 'k');
                    for i = result.anglConstr(:,1)'
                        set(p(i), 'color', rand(1,3), 'LineWidth', 2, 'DisplayName',num2str(i));
                    end
                    x = 0;
                    for iter = 1:(size(result.deform,2)-1)
                        x = x + size(result.deform(iter).interV,2)+0.5;
                        line([x x],[-1.1*pi 1.1*pi], 'Color', [0 0 0])
                    end
                    legend(p(result.anglConstr(:,1)'))
                    saveas(gcf, [nameFolderPlot, nameFilePlot]);
                    savefig([nameFolderPlot,'/',resfilename(1:end-4),'_AnglEv.fig'])
                    close 'all';                    
                    
                    outputResults(extrudedUnitCell,result,opt,resfilename(1:end-4));
                end
                close all;
            end
        end
        fclose('all');
        
end

function [fileMassDist, fileHinge, fileEnergy, fileAngles] = makeFileswHeaders(folderEnergy, folderResults)

if ~exist(folderEnergy, 'dir')
    mkdir(folderEnergy);
end

fileEnergy = strcat(folderEnergy, '/','EnergyData.csv');
if exist(fileEnergy, 'file')
    delete(fileEnergy) % always start with new file
end
headersEnergy = {'Hinge Number'; 'EdgeEnergy';'DiagonalEnergy'; 'FaceEnergy'; ...
    'HingeEnergy'; 'TargetAngleEnergy'; 'Flags'};
writeHeader(fileEnergy, headersEnergy);

fileHinge = strcat(folderEnergy, '/','Hinges.csv');
if exist(fileHinge, 'file')
    delete(fileHinge) % always start with new file
end
headersHinge = {'HingeNumber'; 'Kappa'; 'TargetAngle'; 'NumAngle'; 'NumKappa'};
writeHeader(fileHinge, headersHinge);

fileMassDist = strcat(folderEnergy, '/','PosStad.csv');
if exist(fileMassDist, 'file')
    delete(fileMassDist) % always start with new file
end
headersMassDist = {'Hinge Number';'LowerRadius';'UpperRadius'};
writeHeader(fileMassDist, headersMassDist);

fileAngles = strcat(folderEnergy, '/','Angles.csv');
if exist(fileAngles, 'file')
    delete(fileAngles) % always start with new file
end
headersAngles = {'HingeNumber'; 'All Angles'};
writeHeader(fileAngles, headersAngles);

fileMetadata = strcat(folderEnergy, '/','metadata.txt');
if exist(fileMetadata, 'file')
    delete(fileMetadata) % always start with new file
end
copyfile([folderResults '/metadata.txt'],fileMetadata);

function writeHeader(file, headers)

fid = fopen(file, 'wt') ;                         % Opens file for writing.
for j = 1:size(headers,1)
    fprintf(fid, '%s', headers{j});
    if j ~= size(headers,1)
        fprintf(fid, ',');
    else
        fprintf(fid, '\n');
    end
end
fclose(fid) ;                                          % Closes file.


function [hingeSet, opening] = getHingeSet(fileName)

parsedName = strsplit(fileName(1:end-4), '_');
hingeSetStr = parsedName{1};
% if length(hingeSetStr)>2
%     if strcmp(hingeSetStr(end-1:end),'op')
%         hingeSetStr = hingeSetStr(1:end-2);
        opening = 1;
%     else
%         opening = 0;
%     end
% else
%     opening = 0;
% end
hingeSetStr = strrep(hingeSetStr, '[', '');
hingeSetStr = strrep(hingeSetStr, ']', '');
hingeSetStr = strsplit(hingeSetStr(1:end), ' ');
hinges = str2double(hingeSetStr);
hingeSet = [hinges' zeros(size(hinges))'];

function [lowerR, upperR] = getData(extrudedUnitCell, opt, result)

currIter = length(result.deform);
startPos = extrudedUnitCell.node + result.deform(currIter).interV(1).V;
endPos = extrudedUnitCell.node + result.deform(currIter).interV(end).V;
upperNodes = cellfun(@(v)v(1),extrudedUnitCell.face([extrudedUnitCell.upperFace]));
lowerNodes = cellfun(@(v)v(1),extrudedUnitCell.face([extrudedUnitCell.lowerFace]));
[~,upperR] = sphereFit(startPos(upperNodes,:));
[~,lowerR] = sphereFit(startPos(lowerNodes,:));


function [EhingeIntSum, suminttheta, sumexttheta] = getIntEnergy(result, opt, extrudedUnitCell)
theta=result.theta;
EhingeInt = zeros(size(extrudedUnitCell.innerHinges,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[result.Ve(1:3:end) result.Ve(2:3:end) result.Ve(3:3:end)];

for i= 1:length(extrudedUnitCell.innerHinges)
    EhingeInt(i)=1/2*opt.Khinge*(theta(extrudedUnitCell.innerHinges(i))-...
        extrudedUnitCell.theta(extrudedUnitCell.innerHinges(i)))^2;
end

EhingeIntSum = sum(EhingeInt);
suminttheta = sum(theta(extrudedUnitCell.innerHinges));
theta(extrudedUnitCell.innerHinges) = [];
sumexttheta = sum(theta);

function [maxStrech, minStrech] = getExtremeStreching(u, opt, extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=(L-extrudedUnitCell.edgeL(i))/extrudedUnitCell.edgeL(i);            
end

maxStrech = max(dEdge);
minStrech = min(dEdge);

function normal = getavNormal(nodes,prevnormal) %%%%This can be used to get the normal of the faces

center = mean(nodes,1);
node1 = 1;
node2 = ceil(size(nodes,1)/2);

a=nodes(node1,:)-center;
b=nodes(node2,:)-center;
alpha=acos(sum(a.*b)/(norm(a)*norm(b)));
if imag(alpha) > 0
    alpha = 0;
end
normal=cross(a,b)/(norm(a)*norm(b)*sin(alpha));
if sum(isinf(normal))
    normal = prevnormal;
end

function [Center,Radius] = sphereFit(X)
% this fits a sphere to a collection of data using a closed form for the
% solution (opposed to using an array the size of the data set). 
% Minimizes Sum((x-xc)^2+(y-yc)^2+(z-zc)^2-r^2)^2
% x,y,z are the data, xc,yc,zc are the sphere's center, and r is the radius
% Assumes that points are not in a singular configuration, real numbers, ...
% if you have coplanar data, use a circle fit with svd for determining the
% plane, recommended Circle Fit (Pratt method), by Nikolai Chernov
% http://www.mathworks.com/matlabcentral/fileexchange/22643
% Input:
% X: n x 3 matrix of cartesian data
% Outputs:
% Center: Center of sphere 
% Radius: Radius of sphere
% Author:
% Alan Jennings, University of Dayton
A=[mean(X(:,1).*(X(:,1)-mean(X(:,1)))), ...
    2*mean(X(:,1).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,1).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    mean(X(:,2).*(X(:,2)-mean(X(:,2)))), ...
    2*mean(X(:,2).*(X(:,3)-mean(X(:,3)))); ...
    0, ...
    0, ...
    mean(X(:,3).*(X(:,3)-mean(X(:,3))))];
A=A+A.';
B=[mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,1)-mean(X(:,1))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,2)-mean(X(:,2))));...
    mean((X(:,1).^2+X(:,2).^2+X(:,3).^2).*(X(:,3)-mean(X(:,3))))];
Center=(A\B).';
Radius=sqrt(mean(sum([X(:,1)-Center(1),X(:,2)-Center(2),X(:,3)-Center(3)].^2,2)));
