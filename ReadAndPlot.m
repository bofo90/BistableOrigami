function ReadAndPlot(opt, flag, des)

switch flag
    case{'allDes'}
        for i = des  

            if strcmp(opt.template,'Tessellation')
                fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',i);
            elseif strcmp(opt.template,'SingleVertex')
                fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',i);
            end    

            allFiles = dir(fileContainer);
            for ct = 3:length(allFiles)
                opt.file = strcat(fileContainer,'/',allFiles(ct).name);
                SingleFileRead(opt)
            end
        end
    case{'oneDes'}
        if strcmp(opt.template,'Tessellation')
            fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType);
        elseif strcmp(opt.template,'SingleVertex')
            fileContainer = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType);
        end 
        opt.file = strcat(fileContainer,opt.file);
        SingleFileRead(opt)
end

function SingleFileRead(opt)
file = opt.file;

opt.vertexType = 'non';
if strcmp(opt.template,'Tessellation')
    [opt.xrep, opt.yrep] = getTessell(opt.file);
elseif strcmp(opt.template,'SingleVertex')
    opt.angDesign = getAngles(opt.file);
end
[extrudedUnitCell,opt]=obtainOrigami(opt);

allFiles = dir(file);
for j = 3:length(allFiles)
    if strcmp(allFiles(j).name,'Images')
        continue;
    end
    
    opt.restang = getParam(allFiles(j).name);
    extrudedUnitCell.theta = ones(size(extrudedUnitCell.theta,1),1)*opt.restang;

    allFiles2 = dir(strcat(file, sprintf('/RestAng_%.3f', opt.restang)));
    for k = 3:length(allFiles2)
        opt.Khinge = getParam(allFiles2(k).name);
        opt.file = strcat(file,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));
        temp = char(opt.file);
        fprintf(strcat('Saving data: ', temp(end-65:end), '\n'))
        obtainData(extrudedUnitCell, opt);
    end
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
    
function kappa = getParam(fileName)
parsedName = strsplit(fileName(1:end), '_');
kappa = str2double(parsedName{2});

function obtainData(extrudedUnitCell,opt)  
folderResults = strcat(opt.file,'/mat');
folderEnergy = strcat(opt.file,'/energy');
if ~exist(folderResults, 'dir')
    error('\n----------\nNo data directory found\n----------\n',[])
end 

allFiles = dir(folderResults);

%create folder of data with files
[fMassDist, fHinge, fEnergy, fAngles] = makeFileswHeaders(folderEnergy, folderResults);
Energies = zeros(length(allFiles)-4,7);
PosStad = zeros(length(allFiles)-4,size(extrudedUnitCell.face,2)+1);
Hinges = zeros(length(allFiles)-4,size(extrudedUnitCell.allnodes,2)/4+1);
AllAngles = zeros(length(allFiles)-4,size(extrudedUnitCell.allhinges,2)+1);
dirs = 0;
for ct = 1:length(allFiles)
    if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:end-4), 'metadata')
        % skip all directories and metadata file
        dirs = dirs+1;
        continue;
    end

    % load results from file
    lofile = load(strcat(folderResults,'/', allFiles(ct).name));
%             fprintf('Saving data %d\n', ct);

    sim = getSimulation(allFiles(ct).name);

    curv = getCurvature(extrudedUnitCell, opt, lofile.result);
    areas = getAreas(extrudedUnitCell, lofile.result);
    
    if ~isreal(curv)
        curv = 0;
        lofile.result.exfl(2,end)= -10;
    end
    
    if ~isreal(areas)
        areas = zeros(size(areas));
        lofile.result.exfl(2,end)= -10;
    end
    
    Energies(ct-dirs,:) = [sim, lofile.result.Eedge(2,end),...
        lofile.result.Ediag(2,end), lofile.result.Eface(2,end), lofile.result.Ehinge(2,end),...
        lofile.result.EtargetAngle(2,end), lofile.result.exfl(2,end)];
    PosStad(ct-dirs,:) = [sim, areas'];
    Hinges(ct-dirs,:) = [sim, curv'];%, designAng(1), designAng(2), designAng(3)]];
%             AllAnglesTemp = zeros([size(extrudedUnitCell.theta,1),size(lofile.result.deform,2)]);
%             for iter = 1:size(lofile.result.deform,2)
%                 AllAnglesTemp(:,iter) = lofile.result.deform(iter).theta;
%             end
    AllAngles(ct-dirs,:) = [sim lofile.result.deform(end).theta(extrudedUnitCell.allhinges)'];

    close all;
    clear lofile;
end

dlmwrite(fMassDist, PosStad, 'delimiter', ',', '-append','precision',7);
dlmwrite(fHinge, Hinges, 'delimiter', ',', '-append','precision',7);
dlmwrite(fEnergy, Energies, 'delimiter', ',', '-append','precision',7);
dlmwrite(fAngles, AllAngles, 'delimiter', ',', '-append','precision',7);
fclose('all');
        


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
headersHinge = {'HingeNumber'; 'Curvature'};
writeHeader(fileHinge, headersHinge);

fileMassDist = strcat(folderEnergy, '/','PosStad.csv');
if exist(fileMassDist, 'file')
    delete(fileMassDist) % always start with new file
end
headersMassDist = {'Hinge Number';'Face Areas'};
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
copyfile(strcat(folderResults,'/metadata.txt'),fileMetadata);

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


function hingeSet = getSimulation(fileName)

parsedName = strsplit(fileName(1:end-4), '_');
Angl1 = parsedName{2};
hingeSet = str2double(Angl1);

function [curvature, areaFaces] = getCurvature(extrudedUnitCell, opt, result)
endPos = extrudedUnitCell.node + result.deform(end).interV(end).V;
vertexnodes = reshape(extrudedUnitCell.allnodes,4,[])';
vertexnodes2 = circshift(vertexnodes,1,2);
curvature = zeros(size(vertexnodes,1),1);
for j = 1:size(vertexnodes,1)
    anglesFaces = zeros(opt.numVert,1);
    areaFaces = zeros(opt.numVert,1);
    for i = 1: opt.numVert
        v1 = endPos(vertexnodes(j,i),:)-endPos(extrudedUnitCell.center(j),:);
        v2 = endPos(vertexnodes2(j,i),:)-endPos(extrudedUnitCell.center(j),:);
        l1 = sqrt(v1*v1');
        l2 = sqrt(v2*v2');
        anglesFaces(i) = acos(v1*v2'/l1/l2);
        areaFaces(i) = l1*l2*sin(anglesFaces(i))/2;
    end
    curvature(j) = 3*(2*pi-sum(anglesFaces))/sum(areaFaces);
end

function areas = getAreas(extrudedUnitCell, result)
endPos = extrudedUnitCell.node + result.deform(end).interV(end).V;
areas = zeros(size(extrudedUnitCell.face'));
for i = 1: size(extrudedUnitCell.face,2)
    v1 = endPos(extrudedUnitCell.face{i}(2),:)-endPos(extrudedUnitCell.face{i}(1),:);
    v2 = endPos(extrudedUnitCell.face{i}(3),:)-endPos(extrudedUnitCell.face{i}(1),:);
    l1 = sqrt(v1*v1');
    l2 = sqrt(v2*v2');
    angle = acos(v1*v2'/l1/l2);
    areas(i) = l1*l2*sin(angle)/2;
end

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
