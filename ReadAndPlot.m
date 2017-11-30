function ReadAndPlot(unitCell, extrudedUnitCell, opt)

switch opt.plot
    case 'info'
        result = [];
        outputResults(unitCell,extrudedUnitCell,result,opt);
    case {'result', 'savedata', 'plot'}
        folderResults = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat', opt.saveFile);
        if ~exist(folderResults, 'dir')
            fprintf('No folder with results:' + folderResults + '\n');
        else
            if strcmp(opt.plot,'savedata')
                folderEnergy = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/energy', opt.saveFile);
                [fMassDist, fHinge, fEnergy] = reviewDataFiles(folderEnergy, folderResults);
                
            end
            
            allFiles = dir(folderResults);
            directories = 0;
            succesfullFiles = 0;
            
            for ct = 1:length(allFiles)
                if allFiles(ct).isdir || strcmp(allFiles(ct).name(1:8), 'metadata')
                    % skip all directories and metadata file
                    directories = directories+1;
                    continue;
                end

                % parse the file name to get back hinge set
                hingeSet = getHingeSet(allFiles(ct).name);
                if ~isequal(hingeSet, opt.angleConstrFinal(1).val(:,1)) && strcmp(opt.readAngFile,'off')
                    continue;
                end
                extrudedUnitCell.angleConstr = [hingeSet(:), -pi*0.985 * ones(length(hingeSet), 1)];
                % load results from file
                load(strcat(folderResults,'/', allFiles(ct).name));
                succesfullFiles = succesfullFiles + 1;
                fprintf('Plot of Hinges number %d/%d\n', succesfullFiles, length(allFiles)-directories);
                
                if strcmp(opt.plot, 'savedata')
                    [CM, Radios, Stdev, EhingeInt, maxStrech, minStrech] = getData(extrudedUnitCell, opt, result);
                    Energies = [ones(length(result.E),1)*(ct-directories), result.Eedge,...
                        result.Eface, result.Ehinge, result.EtargetAngle, EhingeInt(1:end-1:end,:), result.exfl];
                    PosStad = cat(3,ones(length(result.E),2,1)*(ct-directories),CM(1:end-1:end,:,:),...
                        Radios(1:end-1:end,:), Stdev(1:end-1:end,:),...
                        maxStrech(1:end-1:end,:), minStrech(1:end-1:end,:));
                    Hinges = [num2str(ct-directories),',',mat2str(hingeSet)];
                    dlmwrite(fMassDist, PosStad, 'delimiter', ',', '-append');
                    dlmwrite(fHinge, Hinges, 'delimiter', '', '-append');
                    dlmwrite(fEnergy, Energies, 'delimiter', ',', '-append');
                end
                if strcmp(opt.createFig, 'on')
                    if strcmp(opt.onlyUnitCell, 'on')
                       [extrudedUnitCell, result] = extrudeInnerPolyhedra(extrudedUnitCell,result,1); 
                    end
                    outputResults(unitCell,extrudedUnitCell,result,opt);
                end
                if strcmp(opt.showFig, 'off')
                    close all;
                end
            end
        end
        fclose('all');
        
end

function [fileMassDist, fileHinge, fileEnergy] = reviewDataFiles(folderEnergy, folderResults)

if ~exist(folderEnergy, 'dir')
    mkdir(folderEnergy);
end

fileEnergy = strcat(folderEnergy, '/','EnergyData.csv');
if exist(fileEnergy, 'file')
    delete(fileEnergy) % always start with new file
end
fileHinge = strcat(folderEnergy, '/','Hinges.csv');
if exist(fileHinge, 'file')
    delete(fileHinge) % always start with new file
end
fileMassDist = strcat(folderEnergy, '/','PosStad.csv');
if exist(fileMassDist, 'file')
    delete(fileMassDist) % always start with new file
end
fileMetadata = strcat(folderEnergy, '/','metadata.txt');
if exist(fileMetadata, 'file')
    delete(fileMetadata) % always start with new file
end
copyfile([folderResults '/metadata.txt'],fileMetadata);


function [hingeSet] = getHingeSet(fileName)

parsedName = strsplit(fileName(1:end-4), '_');
hingeSetStr = parsedName{2};
hingeSetStr = strrep(hingeSetStr, '[', '');
hingeSetStr = strrep(hingeSetStr, ']', '');
hingeSetStr = strsplit(hingeSetStr(1:end), ' ');
hingeSet = str2double(hingeSetStr)';


function [CM, Radios, Stdev, EhingeInt, maxStrech, minStrech] = getData(extrudedUnitCell, opt, result)
CM = zeros(length(result.deform(1).interV),length(result.deform),3);
Radios = zeros(length(result.deform(1).interV),length(result.deform));
Stdev = zeros(length(result.deform(1).interV),length(result.deform));
EhingeInt = zeros(length(result.deform(1).interV),length(result.deform));
maxStrech = zeros(length(result.deform(1).interV),length(result.deform));
minStrech = zeros(length(result.deform(1).interV),length(result.deform));

for iter = 1:length(result.deform)
    for inter = 1:length(result.deform(iter).interV)
        newPos = extrudedUnitCell.node + result.deform(iter).interV(inter).V;
        CM(inter, iter,:) = mean(newPos);
        allRad = sqrt(sum(abs(newPos-CM(inter,iter)).^2,2));
        Radios(inter, iter) = mean(allRad);
        Stdev(inter, iter) = std(allRad);
        EhingeInt(inter, iter) = getIntEnergy(result.deform(iter).interV(inter).Ve, opt, extrudedUnitCell);
        [maxStrech(inter, iter), minStrech(inter, iter)] = ...
            getExtremeStreching(result.deform(iter).interV(inter).Ve, opt, extrudedUnitCell);
    end
end

function [EhingeIntSum] = getIntEnergy(u, opt, extrudedUnitCell)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
EhingeInt = zeros(size(extrudedUnitCell.innerHinges,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    [~,theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end
for i= 1:length(extrudedUnitCell.innerHinges)
    EhingeInt(i)=1/2*opt.Khinge*(theta(extrudedUnitCell.innerHinges(i))-...
        extrudedUnitCell.theta(extrudedUnitCell.innerHinges(i)))^2;
end
EhingeIntSum = sum(EhingeInt);

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

function [extrudedUnitCell, result] = extrudeInnerPolyhedra(extrudedUnitCell,result, extrutionLength)
%%% Assuming only one unit cell
Faces = extrudedUnitCell.face;
Nodes = size(extrudedUnitCell.node,1);
rep = 0;
for i=1:length(Faces)      
    nNo=size(extrudedUnitCell.node,1);
    nodeNum=Faces{i};
    nodeNumNew=nNo+1:nNo+length(nodeNum);
    
%     normal = getNormal(extrudedUnitCell.node, nodeNum(1), nodeNum(2), nodeNum(3),[0 0 0]);
    normal = getavNormal(extrudedUnitCell.node(Faces{i},:), [0 0 0]);
    extrudedUnitCell.node(nodeNumNew,:)=extrudedUnitCell.node(nodeNum,:)+extrutionLength*ones(length(nodeNum),1)*normal;
    
    for j = 1:result.numMode
        for k = 1:length(result.deform(j).interV)
            movedNodes = extrudedUnitCell.node(1:Nodes,:)+result.deform(j).interV(k).V(1:Nodes,:);
%             normal = getNormal(movedNodes, nodeNum(1), nodeNum(2), nodeNum(3),normal);
            normal = getavNormal(movedNodes(Faces{i},:), normal);
            newNodes = movedNodes(nodeNum,:)+extrutionLength*ones(length(nodeNum),1)*normal;
            result.deform(j).interV(k).V(nodeNumNew,:) = newNodes - extrudedUnitCell.node(nodeNumNew,:);
            result.deform(j).interV(k).Ve=result.deform(j).interV(k).V(:,end);
        end
        movedNodes = extrudedUnitCell.node(1:Nodes,:)+result.deform(j).V(1:Nodes,:);
%         normal = getNormal(movedNodes, nodeNum(1), nodeNum(2), nodeNum(3),normal);
        normal = getavNormal(movedNodes(Faces{i},:), normal);
        newNodes = movedNodes(nodeNum,:)+extrutionLength*ones(length(nodeNum),1)*normal;
        result.deform(j).V(nodeNumNew,:)= newNodes - extrudedUnitCell.node(nodeNumNew,:);
        result.deform(j).Ve=result.deform(j).V(:,end);
    end
    index=[1:length(nodeNum) 1];
    for j=1:length(nodeNum)   
        rep=rep+1;
        extrudedUnitCell.face{rep}=[nodeNum(index(j)) nodeNum(index(j+1)) nodeNumNew(index(j+1)) nodeNumNew(index(j))];
    end
end

function normal = getNormal(nodes, node1, node2, node3,prevnormal)

a=nodes(node2,:)-nodes(node1,:);
b=nodes(node3,:)-nodes(node1,:);
alpha=acos(sum(a.*b)/(norm(a)*norm(b)));
if imag(alpha) > 0
    alpha = 0;
end
normal=cross(a,b)/(norm(a)*norm(b)*sin(alpha));
if sum(isinf(normal))
    normal = prevnormal;
end

function normal = getavNormal(nodes,prevnormal)

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

