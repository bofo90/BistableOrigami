function ReadAndPlot(unitCell, extrudedUnitCell, opt)

switch opt.plot
    case 'info'
        result = [];
        outputResults(unitCell,extrudedUnitCell,result,opt);
    case {'result', 'savedata'}
        
        folderResults = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat/', sprintf('kangle%2.4f_khinge%2.4f', opt.KtargetAngle, opt.Khinge));
        if ~exist(folderResults, 'dir')
            fprintf('No folder with results\n');
        else
            if strcmp(opt.plot,'savedata')
                folderEnergy = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/energy/', sprintf('kangle%2.4f_khinge%2.4f', opt.KtargetAngle, opt.Khinge));
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
            end
        
            allFiles = dir(folderResults);
            directories = 0;
            succesfullFiles = 0;
            for ct = 1:length(allFiles)
                if allFiles(ct).isdir
            %         disp('skip all directories...')
                    directories = directories+1;
                    continue;
                end

                % parse the file name to get back hinge set
                fileName = allFiles(ct).name;
                parsedName = strsplit(fileName(1:end-4), '_');
                hingeSet = getHingeSet(parsedName{2});
                if ~isequal(hingeSet, opt.angleConstrFinal(1).val(:,1)) && strcmp(opt.readAngFile,'off')
                    continue;
                end
               
                extrudedUnitCell.angleConstr = [hingeSet(:), -pi*0.985 * ones(length(hingeSet), 1)];
                load(strcat(folderResults,'/', fileName));
                succesfullFiles = succesfullFiles + 1;
                fprintf('Plot of Hinges number %d/%d\n', succesfullFiles, length(allFiles)-directories);
                if strcmp(opt.plot, 'savedata')
                    [CM, Radios, Stdev, EhingeInt, maxStrech, minStrech] = getRadiosStdev(extrudedUnitCell, opt, result);
                    Energies = [ones(length(result.E),1)*(ct-directories), result.Eedge,...
                        result.Eface, result.Ehinge, result.EtargetAngle, EhingeInt(1:end-1:end,:), result.exfl];
                    PosStad = cat(3,ones(length(result.E),2,1)*(ct-directories),CM(1:end-1:end,:,:),...
                        Radios(1:end-1:end,:), Stdev(1:end-1:end,:),...
                        maxStrech(1:end-1:end,:), minStrech(1:end-1:end,:));
                    Hinges = [num2str(ct-directories),',',parsedName{2}];
                    dlmwrite(fileMassDist, PosStad, 'delimiter', ',', '-append');
                    dlmwrite(fileHinge, Hinges, 'delimiter', '', '-append');
                    dlmwrite(fileEnergy, Energies, 'delimiter', ',', '-append');
                end
                if strcmp(opt.createFig, 'on')
                    outputResults(unitCell,extrudedUnitCell,result,opt);
                end
                if strcmp(opt.showFig, 'off')
                    close all;
                end
            end
        end
        fclose('all');
        
end


function [hingeSet] = getHingeSet(Name)

hingeSetStr = Name;
hingeSetStr = strrep(hingeSetStr, '[', '');
hingeSetStr = strrep(hingeSetStr, ']', '');
hingeSetStr = strsplit(hingeSetStr(1:end), ' ');
hingeSet = str2double(hingeSetStr)';


function [CM, Radios, Stdev, EhingeInt, maxStrech, minStrech] = getRadiosStdev(extrudedUnitCell, opt, result)
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
