function ReadAndPlot(unitCell, extrudedUnitCell, opt)

switch opt.plot
    case 'info'
        result = [];
        outputResults(unitCell,extrudedUnitCell,result,opt);
    case {'result', 'savedata'}
        if strcmp(opt.plot,'savedata')
            folderEnergy = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/energy');
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
        
        folderResults = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat');
        if ~exist(folderResults, 'dir')
            fprintf('No folder with results\n');
        else
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
                extrudedUnitCell.angleConstr = [hingeSet(:), -pi*0.985 * ones(length(hingeSet), 1)];
                load(strcat(folderResults,'/', fileName));
                succesfullFiles = succesfullFiles + 1;
                if strcmp(opt.plot, 'savedata')
                    [CM, Radios, Stdev, EhingeInt] = getRadiosStdev(extrudedUnitCell, opt, result);
                    Energies = [ones(length(result.E),1)*(ct-directories), result.Eedge,...
                        result.Eface, result.Ehinge, result.EtargetAngle, EhingeInt];
                    PosStad = cat(3,ones(length(result.E),2,1)*(ct-directories),CM, Radios, Stdev);
                    dlmwrite(fileMassDist, PosStad, 'delimiter', ',', '-append');
                    dlmwrite(fileHinge,[num2str(ct-directories),',',parsedName{2}], 'delimiter', '', '-append');
                    dlmwrite(fileEnergy, Energies, 'delimiter', ',', '-append');
                end
                if strcmp(opt.createFig, 'on')
                    outputResults(unitCell,extrudedUnitCell,result,opt);
                end
                if strcmp(opt.showFig, 'off')
                    close all;
                end
            end
            fprintf('Succesfull result Files read: %d\n', succesfullFiles);
        end
        fclose('all');
        
end


function [hingeSet] = getHingeSet(Name)

hingeSetStr = Name;
hingeSetStr = strrep(hingeSetStr, '[', '');
hingeSetStr = strrep(hingeSetStr, ']', '');
hingeSetStr = strsplit(hingeSetStr(1:end), ' ');
hingeSet = str2double(hingeSetStr)';


function [CM, Radios, Stdev, EhingeInt] = getRadiosStdev(extrudedUnitCell, opt, result)
CM = zeros(length(result.deform(1).interV),length(result.deform),3);
Radios = zeros(length(result.deform(1).interV),length(result.deform));
Stdev = zeros(length(result.deform(1).interV),length(result.deform));
EhingeInt = zeros(length(result.deform(1).interV),length(result.deform));

for iter = 1:length(result.deform)
    for inter = 1:length(result.deform(iter).interV)
        newPos = extrudedUnitCell.node + result.deform(iter).interV(inter).V;
        CM(inter, iter,:) = mean(newPos);
        allRad = sqrt(sum(abs(newPos-CM(inter,iter)).^2,2));
        Radios(inter, iter) = mean(allRad);
        Stdev(inter, iter) = std(allRad);
        EhingeInt(inter, iter) = getIntEnergy(result.deform(iter).interV(inter).Ve, opt, extrudedUnitCell);
    end
end

function [EhingeIntSum] = getIntEnergy(u, opt, extrudedUnitCell)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
innerHinges = [1 2 3 7 8 9 13 14 18 19 23 24 28 32 36 37 41 48]';
EhingeInt = zeros(size(innerHinges,1),1);
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    [~,theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end
for i= 1:length(innerHinges)
    EhingeInt(i)=1/2*opt.Khinge*(theta(innerHinges(i))-extrudedUnitCell.theta(innerHinges(i)))^2;
end
EhingeIntSum = sum(EhingeInt);

