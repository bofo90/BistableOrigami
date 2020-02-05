function extrudedUnitCell = createTessellation(unitCell, opt)

extrudedUnitCell = copyUnitCell(unitCell,opt);
extrudedUnitCell = removeDuplicateNodes(extrudedUnitCell);
extrudedUnitCell = removeDuplicateEdges(extrudedUnitCell,opt);
extrudedUnitCell = removeDuplicateHinges(extrudedUnitCell,opt);
extrudedUnitCell = makeSquareFaces(extrudedUnitCell);
extrudedUnitCell = addDiagonals(extrudedUnitCell);
extrudedUnitCell = calculateLength(extrudedUnitCell);
extrudedUnitCell.theta = ones(size(extrudedUnitCell.nodeHingeEx,1),1)*opt.restang;
if strcmp(opt.periodic,'on')
    extrudedUnitCell = getPeriodicity(extrudedUnitCell,unitCell,opt);
end

function extrudedUnitCell = copyUnitCell(unitCell, opt)

extrudedUnitCell.node = [];
extrudedUnitCell.center = [];
extrudedUnitCell.edge = [];
extrudedUnitCell.face = {};
extrudedUnitCell.nodeHingeEx = [];
% if strcmp(opt.periodic, 'on')
%     extrudedUnitCell.perxedge = zeros(opt.yrep,2);
%     extrudedUnitCell.peryedge = zeros(opt.xrep,2);
%     extrudedUnitCell.perxhinge = zeros(opt.yrep,2);
%     extrudedUnitCell.peryhinge = zeros(opt.xrep,2);
%     
% end
shiftnode = 0;

for i = 0:opt.xrep-1
    for j = 0:opt.yrep-1
        nodePos = unitCell.node + unitCell.latvec(1,:)*i + unitCell.latvec(2,:)*j;
        extrudedUnitCell.node = [extrudedUnitCell.node; nodePos];
        extrudedUnitCell.edge = [extrudedUnitCell.edge; unitCell.edge+shiftnode];
        extrudedUnitCell.nodeHingeEx = [extrudedUnitCell.nodeHingeEx; unitCell.nodeHingeEx+shiftnode];
        extrudedUnitCell.face = [extrudedUnitCell.face addingIndexFace(unitCell.face,shiftnode)];
        extrudedUnitCell.center = [extrudedUnitCell.center unitCell.centerNode+shiftnode];
        shiftnode = size(extrudedUnitCell.node,1);
        
%         if (i == 0) && strcmp(opt.periodic, 'on')
%             extrudedUnitCell.perxedge(j+1,1)= size(extrudedUnitCell.edge,1)-5; 
%             extrudedUnitCell.perxhinge(j+1,1)= size(extrudedUnitCell.nodeHingeEx,1)-1;
%         end
%         if (i == opt.xrep-1) && strcmp(opt.periodic, 'on')
%             extrudedUnitCell.perxedge(j+1,2)= size(extrudedUnitCell.edge,1)-7; 
%             extrudedUnitCell.perxhinge(j+1,2)= size(extrudedUnitCell.nodeHingeEx,1)-3;
%         end
%         if (j == 0) && strcmp(opt.periodic, 'on')
%             extrudedUnitCell.peryedge(i+1,1)= size(extrudedUnitCell.edge,1)-4; 
%             extrudedUnitCell.peryhinge(i+1,1)= size(extrudedUnitCell.nodeHingeEx,1);
%         end
%         if (j == opt.yrep-1) && strcmp(opt.periodic, 'on')
%             extrudedUnitCell.peryedge(i+1,2)= size(extrudedUnitCell.edge,1)-6; 
%             extrudedUnitCell.peryhinge(i+1,2)= size(extrudedUnitCell.nodeHingeEx,1)-2;
%         end
    end
end

extrudedUnitCell.allnodes = 1:size(extrudedUnitCell.node,1);
extrudedUnitCell.allnodes(extrudedUnitCell.center)=[];
extrudedUnitCell.allhinges = 1:size(extrudedUnitCell.nodeHingeEx,1);

function face = addingIndexFace(face,index)
for i = 1:size(face,2)
    face{i} = face{i} + index;
end

function exUC = removeDuplicateNodes(extrudedUnitCell)

exUC = extrudedUnitCell;

[~,b,c] = unique(round(extrudedUnitCell.node,8), 'rows', 'sorted');
exUC.node = extrudedUnitCell.node(b,:);
newnodes = 1:size(b);
newnodes = newnodes(c);

for i = 1:size(extrudedUnitCell.node,1)
    exUC.edge(extrudedUnitCell.edge == i) = newnodes(i);
    exUC.nodeHingeEx(extrudedUnitCell.nodeHingeEx == i) = newnodes(i);
    exUC.center(extrudedUnitCell.center == i) = newnodes(i);
    for f = 1:size(extrudedUnitCell.face,2)
        exUC.face{f}(extrudedUnitCell.face{f}==i) = newnodes(i);
    end
    exUC.allnodes(extrudedUnitCell.allnodes ==i) = newnodes(i);
end

function exUC = removeDuplicateEdges(extrudedUnitCell,opt)

exUC = extrudedUnitCell;

[exUC.edge, ~,newedges] = unique(sort(extrudedUnitCell.edge,2), 'rows');

% if strcmp(opt.periodic,'on')
%     for i = 1:size(extrudedUnitCell.edge,1)
%         exUC.perxedge(extrudedUnitCell.perxedge == i) = newedges(i);
%         exUC.peryedge(extrudedUnitCell.peryedge == i) = newedges(i);
%     end
% end

function extrudedUnitCell = removeDuplicateHinges(extrudedUnitCell,opt)

[a,b,~] = unique((extrudedUnitCell.nodeHingeEx(:,1)==extrudedUnitCell.center) + ...
    (extrudedUnitCell.nodeHingeEx(:,2)==extrudedUnitCell.center),'rows');

hinges = [];
loc1 = [];
loc2 = [];

for i = 1:size(b,1)
    if sum(a(i,:)) == 2
        loc1 = [loc1 b(i)];
        [~,y] = ismember(extrudedUnitCell.nodeHingeEx(b(i),[2,1]),extrudedUnitCell.nodeHingeEx(:,[1,2]), 'rows');
        loc2 = [loc2 y];
        hinges = [hinges; extrudedUnitCell.nodeHingeEx(loc1(end),[1:3]) extrudedUnitCell.nodeHingeEx(loc2(end),3)];
    end
end

numhinges = 1:size(extrudedUnitCell.nodeHingeEx,1);
extrudedUnitCell.nodeHingeEx(loc1,:) = hinges;
extrudedUnitCell.nodeHingeEx(loc2,:) = [];

numhinges(loc2) = loc1;
[a,~,c]=unique(numhinges);
orderedHinges = 1:size(a,2);
orderedHinges = orderedHinges(c);

exUC = extrudedUnitCell;
for i = 1:size(numhinges,2)
    extrudedUnitCell.allhinges(exUC.allhinges==i)= orderedHinges(i);
end

% if strcmp(opt.periodic,'on')
%     for i = 1:size(numhinges,2)
%         extrudedUnitCell.perxhinge(exUC.perxhinge==i)= orderedHinges(i);
%         extrudedUnitCell.peryhinge(exUC.peryhinge==i)= orderedHinges(i);
%     end
% end

function extrudedUnitCell = makeSquareFaces(extrudedUnitCell)

faces = reshape(cell2mat(extrudedUnitCell.face),[3,size(extrudedUnitCell.face,2)])';

squareface = [];
loc1 = [];
loc2 = [];
deledges = [];

for i = 1:size(faces,1)
    [~,y] = ismember(faces(i,[3,2]),faces(i:end,[2,3]),'rows');
    if y ~= 0
        loc1 = [loc1 i];
        loc2 = [loc2 y+i-1];
        squareface = [squareface; faces(i,[1,2]) faces(y+i-1,[1,2])];
        [~,x] = ismember(sort(faces(i,[2,3])),extrudedUnitCell.edge,'rows');
        deledges = [deledges x];
        
    end
end

redfacesloc = 1:size(faces,1);
redfacesloc([loc1,loc2]) = [];
addedges = [];

for i = redfacesloc
    [~,y] = ismember(faces(i,[3,1]), faces(:,[1,2]),'rows');
    if y ~= 0
        loc1 = [loc1 i];
        loc2 = [loc2 y];
        squareface = [squareface; faces(y,[1,2]) faces(i,2) faces(y,3)];
        edges = [faces(i,[2,3]); faces(y,[2,3])];
        addedges = [addedges; faces(i,2) faces(y,3)];
        [~,x] = ismember(sort(edges,2),extrudedUnitCell.edge,'rows');
        deledges = [deledges x'];        
    end
end

[~,y] = unique(sort(squareface,2),'rows');
squareface = squareface(y,:);

extrudedUnitCell.face([loc1,loc2]) = [];
extrudedUnitCell.edge(deledges,:) = [];
extrudedUnitCell.edge = [extrudedUnitCell.edge; addedges];
extrudedUnitCell.face = [extrudedUnitCell.face num2cell(squareface,2)'];

function extrudedUnitCell = calculateLength(extrudedUnitCell)

extrudedUnitCell.edgeL = [];

for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    extrudedUnitCell.edgeL(i)=sqrt(dx*dx');
end

function extrudedUnitCell = addDiagonals(extrudedUnitCell)

extrudedUnitCell.diagonals = [];

for i = 1:length(extrudedUnitCell.face)
    s = length(extrudedUnitCell.face{i});
    if s>3
        for j = 1:round(s/2)
            extrudedUnitCell.edge = [extrudedUnitCell.edge; ...
                extrudedUnitCell.face{i}(j) extrudedUnitCell.face{i}(j+round(length(extrudedUnitCell.face{i})/2))];
            extrudedUnitCell.diagonals = [extrudedUnitCell.diagonals, size(extrudedUnitCell.edge,1)];
        end
    end
end

extrudedUnitCell.edge = sort(extrudedUnitCell.edge,2);

function extrudedUnitCell = getPeriodicity(extrudedUnitCell,unitCell,opt)

extrudedUnitCell.latVec = unitCell.latvec.*[opt.xrep;opt.yrep];
extrudedUnitCell.node = [extrudedUnitCell.node; rand(2,3)*0.01];
extrudedUnitCell.ref = [size(extrudedUnitCell.node,1)-1, size(extrudedUnitCell.node,1)];

extrudedUnitCell.repnumnodes = [];
extrudedUnitCell.repnodes = [];
extrudedUnitCell.repedges = [];
extrudedUnitCell.rephinges = [];

%%%These lattice vectors might change depending on the tessellation
tessvec = [extrudedUnitCell.latVec; sum(extrudedUnitCell.latVec); diff(extrudedUnitCell.latVec)];

for vec = 1:size(tessvec,1)
    repnodes = [];
    for i=1:size(extrudedUnitCell.node,1)-2
        for j=1:size(extrudedUnitCell.node,1)-2
            if norm(extrudedUnitCell.node(i,:)+tessvec(vec,:)-extrudedUnitCell.node(j,:)) < 10^(-8)
                repnodes = [repnodes; i j];
            end
        end
    end
    extrudedUnitCell.repnumnodes = [extrudedUnitCell.repnumnodes; size(repnodes,1)];
    extrudedUnitCell.repnodes = [extrudedUnitCell.repnodes; repnodes];
    extrudedUnitCell.repedges = [extrudedUnitCell.repedges; getLib(extrudedUnitCell.edge, repnodes)];
    extrudedUnitCell.rephinges = [extrudedUnitCell.rephinges; getLib(sort(extrudedUnitCell.nodeHingeEx(:,1:2),2), repnodes)];
end


function lib2 = getLib(hinges, lib)
x = zeros(size(hinges,1),1);
y = zeros(size(hinges,1),1);
for i = 1:size(lib,1)
    x = x + sum(hinges==lib(i,1),2);
    y = y + sum(hinges==lib(i,2),2);
end
numedges = 1:size(hinges,1);
st = numedges(x == 2)';
en = numedges(y == 2)';
a = hinges(st,:);
b = hinges(en,:);
for i = 1:size(lib,1)
    b(b==lib(i,2))=lib(i,1);
end
dela = [];
rel = [];
for i = 1:size(a,1)
    [x,y] = ismember(a(i,:),b,'rows');
    if ~x
        dela = [dela i];
    else
        rel = [rel, y];
    end
end
st(dela,:) = [];
en = en(rel,:);
lib2 = [st,en];




