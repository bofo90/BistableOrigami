function extrudedUnitCell = createTessellation(unitCell, opt)

extrudedUnitCell = copyUnitCell(unitCell,opt);
extrudedUnitCell = removeDuplicateNodes(extrudedUnitCell);
extrudedUnitCell.edge = unique(sort(extrudedUnitCell.edge,2), 'rows');
extrudedUnitCell = removeDuplicateHinges(extrudedUnitCell);
extrudedUnitCell = makeSquareFaces(extrudedUnitCell);
extrudedUnitCell = addDiagonals(extrudedUnitCell);
extrudedUnitCell = calculateLength(extrudedUnitCell);
extrudedUnitCell.theta = ones(size(extrudedUnitCell.nodeHingeEx,1),1)*opt.restang;


function extrudedUnitCell = copyUnitCell(unitCell, opt)

extrudedUnitCell.node = [];
extrudedUnitCell.center = [];
extrudedUnitCell.edge = [];
extrudedUnitCell.face = {};
extrudedUnitCell.nodeHingeEx = [];

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

function extrudedUnitCell = removeDuplicateHinges(extrudedUnitCell)


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

extrudedUnitCell.nodeHingeEx(loc1,:) = hinges;
extrudedUnitCell.nodeHingeEx(loc2,:) = [];

exUC = extrudedUnitCell;
for i = 1:size(loc1,2)
    extrudedUnitCell.allhinges(exUC.allhinges==loc2(i))= loc1(i);
end
[a,~,c] = unique(extrudedUnitCell.allhinges);
orderedHinges = 1:size(a,2);
extrudedUnitCell.allhinges = orderedHinges(c);

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



