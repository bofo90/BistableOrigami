function [extrudedUnitCell, opt] = obtainOrigami(opt)

extrudedUnitCell.node = [];
extrudedUnitCell.edge = [];
extrudedUnitCell.nodeHingeEx = [];
extrudedUnitCell.face = [];

switch opt.template
    case{'snaping'}
        heigth = 0;
        extrudedUnitCell.node = [0 0 0;1 0 -heigth;0 1 heigth;-1 0 -heigth;0 -1 heigth];
        extrudedUnitCell.edge = [1 2;1 3;1 4;1 5;2 3;3 4;4 5;5 1];
        extrudedUnitCell.nodeHingeEx = [1 2 3 5;1 3 4 2; 1 4 5 3;1 5 2 4];
        extrudedUnitCell.face = {[1 2 3] [1 3 4] [1 4 5] [1 5 2]};
        
    case{'SingleVertex'}
        
        if opt.numVert == 4
            opt.angDesign = anglesDesign4Vertex(opt);
        end
                
        radius = opt.Lextrude;
        extrudedUnitCell.node = [0,0,0];
        for phi = opt.angDesign
            extrudedUnitCell.node = [extrudedUnitCell.node; ...
                radius*cos(phi),radius*sin(phi),0];
        end
        
        perimNodes = (1:opt.numVert)+1;
        extrudedUnitCell.edge = [ones(1,opt.numVert);perimNodes]';
        extrudedUnitCell.edge = [extrudedUnitCell.edge; [perimNodes;circshift(perimNodes,-1)]'];
        extrudedUnitCell.face = num2cell([ones(1,opt.numVert);perimNodes;circshift(perimNodes,-1)]',2)';
        
        extrudedUnitCell.nodeHingeEx = [ones(1,opt.numVert);perimNodes;circshift(perimNodes,-1);circshift(perimNodes,1)]';
        extrudedUnitCell.diagonals = [];
        
        extrudedUnitCell = calculateLength(extrudedUnitCell);
        extrudedUnitCell.theta = ones(size(extrudedUnitCell.nodeHingeEx,1),1)*opt.restang;
        
%         extrudedUnitCell.node = extrudedUnitCell.node+rand(size(extrudedUnitCell.node))*0.01;

    case{'Tessellation'}
        
        if opt.numVert ~= 4
            error('\n----------\nFor the moment we only consider tessellations of a 4-Vertex\n----------\n',[])
        end
        
        [unitCell, opt] = unitcell4Vertex(opt);
        
        extrudedUnitCell = createTessellation(unitCell, opt);
        extrudedUnitCell = addDiagonals(extrudedUnitCell);
        extrudedUnitCell = calculateLength(extrudedUnitCell);
        extrudedUnitCell.theta = ones(size(extrudedUnitCell.nodeHingeEx,1),1)*opt.restang;
        
        
    case{'TriangularTiling'}
        tri_node = [-0.5,0,0;0.5,0,0;0,sqrt(3)/2,0];
        extrudedUnitCell.node = [tri_node];

        or.square(1) = createPolygon([1 2 3]);
        or.line = [];
        or.link = [];
        
        %combine all edges, hinges and faces together
        extrudedUnitCell = combineAll(extrudedUnitCell, or);
        %create diagonals at every face
        extrudedUnitCell = addDiagonals(extrudedUnitCell);
        
    case{'SquareTiling2'}
        latDim = opt.ratio/sqrt(2);
        
        unitCell.nodes = [0 0 latDim;1 0 latDim; 1 1 latDim; 0 1 latDim;...
                          1+latDim 1+latDim latDim;2+latDim 1+latDim latDim; 2+latDim 2+latDim latDim; 1+latDim 2+latDim latDim;...
                          1+latDim 0 latDim;2+latDim 0 latDim; 2+latDim 1 latDim; 1+latDim 1 latDim;...
                          0 1+latDim latDim;1 1+latDim latDim; 1 2+latDim latDim; 0 2+latDim latDim];
        unitCell.edge = [1 2; 2 3; 3 4; 4 1; 5 6; 6 7; 7 8; 8 5; 9 10;...
                         10 11; 11 12; 12 9; 13 14; 14 15; 15 16; 16 13;...
                         4 13; 3 14; 2 9;3 12; 8 15; 5 14;5 12; 6 11];
        unitCell.face = {[1 2 3 4];[5 6 7 8]; [9 10 11 12]; [13 14 15 16];...
                         [4 3 14 13];[5 8 15 14];[3 2 9 12];[6 5 12 11]};
        unitCell.nodeHingeEx = [3 4 2 13;13 14 16 3;2 3 1 12;12 9 11 2;...
                                8 5 7 14;14 15 13 8;5 6 8 11;11 12 10 5];            
        
        extrudedUnitCell = unitCell;
        extrudedUnitCell = addDiagonals(extrudedUnitCell);
        
    case{'SquareTiling'}
        layers = opt.layers;
        dimensions = opt.numUnitCell;
        or = struct();
        latDim = opt.ratio/sqrt(2);
        d = 1+latDim;
        
        %Create a checkerboard matrix
        p = mod(1 : dimensions, 2);
        level_squares = ~bsxfun(@xor, p', p);

        %Create the nodes at the squares
        matrix_squares = zeros(dimensions);
        square_num = 1;
        square_node = [];
        initial_square = [0,0,0;1,0,0;1,1,0;0,1,0];
        for i = 0:dimensions-1
            for j = 0:dimensions-1
                if level_squares(i+1,j+1)
                    shift_square = initial_square+[i*d,j*d,level_squares(i+1,j+1)*(latDim)];
                    for k = 1:layers
                        shift_squareDown = initial_square+[i*d,j*d,(level_squares(i+1,j+1)-2*k)*(latDim)];
                        shift_square = [shift_square; shift_squareDown];
                        matrix_squares(i+1,j+1,k) = square_num;
                        square_num = square_num + 1;
                    end
                    square_num = square_num + 1;
                else 
                    shift_square = [];
                    for k = 1:layers
                        shift_squareDown = initial_square+[i*d,j*d,-2*(k-1)*latDim];
                        shift_square = [shift_square; shift_squareDown];
                        matrix_squares(i+1,j+1,k) = square_num;
                        square_num = square_num + 1;
                    end
                end
                square_node = [square_node; shift_square];
            end
        end
        
        %Create the nodes of the flaps to connect squares at the edges
        flapN =[];
        flapW =[];
        flapS =[];
        flapE =[];
        for i = 1:dimensions
            if level_squares(end,i)
                for k = 1:layers
                    flapN = [flapN; [(i-1)*d,dimensions*d,(level_squares(end,i)-1-2*(k-1))*latDim];...
                             [(i-1)*d+1,dimensions*d,(level_squares(end,i)-1-2*(k-1))*latDim]];
                    flapS = [flapS; [(dimensions-i)*d+1,-latDim,(level_squares(1,dimensions-i+1)-1-2*(k-1))*latDim];...
                             [(dimensions-i)*d,-latDim,(level_squares(1,dimensions-i+1)-1-2*(k-1))*latDim]];
                end
            else
                for k = 1:layers-1
                    flapN = [flapN; [(i-1)*d,dimensions*d,(level_squares(end,i)-1-2*(k-1))*latDim];...
                             [(i-1)*d+1,dimensions*d,(level_squares(end,i)-1-2*(k-1))*latDim]];
                    flapS = [flapS; [(dimensions-i)*d+1,-latDim,(level_squares(1,dimensions-i+1)-1-2*(k-1))*latDim];...
                             [(dimensions-i)*d,-latDim,(level_squares(1,dimensions-i+1)-1-2*(k-1))*latDim]];
                end
            end
            if level_squares(i,1)
                for k = 1:layers
                    flapW = [flapW; [-latDim,(i-1)*d,(level_squares(i,1)-1-2*(k-1))*latDim];...
                             [-latDim,(i-1)*d+1,(level_squares(i,1)-1-2*(k-1))*latDim]];
                    flapE = [flapE; [dimensions*d,(dimensions-i)*d+1,(level_squares(dimensions-i+1,end)-1-2*(k-1))*latDim];...
                             [dimensions*d,(dimensions-i)*d,(level_squares(dimensions-i+1,end)-1-2*(k-1))*latDim]];
                end
            else
                for k = 1:layers-1
                    flapW = [flapW; [-latDim,(i-1)*d,(level_squares(i,1)-1-2*(k-1))*latDim];...
                             [-latDim,(i-1)*d+1,(level_squares(i,1)-1-2*(k-1))*latDim]];
                    flapE = [flapE; [dimensions*d,(dimensions-i)*d+1,(level_squares(dimensions-i+1,end)-1-2*(k-1))*latDim];...
                             [dimensions*d,(dimensions-i)*d,(level_squares(dimensions-i+1,end)-1-2*(k-1))*latDim]];
                end                    
            end
        end
        flap_node = [flapW;flapN;flapE;flapS];
        %add all the nodes together
        extrudedUnitCell.node = [square_node; flap_node];
        
        %create faces at each square and edges at each flap
        nodenum = 1;
        for i = 1:size(square_node,1)/4
            or.square(i) = createPolygon([nodenum nodenum+1 nodenum+2 nodenum+3]);
            nodenum = nodenum+4;
        end
        upperFace = matrix_squares(:,:,1);
        lowerFace = matrix_squares(:,:,end);
        extrudedUnitCell.upperFace = upperFace(level_squares);
        extrudedUnitCell.lowerFace = lowerFace(level_squares);
        
        for i = 1:size(flap_node,1)/2
            or.line(i) = createLine([nodenum nodenum+1]);
            nodenum = nodenum+2;
        end
%         or.link = [];

        %connect squares at the edges with their correspondonig flap
        matrix_squares = permute(matrix_squares,[2 1 3]);
        lineW = 1;
        lineN = 1;
        linespEdge = size(flap_node,1)/2/4;
        for i = 1:dimensions
            if level_squares(i,1)
                for k = 1:layers
                    [or.link(lineW), or.line(lineW)] = ...
                    combineS_L(or.square(matrix_squares(i,1,k)),or.square(matrix_squares(i,1,k)+1),or.line(lineW), 4);
                    [or.link(lineW+2*linespEdge), or.line(lineW+2*linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(dimensions-i+1,end,k)),or.square(matrix_squares(dimensions-i+1,end,k)+1),or.line(lineW+2*linespEdge), 2);
                    lineW = lineW + 1; 
                end
            else
                for k = 1:layers-1
                    [or.link(lineW), or.line(lineW)] = ...
                    combineS_L(or.square(matrix_squares(i,1,k)),or.square(matrix_squares(i,1,k)+1),or.line(lineW), 4);
                    [or.link(lineW+2*linespEdge), or.line(lineW+2*linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(dimensions-i+1,end,k)),or.square(matrix_squares(dimensions-i+1,end,k)+1),or.line(lineW+2*linespEdge), 2);
                    lineW = lineW + 1; 
                end
            end
            if level_squares(end,i)
                for k = 1:layers
                    [or.link(lineN+linespEdge), or.line(lineN+linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(end,i,k)),or.square(matrix_squares(end,i,k)+1),or.line(lineN+linespEdge), 3);
                    [or.link(lineN+3*linespEdge), or.line(lineN+3*linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(1,dimensions-i+1,k)),or.square(matrix_squares(1,dimensions-i+1,k)+1),or.line(lineN+3*linespEdge), 1);
                    lineN = lineN + 1;    
                end
            else
                for k = 1:layers-1
                    [or.link(lineN+linespEdge), or.line(lineN+linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(end,i,k)),or.square(matrix_squares(end,i,k)+1),or.line(lineN+linespEdge), 3);
                    [or.link(lineN+3*linespEdge), or.line(lineN+3*linespEdge)] = ...
                    combineS_L(or.square(matrix_squares(1,dimensions-i+1,k)),or.square(matrix_squares(1,dimensions-i+1,k)+1),or.line(lineN+3*linespEdge), 1);
                    lineN = lineN + 1;    
                end
            end
        end

        %connect squares between each other (only the higher ones to all its neighbors
        link_size = size(or.link,2);
        init_SShinge = size(or.line,2)+2*link_size;
        for k = 1:layers
            for i = 1:dimensions
                for j = 1:dimensions
                    if level_squares(i,j)
                        if (i+1) <= dimensions
                            link_size = link_size+1;
                            or.link(link_size) = ...
                            combineS_S(or.square(matrix_squares(i,j,k)),or.square(matrix_squares(i,j,k)+1),...
                                       or.square(matrix_squares(i+1,j,k)), 3, 1);
                        end
                        if (j+1) <= dimensions
                            link_size = link_size+1;
                            or.link(link_size) = ...
                            combineS_S(or.square(matrix_squares(i,j,k)),or.square(matrix_squares(i,j,k)+1),...
                                       or.square(matrix_squares(i,j+1,k)), 2, 4);
                        end
                        if (i-1) > 0
                            link_size = link_size+1;
                            or.link(link_size) = ...
                            combineS_S(or.square(matrix_squares(i,j,k)),or.square(matrix_squares(i,j,k)+1),...
                                       or.square(matrix_squares(i-1,j,k)), 1, 3);
                        end
                        if (j-1) > 0
                            link_size = link_size+1;
                            or.link(link_size) = ...
                            combineS_S(or.square(matrix_squares(i,j,k)),or.square(matrix_squares(i,j,k)+1),...
                                       or.square(matrix_squares(i,j-1,k)), 4, 2);
                        end
                    end
                end
            end
        end
       
        %combine all edges, hinges and faces together
        extrudedUnitCell = combineAll(extrudedUnitCell, or);
        end_SShinge = size(extrudedUnitCell.nodeHingeEx,1);
        extrudedUnitCell.beta = [init_SShinge+2:4:(end_SShinge-init_SShinge)/opt.layers+init_SShinge];
        extrudedUnitCell.alpha = [end_SShinge:-4:end_SShinge-(end_SShinge-init_SShinge)/opt.layers+1];
        %create diagonals at every face
        extrudedUnitCell = addDiagonals(extrudedUnitCell);
        %Determine initial edge length
        extrudedUnitCell = calculateLength(extrudedUnitCell);
        %Determine initial angles
        extrudedUnitCell.theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
        for i=1:size(extrudedUnitCell.nodeHingeEx,1)
            extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
            index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
            index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
            index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
            [~,extrudedUnitCell.theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
        end
    otherwise
        error('Origami not defined\n')
end

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
   
function face = addingIndexFace(face,index)
for i = 1:size(face,2)
    face{i} = face{i} + index;
end

function square = createPolygon(v)

sides = size(v,2);
square.edge = [];
for i = 1:sides
    square.edge =  [square.edge; v(1) v(2)];
    v = circshift(v,-1);
end
square.edge = [square.edge; v(end) v(1)];
square.face = {[v]};

for j = 1:sides
    for i = 1:2
        square.side(j,i).edge = [v(1) 0;v(2) 0];
        
        square.side(j,i).face = {[v(1) v(2) 0 0]};

        square.side(j,i).nodeHingeEx = [v(1) v(2) v(end) 0];
    end
    v = circshift(v,-1);
end

function line = createLine(v)

line.edge = v;
line.nodeHingeEx = [v(1) v(2) 0 0];
for i = 1:2
    line.C(i).edge = [0 v(2);0 v(1)];
    line.C(i).face = {[0 0 v(1) v(2)]};
end

function [link, l] = combineS_L(s1, s2 ,l, dir)

sides = 1:size(s1.face{1},2);
sides = circshift(sides, 1-dir);

link.edge = [s1.side(dir,1).edge + l.C(2).edge;...
             s2.side(dir,2).edge + l.C(1).edge];
link.face = {[s1.side(dir,1).face{1} + l.C(2).face{1}] ...
             [s2.side(dir,2).face{1} + l.C(1).face{1}]};
link.nodeHingeEx = [s1.side(dir,1).nodeHingeEx + [0 0 0 l.edge(1)];...
                    s2.side(dir,2).nodeHingeEx + [0 0 0 l.edge(1)]];
l.nodeHingeEx = l.nodeHingeEx + [0 0 s2.face{1}(sides(2)) s1.face{1}(sides(1))];

function link = combineS_S(s1u, s1d, s2, dir1, dir2)

link.edge = [s1u.side(dir1,1).edge + [[0;0] flipud(s2.side(dir2,2).edge(:,1))];...
             s1d.side(dir1,2).edge + [[0;0] flipud(s2.side(dir2,1).edge(:,1))]];
link.face = {[s1u.side(dir1,1).face{1} + [0 0 s2.side(dir2,2).face{1}(1:2)]]...
             [s1d.side(dir1,2).face{1} + [0 0 s2.side(dir2,1).face{1}(1:2)]]};
         
link.nodeHingeEx = [s1u.side(dir1,1).nodeHingeEx + [0 0 0 s2.face{1}(dir2)];...
                    s2.side(dir2,2).nodeHingeEx + [0 0 0 s1u.face{1}(dir1)];...
                    s1d.side(dir1,2).nodeHingeEx + [0 0 0 s2.face{1}(dir2)];...
                    s2.side(dir2,1).nodeHingeEx + [0 0 0 s1d.face{1}(dir1)]];

function eUC = combineAll(eUC, or)

eUC.edge = [];
eUC.face = {};
eUC.nodeHingeEx = [];

for i = 1:size(or.square,2)
    eUC.edge = [eUC.edge; or.square(i).edge];
    eUC.face = [eUC.face  or.square(i).face];
end

for i = 1:size(or.line,2)
    eUC.edge = [eUC.edge; or.line(i).edge];
    eUC.nodeHingeEx = [eUC.nodeHingeEx; or.line(i).nodeHingeEx];
end

for i = 1:size(or.link,2)
    eUC.edge = [eUC.edge; or.link(i).edge];
    eUC.face = [eUC.face  or.link(i).face];
    eUC.nodeHingeEx = [eUC.nodeHingeEx; or.link(i).nodeHingeEx];
end

