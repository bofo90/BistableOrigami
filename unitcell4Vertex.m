function [unitCell,opt] = unitcell4Vertex(opt)

switch opt.tessellationType
    
    case{'25'}
        if strcmp(opt.vertexType, 'non') 
            if ((opt.angDesign(3)-opt.angDesign(1)) ~= pi || (opt.angDesign(4)-opt.angDesign(2)) ~= pi)
                error('\n----------\nThe angles that you described are incompatible to this tessellation. Please try again or select a compatible 4-Vertex type.\n----------\n',[])
            end
        elseif ~strcmp(opt.vertexType, '2C') && ~strcmp(opt.vertexType, '2CFF')
            error('\n----------\nThe 4-Vertex type is not compatible to this tessellations.Please select either 2C or 2CFF.\n----------\n',[])
        end
        
        opt.angDesign = anglesDesign4Vertex(opt);
        
        unitCell = createOneVertex(opt);
        
        unitCell.latvec(1,:) = unitCell.node(2,:)-unitCell.node(1,:);
        unitCell.latvec(2,:) = unitCell.node(3,:)-unitCell.node(1,:);
            
    otherwise
        error('\n----------\nThe type of 4-Vertex tessellation is non-existent or specified incorrectly\n----------\n',[])

    
    
end

function unitCell = createOneVertex(opt)

radius = opt.Lextrude;

unitCell.centerNode = 1;
unitCell.node = [0,0,0];
for phi = opt.angDesign
    unitCell.node = [unitCell.node; ...
        radius*cos(phi),radius*sin(phi),0];
end

perimNodes = (1:opt.numVert)+1;
unitCell.edge = [ones(1,opt.numVert);perimNodes]';
unitCell.edge = [unitCell.edge; [perimNodes;circshift(perimNodes,-1)]'];
unitCell.face = num2cell([ones(1,opt.numVert);perimNodes;circshift(perimNodes,-1)]',2)';

unitCell.nodeHingeEx = [ones(1,opt.numVert);perimNodes;circshift(perimNodes,-1);circshift(perimNodes,1)]';
