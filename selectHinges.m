function selectHinges(unitCell, extrudedUnitCell, opt)

%This is the function for selecting hinges
if strcmp(opt.plot, 'selectHinges')
    fprintf('Creating graph...\n');
    [G] = buildGraph(unitCell, extrudedUnitCell);
    fprintf('Getting flavours...\n');
    [flavourDict, flavourTypes, flavourNum] = getFlavours(G);
    fprintf('Calculating distances...\n');
    dis = getDistance(G, flavourDict);

    % first generate the list of all hinge sets
    fprintf('Writting hinges...\n');
    getAllHinges(G, dis, flavourTypes, flavourNum, unitCell, extrudedUnitCell, opt);
    fprintf('Done with selecting hinges...\n');
end

% hi all
function [G] = buildGraph(unitCell, extrudedUnitCell)

numEdges = size(unitCell.Polyhedron.edge, 1);
G_adjacency = zeros(numEdges);

hingeNames = cell(numEdges, 1); % corresponds to names in extrudedUnitCell
hingeTypes = cell(numEdges, 1);
isClosed = cell(numEdges, 1);
hingeNodes = cell(numEdges, 1);

center = mean(unitCell.Polyhedron.node, 1);

for ii = 1:numEdges
    n1 = unitCell.Polyhedron.edge(ii,1);
    n2 = unitCell.Polyhedron.edge(ii,2); % endnodes of the edge
    hingeNodes{ii} = sort([n1 n2]);
    % getting the edge numbering from nodeHingeEx (consistent with how
    % they're numbered in outputResults
    [~,idx1]=ismember([n1 n2],extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
    [~,idx2]=ismember([n2 n1],extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
    hingeNames{ii} = num2str(idx1+idx2);
    
    % getting the types of faces associated with current edge
    faceTypes = [0 0];
    for ff = 1:length(unitCell.Polyhedron.face)
        n_face = unitCell.Polyhedron.face{ff};
        if sum(find(n1==n_face)) && sum(find(n2==n_face))
            faceTypes(find(0==faceTypes, 1)) = length(unitCell.Polyhedron.face{ff});
        end
    end
	hingeTypes{ii} = [num2str(min(faceTypes)),'-',num2str(max(faceTypes))];
    isClosed{ii} = '0';
    
    % here starts stuff for adjacency matrix
    for jj = 1:numEdges
        if ii == jj
            G_adjacency(ii,jj) = 0;
        else
            n_edge2 = unitCell.Polyhedron.edge(jj, :);
            
            % see if these two hinges share a common node
            comm_node = [n_edge2(n1==n_edge2), n_edge2(n2==n_edge2)];
            if ~isempty(comm_node)
                % determine the direction between these two hinges
                edgeV1 = - getHingeVector([n1,n2], comm_node, unitCell);
                edgeV2 = getHingeVector(n_edge2, comm_node, unitCell);
                normV = unitCell.Polyhedron.node(comm_node,:) - center;
                direction = dot(normV, cross(edgeV1,edgeV2));
                if direction < 0
                    G_adjacency(ii,jj) = 1;
                elseif direction > 0
                    G_adjacency(jj,ii) = 1;
                end
            end
        end
    end
end

G_adjacency = sparse(G_adjacency);

% finally, make a graph of the hinges
G = digraph(G_adjacency);

% assign names and other properties to the nodes of the graph (edges on the polyherdon)
G.Nodes.Name = hingeNames;
% G.Nodes.Properties.RowNames = hingeNames;
G.Nodes.IsClosed = isClosed;
G.Nodes.Type = hingeTypes;
G.Nodes.HingeNodes = hingeNodes;


function [vec] = getHingeVector(nodes, comm_node, unitCell)
% returns the hinge in the form of a 3D vector, and always use *comm_node*
% as the starting point of the vector
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% nodes     - a 2*1 array, containing the two nodes of a hinge
% comm_node - the node in common between this hinge and the other one 
%             under consideration
% unitCell  - a unit cell
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% vec - a vector of this hinge, pointing from *comm_node* to the other node
%       in *nodes*
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Mar 20, 2017


comm_idx = find(comm_node==nodes);
if comm_idx == 1
    % do nothing
elseif comm_idx == 2
    nodes = flip(nodes);
else
    error('Oops, you''ve got the wrong nodes!')
end

coord1 = unitCell.Polyhedron.node(nodes(1), :);
coord2 = unitCell.Polyhedron.node(nodes(2), :);
vec = coord2 - coord1;

function [flavourDict, flavourTypes, flavourNum] = getFlavours(G)
% Returns the info about flavours of a graph G
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G -  a directed graph
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% flavourTypes - a cell array containing the types of flavours
% flavourNum   - number of different flavours
% flavourDict  - a container map that maps different flavourTypes to
%                different numeric values. E.g., in truncated tetrahedron,
%                flavourDict('3-6') = '1', flavourDict('6-6') = '2'
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun

flavourTypes = unique(G.Nodes.Type);
flavourNum = length(flavourTypes);

flavourDict = containers.Map;
for ii = 1:flavourNum
    flavourDict(flavourTypes{ii}) = num2str(ii);
end

function dis = getDistance(G, flavourDict)
% Returns the ``distance'' matrix with flavours
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - an digraph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% dis - a distance matrix
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun

% initialising...
H = height(G.Nodes);
dis = zeros(H);

%Loop through all the nodes
for ii = 1:H
    [TR,~] = shortestpathtree(G, 'all', ii, 'OutputForm','cell');
    %Loop through all the nodes to get the distance in string
    for jj = 1:H
        currentPath = TR{jj};
        currentPathDisStr = '';
        %loop in each distance to create the string
        for lCt = 1:length(currentPath)
            currentPathDisStr = strcat(currentPathDisStr, ...
                flavourDict(G.Nodes{currentPath(lCt),'Type'}{1}));
        end
        dis(ii,jj) = str2double(currentPathDisStr);
    end
end

%Make it symetric
for ii = 1:H
    for jj = ii:H
        dis(ii,jj) = min(dis(ii,jj), dis(jj,ii));
        dis(jj,ii) = dis(ii,jj);
    end
end

% % loop through everything, because we need to
% for ii = 1:H
%     for jj = ii:H
%         % get all possible paths between two nodes *ii* and *jj* both ways
%         pathsiijj = pathbetweennodes(G.adjacency, ii, jj);
%       
%         % get the ``distance'' from node *ii* to node *jj*
%         % by looping through all minimum length paths
%         pathDis = inf; % initialise with infinity
%         for mCt = 1:length(pathsiijj)
%             currentPath = pathsiijj{mCt};
%             currentPathDisStr = '';
%             
%             % get the string of ``distance''
%             for lCt = 1:length(currentPath)
%             	currentPathDisStr = strcat(currentPathDisStr, ...
%                     flavourDict(G.Nodes{currentPath(lCt),'Type'}{1}));
%             end
%             
%             % convert to number and take the smallest value of the
%             % ``distance'' and the flipped ``distance''
%             currentPathDis = min(str2num(currentPathDisStr),...
%                 str2num(flip(currentPathDisStr)));
%             if currentPathDis < pathDis
%                 pathDis = currentPathDis;
%             end
%         end
% 
%         % update distance matrix
%         dis(ii,jj) = pathDis;
%         dis(jj,ii) = pathDis;
%     end
% end


function pth = pathbetweennodes(adj, src, snk)
%PATHBETWEENNODES Return all paths between two nodes of a graph
%
% pth = pathbetweennodes(adj, src, snk)
% pth = pathbetweennodes(adj, src, snk, vflag)
%
%
% This function returns all simple paths (i.e. no cycles) between two nodes
% in a graph.  Not sure this is the most efficient algorithm, but it seems
% to work quickly for small graphs, and isn't too terrible for graphs with
% ~50 nodes.
%
% Input variables:
%
%   adj:    adjacency matrix
%
%   src:    index of starting node
%
%   snk:    index of target node
%
%   vflag:  logical scalar for verbose mode.  If true, prints paths to
%           screen as it traverses them (can be useful for larger,
%           time-consuming graphs). [false]
%
% Output variables:
%
%   pth:    cell array, with each cell holding the indices of a unique path
%           of nodes from src to snk.

% input validation added by yun, May 02, 2017
if src==snk
    pth = {[src]};
    return;
end
pth = cell(0);
next = cell(size(adj,1),1);
pth{1} = 0;
for in = 1:size(adj,1)
    next{in} = find(adj(in,:));
    pth{1} = [pth{1} in];
end 
%Recursive search of shortest distance between src and snk
[pth,~,~,~] = nextNode(pth, src, snk, next);
%Recursive search of shortest distance between snk and src
[pth,~,~,~] = nextNode(pth, snk, src, next);


function [paths, stack, snk, next] = nextNode(paths, stack, snk, next)

%Check if the distance is longer than the previous found one
if length(stack) > length(paths{end})
    return;
end
%Check if we dont have loops in the stack
if length(unique(stack)) ~= length(stack)
    return;
end
%Check if you already end in the finnal edge
if stack(end) == snk
    %if the distance is shorter, erease the rest, if not just stack them
    if length(stack) < length(paths{end}) 
        paths = {stack};
    else
        paths = [paths; {stack}];
    end
    return;
end
%go through the next edges and call the function again
originalstack = stack;
for j = 1:length(next{originalstack(end)})
    newstack = [originalstack next{originalstack(end)}(j)];
    [paths, stack, snk, next] = nextNode(paths, newstack, snk, next);
end



function getAllHinges(G, dis, flavourTypes, flavourNum, unitCell, extrudedUnitCell, opt)
% get the complete list of hinge sets, and dump them in a csv file
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - a directed graph created from the options
% flavourTypes
% flavourNum
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% A .csv file of the complete list of hinge sets, where each row is one
% such hinge set
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 05, 2017
% yun

% create folder if it doesn't exist
folderName = strcat(pwd, '\Results\hingeList_reduced\');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
fileName = strcat(folderName, opt.template, '.csv');
if exist(fileName, 'file')
    delete(fileName) % always start with new file
end

% name of all hinges converted to an array of numbers
allHinges = zeros(1, height(G.Nodes));
for ct = 1:height(G.Nodes)
    allHinges(ct) = str2num(G.Nodes.Name{ct});
end

% starts generating all hignes
hingeSetsPrev(1).all = [];
for N = 1:height(G.Nodes)-1    %(0.5 * height(G.Nodes))
    fprintf('Calculating distance for %d nodes\n', N);
    if N <= 0.5 * height(G.Nodes)
        hingeSets = chooseHinges(G, dis, flavourTypes, flavourNum, N, hingeSetsPrev(N).all);

        for ct = 1:size(hingeSets, 1)
            % write selections to file
            dlmwrite(fileName, hingeSets(ct,:), 'delimiter', ',', '-append')

            % plot and save figures when specified
            if strcmp(opt.createFig, 'on')
                f = plotSelectedHinges(hingeSets(ct,:), unitCell, extrudedUnitCell, true);
                % set(f, 'visible', 'off');
                if strcmp(opt.saveFig, 'on')
                    folderName = strcat(pwd, '/hingeList_reduced/', opt.template, '/', ...
                                 num2str(N), '/');
                    if ~exist(folderName, 'dir')
                     mkdir(folderName)
                    end
                    saveas(f, strcat(folderName, mat2str(hingeSets(ct,:)), '.png'))
                    close(f)
                end
            end
        end
        % getting ready for the next loop
        hingeSetsPrev(N+1).all = hingeSets;
    end
    
    if N >= 0.5 * height(G.Nodes)
        Nconj = height(G.Nodes)-N;
        cHingeSets = zeros(size(hingeSetsPrev(Nconj+1).all, 1), N);
        for ct = 1:size(cHingeSets, 1)
            % write complement selections to file
            cHingeSets(ct,:) = setdiff(allHinges, hingeSetsPrev(Nconj+1).all(ct,:));
            dlmwrite(fileName, cHingeSets(ct,:), 'delimiter', ',', '-append')

            % plot and save figures when specified
            if strcmp(opt.createFig, 'on')
                f = plotSelectedHinges(cHingeSets(ct,:), unitCell, extrudedUnitCell, true);
                % set(f, 'visible', 'off');
                if strcmp(opt.saveFig, 'on')
                    folderName = strcat(pwd, '/hingeList_reduced/', opt.template, '/', ...
                                num2str(height(G.Nodes)-N), '/');
                    if ~exist(folderName, 'dir')
                        mkdir(folderName)
                    end
                    saveas(f, strcat(folderName,mat2str(cHingeSets(ct,:)), '.png'))
                    close(f)
                end
            end
        end
    end
end

function [hinges] = chooseHinges(G, dis, flavourTypes, flavourNum, N, hingeSetsPrev)
% select *N* number of hinges from graph *G* to actuate, the flavoured
% version
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% G - graph object
% dis - matrix containing the minimum distance between any two nodes
% N - number of hinges to be actuated
% hingeSetsPrev - the complete hinge sets where each set contains (N-1) hinges
% opt - options
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hinges - an x*N array, where each row is a set of hinges to be actuated
%          at the same time
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun

% input validation
if N < 0 || N > size(G.Nodes,1)/2
    hinges = [];
    warning('You are selecting the wrong number of hinges to actuate!')
elseif N == 0
    hinges = [];
elseif N == 1
    % actuate different types of hinges
    hinges = zeros(flavourNum, 1);
    for ct = 1:flavourNum
        idx = find(1==strcmp(G.Nodes.Type, flavourTypes(ct)), 1);
        hinges(ct) = str2num(G.Nodes{idx,'Name'}{1});
    end
    
else% N >= 2
    % initialising...
    hinges = [];
    currentIdx = 1;
%     cache.distance = []; % the matrices representing each selection
    cache.distanceEig = []; % their corresponding eigenvalues    
    % loop through each hinge set
    for ct = 1:size(hingeSetsPrev, 1)
        fprintf('possible hinge slection %d/%d\n', ct, size(hingeSetsPrev, 1));
        hingeSet = hingeSetsPrev(ct, :);
        hingeSetIdx = getHingeIdx(hingeSet, G);
        % first, split the possible candidates into cases of different
        % flavours
        for flavourCt = 1:flavourNum
            candidateIdxs = find(1==strcmp(G.Nodes.Type,flavourTypes{flavourCt}));
            % then start comparing matrices
            for nodeCt = 1:length(candidateIdxs)
                newNodeIdx = candidateIdxs(nodeCt);
                if newNodeIdx <= hingeSetIdx(end)
                    % avoid double counting
                    continue;
                end
                newNode = str2double(G.Nodes{newNodeIdx, 'Name'}{1});
                newSet = [hingeSet, newNode];
                newSetIdx = [hingeSetIdx, newNodeIdx];
                newDisMat = dis(newSetIdx, newSetIdx);
                newDisEig = round(1000*sort(eig(newDisMat))') / 1000; % round off 
                
                % if newDisEig is not in cache, add current solution
                if isempty(cache.distanceEig) || ...
                        ~ismember(newDisEig, cache.distanceEig, 'rows')
                    hinges(currentIdx, :) = newSet;
%                     cache.distance(:,:,currentIdx) = newDisMat;
                    cache.distanceEig(currentIdx, :) = newDisEig;
                    currentIdx = currentIdx + 1;
                end        
            end
        end
    end    
end



function hingeSetIdx = getHingeIdx(hingeSet, G)
% gets the index of hinges in graph, as numbered in *outputResults()*
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingeSet - an array of integers, which represent the names of hinges
% G -  a directed graph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% hingeSetIdx - an array of integers of the same size as hingeSet, which
%               contains the indeces of the hinges in the graph G
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Apr 10, 2017
% yun
hingeSetIdx = zeros(size(hingeSet));
for ct = 1:length(hingeSetIdx)
    hingeSetIdx(ct) = find(strcmp(G.Nodes.Name, num2str(hingeSet(ct))));
end

function f = plotSelectedHinges(hingeList, unitCell, extrudedUnitCell, ...
             orderedNames)
% function plotSelectedHinges(geometry, hingeList)
% highlights the set of selected hinges in a given polyhedron, for better
% visualisation
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% hingelist - a list containing the index of the selected hinges
% unitCell
% extrudedUnitCell
% orderedNames - boolean indicating whether the hinge names are printed as
%                the original order as in 
%                extrudedUnitCell.nodeHingeEx --> orderedNames = false
%                or in ordered version --> orderedNames = true
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% f - a graph object
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on Jun 07, 2017
% yun


if ~exist('orderedNames', 'var')
    orderedNames = true;
end

vertices = unitCell.Polyhedron.node;

% plot transparent polyhedra with patch()
faces = unitCell.Polyhedron.face;
figure; hold on;
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]); view([20 10]); axis('equal')
axis('off')
xlabel('x'); ylabel('y'); zlabel('z')
if 1==length(hingeList)
    hingeStr = num2str(hingeList);
else
    hingeStr = mat2str(hingeList);
end

title(strcat('Hinges selected: ',hingeStr(2:(end-1))))
for ct = 1:length(faces)
    patch('Vertices',vertices, 'Faces',faces{ct},...
          'FaceColor',[.73 .86 .95], 'FaceAlpha',0.3, ...
          'EdgeColor',[.2 .2 .2], 'lineWidth',0.5);
end
light('Position',[1 -1 1],'Style','local')


% find corresponding nodes on the hingeList and highlight those edges
for ed = 1:length(hingeList)
    endNodes = extrudedUnitCell.nodeHingeEx(hingeList(ed), 1:2);
    X = vertices(endNodes, 1);
    Y = vertices(endNodes, 2);
    Z = vertices(endNodes, 3);
    line(X, Y, Z, 'lineStyle',':', 'color', [.01 .23 .45], 'LineWidth', 8)
end


% add lables for the edges
for ct = 1:size(unitCell.Polyhedron.edge, 1)
    endNodes = unitCell.Polyhedron.edge(ct, :);
    
    if orderedNames
        hingeName = num2str(ct);
    else
        [~, hinge] = ismember(endNodes,...
                     extrudedUnitCell.nodeHingeEx(:,1:2),'rows');
        hingeName = num2str(hinge);
    end
    
    textPos = 0.5 * (vertices(endNodes(1), :) + vertices(endNodes(2), :));
    text(textPos(1)+.01, textPos(2)+.01, textPos(3)-.04, hingeName, ...
        'FontSize', 26, 'FontWeight', 'bold', 'color', [.2 .2 .2]);
end

f = gcf;