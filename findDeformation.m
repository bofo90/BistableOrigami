function findDeformation(unitCell,extrudedUnitCell,opt)

%Show details geometries (if requested)
if strcmp(opt.plot,'result')
    if strcmp(opt.readAngFile,'off')
        nonlinearFolding(unitCell,extrudedUnitCell,opt);
    else
        opt.angleConstrFinal = [];
        fileHinges = strcat(pwd, '\Results\hingeList_reduced\', opt.template, '.csv');
        hingeList = dlmread(fileHinges);
        for i = 1:size(hingeList, 1)
            row = hingeList(i, :);
            hinges = row(0~=row);
            opt.angleConstrFinal(1).val = [hinges(:), -pi*0.985 * ones(length(hinges), 1)];
            fprintf('Hinge selection number %d/%d. ', i, size(hingeList, 1))
            nonlinearFolding(unitCell,extrudedUnitCell,opt);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NON-LINEAR ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonlinearFolding(unitCell,extrudedUnitCell,opt)

%INITIALIZE LINEAR CONSTRAINTS
[Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt);

%FIND DEFORMATION
theta0=extrudedUnitCell.theta;
max_iter = 100000;
extrudedUnitCell.angleConstr=[];

folderName = strcat(pwd, '/Results/', opt.template,'/',opt.relAlgor,'/mat');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

for iter=1:length(opt.angleConstrFinal)
    clearvars V;
    u0=zeros(3*size(extrudedUnitCell.node,1),1);
    V(:,1)=u0;
    [E(1,1),~,Eedge(1,1),Eface(1,1),Ehinge(1,1),EtargetAngle(1,1), ~]=Energy(u0,extrudedUnitCell,opt);
    exfl(1,1) = 1;
    fprintf(['Angle contrain:', mat2str(opt.angleConstrFinal(iter).val(:,1)') ,'\n']);
    extrudedUnitCell.angleConstr=opt.angleConstrFinal(iter).val;
    for inter=1:opt.interval
        fprintf('load: %d, step: %1.2f',2*iter-1,inter/opt.interval);
        t1 = toc;
        extrudedUnitCell.angleConstr(:,2)=theta0(extrudedUnitCell.angleConstr(:,1))+inter*(-theta0(extrudedUnitCell.angleConstr(:,1))+opt.angleConstrFinal(iter).val(:,2))/opt.interval;
        %Determine new equilibrium
        [V(:,inter+1),~,exfl(inter+1,1)]=fmincon(@(u) Energy(u,extrudedUnitCell,opt),u0,[],[],Aeq,Beq,[],[],@(u) nonlinearConstr(u,extrudedUnitCell,opt),opt.options);
        u0 = V(:,inter+1);
        %Determine energy of that equilibrium
        [E(inter+1,1),~,Eedge(inter+1,1),Eface(inter+1,1),Ehinge(inter+1,1),EtargetAngle(inter+1,1), ~]=Energy(u0,extrudedUnitCell,opt);
        t2 = toc;
        fprintf(', time: %1.2f, exitflag: %d\n',t2-t1,exfl(inter+1,1));
    end
    
    result.deform(1).V=[V(1:3:end,end) V(2:3:end,end) V(3:3:end,end)];
    result.deform(1).Ve=V(:,end);
    for j=1:size(V,2)
        result.deform(1).interV(j).V=[V(1:3:end,j) V(2:3:end,j) V(3:3:end,j)];
        result.deform(1).interV(j).Ve=V(:,j);
    end
    
    clearvars V;
    V(:,1)=u0;

    switch opt.relAlgor
        case 'gradDesc'
            fprintf('load: %d, step: 1.00',2*iter)
            t1 = toc;
            extrudedUnitCell.angleConstr=[];
            [E(1,2),~,Eedge(1,2),Eface(1,2),Ehinge(1,2),EtargetAngle(1,2), ~]=Energy(u0,extrudedUnitCell,opt);
            exfl(1,2) = 1;
            [V, exfl(2,2)] = gradDescent(extrudedUnitCell, opt, u0, max_iter);
            u0=V(:,end);
            [E(2,2),~,Eedge(2,2),Eface(2,2),Ehinge(2,2),EtargetAngle(2,2), ~]=Energy(u0,extrudedUnitCell,opt);

            t2 = toc;
            fprintf(', time: %1.2f, exitflag: %d\n',t2-t1,exflexfl(2,2))
        case {'interior-point', 'sqp'}
            [E(1,2),~,Eedge(1,2),Eface(1,2),Ehinge(1,2),EtargetAngle(1,2), ~]=Energy(u0,extrudedUnitCell,opt);
            exfl(1,2) = 1;
            prevKtargetAngle = opt.KtargetAngle;
            opt.options.Algorithm = opt.relAlgor;
            for inter = 1:opt.relInterval
                fprintf('load: %d, step: %1.2f',2*iter,inter/opt.relInterval)
                t1 = toc;
                opt.KtargetAngle = prevKtargetAngle - (prevKtargetAngle/opt.relInterval)*inter;
                [V(:,inter+1),~,exfl(inter+1,2)]=fmincon(@(u) Energy(u,extrudedUnitCell,opt),u0,[],[],Aeq,Beq,[],[],@(u) nonlinearConstr(u,extrudedUnitCell,opt),opt.options);
                u0 = V(:,inter+1);
                [E(inter+1,2),~,Eedge(inter+1,2),Eface(inter+1,2),Ehinge(inter+1,2),EtargetAngle(inter+1,2), ~]=Energy(u0,extrudedUnitCell,opt);
                t2 = toc;
                fprintf(', time: %1.2f, exitflag: %d\n',t2-t1,exfl(inter+1,2))
            end
            opt.options.Algorithm = opt.folAlgor;
            opt.KtargetAngle = prevKtargetAngle;
            extrudedUnitCell.angleConstr=[];
    end
    
    result.deform(2).V=[V(1:3:end,end) V(2:3:end,end) V(3:3:end,end)];
    result.deform(2).Ve=V(:,end);
    for j=1:size(V,2)
        result.deform(2).interV(j).V=[V(1:3:end,j) V(2:3:end,j) V(3:3:end,j)];
        result.deform(2).interV(j).Ve=V(:,j);
    end
    
    result.E=E;
    result.Eedge=Eedge;
    result.Eface=Eface;
    result.Ehinge=Ehinge;
    result.EtargetAngle=EtargetAngle;    
    result.exfl = exfl;
    result.numMode=length(result.deform);
    
    fileName = strcat(folderName,'/',opt.template,'_',...
        mat2str(opt.angleConstrFinal(iter).val(:,1)'),'_','.mat');
    save(fileName, 'result');
    
    clearvars result E Eedge Eface Ehinge EtargetAngle exfl;
    fclose('all');
    
    %update angle starting point
%     [~,~,~,~,~,~,theta0]=Energy(u0,extrudedUnitCell,opt);
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENERGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, dE,Eedge,Eface,Ehinge,EtargetAngle,theta]=Energy(u,extrudedUnitCell,opt)

E=0; dE=zeros(3*size(extrudedUnitCell.node,1),1);
Eedge=0;
Eface=0;
Ehinge=0;
EtargetAngle=0;

%APPLY DEFORMATION NODES
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%ENERGY ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'off')
    [dEdge, Jedge]=getEdge(extrudedUnitCell);
    Eedge=1/2*opt.Kedge*sum(dEdge.^2);
    dE=dE+opt.Kedge*Jedge'*dEdge;
end

%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'off')
    [dFace, Jface]=getFace(extrudedUnitCell);
    Eface=opt.Kface/2*sum(dFace.^2);
    dE=dE+opt.Kface*Jface'*dFace;
end

%ENERGY ASSOCIATED TO HINGE BENDING
[theta, Jhinge]=getHinge(extrudedUnitCell);
Ehinge=1/2*opt.Khinge*sum((theta-extrudedUnitCell.theta).^2);
dE=dE+opt.Khinge*(Jhinge'*(theta-extrudedUnitCell.theta));

%ENERGY ASSOCIATED TO TARGET HINGE ANGLES
if size(extrudedUnitCell.angleConstr,1)==0
    dtheta=[];
else
    dtheta=theta(extrudedUnitCell.angleConstr(:,1))-extrudedUnitCell.angleConstr(:,2);
    dE=dE+opt.KtargetAngle*Jhinge(extrudedUnitCell.angleConstr(:,1),:)'*dtheta;
end
EtargetAngle=1/2*opt.KtargetAngle*sum(dtheta.^2);

%TOTAL ENERGY
E=Eedge+Eface+Ehinge+EtargetAngle;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NONLINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,Ceq,DC,DCeq]=nonlinearConstr(u,extrudedUnitCell,opt)

C1=[]; C2=[];
DC1=[]; DC2=[];
Ceq1=[]; Ceq2=[]; Ceq3=[];
DCeq1=[]; DCeq2=[]; DCeq3=[];

%APPLY DEFORMATION NODES
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%%%%%%%%%%%%%%%%%%%%%%
%INEQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM AND MINIMUM ANGLES
[angles, Dangles]=getHinge(extrudedUnitCell);
C1 = [-angles-pi*opt.constAnglePerc; angles-pi*opt.constAnglePerc];
DC1 = [-Dangles; Dangles];
%MAXIMUM AND MINIMUM EDGE STRECHING
if strcmp(opt.constrEdge,'off') && ~isnan(opt.maxStretch)
    [normStrech, DnormStrech]=getEdgeNorm(extrudedUnitCell);
    C2 = [-normStrech-opt.maxStretch; normStrech-opt.maxStretch];
    DC2 = [-DnormStrech; DnormStrech];
end
C = [C1; C2];
DC = [DC1; DC2]';

%%%%%%%%%%%%%%%%%%%%%%
%EQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%CONSTRAINT ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'on')
    [Ceq1, DCeq1]=getEdge(extrudedUnitCell);
end
%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'on')
    [Ceq2, DCeq2]=getFace(extrudedUnitCell);
%     [Ceq3, DCeq3]=getConvFace(extrudedUnitCell);
end
%FINAL CONSTRAINTS
Ceq=[Ceq1; Ceq2; Ceq3];
DCeq=[DCeq1; DCeq2; DCeq3]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aeq, Beq]=linearConstr(unitCell,extrudedUnitCell,opt)

%FIX NODE CONSTRAINTS
%IMPROVE FOLLOWING - AUTOMATIC DEPENDING ON NODES OF FACE 1
nodeFix=extrudedUnitCell.face{1};
e1=extrudedUnitCell.node(nodeFix(2),:)-extrudedUnitCell.node(nodeFix(1),:);
e2=extrudedUnitCell.node(nodeFix(3),:)-extrudedUnitCell.node(nodeFix(2),:);
e1=e1/norm(e1);
e2=e2/norm(e2);
e3=cross(e2,e1);
e3=e3/norm(e3);

Aeq=zeros(6,3*size(extrudedUnitCell.node,1));
Beq=zeros(6,1);
Aeq(1,3*nodeFix(2)-2)=1;
Aeq(2,3*nodeFix(2)-1)=1;
Aeq(3,3*nodeFix(2))=1;
Aeq(4,3*nodeFix(1)-2:3*nodeFix(1))=e3;
Aeq(5,3*nodeFix(1)-2:3*nodeFix(1))=e2;
Aeq(6,3*nodeFix(3)-2:3*nodeFix(3))=e3;

%MERGE NODES AT INITIALLY SAME LOCATION
rep=size(Aeq,1);
for i=1:size(unitCell.internalFacePairs,1)   
    for j=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)})
        index1=unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)}(j);
        for k=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)})
            index2=unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)}(k);
            if norm(extrudedUnitCell.node(index2,:)'-extrudedUnitCell.node(index1,:)')<opt.Lextrude/1e6
                rep=rep+1;
                Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                Aeq(3*rep-2:3*rep,3*index1-2:3*index1)=[1 0 0; 0 1 0; 0 0 1];
                Aeq(3*rep-2:3*rep,3*index2-2:3*index2)=[-1 0 0; 0 -1 0; 0 0 -1];
                Beq(3*rep-2:3*rep)=0;
            end
        end
     end
end

%PERIODIC NODAL CONSTRAINTS
if strcmp(opt.periodic,'on')
    nref=length(extrudedUnitCell.ref);
    rep=size(Aeq,1);
    for i=1:size(unitCell.possibleAlpha,1)
        for j=1:size(extrudedUnitCell.node,1)-nref
            coor1=extrudedUnitCell.node(j,:)';
            for k=1:size(extrudedUnitCell.node,1)-nref
                coor2=extrudedUnitCell.node(k,:)';
                if norm(coor2-coor1-unitCell.l'*unitCell.possibleAlpha(i,:)')<1e-6
                    rep=rep+1;
                    %sprintf('%d, node 1 = %d, node 2 =%d',[rep,j,k])
                    Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
                    Aeq(3*rep-2:3*rep,3*j-2:3*j)=[1 0 0; 0 1 0; 0 0 1];
                    Aeq(3*rep-2:3*rep,3*k-2:3*k)=[-1 0 0; 0 -1 0; 0 0 -1];
                    for l=1:nref
                        Aeq(3*rep-2:3*rep,3*extrudedUnitCell.ref(l)-2:3*extrudedUnitCell.ref(l))=unitCell.possibleAlpha(i,l)*[-1 0 0; 0 -1 0; 0 0 -1];
                    end
                    Beq(3*rep-2:3*rep)=0;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDGE LENGTH AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dEdge, Jedge]=getEdge(extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
Jedge=zeros(length(extrudedUnitCell.edge),size(extrudedUnitCell.node,1)*3);   
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=L-extrudedUnitCell.edgeL(i);            
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))=-dx/L;
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))=dx/L;
end

function [dEdge, Jedge]=getEdgeNorm(extrudedUnitCell)
dEdge=zeros(size(extrudedUnitCell.edge,1),1);
Jedge=zeros(length(extrudedUnitCell.edge),size(extrudedUnitCell.node,1)*3);   
for i=1:size(extrudedUnitCell.edge,1)
    coor1=extrudedUnitCell.node(extrudedUnitCell.edge(i,1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.edge(i,2),:);
    dx=coor2-coor1;
    L=sqrt(dx*dx');
    dEdge(i)=(L-extrudedUnitCell.edgeL(i))/extrudedUnitCell.edgeL(i);            
    Jedge(i,3*extrudedUnitCell.edge(i,1)-2:3*extrudedUnitCell.edge(i,1))=-dx/L/extrudedUnitCell.edgeL(i);
    Jedge(i,3*extrudedUnitCell.edge(i,2)-2:3*extrudedUnitCell.edge(i,2))=dx/L/extrudedUnitCell.edgeL(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FACE OUT OF PLANE DEFORMATION AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dFace, Jface]=getFace(extrudedUnitCell)
rep=0;
tnf=0;
for i=1:length(extrudedUnitCell.face)
    tnf=tnf+length(extrudedUnitCell.face{i})-3;
end
dFace=zeros(tnf,1);
Jface=zeros(tnf,size(extrudedUnitCell.node,1)*3);
for i=1:length(extrudedUnitCell.face)
    coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(2),:);
    coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(3),:);
    a=cross(coor2-coor1,coor3-coor1);
    for j=1:length(extrudedUnitCell.face{i})-3
        rep=rep+1;
        coor4=extrudedUnitCell.node(extrudedUnitCell.face{i}(3+j),:);
        Jface(rep,3*extrudedUnitCell.face{i}(1)-2:3*extrudedUnitCell.face{i}(1))=cross((coor3-coor2),(coor3-coor4));
        Jface(rep,3*extrudedUnitCell.face{i}(2)-2:3*extrudedUnitCell.face{i}(2))=cross((coor3-coor1),(coor4-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3)-2:3*extrudedUnitCell.face{i}(3))=cross((coor4-coor1),(coor2-coor1));
        Jface(rep,3*extrudedUnitCell.face{i}(3+j)-2:3*extrudedUnitCell.face{i}(3+j))=cross((coor2-coor1),(coor3-coor1));
        dFace(rep)=(coor4-coor1)*a';
    end
end

function [dFace, Jface]=getConvFace(extrudedUnitCell)
rep=0;
tnf=0;
for i=1:length(extrudedUnitCell.face)
    tnf=tnf+length(extrudedUnitCell.face{i})-1;
end
dFace=zeros(tnf,1);
Jface=zeros(tnf,size(extrudedUnitCell.node,1)*3);
for i=1:length(extrudedUnitCell.face)
    endpoint = length(extrudedUnitCell.face{i});
    coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(2),:);
    coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(3),:);
    a=cross(coor2-coor1,coor3-coor2);
    a = a/norm(a);
    for j=1:endpoint-1
        rep=rep+1;
        coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(mod(j,endpoint)+1),:);
        coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(mod(j+1,endpoint)+1),:);
        coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(mod(j+2,endpoint)+1),:);
        b = cross(coor2-coor1,coor3-coor2);
        b = b/norm(b);
        dFace(rep) = a*b'-1;
        Jface(rep,3*extrudedUnitCell.face{i}(1)-2:3*extrudedUnitCell.face{i}(1))=b;
        Jface(rep,3*extrudedUnitCell.face{i}(j+1)-2:3*extrudedUnitCell.face{i}(j+1))=a;
        
        
%         coor4=extrudedUnitCell.node(extrudedUnitCell.face{i}(3+j),:);
%         Jface(rep,3*extrudedUnitCell.face{i}(1)-2:3*extrudedUnitCell.face{i}(1))=cross((coor3-coor2),(coor3-coor4));
%         Jface(rep,3*extrudedUnitCell.face{i}(2)-2:3*extrudedUnitCell.face{i}(2))=cross((coor3-coor1),(coor4-coor1));
%         Jface(rep,3*extrudedUnitCell.face{i}(3)-2:3*extrudedUnitCell.face{i}(3))=cross((coor4-coor1),(coor2-coor1));
%         Jface(rep,3*extrudedUnitCell.face{i}(3+j)-2:3*extrudedUnitCell.face{i}(3+j))=cross((coor2-coor1),(coor3-coor1));
%         dFace(rep)=(coor4-coor1)*a';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HINGE ANGLE AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, Jhinge]=getHinge(extrudedUnitCell)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
Jhinge=zeros(size(extrudedUnitCell.nodeHingeEx,1),size(extrudedUnitCell.node,1)*3);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [Jhinge(i,index),theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%         YUN's addition                                        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [deformationV, exitFlag] = ...
    gradDescent(original_extrudedUnitCell, opt, u0, max_iter)
% gradient descent method with a penalty imposed for preventing faces from
% corssing over each other.
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% INPUT
% original_extrudedUnitCell - extruded unit cell without applying any
%                             deformations
% opt      - options
% u0       - coordinates of all nodes
% max_iter - maximum number of iterations
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% OUTPUT
% deformationV - an array containing the deformation in each iteration
% exitFlag     - if exitFlag = 1, convergence occured within max_iter
%                if exitFlag = 0, max_iter is reached without convergence
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun


original_extrudedUnitCell.angleConstr = [];
exitFlag = 2;

% getting ready for the loop
u = u0(:);
deformationV(:,1) = u0(:);
ct = 0;
difference = Inf;
e_prev = Inf;
% fprintf('\n_ %d _',difference)
while ct <= max_iter && difference > opt.gradDescTol
    % energy of current u
    [e, de, ~, ~, ~, ~, theta] = Energy(u, original_extrudedUnitCell, opt);
    % uses penalty to prevent faces from crossing over
    e = e + penalty(theta);
   
    
    % get the new solution from jacobian
    u_new = u - opt.gradDescStep .* de(:);
    
    % calculate difference between two consecutive solutions
    difference = abs(e - e_prev) / ...
        size(original_extrudedUnitCell.edge,1); % normalised to the number of hinges
    
    % prepare for next step
    ct = ct + 1;
    if isequal(mod(ct,1000), 0)
        deformationV(:,size(deformationV,2)+1) = u_new(:);
    end
    u = u_new;
    e_prev = e;
%     fprintf(' %d_\n', e)
end

if ct > max_iter
    exitFlag = 0;
    deformationV(:,size(deformationV,2)+1) = u(:);
end

if difference <= opt.gradDescTol
    exitFlag = 1;
    deformationV(:,size(deformationV,2)+1) = u(:);
end

function ePen = penalty(theta)
% A penalty function which punished theta value outside the range of
% [-0.95*pi, 0.95*pi].
% The range is not [-pi, pi]
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% last modified on May 22, 2017
% yun

ePen = sum(max(0, 0.5 .* (theta - 0.95*pi) .* (theta + 0.95*pi)));