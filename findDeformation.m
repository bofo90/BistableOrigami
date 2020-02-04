function findDeformation(opt, des, xrep, yrep, ang, kap)
file = opt.file;
for d = des
    opt.vertexType = d;
    
    if strcmp(opt.template, 'SingleVertex')
        xrep = 1;
        yrep = 1;
    end
    
    for x = xrep
        opt.xrep = x;
        for y = yrep
            opt.yrep = y;

            [extrudedUnitCell,opt]=obtainOrigami(opt);
            if strcmp(opt.template, 'Tessellation')
                opt.file = strcat(file,num2str([opt.xrep,opt.yrep],'%d_'));
            elseif strcmp(opt.template, 'SingleVertex')
                opt.file = strcat(file,num2str(opt.angDesign*180/pi,'%.2f_'));
            end
            
            ChangeParam(extrudedUnitCell,opt, ang, kap);
            
        end
    end
end


function ChangeParam(extrudedUnitCell,opt, ang, kap)
file = opt.file;
for a = ang
    opt.restang = a;
    extrudedUnitCell.theta = ones(size(extrudedUnitCell.theta,1),1)*opt.restang;
    for k = kap
        opt.Khinge = k;

        if strcmp(opt.template,'Tessellation')
            opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.tessellationType,'/',opt.vertexType,file,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));
        elseif strcmp(opt.template,'SingleVertex')
            opt.file = strcat(pwd,'/Results/',opt.template,num2str(opt.numVert),'/',opt.vertexType,file,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge));
        end
        
        fprintf('Start folding...\n');
        metadataFile(opt, extrudedUnitCell);
        switch opt.analysisType
            case 'single'
                nonlinearFoldingOne(extrudedUnitCell, opt, opt.angleConstrFinal(1).val);
            case 'multiple'
                nonlinearFoldingMulti(extrudedUnitCell, opt, opt.angleConstrFinal(1).val);
            case 'randomPert'
                nonlinearFoldingRand(extrudedUnitCell, opt, opt.RandstDev);
        end
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NON-LINEAR ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function nonlinearFoldingMulti(extrudedUnitCell,opt,angtemp)

savefile = opt.file;

%INITIALIZE LINEAR CONSTRAINTS
[Aeq, Beq]=linearConstr(extrudedUnitCell,opt);

foldAngl = 6;
angles1_1 = 0:(foldAngl*pi/180):(0.985*pi);
angles1_2 = -angles1_1;
angles1 = [angles1_1 angles1_2(2:end)];
kappas = [logspace(-3.875,-2.625,6) logspace(0.375,0.5,2)];
    
%%%%%% Folding part %%%%%%
%Run the Folding of the structure
for kappa = 1:size(kappas,2)
    extrudedUnitCell.angleConstr=[];
    opt.angleConstrFinal = [];
    result = [];
    E=[];
    exfl= [];
    
    u0=zeros(3*size(extrudedUnitCell.node,1),1);
    theta0=zeros(size(extrudedUnitCell.theta));
    flag1 = false;
    
    opt.Khinge = kappas(kappa);
    %Create file for saving the results
    opt.file = strcat(savefile,sprintf('/RestAng_%.3f/kappa_%2.5f', opt.restang, opt.Khinge)); 
    folderName = strcat(opt.file,'/mat');
    if ~exist(folderName, 'dir')
        mkdir(folderName);
    end
    
    opt.angleConstrFinal(1).val = [angtemp;[2 0; 4 0]];
    initialiseGlobalx(u0, theta0);
    [V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, 1, 1, Aeq, Beq);
    [result, thetainit, uinit] = SaveResultPos(result, opt, V, output, 1);
    
    u0 = uinit;
    theta0 = thetainit;
    
    for ang1 = 1:size(angles1,2)

        if angles1(ang1) < 0 && ~flag1
            u0 = uinit;
            theta0 = thetainit;
            result.deform(2) = [];
            flag1 = true;
        end 
        
        opt.angleConstrFinal(2).val = [angtemp; [2 angles1(ang1); 4 0]];
        initialiseGlobalx(u0, theta0);
        [V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, 2, 1, Aeq, Beq);
        [result, theta0, u0] = SaveResultPos(result, opt, V, output, 2);

        u1 = u0;
        theta1 = theta0;  
        flag2 = false;
        
        for ang2 = 1:size(angles1,2)

            if angles1(ang2) < 0 && ~flag1
                u1 = u0;
                theta1 = theta0;
                result.deform(3) = [];
                flag2 = true;
            end 

            opt.angleConstrFinal(3).val = [angtemp; [2 angles1(ang1); 4 angles1(ang2)]];
            initialiseGlobalx(u1, theta1);
            [V, exfl, output, E] = FoldStructure(u1, E, exfl, extrudedUnitCell, opt, 3, 1, Aeq, Beq);
            [result, theta1, u1] = SaveResultPos(result, opt, V, output, 3);

            %%%%%% Releasing part %%%%%%
            %change algorithm for releasing
            opt.angleConstrFinal(4).val = [];
            initialiseGlobalx(u1, theta1);
            [V, exfl, output, E] = FoldStructure(u1, E, exfl, extrudedUnitCell, opt, 4, 1, Aeq, Beq);
            [result, ~, ~] = SaveResultPos(result, opt, V, output, 4);

            result = SaveResultEnergy(result, E, exfl, opt);
            result.angNum = [kappa ang1 ang2];
            result.angVal = [angles1(ang1) angles1(ang2)];
            result.kappa = opt.Khinge;

            %Save the result in a file
            fileName = strcat(folderName,'/','_Ang1_',int2str(ang1),'_Ang2_',int2str(ang2),'.mat');
            save(fileName, 'result');

            result.deform(4) = [];
            fclose('all');
        end
        
        result.deform(3) = [];
    end
    
    result.deform(2) = [];
end

function nonlinearFoldingOne(extrudedUnitCell,opt,angtemp)

%INITIALIZE LINEAR CONSTRAINTS
[Aeq, Beq]=linearConstr(extrudedUnitCell,opt);

%Save some variables
opt.angleConstrFinal(1).val = angtemp;
extrudedUnitCell.angleConstr=[];
result = [];
E=[];
exfl= [];
u0=zeros(3*size(extrudedUnitCell.node,1),1);
theta0=zeros(size(extrudedUnitCell.theta));

%Create file for saving the results
folderName = strcat(opt.file,'/mat');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
    
%%%%%% Folding part %%%%%%
%Run the Folding of the structure
initialiseGlobalx(u0, theta0);
[V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, 1, opt.steps, Aeq, Beq);
[result, theta1,u1] = SaveResultPos(result, opt, V, output, 1);

%%%%%% Releasing part %%%%%%
%change algorithm for releasing
opt.options.Algorithm = opt.relAlgor;
opt.angleConstrFinal(2).val = [];
opt.KtargetAngle = 0;

initialiseGlobalx(u1, theta1);
[V, exfl, output, E] = FoldStructure(u1, E, exfl, extrudedUnitCell, opt, 2, 1, Aeq, Beq);
[result, ~,~] = SaveResultPos(result, opt, V, output, 2);

%Return to original options
opt.options.Algorithm = opt.folAlgor;
result = SaveResultEnergy(result, E, exfl, opt);

%Save the result in a file
fileName = strcat(folderName,'/',mat2str(sign(opt.angleConstrFinal(1).val(:,2))'),'.mat');
save(fileName, 'result');


%Clear variables for next fold
clearvars result E exfl output;
fclose('all');

function nonlinearFoldingRand(extrudedUnitCell,opt,stdev)

%INITIALIZE LINEAR CONSTRAINTS
[Aeq, Beq]=linearConstr(extrudedUnitCell,opt);

%Save some variables
opt.angleConstrFinal(1).val = [];
extrudedUnitCell.angleConstr=[];
rng('shuffle');

%Create file for saving the results
folderName = strcat(opt.file,'/mat');
if ~exist(folderName, 'dir')
    mkdir(folderName);
end
    
opt.KtargetAngle = 0;
parfor i = 1:opt.numIterations
    
    %%%%%% Folding part %%%%%%
    %Perturb the structure
    result = [];
    E=[];
    exfl= [];
    u0 = randn(3*size(extrudedUnitCell.node,1),1)*stdev^2;
    theta0 = getSimpleAngle(u0, extrudedUnitCell);
    
    %%%%%% Releasing part %%%%%%
    %change algorithm for releasing
    initialiseGlobalx(u0, theta0);
    [V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, 1, 1, Aeq, Beq);
    [result, ~, ~] = SaveResultPos(result, opt, V, output, 1);

    %Return to original options
    result = SaveResultEnergy(result, E, exfl, opt);

    %Save the result in a file
    fileName = strcat(folderName,'/iter_',mat2str(i),'.mat');
    parsave(fileName, result);

end
%Clear variables for next fold
clearvars result E exfl output;
fclose('all');

function parsave(filename,result)

save(filename, 'result')

function [V, exfl, output, E] = FoldStructure(u0, E, exfl, extrudedUnitCell, opt, iter, steps, Aeq, Beq)

V = [];
V(:,1)=u0;
[E.E(1,iter),~,E.Eedge(1,iter),E.Ediag(1,iter),E.Eface(1,iter),E.Ehinge(1,iter),E.EtargetAngle(1,iter), theta0]=Energy(u0,extrudedUnitCell,opt);
exfl(1,iter) = 1;

exfl(2,iter) = 1;

%Run the Folding of the structure
if isempty(opt.angleConstrFinal(iter).val)
    fprintf('Angle contrain: None\n');
    steps = 1;
else
    fprintf(['Angle contrain:', mat2str(opt.angleConstrFinal(iter).val(:,1)') ,'\n']);
    angles = theta0(opt.angleConstrFinal(iter).val(:,1)) + linspace(0,1,steps+1).*(opt.angleConstrFinal(iter).val(:,2) - theta0(opt.angleConstrFinal(1).val(:,1)));
end

for anglestep = 2:steps+1
    if ~isempty(opt.angleConstrFinal(iter).val)
        opt.angleConstrFinal(iter).val(:,2) = angles(:,anglestep);
    end
    extrudedUnitCell.angleConstr=opt.angleConstrFinal(iter).val;
    % fprintf('Folding:\t');
    % t1 = toc;
    %Determine new equilibrium
    [V(:,2),~,exfltemp,output]=fmincon(@(u) Energy(u,extrudedUnitCell,opt),u0,[],[],Aeq,Beq,[],[],@(u) nonlinearConstr(u,extrudedUnitCell,opt),opt.options);
    u0 = V(:,2);
    if exfltemp ~= 1
        exfl(2,1) = exfltemp;
    end 
end
    
%Determine energy of that equilibrium
[E.E(2,iter),~,E.Eedge(2,iter),E.Ediag(2,iter),E.Eface(2,iter),E.Ehinge(2,iter),E.EtargetAngle(2,iter), ~]=Energy(u0,extrudedUnitCell,opt);
% t2 = toc;
% fprintf('time %1.2f, exitflag %d\n',t2-t1,exfl(2,iter));


function [result, lastAngle, lastPosition] = SaveResultPos(result, opt, Positions, minimizationOuput, state)

%Get Angles from global variable
angles = getGlobalAngles;     
switch opt.gethistory
    case 'on'
        Positions = getGlobalAllx;
    case 'off'
        angles = [angles(:,1) angles(:,end)];
end

%Arrange results in a better way to save them
result.deform(state).V=[Positions(1:3:end,end) Positions(2:3:end,end) Positions(3:3:end,end)];
result.deform(state).Ve=Positions(:,end);
result.deform(state).theta = angles(:,end);
result.deform(state).output = minimizationOuput;

if ~isfield(result.deform(state), 'interV')
    result.deform(state).interV = [];
end

for j=1:size(Positions,2)
    result.deform(state).interV(end+1).V=[Positions(1:3:end,j) Positions(2:3:end,j) Positions(3:3:end,j)];
    result.deform(state).interV(end).Ve=Positions(:,j);
    result.deform(state).interV(end).theta = angles(:,j);
end

lastAngle = angles(:,end);
lastPosition = Positions(:,end);

function result = SaveResultEnergy(result, E, exfl, opt)

result.E=E.E;
result.Eedge=E.Eedge;
result.Ediag=E.Ediag;
result.Eface=E.Eface;
result.Ehinge=E.Ehinge;
result.EtargetAngle=E.EtargetAngle;    
result.exfl = exfl;
result.numMode=length(result.deform);
result.anglConstr = opt.angleConstrFinal(1).val;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENERGY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, dE,Eedge,Ediag,Eface,Ehinge,EtargetAngle,theta]=Energy(u,extrudedUnitCell,opt)

E=0; dE=zeros(3*size(extrudedUnitCell.node,1),1);
Eedge=0;
Ediag = 0;
Eface=0;
Ehinge=0;
EtargetAngle=0;

%APPLY DEFORMATION NODES
extrudedUnitCellPrev = extrudedUnitCell;
extrudedUnitCell.node = extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%ENERGY ASSOCIATED TO EDGE STRETCHING
if strcmp(opt.constrEdge,'off')
    [dEdge, Jedge]=getEdgeNorm(extrudedUnitCell);
%     %shearing energy
%     Ediag=1/2*opt.Kdiag*sum(dEdge(extrudedUnitCell.diagonals).^2);
%     dE=dE+opt.Kdiag*Jedge(extrudedUnitCell.diagonals,:)'*dEdge(extrudedUnitCell.diagonals);
%     %streching energy
%     notdiagonal = 1:size(extrudedUnitCell.edge,1);
%     notdiagonal(extrudedUnitCell.diagonals) = [];
    Eedge=1/2*opt.Kedge*sum(dEdge.^2);
    dE=dE+opt.Kedge*(Jedge'*dEdge);
end

%ENERGY ASSOCIATED TO FACE BENDING
if strcmp(opt.constrFace,'off')
    [dFace, Jface]=getFace(extrudedUnitCell);
    Eface=opt.Kface/2*sum(dFace.^2);
    dE=dE+opt.Kface*Jface'*dFace;
end

%ENERGY ASSOCIATED TO HINGE BENDING
[theta, Jhinge]=getHinge(extrudedUnitCell, extrudedUnitCellPrev);
ThetaVar = theta-extrudedUnitCell.theta;
Ehinge=4*opt.Khinge/(opt.restang^4)*sum(0.25*ThetaVar.^4+...
    ThetaVar.^3.*extrudedUnitCell.theta+...
    ThetaVar.^2.*extrudedUnitCell.theta.^2);
dE=dE+4*opt.Khinge/(opt.restang^4)*(Jhinge'*(ThetaVar.^3+...
    3*ThetaVar.^2.*extrudedUnitCell.theta+...
    2*ThetaVar.*extrudedUnitCell.theta.^2));

%ENERGY ASSOCIATED TO TARGET HINGE ANGLES
if size(extrudedUnitCell.angleConstr,1)==0
    dtheta=[];
else
    dtheta=theta(extrudedUnitCell.angleConstr(:,1))-extrudedUnitCell.angleConstr(:,2);
    dE=dE+opt.KtargetAngle*Jhinge(extrudedUnitCell.angleConstr(:,1),:)'*dtheta;
end
EtargetAngle=1/2*opt.KtargetAngle*sum(dtheta.^2);

%TOTAL ENERGY
E=Eedge+Ediag+Eface+Ehinge+EtargetAngle;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NONLINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,Ceq,DC,DCeq]=nonlinearConstr(u,extrudedUnitCell,opt)

C1=[]; C2=[]; C3=[];
DC1=[]; DC2=[];
Ceq1=[]; Ceq2=[]; Ceq3=[];
DCeq1=[]; DCeq2=[]; DCeq3=[];

%APPLY DEFORMATION NODES
extrudedUnitCellPrev = extrudedUnitCell;
extrudedUnitCell.node=extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];

%%%%%%%%%%%%%%%%%%%%%%
%INEQUALITY CONSTRAINS
%%%%%%%%%%%%%%%%%%%%%%
%MAXIMUM AND MINIMUM ANGLES
if ~isnan(opt.constAnglePerc)
    [angles, Dangles]=getHinge(extrudedUnitCell, extrudedUnitCellPrev);
    C1 = [-angles-pi*opt.constAnglePerc; angles-pi*opt.constAnglePerc];
    DC1 = [-Dangles; Dangles];
end
%MAXIMUM AND MINIMUM EDGE STRECHING
if strcmp(opt.constrEdge,'off') && ~isnan(opt.maxStretch)
    [normStrech, DnormStrech]=getEdgeNorm(extrudedUnitCell);
    C2 = [-normStrech-opt.maxStretch; normStrech-opt.maxStretch];
    DC2 = [-DnormStrech; DnormStrech];
end
%MINIMUM FACE AREA
if ~isnan(opt.maxArea)
    dFace = getFace2(extrudedUnitCell);
    C3 = [opt.maxArea-dFace];
end

C = [C1; C2; C3];
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
end
%FINAL CONSTRAINTS
Ceq=[Ceq1; Ceq2; Ceq3];
DCeq=[DCeq1; DCeq2; DCeq3]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LINEAR CONSTRAINTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Aeq, Beq]=linearConstr(extrudedUnitCell, opt)

%FIX NODE CONSTRAINTS
%IMPROVE FOLLOWING - AUTOMATIC DEPENDING ON NODES OF FACE 1
nodeFix=extrudedUnitCell.face{1};
e1=extrudedUnitCell.node(nodeFix(1),:)-extrudedUnitCell.node(nodeFix(2),:);
e2=extrudedUnitCell.node(nodeFix(3),:)-extrudedUnitCell.node(nodeFix(2),:);
e1=e1/norm(e1);
e2=e2/norm(e2);
e3=cross(e2,e1);
e3=e3/norm(e3);

Aeq=zeros(6,3*size(extrudedUnitCell.node,1));
Beq=zeros(6,1);
Aeq(1,3*nodeFix(1)-2)=1;
Aeq(2,3*nodeFix(1)-1)=1;
Aeq(3,3*nodeFix(1))=1;
Aeq(4,3*nodeFix(2)-2:3*nodeFix(2))=e3;
Aeq(5,3*nodeFix(2)-2:3*nodeFix(2))=e2;
Aeq(6,3*nodeFix(3)-2:3*nodeFix(3))=e3;

%MERGE NODES AT INITIALLY SAME LOCATION
% rep=size(Aeq,1);
% for i=1:size(unitCell.internalFacePairs,1)   
%     for j=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)})
%         index1=unitCell.Polyhedron(unitCell.internalFacePairs(i,2)).faceNodeExtrude{unitCell.internalFacePairs(i,1)}(j);
%         for k=1:length(unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)})
%             index2=unitCell.Polyhedron(unitCell.internalFacePairs(i,4)).faceNodeExtrude{unitCell.internalFacePairs(i,3)}(k);
%             if norm(extrudedUnitCell.node(index2,:)'-extrudedUnitCell.node(index1,:)')<opt.Lextrude/1e6
%                 rep=rep+1;
%                 Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
%                 Aeq(3*rep-2:3*rep,3*index1-2:3*index1)=[1 0 0; 0 1 0; 0 0 1];
%                 Aeq(3*rep-2:3*rep,3*index2-2:3*index2)=[-1 0 0; 0 -1 0; 0 0 -1];
%                 Beq(3*rep-2:3*rep)=0;
%             end
%         end
%      end
% end

%PERIODIC NODAL CONSTRAINTS
% if strcmp(opt.periodic,'on')
%     nref=length(extrudedUnitCell.ref);
%     rep=size(Aeq,1);
%     for i=1:size(unitCell.possibleAlpha,1)
%         for j=1:size(extrudedUnitCell.node,1)-nref
%             coor1=extrudedUnitCell.node(j,:)';
%             for k=1:size(extrudedUnitCell.node,1)-nref
%                 coor2=extrudedUnitCell.node(k,:)';
%                 if norm(coor2-coor1-unitCell.l'*unitCell.possibleAlpha(i,:)')<1e-6
%                     rep=rep+1;
%                     %sprintf('%d, node 1 = %d, node 2 =%d',[rep,j,k])
%                     Aeq(3*rep-2:3*rep,:)=zeros(3,size(extrudedUnitCell.node,1)*3);
%                     Aeq(3*rep-2:3*rep,3*j-2:3*j)=[1 0 0; 0 1 0; 0 0 1];
%                     Aeq(3*rep-2:3*rep,3*k-2:3*k)=[-1 0 0; 0 -1 0; 0 0 -1];
%                     for l=1:nref
%                         Aeq(3*rep-2:3*rep,3*extrudedUnitCell.ref(l)-2:3*extrudedUnitCell.ref(l))=unitCell.possibleAlpha(i,l)*[-1 0 0; 0 -1 0; 0 0 -1];
%                     end
%                     Beq(3*rep-2:3*rep)=0;
%                 end
%             end
%         end
%     end
% end


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

function dFace=getFace2(extrudedUnitCell)

dFace=zeros(length(extrudedUnitCell.face),1);
for i=1:length(extrudedUnitCell.face)
    coor1=extrudedUnitCell.node(extrudedUnitCell.face{i}(1),:);
    coor2=extrudedUnitCell.node(extrudedUnitCell.face{i}(2),:);
    coor3=extrudedUnitCell.node(extrudedUnitCell.face{i}(3),:);
    a=cross(coor2-coor1,coor3-coor1);
    dFace(i) = 1/2*sqrt(a*a');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HINGE ANGLE AND JACOBIAN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta, Jhinge]=getHinge(extrudedUnitCell, extrudedUnitCellPrev)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
thetaPrev = getGlobalx(extrudedUnitCellPrev);
Jhinge=zeros(size(extrudedUnitCell.nodeHingeEx,1),size(extrudedUnitCell.node,1)*3);
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [Jhinge(i,index),theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
    if abs(thetaPrev(i)-theta(i)) >= 0.50*2*pi
        theta(i) = theta(i)+sign(thetaPrev(i))*2*pi;
%         fprintf('angle %d energy\n', i);
    end
end

function theta = getSimpleAngle(u, extrudedUnitCell)
theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
extrudedUnitCell.node = extrudedUnitCell.node+[u(1:3:end) u(2:3:end) u(3:3:end)];
for i=1:size(extrudedUnitCell.nodeHingeEx,1)
    extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:);
    index(1:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-2;
    index(2:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:)-1;
    index(3:3:12)=3*extrudedUnitCell.nodeHingeEx(i,:);
    [~,theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions to extract information from Global Variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function initialiseGlobalx(u0, theta)
global poop
poop.x = [];
poop.theta = [];
poop.x = [poop.x u0];
poop.theta = [poop.theta theta];
poop.flag = 0;

function theta = getGlobalx(extrudedUnitCell)
global poop
if poop.flag
    theta=zeros(size(extrudedUnitCell.nodeHingeEx,1),1);
    extrudedUnitCell.node=extrudedUnitCell.node+[poop.x(1:3:end,end) poop.x(2:3:end,end) poop.x(3:3:end,end)];
    for i=1:size(extrudedUnitCell.nodeHingeEx,1)
        [~,theta(i)]=JacobianHinge(extrudedUnitCell.node(extrudedUnitCell.nodeHingeEx(i,:),:));
        if abs(poop.theta(i,end)-theta(i)) >= 0.5*2*pi
            theta(i) = theta(i)+sign(poop.theta(i,end))*2*pi;
%             fprintf('angle %d global\n', i);
        end
    end
    poop.theta = [poop.theta theta];
    %turn the flag off after analising the angles just after an iteration,
    %so it doesnt analyse it every time you get the energy
    poop.flag = 0;
else
    theta = poop.theta(:,end);
end

function r = getGlobalAllx
global poop
r = poop.x;

function theta = getGlobalAngles
global poop
theta = poop.theta;

