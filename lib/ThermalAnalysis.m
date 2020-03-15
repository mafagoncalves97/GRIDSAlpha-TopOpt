function [Q,T1]=ThermalAnalysis(elNodeCoord,elNodes,connectInfo,numNodes,nodeShapeFun,intPts,k,t,nodeCoord)
%% Degrees of freedom
dOf=1;
%%Condições-fronteira
ymax=max(nodeCoord(:,3));
xmax=max(nodeCoord(:,2));
ymin=min(nodeCoord(:,3));
xmin=min(nodeCoord(:,2));
nEly=nodeCoord(:,3)==max(nodeCoord(:,3));%ymax
nElx=nodeCoord(:,2)==max(nodeCoord(:,2)); %xmax
nElymin=nodeCoord(:,3)==min(nodeCoord(:,3));%ymin
nElxmin=nodeCoord(:,2)==min(nodeCoord(:,2)); %xmin

%%BC-MO
% fixedNodesth=find(nodeCoord(:,2)==xmax&nodeCoord(:,3)==ymin);
% sSnodesth.BC=find(nodeCoord(:,2)==xmin|nodeCoord(:,3)==ymax);

%%BC-termico
% fixedNodesth=find(nodeCoord(:,2)==xmax);
q=zeros(numNodes,1);
sSnodesth.BC=(find(nElx==1|nElymin==1|nElxmin==1));
nodesFlux=find(nEly==1);

nodes=nodeCoord(:,1)';
activeDoFth=nodes;

q(nodesFlux)=0.5;
%activeDoFth(fixedNodesth)=[];
%% Input material properties
C=[k 0;0 k];
tP=zeros(1,9);
%% Temperatura prescrita
if ~isempty(sSnodesth.BC)
    PrescTemp=30;
    nodesPrescrT=sSnodesth.BC;
    tPrescrita(nodesPrescrT)=PrescTemp*ones(1,size(nodesPrescrT,2));
    tP=tPrescrita==30;
end
%%Fluxos
 if sum(tP)>0
     q(tP==1)=tPrescrita(tP==1)*10^6;
 end
 
% flux=5;
% for i=1:size(nodesFlux,2)
%     %Node exterior
%     if (sum(nodeCoord(nodesFlux,2)==xmax)==size(nodesFlux,2)|sum(nodeCoord(nodesFlux,2)==xmin)==size(nodesFlux,2))
%         if (nodeCoord(i,2)==xmax|nodeCoord(i,2)==xmin)
%             q(nodesFlux(i))=flux/(round(sum(nElx)/2)*2);
%         else
%             q(nodesFlux(i))=flux/(round(sum(nElx)/2));
%         end
%     %Node exterior
%     elseif (sum(nodeCoord(nodesFlux,2)==ymax)==size(nodesFlux,2)|sum(nodeCoord(nodesFlux,2)==ymin)==size(nodesFlux,2))
%         if (nodeCoord(i,2)==ymax|nodeCoord(i,2)==ymin)
%             q(nodesFlux(i))=flux/(round(sum(nEly)/2)*2);
%         else
%             q(nodesFlux(i))=flux/(round(sum(nEly)/2));
%         end
%     else %interior
%         q(nodesFlux(i))=0;
%     end
% end
%% Element stiffness matrix
disp('> Assembling global stiffness matrix...')
%Assembly of the global stiffness matrix
K = assemGlobalStiffT(nodeShapeFun,intPts,elNodeCoord,dOf,...
                     elNodes,numNodes,connectInfo,C,t,tP);
%% Global system of equations
T=K(activeDoFth,activeDoFth)\q(activeDoFth);
T1=zeros(dOf*numNodes,1);
T1(activeDoFth)=T;
%% Reaction Forces
Q=K*T1;
end

