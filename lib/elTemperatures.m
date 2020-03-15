%Deslocamentos por elemento, cada linha apresenta com base na tabela de
%connectividades, os deslocamentos de cada nó num elemento
function [elT]=elTemperatures(T1,dOf,connectInfo,numElements)

indexB=getGlobalDoF(connectInfo, dOf);

indexB(:,1)=[];
% indexB(:,1)=[];
elT= [];

for i=1:1:numElements
    
   
    elT=[elT ; T1(sort(indexB(i,:)))'];
    
end 
elT=elT';
end 
