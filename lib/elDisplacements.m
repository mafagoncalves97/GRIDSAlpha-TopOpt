%Deslocamentos por elemento, cada linha apresenta com base na tabela de
%connectividades, os deslocamentos de cada nó num elemento

function [elU]=elDisplacements(U1,dOf,connectInfo,numElements)

indexB=getGlobalDoF(connectInfo, dOf);

indexB(:,1)=[];
indexB(:,1)=[];
elU= [];

for i=1:1:numElements
    
   
    elU=[elU ; U1(sort(indexB(i,:)))'];
    
end 
elU=elU';
end 
