
function [ elStiffness ] = Globalteste(elementK,numElements,connectInfo,dOf,numNodes,x,penal)
%% Initiating variables


K=sparse(dOf*numNodes,dOf*numNodes);  
indexB=getGlobalDoF(connectInfo,dOf);
indexB(:,1)=[];
indexB(:,1)=[];
for i=1:1:numElements
    
    
        K(sort(indexB(i,:)),sort(indexB(i,:)))=K(sort(indexB(i,:)),sort(indexB(i,:)))+(x(i)^penal)*elementK(:,:,i);
        
end
%% Output

elStiffness=K;


end