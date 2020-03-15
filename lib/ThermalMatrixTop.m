%======================================================================
%>@brief Computes the global stiffness matrix of the problem.
%>
%>@param nodeShapeFun (@b sym @b matrix) Matrix of nodal shape functions
%>@param intPts (@b matrix) Matrix of integration points and weighting
%>factors
%>@param elNodeCoord (@b matrix) Multidimensional array of nodal coordinates
%>@param dOf (@b integer) Number of degrees of freedom
%>@param elNodes (@b integer) Number of nodes per element
%>@param connectInfo (@b matrix) Element connectivity matrix
%>@param D (@b matrix) Elasticity matrix
%>@param elThick (@b double) Element thickness
%>
%>@retval elStiffness (@b matrix) Global stiffness matrix
%>
%>@details
%>For each element, computes the corresponding stiffness matrix and
%>proceeds to assemble it into the global stiffness matrix, based on the
%>global degrees of freedom associated with the element.
%>@todo
%>Devise a way to reduce the number of parameters needed to call the
%>function.
%>@todo
%>The use of symbolic variables may be affecting the performance for models
%>with higher number of elements. Using hardcodded shape functions may
%>improve performance in those situations. Additional code changes to this 
%>function need be considered to evaluate this. The algorithm used to generate 
%>the shape functions also needs to be changed.
%======================================================================
function [ elStiffness] = ThermalMatrixTop(dOf,elNodes, numNodes, connectInfo,x,penal,tP,elementK)

%% Initiating variables

connectInfo(:,1)=[];    %Cleaning element labels

dimen=dOf*elNodes;      %Setting global dimensions

%TODO-verify condition for other types of analysis
elK=sparse(dOf*elNodes,dOf*elNodes);
K=sparse(dOf*numNodes,dOf*numNodes);  %Global stiffness matrix

indexB=getGlobalDoF(connectInfo,dOf);

for i=1:size(connectInfo,1)

        %Direct assembly of the element stiffness matrix into global stiffness matrix
        K(sort(indexB(i,:)),sort(indexB(i,:)))=K(sort(indexB(i,:)),sort(indexB(i,:)))+(x(i)^penal).*elementK(:,:,i);
   
  
    
ind=find(tP==1);
for i=1:size(ind,2)
    K(ind(i),ind(i))=K(ind(i),ind(i))+10^6;
end

end
%% Output

elStiffness=K;


end

%======================================================================
%>@file elementStiffnessMatrix.m
%>@brief Algorithm to define the global stiffness matrix.
%>@details
%>
%>@author Rúben Lourenço
%>@dat
%28-Oct-2018
%======================================================================


