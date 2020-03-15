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
function [ elStiffness ] = assemGlobalStiff( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                             elNodes, numNodes, connectInfo, D, elThick)
%% Symbolic variables

syms csi eta

%% Initiating variables

connectInfo(:,1)=[];    %Cleaning element labels

dimen=dOf*elNodes;      %Setting global dimensions

dNCsiEta=zeros(2,elNodes);  %Derivative matrix for jacobian computation
dNCsiEta=sym(dNCsiEta);

dNCsiEta=[diff(nodeShapeFun,csi).';diff(nodeShapeFun,eta).'];

%TODO-verify condition for other types of analysis
B=zeros(3,dimen);   %Shape function derivative matrix

K=zeros(dOf*numNodes,dOf*numNodes);  %Global stiffness matrix

l=1:dOf:dimen;    %Index for derivative placement
m=0:dOf:dimen;    %Index for derivative placement
m(1)=[];

indexB=getGlobalDoF(connectInfo,dOf);

for i=1:size(connectInfo,1)

    for j=1:size(intPts,1)
        
        dN=dNCsiEta;
        dN=subs(dN,[csi,eta],[intPts(j,1),intPts(j,2)]);  %Calculating derivatives in the integration point
        dN=double(dN);                                    %Converting from symbolic to double
        J=dN*elNodeCoord(:,:,i);                          %Jacobian
        
        dNdXdY=J\dN;           %Transforming to global coordinates
        
        B(1,l)=dNdXdY(1,:);    %Derivative matrix on global coordinates
        B(2,m)=dNdXdY(2,:);
             
        B(3,l)=dNdXdY(2,:);
        B(3,m)=dNdXdY(1,:);        
        
        wCsi=intPts(j,3);
        wEta=intPts(j,4);
        
        %Direct assembly of the element stiffness matrix into global stiffness matrix
        K(indexB(i,:),indexB(i,:))=K(indexB(i,:),indexB(i,:))+B'*D{1}*B*wCsi*wEta*det(J)*elThick;

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
%>@date 28-Oct-2018
%======================================================================
