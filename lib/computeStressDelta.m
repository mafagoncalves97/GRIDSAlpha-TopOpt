%======================================================================
%>@brief Performs the stress recovery
%>
%>@param nodeShapeFun (@b sym @b matrix) Matrix of nodal shape functions
%>@param intPts (@b matrix) Matrix of integration points and weighting
%>factors
%>@param elNodeCoord (@b matrix) Multidimensional array of nodal coordinates
%>@param dOf (@b integer) Number of degrees of freedom
%>@param elNodes (@b integer) Number of nodes per element
%>@param numNodes (@b integer) Total number of nodes
%>@param numElements (@b integer) Total number of elements
%>@param connectInfo (@b matrix) Element connectivity matrix
%>@param D (@b matrix) Elasticity matrix
%>@param U (@b matrix) Displacement field
%>
%>@retval stressField (@b matrix) Matrix of stress fields
%>
%>@details
%>For each element, computes the stress fields at the integration points
%>and extrapolates the values to the corner nodes. At a node-level the
%>stresses are discontinuous, so a stress smoothing is performed for each
%>node, based on the inverse connectivity.
%>@note
%> <strong>`stressField`</strong> has dimensions:
%> - (numNodes x 4)
%>
%>@note
%>Each column of the matrix <strong>`stressField`</strong> corresponds to a
%>stress field, such that:
%> - **1st column**: Sx
%> - **2nd column**: Sy
%> - **3rd column**: Sxy
%> - **4th column**: Svm (von Mises)
%>@see
%> Theory manual<br>
%> [Chapter 4, Subsection 4.5.1 - Stress recovery](http://grids.web.ua.pt/wp-content/uploads/2018/09/DISSERTACAO_Ruben_Lourenco_MIEM_1718_final.pdf#page=59)
%======================================================================

function [ stressField,stressInteressante ] = computeStressDelta( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                          elNodes, numNodes, connectInfo, D, U,deltaT,alpha)
%% Symbolic variables
% tic
syms csi eta

%% Initiating variables

connectInfo(:,1)=[];    %Cleaning element labels

dimen=dOf*elNodes;      %Setting global dimensions

dNCsiEta=zeros(2,elNodes);  %Derivative matrix for jacobian computation
dNCsiEta=sym(dNCsiEta);

dNCsiEta=[diff(nodeShapeFun,csi).';diff(nodeShapeFun,eta).'];

%TODO-verify condition for other types of analysis
B=zeros(3,dimen);   %Shape function derivative matrix

% S=zeros(numElements,elNodes,3);

l=1:dOf:dimen;    %Index for derivative placement
m=0:dOf:dimen;    %Index for derivative placement
m(1)=[];

%Reordering integration poinst counterclockwise (needs revision for
%higher-order elements

intPts=intPts([3 4 2 1],:);

indexB=getGlobalDoF(connectInfo,dOf);
beta=D{1}*alpha;
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
                
        e=B*U(indexB(i,:));         %Computing the strains
        
        nodeShapeFunINT=subs(nodeShapeFun,[csi,eta],[intPts(j,1),intPts(j,2)]);
        nodeShapeFunINT=double(nodeShapeFunINT);
        nodeShapeFunTemp=deltaT*nodeShapeFunINT;
        
        if size(D,2)==2
            S(i,j,:)=(D{1}+D{2})*e-nodeShapeFunTemp*beta;     %Computing the stresses
        else
            S(i,j,:)=D{1}*e-nodeShapeFunTemp*beta; 
        end
    end   

end

%Swapping dimensions
S=permute(S,[2 1 3]);

%Extrapolating stresses to the element's corner nodes
for i=1:size(S,3)
    
   S(:,:,i)=[1+0.5*3^0.5 -0.5 1-0.5*3^0.5 -0.5;
            -0.5 1+0.5*3^0.5 -0.5 1-0.5*3^0.5;
            1-0.5*3^0.5 -0.5 1+0.5*3^0.5 -0.5;
            -0.5 1-0.5*3^0.5 -0.5 1+0.5*3^0.5]*S(:,:,i);

end

Sx=S(:,:,1)';
Sy=S(:,:,2)';
Sxy=S(:,:,3)';
Svm=(Sx.^2+Sy.^2-Sx.*Sy+3.*Sxy.^2).^0.5;

meanStress=zeros(numNodes,4);
%Stress smoothing at node level
for i=1:numNodes
    
    mask=ismember(connectInfo,i);
    meanStress(i,:)=mean([Sx(mask) Sy(mask) Sxy(mask) Svm(mask)],1);

end

%% Output

stressField=meanStress;
% toc
end

%======================================================================
%>@file computeStress.m
%>@brief Algorithm to perform stress recovery.
%>@details
%>
%>@author Rúben Lourenço
%>@date 2-Nov-2018
%======================================================================