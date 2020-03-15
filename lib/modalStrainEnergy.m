function [ strainEnergy ] = modalStrainEnergy( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                   elNodes, numNodes, numElements, connectInfo, D, U)
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

strainEnergy=zeros(size(intPts,1),size(U,2));

l=1:dOf:dimen;    %Index for derivative placement
m=0:dOf:dimen;    %Index for derivative placement
m(1)=[];

%Reordering integration poinst counterclockwise (needs revision for
%higher-order elements. 
if size(intPts,1)>1
    intPts=intPts([3 4 2 1],:);
end

indexB=getGlobalDoF(connectInfo,dOf);
for i=1:numElements 
   
    for k=1:size(U,2)
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

            e=B*U(indexB(i,:),k);         %Computing the strains
            strainEnergy(j,k)=0.5*e'*D*e;

        end   
    end
end

end
