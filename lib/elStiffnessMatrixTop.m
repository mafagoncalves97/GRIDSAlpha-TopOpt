function [elementK ] = elStiffnessMatrixTop( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                   elNodes, connectInfo, D,elThick)
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
elementK=[];
elK=sparse(dOf*elNodes,dOf*elNodes);
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
        
      
        z=[1:4,7,8,5,6];
        elK(z,z) = elK(z,z)+ B'*D{1}*B*wCsi*wEta*det(J)*elThick;
      
        
    end
    elementK(:,:,i)=elK;
    elK=sparse(dOf*elNodes,dOf*elNodes);

end

end