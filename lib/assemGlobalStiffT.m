function [ elStiffness ] = assemGlobalStiffT( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                             elNodes, numNodes, connectInfo, C, elThick,tP)
%% Symbolic variables
syms csi eta

%% Initiating variables

connectInfo(:,1)=[];    %Cleaning element labels

dimen=dOf*elNodes;      %Setting global dimensions

dNCsiEta=zeros(2,elNodes);  %Derivative matrix for jacobian computation
dNCsiEta=sym(dNCsiEta);

dNCsiEta=[diff(nodeShapeFun,csi).';diff(nodeShapeFun,eta).'];

%TODO-verify condition for other types of analysis
B=zeros(2,dimen);   %Shape function derivative matrix

K=zeros(dOf*numNodes,dOf*numNodes);  %Global stiffness matrix

l=1:dOf:dimen;    %Index for derivative placement


indexB=getGlobalDoF(connectInfo,dOf);

for i=1:size(connectInfo,1)

    for j=1:size(intPts,1)
        
        dN=dNCsiEta;
        dN=subs(dN,[csi,eta],[intPts(j,1),intPts(j,2)]);  %Calculating derivatives in the integration point
        dN=double(dN);                                    %Converting from symbolic to double
        J=dN*elNodeCoord(:,:,i);                          %Jacobian
        
        dNdXdY=J\dN;           %Transforming to global coordinates
        
        B(1,l)=dNdXdY(1,:);    %Derivative matrix on global coordinates
        B(2,l)=dNdXdY(2,:);       
        
        wCsi=intPts(j,3);
        wEta=intPts(j,4);
        
      
            %Direct assembly of the element stiffness matrix into global stiffness matrix
            K(indexB(i,:),indexB(i,:))=K(indexB(i,:),indexB(i,:))+B'*C*B*wCsi*wEta*det(J)*elThick;
            
            
    end   
            
end

ind=find(tP==1);
for i=1:size(ind,2)
    K(ind(i),ind(i))=K(ind(i),ind(i))+10^6;
end
%% Output

elStiffness=K;

end
