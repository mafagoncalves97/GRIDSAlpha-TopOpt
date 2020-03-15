function [ elStiffness ] = assemGlobalStiffSRI( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                elNodes, numNodes, connectInfo, D, t)
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
B2=zeros(3,dimen);

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
        
        K(indexB(i,:),indexB(i,:))=K(indexB(i,:),indexB(i,:))+B'*D{2}*B*wCsi*wEta*det(J)*t;

    end 
    
    dN2=subs(dN,[csi,eta],[0,0]);
    dN2=double(dN2);
    J2=dN2*elNodeCoord(:,:,i);
    dNdXdY2=J2\dN2;
    
    B2(1,l)=dNdXdY2(1,:);    %Derivative matrix on global coordinates
    B2(2,m)=dNdXdY2(2,:);
            
    B2(3,l)=dNdXdY2(2,:);
    B2(3,m)=dNdXdY2(1,:);        
        
    wCsi2=2;
    wEta2=2;
    
    K(indexB(i,:),indexB(i,:))=K(indexB(i,:),indexB(i,:))+B2'*D{1}*B2*wCsi2*wEta2*det(J2)*t;  

end

%% Output

elStiffness=K;

end

