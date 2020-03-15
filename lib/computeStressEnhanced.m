function [ stressField ] = computeStressEnhanced( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                  elNodes, numNodes, connectInfo, D, U, alphas, choiceElement)
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

l=1:dOf:dimen;    %Index for derivative placement
m=0:dOf:dimen;    %Index for derivative placement
m(1)=[];

S=zeros(size(connectInfo,1),size(intPts,1),3);

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
        
        switch choiceElement
        
            case 'PE/FI/IM - Q4/6I'
                
                dNAlfa=zeros(2,4);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-csi*eta*(eta-1) 0.5*(2*csi+1)*(1-eta^2) -csi*eta*(eta+1) 0.5*(2*csi-1)*(1-eta^2);
                        0.5*(2*eta-1)*(1-csi^2) -csi*eta*(csi+1) 0.5*(2*eta+1)*(1-csi^2) -csi*eta*(csi-1)];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates
                
                Balfa=[dNAlfaDxDY(1,1) 0 dNAlfaDxDY(1,2)+dNAlfaDxDY(1,4) 0 dNAlfaDxDY(1,3) 0;...
                       0 dNAlfaDxDY(2,1)+dNAlfaDxDY(2,3) 0 dNAlfaDxDY(2,2) 0 dNAlfaDxDY(2,4);...
                       dNAlfaDxDY(2,1)+dNAlfaDxDY(2,3) dNAlfaDxDY(1,1) dNAlfaDxDY(2,2) dNAlfaDxDY(1,2)+dNAlfaDxDY(1,4) dNAlfaDxDY(2,4) dNAlfaDxDY(1,3)];
                   
            case 'PE/FI/IM - Q4/4I'
                
                dNAlfa=zeros(2,4);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-csi*eta*(eta-1) 0.5*(2*csi+1)*(1-eta^2) -csi*eta*(eta+1) 0.5*(2*csi-1)*(1-eta^2);...
                              0.5*(2*eta-1)*(1-csi^2) -csi*eta*(csi+1) 0.5*(2*eta+1)*(1-csi^2) -csi*eta*(csi-1)];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                          
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates
                
                Balfa=[dNAlfaDxDY(1,1)+dNAlfaDxDY(1,3) 0 dNAlfaDxDY(1,2)+dNAlfaDxDY(1,4) 0;...
                       0 dNAlfaDxDY(2,1)+dNAlfaDxDY(2,3) 0 dNAlfaDxDY(2,2)+dNAlfaDxDY(2,4);...
                       dNAlfaDxDY(2,1)+dNAlfaDxDY(2,3) dNAlfaDxDY(1,1)+dNAlfaDxDY(1,3) dNAlfaDxDY(2,2)+dNAlfaDxDY(2,4) dNAlfaDxDY(1,2)+dNAlfaDxDY(1,4)];
            
            case 'PE/FI/IM - Q6'
                
                dNAlfa=zeros(2,2);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-2*csi 0;
                        0 -2*eta];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                          
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates
                
                Balfa=[dNAlfaDxDY(1,1) 0 dNAlfaDxDY(1,2) 0;
                       0 dNAlfaDxDY(2,1) 0 dNAlfaDxDY(2,2);
                       dNAlfaDxDY(2,1) dNAlfaDxDY(1,1) dNAlfaDxDY(2,2) dNAlfaDxDY(1,2)];
                   
             case 'PE/FI/IM - QM6'
                
                dNAlfa=zeros(2,2);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-2*csi 0;
                        0 -2*eta];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                          
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates           
                   
                DN0=subs(dN,[csi,eta],[0,0]);
                DN0=double(DN0);
                J0=DN0*elNodeCoord(:,:,i);
                             
                F0=[J0(1,1)^2 J0(2,1)*J0(1,2) 2*J0(1,1)*J0(1,2);
                    J0(2,1)*J0(1,2) J0(2,2)^2 2*J0(2,1)*J0(2,2);
                    2*J0(1,1)*J0(2,1) 2*J0(1,2)*J0(2,2) (J0(1,1)*J0(2,2)+J0(1,2)*J0(2,1))];
                
                BalfaI=[dNAlfaDxDY(1,1) 0 dNAlfaDxDY(1,2) 0;
                       0 dNAlfaDxDY(2,1) 0 dNAlfaDxDY(2,2);
                       dNAlfaDxDY(2,1) dNAlfaDxDY(1,1) dNAlfaDxDY(2,2) dNAlfaDxDY(1,2)];
                Balfa=(det(J0)/det(J))*inv(F0)'*BalfaI;
                   
            case 'PE/FI/CM - Qi5'
                
                dNAlfa=zeros(2,1);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-2*csi*(1-eta^2);-2*eta*(1-csi^2)];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                          
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates
                
                Balfa=[dNAlfaDxDY(1) 0;
                       0 dNAlfaDxDY(2);
                       dNAlfaDxDY(2) dNAlfaDxDY(1)];
                   
            case 'PE/FI/CM - Qi6'
                
                dNAlfa=zeros(2,1);
                dNAlfa=sym(dNAlfa);
                
                dNAlfa=[-2*csi*(1-eta^2);-2*eta*(1-csi^2)];
                          
                dN2=subs(dNAlfa,[csi,eta],[intPts(j,1),intPts(j,2)]); %Calculating derivatives in the integration point
                dN2=double(dN2);
                          
                dNAlfaDxDY=J\dN2; %Transforming to global coordinates
                
                Balfa=[dNAlfaDxDY(1) 0 0 0;...
                       0 dNAlfaDxDY(2) 0 0;...
                       0 0 dNAlfaDxDY(2) dNAlfaDxDY(1)];
                
        
        end
        
        e=B*U(indexB(i,:))+Balfa*alphas(i,:)';
        S(i,j,:)=D{1}*e;

    end

end

S=permute(S,[2 1 3]);

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

for i=1:numNodes
    
    mask=ismember(connectInfo,i);
    meanStress(i,:)=mean([Sx(mask) Sy(mask) Sxy(mask) Svm(mask)],1);

end

%% Output
stressField=meanStress;
end

