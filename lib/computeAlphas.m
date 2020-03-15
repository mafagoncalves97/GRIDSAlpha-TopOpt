function alphas = computeAlphas( dOf, connectInfo, U, choiceElement, kAlpha, kAlphaD )
%% Initiating variables

connectInfo(:,1)=[];    %Cleaning element labels

switch choiceElement
    
    case 'PE/FI/IM - Q4/6I'
        
        genDisp=zeros(size(connectInfo,1),6);
        
    case 'PE/FI/IM - Q4/4I'
        
        genDisp=zeros(size(connectInfo,1),4);
    
    case 'PE/FI/IM - Q6'
        
        genDisp=zeros(size(connectInfo,1),4);
    
    case 'PE/FI/IM - QM6'
    
        genDisp=zeros(size(connectInfo,1),4);
    
    case 'PE/FI/CM - Qi5'
        
        genDisp=zeros(size(connectInfo,1),2);
        
    case 'PE/FI/CM - Qi6'
    
        genDisp=zeros(size(connectInfo,1),4);
        
end

indexB=getGlobalDoF(connectInfo,dOf);

for i=1:size(connectInfo,1)

      genDisp(i,:)=(-inv(kAlpha(:,:,i))*kAlphaD(:,:,i))*U(indexB(i,:));

end

%% Output

alphas=genDisp;
% toc
end



