function TopOptMecDelta(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,D,activeDoF,nodeShapeFun,intPts,F)
%% Graus de liberdade por nó
dOf=2;
FthConst=1;
%% Definição das variáveis 
volfrac=0.5;
xnew(1:size(connectInfo,1))=volfrac;
penal=3;
rmin=1.5*(nodeCoord(2,2)-nodeCoord(1,2));
%% Forças térmicas
if FthConst==1
    %% Forças térmicas
    deltaT=0*ones(1,4);
    alpha=[11*10^(-6);11*10^(-6);0];
    %Creates the global thermal load vector
    [ elF ] = assemGlobalF( nodeShapeFun, intPts, elNodeCoord, dOf,elNodes, numNodes, connectInfo, D,alpha,deltaT, t);
    %Creates the global vector
    f=F+elF;
end
%% Optimização
loop=0;
change=1;
[elementK ] = elStiffnessMatrixTop( nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                 elNodes, connectInfo, D,t);
c0m=0;  
while change>0.001
loop=loop+1;
x=xnew;
%% Matriz rigidez global e Matrizes elementares (Multidimensional array)
[K ] = StiffnessMatrixTop(dOf,numNodes, connectInfo,x,penal,elementK);

if FthConst==0
    %% Forças térmicas
    deltaT=1*ones(1,4);
    %Creates the global thermal load vector
    [ elF ] = assemGlobalF( nodeShapeFun, intPts, elNodeCoord, dOf,elNodes, numNodes, connectInfo, D,alpha,deltaT, t);
    %Creates the global vector
    f=F+elF;
end
%% Sistema global 
U=K(activeDoF,activeDoF)\f(activeDoF);
U1=zeros(dOf*numNodes,1);
U1(activeDoF)=U;
%% Deslocamentos elementares
[elU]=elDisplacements(U1,dOf,connectInfo,numElements);
%% Função objetivo 
dc=sparse(1,numElements);
if loop==1
for i=1:numElements
    c0m=c0m+((x(i)^penal)*elU(:,i)'*elementK(:,:,i)*elU(:,i));
end
end
c=0;
for i=1:1:numElements
    c=c+((x(i)^penal)*elU(:,i)'*elementK(:,:,i)*elU(:,i))/c0m;
    if FthConst==0
        dc(i)=-penal*(x(i)^(penal-1))*elU(:,i)'*elementK(:,:,i)*elU(:,i)+2*penal*(x(i)^(penal-1))*elU(:,i)'*f(i);
    elseif FthConst==1
        dc(i)=-penal*(x(i)^(penal-1))*elU(:,i)'*elementK(:,:,i)*elU(:,i);
    end
end 
%% Filtro de sensibilidades (Novas sensibilidades)
[dcNew]=newDc(numElements,nodeCoord,connectInfo,rmin,x,dc);
% dcNew=dc;  %Sem filtro
%% Critério de otimização
l1 = 0; move = 0.2; l2=100000;
while (l2-l1 > 1e-4)
  lmid = 0.5*(l2+l1);
  xnew = max(0.001,max(x-move,min(1.,min(x+move,real(x.*sqrt(-dcNew./lmid))))));
  if (sum(xnew))-volfrac*numElements > 0
     l1 = lmid;
  else
     l2 = lmid;
  end
end
%%
 change = max(abs(xnew-x));

%% PLOT DENSITIES  
matcol=nodeCoord(nodeCoord(:,3)~=nodeCoord(1,3));
mat = vec2mat(xnew,matcol(1)-2);

colormap(gray); imagesc(-mat);axis equal; axis tight; axis off; pause(1e-6);
disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(full(xnew))/(numElements)) ...
       ' ch.: ' sprintf('%6.3f',full(change) )])

   
if loop > 100
    change=0;
end
FO(loop)=c;
end
end