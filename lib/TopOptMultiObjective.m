function TopOptMultiObjective(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,D,k,activeDoF,nodeShapeFun,intPts,f)
%% Definição das variáveis 
volfrac=0.4;
xnew(1:size(connectInfo,1))=volfrac;
penal=3;
rmin=1.25*(nodeCoord(2,2)-nodeCoord(1,2));
%% Optimização
loop=0;
change=1;
c0m=0; c0t=0;
wt=1;
dOfm=2;
dOft=1;
C=[k 0;0 k];
%% Matrizes elementares
[elementKth]=elStiffnessTopOptThermal(nodeShapeFun, intPts, elNodeCoord, dOft,...
                                                 elNodes,connectInfo, C,t);

[elementKm ] = elStiffnessMatrixTop( nodeShapeFun, intPts, elNodeCoord, dOfm,...
                                                 elNodes, connectInfo, D,t);
                                             
while change>0.001
loop=loop+1;
x=xnew;
%% Análise mecânica
    %%Matriz rigidez global e Matrizes elementares (Multidimensional array)
    [Km] = StiffnessMatrixTop(dOfm,numNodes, connectInfo,x,penal,elementKm);
    %%Sistema global 
    U=Km(activeDoF,activeDoF)\f(activeDoF);
    U1=zeros(dOfm*numNodes,1);
    U1(activeDoF)=U;

    %%Deslocamentos elementares
    [elU]=elDisplacements(U1,dOfm,connectInfo,numElements);

%% Análise térmica
    %%BC-termico
    nodes=nodeCoord(:,1)';
%     fixedNodesth=[211];
    fixedNodesth=[3241];
    nodesFlux=nodes;
    nodesFlux(fixedNodesth)=[];
    q=zeros(numNodes,1);
    q(nodesFlux)=0.5;

    activeDoFth=nodes;
    activeDoFth(fixedNodesth)=[];
    tP=zeros(1,9);

    %%Sistema global 
    [Kth] = ThermalMatrixTop(dOft,elNodes, numNodes, connectInfo,x,penal,tP,elementKth);
    T=Kth(activeDoFth,activeDoFth)\q(activeDoFth);
    T1=zeros(dOft*numNodes,1);
    T1(activeDoFth)=T;
    [elT]=elTemperatures(T1,dOft,connectInfo,numElements);
%% Função objetivo 
dc=sparse(1,numElements);
c=0;cm=0;ct=0;CCm=0;CCt=0;
if loop==1
    for i=1:numElements
        c0m=c0m+(x(i)^penal)*elU(:,i)'*elementKm(:,:,i)*elU(:,i);
        c0t=c0t+((x(i)^penal).*elT(:,i)'*elementKth(:,:,i)*elT(:,i));
    end
end
for i=1:1:numElements
    cm=cm+(1-wt)*((x(i)^penal)*elU(:,i)'*elementKm(:,:,i)*elU(:,i))/(c0m);
    ct=ct+wt*(((x(i)^penal).*elT(:,i)'*elementKth(:,:,i)*elT(:,i))/c0t);
    c=ct+cm;
    CCm=CCm+((x(i)^penal)*elU(:,i)'*elementKm(:,:,i)*elU(:,i))/(c0m);
    CCt=CCt+(((x(i)^penal).*elT(:,i)'*elementKth(:,:,i)*elT(:,i))/c0t);
    dc(i)=((1-wt)*(-penal*(x(i)^(penal-1))*elU(:,i)'*elementKm(:,:,i)*elU(:,i)))/c0m+(wt*(-penal*(x(i)^(penal-1))*elT(:,i)'*elementKth(:,:,i)*elT(:,i)))/c0t;
end 
%% Filtro de sensibilidades (Novas sensibilidades)
[dcNew]=newDc(numElements,nodeCoord,connectInfo,rmin,x,dc);
%% Critério de otimização
l1 = 0; l2 = 100000; move = 0.2;
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
matcol=nodeCoord(find(nodeCoord(:,3)~=nodeCoord(1,3)));
mat = vec2mat(xnew,matcol(1)-2);

colormap(gray); imagesc(-mat);axis equal; axis tight; axis off; pause(1e-6);
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c)  ' ObjM.: ' sprintf('%10.4f',CCm) ' ObjT.: ' sprintf('%10.4f',CCt) ...
       ' Vol.: ' sprintf('%6.3f',sum(full(xnew))/(numElements))])
   
if loop > 100
    change=0;
end
FO(loop)=c;
end 
end