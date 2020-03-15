function TopOptThermal(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,k,nodeShapeFun,intPts)

%% Graus de liberdade por nó
dOf=1;
%% Definição das variáveis 
volfrac=0.4;
xnew(1:size(connectInfo,1))=volfrac;
penal=3;
rmin=1.25*(nodeCoord(2,2)-nodeCoord(1,2));
c0t=0;
wt=0.2;
%% Optimização
loop=0;
change=1;
%%Condições-fronteira
ymax=max(nodeCoord(:,3));
xmax=max(nodeCoord(:,2));
ymin=min(nodeCoord(:,3));
xmin=min(nodeCoord(:,2));
nEly=nodeCoord(:,3)==max(nodeCoord(:,3));
nElx=nodeCoord(:,2)==max(nodeCoord(:,2));

%%BC-MO
% fixedNodesth=find(nodeCoord(:,2)==xmax&nodeCoord(:,3)==ymin);
% sSnodesth.BC=find(nodeCoord(:,2)==xmin|nodeCoord(:,3)==ymax);

%%BC-termico

nodes=nodeCoord(:,1)';
fixedNodesth=[3241];
%fixedNodesth=[211];
nodesFlux=nodes;
nodesFlux(fixedNodesth)=[];
q=zeros(numNodes,1);
q(nodesFlux)=0.1;

activeDoF=nodes;
activeDoF(fixedNodesth)=[];
%% Input material properties
C=[k 0;0 k];
tP=zeros(1,9);
%% Temperatura prescrita
    [elementK]=elStiffnessTopOptThermal(nodeShapeFun, intPts, elNodeCoord, dOf,...
                                                   elNodes,connectInfo, C,t);
while change>0.001
    loop=loop+1;
    x=xnew;
    %% Matriz rigidez global e Matrizes elementares (Multidimensional array)
    [K] = ThermalMatrixTop(dOf,elNodes, numNodes, connectInfo,x,penal,tP,elementK);
    
    %% Sistema global 
    T=K(activeDoF,activeDoF)\q(activeDoF);
    T1=zeros(dOf*numNodes,1);
    T1(activeDoF)=T;

    %% Deslocamentos elementares
    [elT]=elTemperatures(T1,dOf,connectInfo,numElements);

    %% Função objetivo 
    dc=sparse(1,numElements);
    c=0;
    if loop==1
    for i=1:numElements
        c0t=c0t+((x(i)^penal).*elT(:,i)'*elementK(:,:,i)*elT(:,i));
    end
    end
    for i=1:1:numElements
        c=c+(((x(i)^penal)*elT(:,i)'*elementK(:,:,i)*elT(:,i))/c0t)*wt;
        dc(i)=-penal*(x(i)^(penal-1))*elT(:,i)'*elementK(:,:,i)*elT(:,i);
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
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(full(xnew))/(numElements)) ...
       ' ch.: ' sprintf('%6.3f',full(change) )])

        if loop > 100
            change=0;
        end
        FO(loop)=c;
end 
end