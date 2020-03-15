 function [dcNew]=newDc(numElements,nodeCoord,connectInfo,rmin,x,dc)
%Matriz das coordenadas dos centros dos elementos(1coluna - x 2coluna -y)
elCoord=[];
for i=1:1:numElements 

        centro= (nodeCoord(connectInfo(i,3),1)-nodeCoord(connectInfo(i,2),1))/2;
        elCoord=[elCoord ;nodeCoord(connectInfo(i,2),2)+centro,nodeCoord(connectInfo(i,2),3)+centro];
end

%Matriz com as distâncias entre os centros dos elementos
for i=1:numElements
    for j=1:numElements 
        dist(i,j)= sqrt((elCoord(i,1)-elCoord(j,1))^2+(elCoord(i,2)-elCoord(j,2))^2);
        fac(i,j) = rmin-dist(i,j);
    end
end 

%Apenas interessam os elementos que se encontram dentro da circunferência
%de raio rmin 
 for i=1:numElements
   for j=1:numElements
       if i==j || fac(i,j)<0
           fac(i,j)=0;
       end
   end
 end
 
%Cálculo das novas sensibilidades tendo em conta apenas os elementos no
%interior da área de raio rmin
sumf=sparse(numElements,1);
dcn=sparse(numElements,1);
for i=1:numElements
    for j=1:numElements
   
            sumf(i)= sumf(i)+ fac(i,j);
            dcn(i)=dcn(i)+fac(i,j)*x(j)*dc(j);
        
    end
    dcNew(i)=1/(x(i)*sumf(i))*dcn(i);
end
 end
 