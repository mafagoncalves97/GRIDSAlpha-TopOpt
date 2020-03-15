%======================================================================
%>@brief Defines the vector of global degrees of freedom.
%>
%>@param connectInfo (@b matrix) Matrix with connectivity information
%>@param dOf (@b integer) Degrees of freedom of the problem
%>
%>@retval indexB (@b matrix) Matrix of global degrees of freedom
%>
%>@details
%>Based on the connectivity matrix, calculates the global degrees 
%>of freedom associated with each element. 
%>
%>@note
%> <strong>`indexB`</strong> has dimensions:
%> - (nElements x (nNodes*dOf)),
%>@note
%>where:
%> - <strong>nElements</strong>: total number of elements
%> - <strong>elNodes</strong>: number of nodes per element
%======================================================================
 
function [ indexB ] = getGlobalDoF( connectInfo, dOf )

globalDoF=zeros(size(connectInfo,1),size(connectInfo,2)*dOf);

for i=1:size(connectInfo,1)
   
    idx=1;
    for x=1:length(connectInfo(i,:))
        for y=1:dOf                
              %Index vector of global degrees of freedom             
              globalDoF(i,idx)=(connectInfo(i,x)-1)*dOf+y;
              idx=idx+1;
        end        
    end   
    
end

indexB=globalDoF;

end

%======================================================================
%>@file getGlobalDoF.m
%>@brief Algorithm to define the vector of global degrees of freedom
%>@details
%>
%>@author Rúben Lourenço
%>@date 28-Oct-2018
%======================================================================
