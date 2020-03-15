%======================================================================
%>@brief Defines the vector of active degrees of freedom.
%>
%>@param fixedNodes (@b matrix) Matrix identifying the fixed nodes
%>@param sSnodes (@b struct) Structure array containing information about
%>the simply supported nodes
%>@param numNodes (@b integer) Total number of nodes
%>@param dOf (@b integer) Degrees of freedom of the problem
%>
%>@retval activeDoF (@b matrix) Vector of active degrees of freedom
%>
%>@details
%>Given the information about the fixed and simply supported nodes, defines
%>a vector of labels indicating the active degrees of freedom. This is
%>later used as a mask to solve the reduced system of equations.
%>
%======================================================================
 
function [ activeDoF ] = setFreeNodes( fixedNodes, sSnodes, numNodes, dOf)

nodeOrder=1:numNodes;
dofs=[];
for k=1:dOf    
    nodes=nodeOrder(fixedNodes);
    dofs=[dofs,(nodes-1).*dOf+k];
end

for i=1:numel(sSnodes.BC)
    
    dir=sSnodes.BC(i).FixedDir;
    nodes=nodeOrder(sSnodes.BC(i).nodes);
    switch dir
        case 'x'
            k=1;            
            
        case 'y'
            k=2;            
    end
    dofs=[dofs,(nodes-1).*dOf+k];
end

globalDofs=1:dOf*numNodes;
activeDoF=setdiff(globalDofs,dofs);

end

%======================================================================
%>@file setFreeNodes.m
%>@brief Functions to define the active degrees of freedom
%>@details
%>
%>@author Rúben Lourenço
%>@date 28-Oct-2018
%======================================================================