%======================================================================
%>@brief Creates the global nodal load vector.
%>
%>@param loads (@b struct) Structure array containing the applied loads
%>@param numNodes (@b integer) Total number of nodes
%>@param dOf (@b integer) Degrees of freedom of the problem
%>
%>@retval nodalForces (@b matrix) Column vector of global loads
%>
%>@details
%>Cycles through the structure array <strong>`loads`</strong>, which
%>contains information about the the magnitude, direction and point of
%>application of all the loads acting on the model. The structure array
%> <strong>`loads`</strong> has the following format:
%>\htmlonly <style>div.image img[src="struct3.png"]{width:20%;height=20%}</style> \endhtmlonly 
%>@image html struct3.png
%>@image latex struct3.pdf width=4cm
%> <br>
%>
%======================================================================
 
function [ nodalForces ] = setNodalForces( loads, numNodes, dOf)

f=zeros(numNodes*dOf,1);
nodeOrder=1:numNodes;

for k=1:dOf
    
    switch k
        
        case 1
            
            for i=1:numel(loads.Xdirection)
                node=nodeOrder(loads.Xdirection(i).nodes);
                index=(node-1)*dOf+k;
                f(index)=f(index)+loads.Xdirection(i).value;
    
            end
            
            
        case 2
            
            for i=1:numel(loads.Ydirection)
                node=nodeOrder(loads.Ydirection(i).nodes);
                index=(node-1)*dOf+k;
                f(index)=f(index)+loads.Ydirection(i).value;
    
            end
    
    end


end

nodalForces=f;

end

%======================================================================
%>@file setNodalForces.m
%>@brief Functions to define the global load vector
%>@details
%>
%>@author Rúben Lourenço
%>@date 28-Oct-2018
%======================================================================
