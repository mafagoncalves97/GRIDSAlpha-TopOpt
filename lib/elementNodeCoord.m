%======================================================================
%>@brief Organizes element node coordinates based on connectivity
%>information
%>
%>@param nodeCoord (@b matrix) Nodal coordinates matrix
%>@param connectInfo (@b matrix) Connectivity matrix
%>@param elNodes (@b int) Number of nodes per element
%>@param numElements (@b int) Total number of elements
%>
%>@retval elNodeCoord (@b matrix) Multidimensional array with nodal coordinates
%>organized on a per element basis
%>
%>@details Organizes the node coordinates on a per element basis. For
%>each element, a matrix containing the respective set of node coordinates, 
%>is generated, where the first column is the x-coordinate and the
%>second column is the y-coordinate. This output is stacked in a
%>multi-dimensional array <strong>`elNodeCoord`</strong>.
%>
%>@note
%> <strong>`elNodeCoord`</strong> has dimensions:
%> - (<strong>`elNodes`</strong> x 2 x <strong>`numElements`</strong>) for
%>2D problems;
%> - (<strong>`elNodes`</strong> x 3 x <strong>`numElements`</strong>) for
%>3D problems.
%>
%>@see
%> MATLAB documentation<br>
%> [Multidimensional arrays](https://www.mathworks.com/help/matlab/math/multidimensional-arrays.html)
%======================================================================
function [elNodeCoord]=elementNodeCoord(nodeCoord,connectInfo,elNodes,numElements)
    %Deleting node and element labels
    nodeCoord(:,1)=[];
    connectInfo(:,1)=[];
    
    %Initializing matrices
    element=zeros(elNodes,size(nodeCoord,2),numElements);
    nodalCoord=zeros(elNodes,size(nodeCoord,2));

    for i=1:numElements

        for j=1:size(connectInfo,2)
            %Extracts connectivity info
            node=connectInfo(i,j);
            %Extracts nodal coordinates based on connectivity info
            nodalCoord(j,:)=nodeCoord(node,:);

        end
        %Stacks the individual matrices into a multidimensional array
        element(:,:,i)=nodalCoord;

    end

    elNodeCoord=element;

end
%======================================================================
%>@file elementNodeCoord.m
%>@brief Functions to organize mesh data
%>@details
%>
%>@author Rúben Lourenço
%>@date 15-Oct-2018
%======================================================================