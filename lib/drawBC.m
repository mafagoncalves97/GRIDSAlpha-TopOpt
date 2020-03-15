%======================================================================
%>@brief Draws the appropriate boundary condition symbol for a given set of
%>selected nodes.
%>
%>@param vertexData (@b matrix) Matrix containing the selected nodes'
%coordinates
%>@param bcType (@b string) String indicating the BC type
%>
%>@retval none
%>
%>@details
%>For a set previously selected nodes, draws the appropriate boundary
%>condition symbols at the corresponding nodal coordinates.
%>@note
%>The boundary condition symbols have the following meanings: 
%>\htmlonly <style>div.image img[src="fixed.png"]{width:25%;height=25%}</style> \endhtmlonly
%>\htmlonly <style>div.image img[src="ss1.png"]{width:25%;height=25%}</style> \endhtmlonly
%>\htmlonly <style>div.image img[src="ss2.png"]{width:25%;height=25%}</style> \endhtmlonly
%>\htmlonly <style>table{table-layout:fixed;}</style> \endhtmlonly
%> <br><br>
%> <table border="0">
%> <tr align="center">
%> <td>@image html fixed.png</td>
%> <td>@image html ss1.png</td>
%> <td>@image html ss2.png</td>
%> </tr><tr align="center"><td><b>Fixed node</b></td><td><b>Simply 
%>supported, free along x-axis<b></td><td><b>Simply supported, free along 
%>y-axis</b></td></tr></table>
%>
%======================================================================
function drawBC( vertexData, bcType )

v=[vertexData(1,:)' vertexData(2,:)'];

switch bcType
   
    case 'fixedNodes'
        
        col=[1 0 0; 1 0 0];
        
    case 'fixedX'
        
        col=[0 1 0; 1 0 0];
        
    case 'fixedY'
        
        col=[1 0 0; 0 1 0];
    
end

for i=1:size(v,1)

drawPatch(v(i,:),col,i);

end

end

%======================================================================
%>@brief Draws a patch object for a given set of vertices.
%>
%>@param v (@b matrix) Matrix containing the vertices' coordinates
%>@param col (@b matrix) Matrix of RGB color triplets defining the color 
%>for each face 
%>@param i (@b integer) Number to be used as a tag (can be ignored)
%>@retval none
%>
%>@details
%>Draws a patch object with two triangluar shapes, based on the set of
%>corresponding vertex coordinates. Uses a matrix of RGB color triplets to
%>define the color of each face.
%======================================================================
function drawPatch(v,col,i)
    x=v(1);
    y=v(2);
    
    p=findobj(findobj(gcf,'Type','line'),'Tag','nodes');
    maxX=max(get(p,'XData'));
    maxY=max(get(p,'YData'));
    c=(maxX^2+maxY^2)^0.5;
    
    vert=[x y; x-c*0.0045 y+c*0.009; x+c*0.0045 y+c*0.009; x+c*0.009 y-c*0.0045; x+c*0.009 y+c*0.0045];
    patch(gca,'Faces',[1 2 3; 1 4 5],'Vertices',vert,'FaceVertexCData',col,'FaceColor','flat','Tag',num2str(i),'FaceAlpha',0.65);

end

%======================================================================
%>@file drawBC.m
%>@brief Functions to draw boundary condition symbols
%>@details
%>
%>@author Rúben Lourenço
%>@date 28-Oct-2018
%======================================================================