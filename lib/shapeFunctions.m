%======================================================================
%>@brief Provides the nodal shape functions for quadrilateral and
%>triangular elements
%>
%>@param elementType (@b string) The finite element type
%>@param nNodes (@b int) The number of nodes of the finite element
%>
%>@retval shapeFun (@b sym @b matrix) Nodal shape functions of the finite 
%>element
%>
%>@details
%>Provides the algorithm to generate the nodal shape functions of a finite
%>element, based on its type and number of nodes.
%>@note
%> <strong>`shapeFun`</strong> is written in symbolic form and has
%>dimensions:
%> - (<strong>`nNodes`</strong> x 1).
%>
%>@todo
%>Add nodal shape functions generation for 3D elements: tetrahedrons and
%>hexahedrons
%>@todo
%>Check correctness of solutions for non-linear 2D finite elements
%>@see
%> Theory manual<br>
%> [Chapter 4, Subsection 4.4.1 - Shape functions generation](http://grids.web.ua.pt/wp-content/uploads/2018/09/DISSERTACAO_Ruben_Lourenco_MIEM_1718_final.pdf#page=55)
%> <br>MATLAB documentation<br>
%> [Symbolic variables](https://www.mathworks.com/help/symbolic/create-symbolic-numbers-variables-and-expressions.html)
%====================================================================== 
function [shapeFun]=shapeFunctions(elementType,nNodes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function shapeFunctions provides the nodal shape functions for      %
%   quadrilateral and triangular elements                               %                  %
%                                                                       %
%   Input: elementType - the element type ('quad' or 'tri')             %
%          nNodes - the number of nodes of the element                  %
%                                                                       %
%   Output: points - (N*1) matrix:                                      %
%                  1st column gives nodal shape function                %
%                                                                       %
%                  N is the number of nodes of the element              %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %Initializing symbolic variables
    syms csi eta N L1 L2 L3
    
    %Initializing shape function matrix
    nodeFun=zeros(nNodes,1);
    nodeFun=sym(nodeFun);
       
    switch elementType
        
        %Triangular elements
        case 'tri'
           %Shape functions for 3-node triangle
           if nNodes==3
           
                nodeFun=[L1
                         L2
                         L3];
                     
           %Shape functions for 6-node triangle     
           elseif nNodes==6
               
               nodeFun=[L1*(2*L1-1)
                        4*L1*L2
                        L2*(2*L2-1)
                        4*L2*L3
                        L3*(2*L3-1)
                        4*L1*L3];
                    
           %Shape functions for 10-node triangle         
           elseif nNodes==10
               
               nodeFun=[0.5*L1*(3*L1-1)*(3*L1-2)
                        0.5*9*L2*L1*(3*L1-1)
                        0.5*9*L1*L2*(3*L2-1)
                        0.5*L2*(3*L2-1)*(3*L2-2)
                        0.5*9*L2*L3*(3*L2-1)
                        0.5*9*L2*L3*(3*L3-1)
                        0.5*L3*(3*L3-1)*(3*L3-2)
                        0.5*9*L1*L3*(3*L3-1)
                        0.5*9*L1*L3*(3*L1-1)
                        27*L1*L2*L3];
           else
               %Returns message when trying to analyze different elements
               msgbox('Predefined elements only have 3, 6 or 10 nodes...')
           
           end        
            
        %Quadrilateral elements
        case 'quad'            
           
           %Setting matrix dimensions based on node number
           if floor(nNodes^0.5)==nNodes^0.5
               dim=nNodes^0.5;
           elseif (nNodes^0.5-floor(nNodes^0.5))>0.5 || (nNodes^0.5-floor(nNodes^0.5))<0.5
               dim=floor(nNodes^0.5)+1;
           end
           
           %Initializing variables
           nCsiI=zeros(dim,1);
           nEtaI=zeros(dim,1);

           nCsiI=sym(nCsiI);
           nEtaI=sym(nEtaI);
       
           nodeCoord=linspace(-1,1,dim);
          
           %Constructing one-dimensional shape functions over csi and eta
           for i=1:length(nodeCoord)
              N=1;  
              csiI=setdiff(nodeCoord,nodeCoord(i));
               
              for j=1:length(csiI)
               
                   N=N*(csi-csiI(j))/(nodeCoord(i)-csiI(j));
                   
              end
              
              nCsiI(i,1)=N;
              nEtaI=subs(nCsiI,csi,eta);
              
           end
           
           %Constructing nodal bi-dimensional shape functions (csi,eta)
           k=1; 
 
           for i=1:length(nodeCoord)-1               
               nodeFun(k,1)=nCsiI(i,1)*nEtaI(1,1);
               k=k+1;
           end
           
           for i=1:length(nodeCoord)-1               
               nodeFun(k,1)=nCsiI(end,1)*nEtaI(i,1);
               k=k+1;
           end
           
           for i=length(nodeCoord):-1:2
               nodeFun(k,1)=nCsiI(i,1)*nEtaI(end,1);
               k=k+1;
           end
           
           for i=length(nodeCoord):-1:2
               nodeFun(k,1)=nCsiI(1,1)*nEtaI(i,1);
               k=k+1;
           end
           
           if nNodes==9
               
               %Central node
               nodeFun(k,1)=nCsiI(2,1)*nEtaI(2,1);
               
           elseif nNodes==16
               
               %Central nodes
               for i=2:length(nodeCoord)-1
                   nodeFun(k,1)=nCsiI(i,1)*nEtaI(2,1);
                   k=k+1;
               end
               
               for i=length(nodeCoord)-1:-1:2
                   nodeFun(k,1)=nCsiI(i,1)*nEtaI(3,1);
                   k=k+1;
               end
               
           end
    
    end
    
    %Setting function output
    shapeFun=nodeFun;
end
%======================================================================
%>@file shapeFunctions.m
%>@brief Functions to generate element nodal shape functions
%>@details
%>
%>@author Rúben Lourenço
%>@date 16-Oct-2018
%======================================================================