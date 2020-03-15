%======================================================================
%>@brief Provides the (&xi;,&eta;) coordinates of the Gauss points and
%>respective weights for quadrilateral and triangular elements
%>
%>@param elementType (@b string) The element type
%>@param nNodes (@b int) The number of nodes per element
%>
%>@retval points (@b matrix) The coordinates of the Gaussian points and
%>associated weighting factors
%>
%>Based on the quadrature order and element type, provides the (&xi;,&eta;) 
%>coordinates and respective weighting factors in the appropriate format.
%>
%>@note
%> <strong>`points`</strong> has dimensions:
%> - (N<sup>2</sup> x 4).
%>
%>@note
%>Columns are organized such that:
%> - __First column__: &xi;-coordinate;
%> - __Second column__: &eta;-coordinate;
%> - __Third column__: &xi;-weight;
%> - __Fourth column__: &eta;-weight.
%>
%>@todo
%>Add Gauss points generation for 3D elements: tetrahedrons and
%>hexahedrons.
%>@todo
%>Include necessary code logic to deal with reduced integration schemes. 
%>If a reduced integration scheme is used, the integration points will need
%>to be calculated accordingly.
%>
%>@see
%>Theory manual<br>
%>[Chapter 2, Subsection 2.2.3 - Isoparametric four-node quadrilateral
%>element](http://grids.web.ua.pt/wp-content/uploads/2018/09/DISSERTACAO_Ruben_Lourenco_MIEM_1718_final.pdf#page=30)
%======================================================================
function [ points ] = intPoints( elementType, nNodes )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function intPoints provides the (csi, eta) coordinates of the       %
%   Gauss points and weights for quadrilateral and triangular        %
%   elements up to the 3rd order                                        %                  
%                                                                       %
%   Input: elementType - the element type ('quad' or 'tri')             %
%          nNodes - the number of nodes of the element                  %
%                                                                       %
%   Output: points - (N*N by 4) matrix:                                 %
%                  1st column gives the csi-coordinates                 %
%                  2nd column gives the eta-coordinates                 %
%                  3rd column gives the csi-weights                     %
%                  4th column gives the eta-weights                     %
%                  N is the order of quadrature                         %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch elementType
    
    %Triangular elements    
    case 'tri'
        
      %Order of quadrature
      quadOrder=round(nNodes/3+1);
      
      if (quadOrder == 1)
         points=[0.33333333333333    0.33333333333333    1.00000000000000];

      elseif (quadOrder == 2)
         points=[0.16666666666667    0.16666666666667    0.33333333333333
                 0.16666666666667    0.66666666666667    0.33333333333333
                 0.66666666666667    0.16666666666667    0.33333333333333];

      elseif (quadOrder == 3)
         points=[0.33333333333333    0.33333333333333   -0.56250000000000
                 0.20000000000000    0.20000000000000    0.52083333333333
                 0.20000000000000    0.60000000000000    0.52083333333333
                 0.60000000000000    0.20000000000000    0.52083333333333];

      elseif (quadOrder == 4)
         points=[0.44594849091597    0.44594849091597    0.22338158967801
                 0.44594849091597    0.10810301816807    0.22338158967801
                 0.10810301816807    0.44594849091597    0.22338158967801
                 0.09157621350977    0.09157621350977    0.10995174365532
                 0.09157621350977    0.81684757298046    0.10995174365532
                 0.81684757298046    0.09157621350977    0.10995174365532];

      elseif (quadOrder == 5)
         points=[0.33333333333333    0.33333333333333    0.22500000000000
                 0.47014206410511    0.47014206410511    0.13239415278851
                 0.47014206410511    0.05971587178977    0.13239415278851
                 0.05971587178977    0.47014206410511    0.13239415278851
                 0.10128650732346    0.10128650732346    0.12593918054483
                 0.10128650732346    0.79742698535309    0.12593918054483
                 0.79742698535309    0.10128650732346    0.12593918054483];

      elseif (quadOrder == 6)
         points=[0.24928674517091    0.24928674517091    0.11678627572638
                 0.24928674517091    0.50142650965818    0.11678627572638
                 0.50142650965818    0.24928674517091    0.11678627572638
                 0.06308901449150    0.06308901449150    0.05084490637021
                 0.06308901449150    0.87382197101700    0.05084490637021
                 0.87382197101700    0.06308901449150    0.05084490637021
                 0.31035245103378    0.63650249912140    0.08285107561837
                 0.63650249912140    0.05314504984482    0.08285107561837
                 0.05314504984482    0.31035245103378    0.08285107561837
                 0.63650249912140    0.31035245103378    0.08285107561837
                 0.31035245103378    0.05314504984482    0.08285107561837
                 0.05314504984482    0.63650249912140    0.08285107561837];

      elseif (quadOrder == 7)
         points=[0.33333333333333    0.33333333333333   -0.14957004446768
                 0.26034596607904    0.26034596607904    0.17561525743321
                 0.26034596607904    0.47930806784192    0.17561525743321
                 0.47930806784192    0.26034596607904    0.17561525743321
                 0.06513010290222    0.06513010290222    0.05334723560884
                 0.06513010290222    0.86973979419557    0.05334723560884
                 0.86973979419557    0.06513010290222    0.05334723560884
                 0.31286549600487    0.63844418856981    0.07711376089026
                 0.63844418856981    0.04869031542532    0.07711376089026
                 0.04869031542532    0.31286549600487    0.07711376089026
                 0.63844418856981    0.31286549600487    0.07711376089026
                 0.31286549600487    0.04869031542532    0.07711376089026
                 0.04869031542532    0.63844418856981    0.07711376089026];

      elseif (quadOrder == 8)
         points=[0.33333333333333    0.33333333333333    0.14431560767779
                 0.45929258829272    0.45929258829272    0.09509163426728
                 0.45929258829272    0.08141482341455    0.09509163426728
                 0.08141482341455    0.45929258829272    0.09509163426728
                 0.17056930775176    0.17056930775176    0.10321737053472
                 0.17056930775176    0.65886138449648    0.10321737053472
                 0.65886138449648    0.17056930775176    0.10321737053472
                 0.05054722831703    0.05054722831703    0.03245849762320
                 0.05054722831703    0.89890554336594    0.03245849762320
                 0.89890554336594    0.05054722831703    0.03245849762320
                 0.26311282963464    0.72849239295540    0.02723031417443
                 0.72849239295540    0.00839477740996    0.02723031417443
                 0.00839477740996    0.26311282963464    0.02723031417443
                 0.72849239295540    0.26311282963464    0.02723031417443
                 0.26311282963464    0.00839477740996    0.02723031417443
                 0.00839477740996    0.72849239295540    0.02723031417443];

      elseif (quadOrder == 9)
         points=[0.33333333333333    0.33333333333333    0.09713579628280
                 0.48968251919874    0.48968251919874    0.03133470022714
                 0.48968251919874    0.02063496160252    0.03133470022714
                 0.02063496160252    0.48968251919874    0.03133470022714
                 0.43708959149294    0.43708959149294    0.07782754100477
                 0.43708959149294    0.12582081701413    0.07782754100477
                 0.12582081701413    0.43708959149294    0.07782754100477
                 0.18820353561903    0.18820353561903    0.07964773892721
                 0.18820353561903    0.62359292876193    0.07964773892721
                 0.62359292876193    0.18820353561903    0.07964773892721
                 0.04472951339445    0.04472951339445    0.02557767565870
                 0.04472951339445    0.91054097321109    0.02557767565870
                 0.91054097321109    0.04472951339445    0.02557767565870
                 0.22196298916077    0.74119859878450    0.04328353937729
                 0.74119859878450    0.03683841205474    0.04328353937729
                 0.03683841205474    0.22196298916077    0.04328353937729
                 0.74119859878450    0.22196298916077    0.04328353937729
                 0.22196298916077    0.03683841205474    0.04328353937729
                 0.03683841205474    0.74119859878450    0.04328353937729];

      elseif (quadOrder == 10)
         points=[0.33333333333333    0.33333333333333    0.09081799038275
                 0.48557763338366    0.48557763338366    0.03672595775647
                 0.48557763338366    0.02884473323269    0.03672595775647
                 0.02884473323269    0.48557763338366    0.03672595775647
                 0.10948157548504    0.10948157548504    0.04532105943553
                 0.10948157548504    0.78103684902993    0.04532105943553
                 0.78103684902993    0.10948157548504    0.04532105943553
                 0.30793983876412    0.55035294182100    0.07275791684542
                 0.55035294182100    0.14170721941488    0.07275791684542
                 0.14170721941488    0.30793983876412    0.07275791684542
                 0.55035294182100    0.30793983876412    0.07275791684542
                 0.30793983876412    0.14170721941488    0.07275791684542
                 0.14170721941488    0.55035294182100    0.07275791684542
                 0.24667256063990    0.72832390459741    0.02832724253106
                 0.72832390459741    0.02500353476269    0.02832724253106
                 0.02500353476269    0.24667256063990    0.02832724253106
                 0.72832390459741    0.24667256063990    0.02832724253106
                 0.24667256063990    0.02500353476269    0.02832724253106
                 0.02500353476269    0.72832390459741    0.02832724253106
                 0.06680325101220    0.92365593358750    0.00942166696373
                 0.92365593358750    0.00954081540030    0.00942166696373
                 0.00954081540030    0.06680325101220    0.00942166696373
                 0.92365593358750    0.06680325101220    0.00942166696373
                 0.06680325101220    0.00954081540030    0.00942166696373
                 0.00954081540030    0.92365593358750    0.00942166696373];

      elseif (quadOrder == 11)
         points=[0.53461104827076    0.53461104827076    0.00092700632896
                 0.53461104827076   -0.06922209654152    0.00092700632896
                -0.06922209654152    0.53461104827076    0.00092700632896
                 0.39896930296585    0.39896930296585    0.07714953491481
                 0.39896930296585    0.20206139406829    0.07714953491481
                 0.20206139406829    0.39896930296585    0.07714953491481
                 0.20330990043128    0.20330990043128    0.05932297738077
                 0.20330990043128    0.59338019913744    0.05932297738077
                 0.59338019913744    0.20330990043128    0.05932297738077
                 0.11935091228258    0.11935091228258    0.03618454050342
                 0.11935091228258    0.76129817543484    0.03618454050342
                 0.76129817543484    0.11935091228258    0.03618454050342
                 0.03236494811128    0.03236494811128    0.01365973100268
                 0.03236494811128    0.93527010377745    0.01365973100268
                 0.93527010377745    0.03236494811128    0.01365973100268
                 0.35662064826129    0.59320121342821    0.05233711196220
                 0.59320121342821    0.05017813831050    0.05233711196220
                 0.05017813831050    0.35662064826129    0.05233711196220
                 0.59320121342821    0.35662064826129    0.05233711196220
                 0.35662064826129    0.05017813831050    0.05233711196220
                 0.05017813831050    0.59320121342821    0.05233711196220
                 0.17148898030404    0.80748900315979    0.02070765963914
                 0.80748900315979    0.02102201653617    0.02070765963914
                 0.02102201653617    0.17148898030404    0.02070765963914
                 0.80748900315979    0.17148898030404    0.02070765963914
                 0.17148898030404    0.02102201653617    0.02070765963914
                 0.02102201653617    0.80748900315979    0.02070765963914];

      elseif (quadOrder == 12)
         points=[0.48821738977381    0.48821738977381    0.02573106644045
                 0.48821738977381    0.02356522045239    0.02573106644045
                 0.02356522045239    0.48821738977381    0.02573106644045
                 0.43972439229446    0.43972439229446    0.04369254453804
                 0.43972439229446    0.12055121541108    0.04369254453804
                 0.12055121541108    0.43972439229446    0.04369254453804
                 0.27121038501212    0.27121038501212    0.06285822421789
                 0.27121038501212    0.45757922997577    0.06285822421789
                 0.45757922997577    0.27121038501212    0.06285822421789
                 0.12757614554159    0.12757614554159    0.03479611293071
                 0.12757614554159    0.74484770891683    0.03479611293071
                 0.74484770891683    0.12757614554159    0.03479611293071
                 0.02131735045321    0.02131735045321    0.00616626105156
                 0.02131735045321    0.95736529909358    0.00616626105156
                 0.95736529909358    0.02131735045321    0.00616626105156
                 0.27571326968551    0.60894323577979    0.04037155776638
                 0.60894323577979    0.11534349453470    0.04037155776638
                 0.11534349453470    0.27571326968551    0.04037155776638
                 0.60894323577979    0.27571326968551    0.04037155776638
                 0.27571326968551    0.11534349453470    0.04037155776638
                 0.11534349453470    0.60894323577979    0.04037155776638
                 0.28132558098994    0.69583608678780    0.02235677320230
                 0.69583608678780    0.02283833222226    0.02235677320230
                 0.02283833222226    0.28132558098994    0.02235677320230
                 0.69583608678780    0.28132558098994    0.02235677320230
                 0.28132558098994    0.02283833222226    0.02235677320230
                 0.02283833222226    0.69583608678780    0.02235677320230
                 0.11625191590760    0.85801403354407    0.01731623110866
                 0.85801403354407    0.02573405054833    0.01731623110866
                 0.02573405054833    0.11625191590760    0.01731623110866
                 0.85801403354407    0.11625191590760    0.01731623110866
                 0.11625191590760    0.02573405054833    0.01731623110866
                 0.02573405054833    0.85801403354407    0.01731623110866];

      elseif quadOrder>12
          
         msgbox('Highest possible order of quadrature is 12.');
         
      end

    %Quadrilateral elements
    case 'quad'
        
        %Order of quadrature
        quadOrder=round(nNodes^0.5);
        %Getting quadrature points from @gaussLegQ function        
        quadrature=gaussLegQ(quadOrder);

        %Number of integration points
        numPoints=quadOrder^2;
        
        %Initializing variables
        coord=ones(numPoints,4);
        csi=repmat(quadrature(:,1),quadOrder,1);
        wCsi=repmat(quadrature(:,2),quadOrder,1);
        coord(:,1)=csi;
        coord(:,3)=wCsi;
        
        %Setting integration points coordinates
        k=1:size(quadrature,1):numPoints;
        for i=1:size(quadrature,1)
            coord(k(i):k(i)+size(quadrature,1)-1,2)=repmat(-quadrature(i,1),size(quadrature,1),1);
            coord(k(i):k(i)+size(quadrature,1)-1,4)=repmat(quadrature(i,2),size(quadrature,1),1);
        end
                  
        points=coord;
        
    end
    
    
end

%======================================================================
%>@brief Provides the Gaussian points and weights for the n-th order
%>quadrature over a quadrilateral domain
%>
%>@param N (@b int) Order of quadrature
%>
%>@retval q (@b matrix) Gauss points and respective weights
%>
%>@details
%>Calculates the roots of the n-th order Legendre polynomials, defined 
%>recursively. The process has no analytical solution, thus the roots are 
%>numerically approximated using the Newton-Raphson method.
%>
%>@note
%> <strong>`q`</strong> has dimensions:
%> - (<strong>`N`</strong> x 2).
%>
%>@note
%>Columns are organized such that:
%> - __First column__: roots of the Legendre polynomials;
%> - __Second column__: weighting factors.
%>
%>@see
%>Rosetta code<br>
%>[Numerical integration/Gauss-Legendre quadrature](https://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature)
%> <br>Theory manual<br>
%>[Chapter 4, Subsection 4.4.2 - Gauss points generation](http://grids.web.ua.pt/wp-content/uploads/2018/09/DISSERTACAO_Ruben_Lourenco_MIEM_1718_final.pdf#page=57)
%======================================================================
function [ q ] = gaussLegQ( N )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function gaussLegQ provides the Gaussian points and weights         %
%   for the quadrature of order n over the quadrilateral domain.        %
%                                                                       %
%   Input: N   - the order of the Gaussian quadrature                   %
%                                                                       %
%   Output: q - (N by 2) matrix:                                        %
%               1st column gives the x-coordinates of points            %
%               2nd column gives the weights                            %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

   % Initial guess
   x=cos(pi*((1:N)'-0.25)/(N+0.5));  
   % Legendre-Gauss Matrix
   L=zeros(N,N+1);                     

   x0=2;
   while max(abs(x-x0))>eps 
        L(:,1)=1; 
        L(:,2)=x;
   
        for k=2:N
            % Legendre polynomial Pn defined by recursive rule
            L(:,k+1)=((2*k-1)*x.*L(:,k)-(k-1)*L(:,k-1))/k; 
        end
        % Derivative Pn'
        Lp=N*(x.*L(:,N+1)-L(:,N))./(x.^2-1);       
        % Newton-Raphson method for the roots
        x0=x;
        x=x0-L(:,N+1)./Lp;         
   end
   %Computing weights
   w=2./((1-x.^2).*Lp.^2); 
   q=[-x w];

end

%======================================================================
%>@file intPoints.m
%>@brief Functions for Gaussian points generation
%>@details
%>
%>@author R�ben Louren�o
%>@date 16-Oct-2018
%======================================================================