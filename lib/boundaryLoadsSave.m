function [answer]=boundaryLoadsSave()

 prompt={'Use the same conditions than last time:'};
 title='Boundary conditions and loads';
 dims=[1 35];
 definput={'Yes'};
 answer=inputdlg(prompt,title,dims,definput);
 
    
end 