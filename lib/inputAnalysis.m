function [D, choiceElement, intMode] = inputAnalysis(E, v)
    
    D=cell(1,2);
    
    list = {
            'PS/FI - Q4',...
            'PE/FI - Q4',...
            'PE/SRI - Q4',...
            'PE/FI/IM - Q6',...
            'PE/FI/IM - QM6',...
            'PE/FI/IM - Q4/6I',...
            'PE/FI/IM - Q4/4I',...
            'PE/FI/CM - Qi5',...
            'PE/FI/CM - Qi6'
            };

    WAIT_USER_INPUT=true;    
    while WAIT_USER_INPUT

        [indx,tf] = listdlg('ListString',list,'SelectionMode','single',...
                            'ListSize',[350,150],'Name','Select finite element type...');    

        if tf==1
            splitResult=strsplit(list{indx},' - ');
            splitAnalysis=strsplit(splitResult{1},'/');
            analysis=splitAnalysis{1};
            intMode=splitAnalysis{2};
            choiceElement=list{indx};
            D = getElasticMat(analysis, intMode,E,v);
            WAIT_USER_INPUT=false;
        else
            quest=sprintf('The analysis could not proceed, unless a finite element formulation is selected.\n\nClick Continue to return to the element selection dialog.\nClick Quit to stop the program.');
            answer = questdlg(quest,'Warning','Continue','Quit','Continue');
            if strcmp(answer,'Quit')
               error('Program execution stopped by user');
            end            
        end

    end
%/////////////////////////////////////////////////////////////////////////   
%     function D = getElasticMat(analysis, intMode)
% 
%         switch analysis
% 
%             case 'PS'
%                 D{1}=(E/(1-v^2))*[1 v 0; v 1 0; 0 0 (1-v)/2];
%             case 'PE'
%                 if strcmp(intMode,'SRI')
%                     Dp=(E/3)*[1/(1-2*v) 1/(1-2*v)      0;
%                           1/(1-2*v) 1/(1-2*v)      0;
%                           0         0         3/(2*(1+v))];
% 
%                     Ds=E/(3*(1+v))*[2 -1 0;
%                                    -1  2 0;
%                                     0  0 0];
%                     D{1}=Dp;
%                     D{2}=Ds;
%                 else
%                     D{1}=(E/((1+v)*(1-2*v)))*[1-v v 0; v 1-v 0; 0 0 (1-2*v)/2];
%                 end              
% 
%         end
% 
% 
%     end
%/////////////////////////////////////////////////////////////////////////
end
