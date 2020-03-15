function [ProblemType] = inputProblemType()
    
    list = {
            'Elástico',
            'Térmico',
            'Termoelástico'
            };

    WAIT_USER_INPUT=true;    
    while WAIT_USER_INPUT

        [indx,tf] = listdlg('ListString',list,'SelectionMode','single',...
                            'ListSize',[350,150],'Name','Select problem type...');    
            if tf==1
                ProblemType=list{indx};
                WAIT_USER_INPUT=false;
        else
            quest=sprintf('The analysis could not proceed, unless a finite element formulation is selected.\n\nClick Continue to return to the element selection dialog.\nClick Quit to stop the program.');
            answer = questdlg(quest,'Warning','Continue','Quit','Continue');
            if strcmp(answer,'Quit')
               error('Program execution stopped by user');
            end            
        end
    
    
    end
    
end

