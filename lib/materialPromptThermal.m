%======================================================================
%>@brief Launches an input dialog
%>
%>@retval materialProps (@b cell) Cell array of character vectors 
%>containing one input per edit field.
%>
%>@details
%>Provides and input dialog with the necesary edit field for the user to
%>input the material properties and thickness of the structure.
%>
%>@todo
%>Review code to handle multiple-part models which might have different
%>material properties.
%>@todo
%>This could be transferred to a separate file and modified to serve more
%>general purposes related with user inputs.
%>
%>@see
%>MATLAB documentation<br>
%>[inputdlg](https://www.mathworks.com/help/matlab/ref/inputdlg.html#d120e619809)
%======================================================================
function [materialPropsThermal]=materialPromptThermal()
    
    promptMaterial = ({'Thermal conductivity:','MafaLinda:','Thickness:'});
    definput={'50.2','0.3','0.001'};
    
    WAIT_USER_INPUT=true;
    while WAIT_USER_INPUT
        materialPropsThermal = inputdlg(promptMaterial,'Material properties',[1 35],definput);
        if ~isempty(materialPropsThermal)
            checkAnswer(materialPropsThermal);
        else
            quest=sprintf('The analysis could not proceed, unless the material properties are defined.\n\nClick Continue to return to the element selection dialog.\nClick Quit to stop the program.');
            answer = questdlg(quest,'Warning','Continue','Quit','Continue');
            if strcmp(answer,'Quit')
               error('Program execution stopped by user');
            end
        end
    end   

%======================================================================
%>@brief Checks the validity of the answer coming from an input dialog
%>
%>@param materialProps (@b cell) Cell array of character vectors 
%>containing one input per edit field.
%>
%>@retval answer (@b int) Value corresponding to the validity of the
%>answer.
%>
%>@details
%>Checks if the edit fields of an input dialog are correctly filled by the
%>user by evaluating the answer coming from the input dialog. A message box
%>is launched if some fields are left blank.
%>
%>@todo
%>Define this function in a separate .m file, where future function 
%>definitions aimed at error handling would be created.
%>@todo
%>Add code logic to check if the values coming from the edit fields are in
%>the correct format. Code should be dynamic enough to be able to handle the 
%>results from multiple edit fields.
%======================================================================
function checkAnswer(materialProps)
  
    if isempty(materialProps{1}) || isempty(materialProps{2}) || isempty(materialProps{3})
        message = sprintf('One or more fields are not defined.\n');
        uiwait(msgbox(message));    
    else
        WAIT_USER_INPUT=false;
    end

end
    
end

