
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
function [answer]=checkAnswer(materialProps)
  
    if isempty(materialProps{1}) || isempty(materialProps{2})
        message = sprintf('One or more fields are not defined.\n');
        uiwait(msgbox(message)); 
        answer=1;   
    else
        answer=0;
    end

end