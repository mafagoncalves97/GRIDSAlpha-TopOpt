%======================================================================
%>@file
%>@brief Functions for input file reading and parsing
%>@details
%>Contains function definitions related to input file importing 
%>to MATLAB and subsequent data reading and parsing into appropriate
%>variables.
%>
%> ### What is an input file?
%>
%>An input file stores the necessary information to define a 
%>given Finite Element model. The items listed in such file may include, 
%>for example:
%>
%> - nodal coordinates;
%> - connectivity data;
%> - sets of nodes/elements with applied loads and boundary conditions;
%> - contact conditions; 
%> - material properties;
%> - etc.
%> 
%>On GRIDS Alpha, the material properties, applied loads and boundary 
%>conditions are defined via the graphical user interface. Thus, only the 
%>nodal coordinates and connectivity data need be listed initially in the 
%>input file, in order to define the mesh.
%>
%> ### Input file structure
%>
%>GRIDS Alpha expects an Abaqus-style input file with extension
%>`.inp` or `.txt`. In this format, data is organized into blocks inside 
%>labels which are identified by certain keywords. For a single-part 
%>two-dimensional model is structured as in the example below. For a 
%>multiple-part model, the blocks of code may be replicated analogously.</p>
%>
%>@par Example
%>@code
%>*Part, name=Part-1
%>*Node
%>      1,           0.,           0.
%>      2,          50.,           0.
%>      3,         100.,           0.
%>      4,           0.,          50.
%>      5,          50.,          50.
%>      6,         100.,          50.
%>      7,           0.,         100.
%>      8,          50.,         100.
%>      9,         100.,         100.
%>*Element, type=CPE4H
%>1, 1, 2, 5, 4
%>2, 2, 3, 6, 5
%>3, 4, 5, 8, 7
%>4, 5, 6, 9, 8
%>*End Part
%>@endcode
%>The input data in the example above refers to a mesh of two-dimensional
%>quadrilateral elements, organized as depicted in the following figure:
%>\htmlonly <style>div.image img[src="mesh.png"]{width:20%;height=20%}</style> \endhtmlonly 
%>@image html mesh.png
%>@image latex mesh.pdf width=6cm
%>
%>@warning
%>Even if the model is composed by only one part, the nodal coordinates 
%>and connectivity data still have to be listed between the keywords 
%> <strong>`*Part, name=`</strong> and <strong>`*End Part`</strong>. <br>
%>@warning
%>For multiple-part models, each part has to be given a unique name.
%>
%>@note
%>Dots (<strong>`.`</strong>) are used as decimal separators, while commas 
%>(<strong>`,`</strong>) are used as data separators.
%>@note
%>Nodal coordinates matrix is written such that each line corresponds to 
%>one node. Columns are organized such that:
%> - First column: node label;
%> - Second column: x coordinate;
%> - Third column: y coordinate.
%>
%>@note
%>An element is identified by an element label and defined by the set of 
%>it comprises, by means of the respective node labels. It is assumed that 
%>the nodes are numbered in ascending order on the counterclockwise
%>direction.
%>
%>@see
%> Theory manual<br>
%> [Chapter 4, Subsection 4.3.1 - Model definition](http://grids.web.ua.pt/wp-content/uploads/2018/09/DISSERTACAO_Ruben_Lourenco_MIEM_1718_final.pdf#page=51)
%>
%>@author Rúben Lourenço
%>@date 12-Oct-2018
%>
%>@fn readInput()
%>@brief Reads and parses the information stored in an input file.
%>@retval model (@b struct) Structure array containing the Finite Element 
%>model information
%>
%>@details
%>
%> <p>Provides the necessary algorithm to access the input file and read its
%>contents in order to extract the above information. When the function is
%>called the user is prompted to select the desired file to import.</p>
%>
%> <p>For each part of the model, the algorithm extracts the nodal 
%>coordinates and connectivity information. The data is stored in a 
%>struct array with the name <strong>`model`</strong> in a suitable format, for example:</p>
%>
%>\htmlonly <style>div.image img[src="struct.png"]{width:25%;height=25%}</style> \endhtmlonly 
%>@image html struct.png
%>@image latex struct.pdf width=4cm
%>@see
%> MATLAB Documentation<br>
%> [uigetfile](https://www.mathworks.com/help/matlab/ref/uigetfile.html) | 
%> [fopen](https://www.mathworks.com/help/matlab/ref/fopen.html) | 
%> [textscan](https://www.mathworks.com/help/matlab/ref/textscan.html) | 
%> [find](https://www.mathworks.com/help/matlab/ref/find.html) | 
%> [regexp](https://www.mathworks.com/help/matlab/ref/regexp.html) |
%> [regexprep](https://www.mathworks.com/help/matlab/ref/regexprep.html) | 
%> [structure array](https://www.mathworks.com/help/matlab/ref/struct.html)
%======================================================================
function [model]=readInput()

    disp('-----------------------------------------------------')
    disp('|              INPUT FILE READING MODULE            |')
    disp('-----------------------------------------------------')
    
    disp('>> Waiting for user file selection...')

    %Asks for user to select the input file
    [fileName,pathName,~]=uigetfile({'*.inp';'*.txt'},'Select input file');

    %Opens input file for reading
    fid=fopen(strcat(pathName,'/',fileName),'r');

    disp('>> Reading file contents:')
    disp(['          ' fileName])

    %>Reads the input file contents
    fileContent=textscan(fid,'%[^\n]');

    %>Closes the file
    fclose(fid);

    %Stores the input file lines
    fileLines=fileContent{1,1};

    %Getting indices of lines containing keywords '*Part, name=', '*Node' and 
    %'*Element, type='
    idxParts=find(contains(fileLines,'*Part, name='));
    idxEndParts=find(contains(fileLines,'*End Part'));
    idxNode=find(contains(fileLines,'*Node'));
    idxNode=idxNode(1:size(idxParts,1));
    idxElement=find(contains(fileLines,'*Element, type='));

    %Cell arrays to hold part names and element types
    partNames=cell(size(idxParts,1),1);
    elementTypes=cell(size(idxElement,1),1);

    %Extracting part names
    disp('>> Detecting part names and element types...')
    for i=1:length(partNames)
        splitName=regexp(fileLines(idxParts(i),1),'=','split');
        splitName = regexprep(splitName{1}(1,2), '.$', '', 'lineanchors');
        partNames{i,1}=splitName{1};    
    end

    %Extracting element types
    for i=1:length(elementTypes)
        splitElementType=regexp(fileLines(idxElement(i),1),'=','split');
        splitElementType=regexprep(splitElementType{1}(1,2), '.$', '', 'lineanchors');
        elementTypes{i,1}=splitElementType{1};
    end

    %Initializing struct array to hold model information
    model=struct;

    disp('>> Parsing nodal and connectivity data...')
    for i=1:length(partNames)
      %Storing part names into struct array
      model.Part(i).Name=partNames{i}; 
      %Limits to use in the loops ensuring we're working inside each PART  
      if i==1
          limits=find(idxElement<idxEndParts(i));
      else
          limits=find(idxElement<idxEndParts(i)&idxElement>idxEndParts(i-1));
      end
      j=idxNode(i)+1;
      %Selects lines with nodal coordinates
      while ~contains(fileLines(j,1),'*')
        j=j+1;
      end
      %Extracts nodal coordinates
      nodeRead=fileLines(idxNode(i)+1:j-1,1);
      %Parses the comma separated strings
      parseNodes=regexp(nodeRead, ',', 'split');
      %--------DEBUG--------
      %Displays de node matrix
      %sprintf('PART: %s',partNames{i})
      %------------------------------------------
      nodes = str2double(vertcat(parseNodes{:}));
      %Storing nodal coortinates into struct array
      model.Part(i).Node=nodes;
      z=1;
      for l=min(limits):length(elementTypes)
        k=idxElement(l)+1;
        %Checking if the loop is working inside the corresponding PART
        if k>idxEndParts(i)
            z=1;
            break;
        else
            %Selects the lines with connectivity information
            while ~contains(fileLines(k,1),'*')
                k=k+1;
            end
            %Extracts connectivity information
            elementRead=fileLines(idxElement(l)+1:k-1,1);
            %Parses the comma separated strings
            parseElements=regexp(elementRead, ',', 'split');
            %--------DEBUG--------
            %sprintf('Element type: %s',elementTypes{l})
            %Displays de connectivity matrix
            %------------------------------------------
            elements=str2double(vertcat(parseElements{:}));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nNodes=size(elements,2)-1;
            if (nNodes==3 || nNodes==6 || nNodes==10)
                type='tri';
            elseif (nNodes==4 || nNodes==8 || nNodes==9 || nNodes==16)
                type='quad';
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end 
        %Storing element types and connectvity data into struct array
        model.Part(i).Connect(z).Type=type;
        model.Part(i).Connect(z).Nodes=size(elements,2)-1;
        model.Part(i).Connect(z).Elements=elements;
        z=z+1;
      end

    end
    disp('>> Input file data read successfully.')
    disp('>> Terminating...')
    disp('-----------------------------------------------------')
    disp(' ')
end