%======================================================================
%>@brief Launches the GUI window to define loads and boundary conditions
%>
%>@param nodeCoord (@b matrix) Nodal coordinates matrix.
%>@param connecInfo (@b matrix) Connectivity matrix.
%>
%>@retval fixedNodes (@b matrix) Labels identifying the fixed nodes.
%>@retval sSnodes (@b struct) Structure array with labels identifying the simply
%>supported nodes. Also indicates which direction is fixed.
%>@retval loads (@b struct) Structure array containing the magnitude,
%>direction and point of application of the loads.
%>
%>@details
%>Launches a figure window with a representation of the meshed structure,
%>allowing the user to select the nodes with applied loads and boundary
%>conditions.<br><br>
%>Structure array <strong>`sSnodes`</strong> is organized in the following
%>format:
%>\htmlonly <style>div.image img[src="struct2.png"]{width:20%;height=20%}</style> \endhtmlonly 
%>@image html struct2.png
%>@image latex struct2.pdf width=4cm
%> <br><br>
%>%>Structure array <strong>`loads`</strong> is organized in the following
%>format:
%>\htmlonly <style>div.image img[src="struct3.png"]{width:20%;height=20%}</style> \endhtmlonly 
%>@image html struct3.png
%>@image latex struct3.pdf width=4cm
%>
%>@todo
%>Review code logic to maintain the correct geometry ratio of the
%>structure, given the plot window size. The current code used for this is a 
%>complete mess. Should be converted into a function to be used in more 
%>general cases.
%>@bug
%>There's a bug when the user closes the figure window without defining any
%>BC's. The program launches a message box warning that the program could
%>not proceed and will quit, but the main routine continues normally
%>afterwards.
%======================================================================
function [fixedNodes, sSnodes, loads]=nodeSelectionThermal(nodeCoord, connectInfo)

%Separating nodal coordinates into vectors
xData=nodeCoord(:,2);
yData=nodeCoord(:,3);
%Array with different window titles
strTitle={'Select nodes with 0ºC...','Select node with prescribed temperature...',...
          'Select faces with heat fluxes...'};
%Keywords for BC definition
bcType={'zeroTemp','tempPresc','fluxes'};

%Initializing variables
fixed=[];
simplySup=struct('BC',[]);
load=struct('Xdirection',[],'Ydirection',[]);

%Launching figure window
fig=figure('Name','Visualization','NumberTitle','off',...
            'Position', [1 41 1440 700],'Visible','off','color','w');

%Initializing axes handle and defining its position
ax=axes('Position',[0.05 0.07 0.75 0.86786]);

%Creating button
btnNext=createUIControl(fig,'Next...','pushbutton',[0.815 0.03 0.0893 0.0476]);

%Representing the meshed structure
patch(ax,'faces',connectInfo(:,2:end),'vertices',[xData(:) yData(:)],'facecolor',[0.9 0.9 0.9],'FaceAlpha',0.65);
hold on
%Overlaying the mesh nodes
nodes=plot(ax,xData,yData,'.','MarkerSize',8,'Tag','nodes');

%------------------------------------------------------------
% Trying to define the correct geometry ratio
%------------------------------------------------------------
ratio=max(xData)/max(yData);
axis equal;
% xlim(ax,[min(xData)-max(xData)*0.2 max(xData)*1.2])
% ylim(ax,[min(yData)-max(yData)*0.2 max(yData)*1.2])

if ratio>1
    ax.XLim=ax.XLim+[-ax.XLim(2)*0.05 ax.XLim(2)*0.05];
%     ax.DataAspectRatio=[1 0.75 1];
elseif ratio<1
    ax.YLim=ax.YLim+[-ax.YLim(2)*0.05 ax.YLim(2)*0.05];
%     ax.DataAspectRatio=[0.75 1 1];
else
    ax.YLim=ax.YLim+[-ax.XLim(2)*0.05 ax.XLim(2)*0.05];    
end
%-----------------------------------------------------------

%Turning axes off
axis off
%Initializing brush tool
b=brush(fig);

%-------------------------------------------------------------------------
% Code for going through the multiple steps of defining the BC's and loads
%--------------------------------------------------------------------------
i=1;
while i<=length(bcType)
    
    if isempty(findobj('Type','figure'))
        message=sprintf('Window closed by user.\nProgram cannot continue without boundary conditions and will now quit.');
        uiwait(msgbox(message));
        break;
    end
    
    switch bcType{i}

        case 'zeroTemp'
            b.Color=[1 0 0];
            set(b,'Enable','on','ActionPostCallBack',@registerValsFixed);
        case 'tempPresc'
            b.Color=[0 1 0];
            btnRegister=createUIControl(fig,'Register','pushbutton',[0.82 0.6 0.0893 0.0476]);

            %Creating radiobutton group
%             bg = uibuttongroup(fig,...
%                       'Title','Select fixed direction...',...
%                       'Position', [0.82 0.69 0.15 0.1],'Tag','radioGroup');
% 
% 
%             % Create two radio buttons in the button group.
%             r1 = createUIControl(bg,'x','radiobutton',[30 5 30 20]);
%             r2 = createUIControl(bg,'y','radiobutton',[70 5 30 20]);  
%             
%             set(btnRegister,'Callback',@registerValsSS)
            set(b,'Enable','on','ActionPostCallBack',@onBrushData);
            
        case 'fluxes'
            b.Color=[0 0 1];
%             delete([r1 r2 bg])
            set(btnNext,'String','Finish');
            set(btnRegister,'Position',[0.85 0.55 0.0893 0.0476]);
            loadX=createUIControl(fig,'none','edit',[0.85 0.69 0.1 0.0476]);
            createUIControl(fig,'Fx','text',[0.8 0.69 0.04 0.0476]);
            loadY=createUIControl(fig,'none','edit',[0.85 0.6 0.1 0.0476]);
            createUIControl(fig,'Fy','text',[0.8 0.6 0.04 0.0476]);
            
%             set(btnRegister,'Callback',@registerValLoads);
            set(b,'Enable','on','ActionPostCallBack',@onBrushData);

    end     
    
    if i==1
        fig.Visible='on';
    end
    
    title(ax,strTitle{i});

    set(btnNext,'Callback',@goNext);
    
    drawnow
    uiwait(gcf);
    i=i+1;    

end
%--------------------------------------------------------------------------

%Output variables
fixedNodes=fixed;
sSnodes=simplySup;
loads=load;

%Closing any figure windows that may exist
if ~isempty(findobj('Type','figure'))
    close(fig)
end

%----------------------------------------------------------------
%                  NESTED CALLBACK FUNCTIONS             
%----------------------------------------------------------------                                                   

    %======================================================================
    %>@brief Registers the labels of the fixed nodes
    %>
    %>@param none
    %>@retval none
    %>
    %>@details
    %>Registers the labels of the fixed nodes by converting the indices
    %>returned by the brush tool.
    %>@note
    %>Parameters represented by a tilde (`~`) are ignored.
    %>
    %======================================================================
    function registerValsFixed(~,~)
 
          fixed=find(nodes.BrushData);
                   
    end

    %======================================================================
    %>@brief Registers the labels of the simply supported nodes and the
    %>corresponding fixed direction.
    %>
    %>@param hObject (@b object) Handle to the object that has provided the
    %>interaction
    %>@retval none
    %>
    %>@details
    %>Registers the labels of the simply supported nodes and corresponding 
    %>fixed direction by converting the indices returned by the brush tool 
    %>and inspecting the state of the radiobuttons.
    %>@note
    %>Parameters represented by a tilde (`~`) are ignored.
    %======================================================================
    function registerValsSS(hObject,~)

        persistent k

        if isempty(k)
            k=1;
        end
        
        bc=struct('FixedDir',[],'nodes',[]);
        
        if isempty(simplySup.BC)
            simplySup.BC=bc;
        end
        
        if ~isempty(nodes.BrushData)
%             set(bg.Children,'Enable','off');
%             bc.FixedDir=get(bg.SelectedObject,'String');
%             bc.nodes=find(nodes.BrushData);
            
            simplySup.BC(k)=bc;
            
            k=k+1;
        end

        hBrushHandles = nodes.BrushHandles;
        hBrushChildrenHandles = hBrushHandles.Children;
        if get(bg.SelectedObject,'String')=='x'
             drawBC(hBrushChildrenHandles(1).VertexData, 'fixedX');
        else
             drawBC(hBrushChildrenHandles(1).VertexData, 'fixedY');
        end
        hBrushChildrenHandles(1).VertexData=[];
        
        set(hObject, 'Enable', 'off');
    
    end

    %======================================================================
    %>@brief Registers the labels of the nodes with applied loads, the
    %>direction of the load and the corresponding values.
    %>
    %>@param hObject (@b object) Handle to the object that has provided the
    %>interaction
    %>@retval none
    %>
    %>@details
    %>Registers the labels of the nodes with applied loads, the direction
    %>of the load and the corresponding values, by converting the indices
    %>of returned by the brush tool and inspecting the values of the edit
    %>fields.
    %>@note
    %>Parameters represented by a tilde (`~`) are ignored.
    %======================================================================
    function registerValLoads(hObject,~)
        
        persistent k
        
        if isempty(k)
            k=1;
        end
        
        fX=struct('value',[],'nodes',[]);
        fY=struct('value',[],'nodes',[]);
        
        if isempty(load.Xdirection) || isempty(load.Ydirection)
            load.Xdirection=fX;
            load.Ydirection=fY;
        end
        
        if ~isempty(nodes.BrushData)
          
            if ~isempty(get(loadX,'String'))
                fX.value=str2double(get(loadX,'String'));
            else
                fX.value=0;
            end
            
            fX.nodes=find(nodes.BrushData);
        
            if ~isempty(get(loadY,'String'))
                fY.value=str2double(get(loadY,'String'));
            else
                fY.value=0;
            end
            
            fY.nodes=find(nodes.BrushData);
            
            load.Xdirection(k)=fX;
            load.Ydirection(k)=fY;
            
            k=k+1;
            
        end 
        
        hBrushHandles = nodes.BrushHandles;
        hBrushChildrenHandles = hBrushHandles.Children;
        hBrushChildrenHandles(1).VertexData=[];
        set(findobj('Style','edit'),'String','','Enable','off');
        
        set(hObject, 'Enable', 'off');
    end

    %======================================================================
    %>@brief Provides the necessary user-interface objects on brush
    %>selection.
    %>
    %>@param none
    %>@retval none
    %>
    %>@details
    %>Upon node selection with the brush tool, enables the radiobuttons and
    %>edit fields inside the radiogroup. Also enables a button to register
    %>the values from these objects.
    %>@note
    %>Parameters represented by a tilde (`~`) are ignored.
    %======================================================================
    function onBrushData(~,~)
       
        if ~isempty(nodes.BrushData)
            
            if ~isempty(findobj('Tag','radioGroup'))
                set(bg.Children, 'Enable', 'on');
            end
            if ~isempty(findobj('Style','edit'))
                set(findobj('Style','edit'),'Enable','on');
            end
            set(btnRegister, 'Enable', 'on');            
            
        end
               
    end

    %======================================================================
    %>@brief Performs back-end operations when the user clicks "Next".
    %>
    %>@param none
    %>@retval none
    %>
    %>@details
    %>When the user hits the "Next" button reinitializes some variables.
    %>@note
    %>Parameters represented by a tilde (`~`) are ignored.
    %======================================================================
    function goNext(~,~)
       
        try
            hBrushHandles = nodes.BrushHandles;
            hBrushChildrenHandles = hBrushHandles.Children;
            if bcType{i}=='fixedNodes'
                drawBC(hBrushChildrenHandles(1).VertexData, bcType{i});
            end
            hBrushChildrenHandles(1).VertexData=[];
            uiresume(gcbf)
        catch
            uiresume(gcbf)   
        end

    end

%----------------------------------------------------------------

end

%======================================================================
%>@brief Simple function to make it easier to create user-interface
%>objects.
%>
%>@param objHandler (@b obj) Handler to the where the object is to be
%>created.
%>@param string (@b string) String that the object should contain.
%>@param style (@b string) Keyword to the object type being created.
%>@param position (@b array) 4x1 array of values setting the position of
%>the object.
%>@retval uiObject (@b obj) Handle to the created object.
%>
%>@details
%>Registers the labels of the nodes with applied loads, the direction
%>of the load and the corresponding values, by converting the indices
%>of returned by the brush tool and inspecting the values of the edit
%>fields.
%>@note
%> <strong>`position`</strong> should be defined in normalized units and
%>has the following structure:
%>  - [x y width height]
%>@todo
%>Code is replicated in other files, this should be transferred to a
%>separate file and taylored to be suitable for more general cases.
%>@see
%>[uicontrol] (https://www.mathworks.com/help/matlab/ref/uicontrol.html)
%>| [UIControl Properties](https://www.mathworks.com/help/matlab/ref/matlab.ui.control.uicontrol-properties.html)
%======================================================================
function [uiObject]=createUIControl(objHandler,string,style,position)
    
    switch style
        case 'radiobutton'
            uiObject= uicontrol(objHandler,'Style', style, 'String', string,...
            'Position', position,'Units','normal','Enable','off');
        case 'edit'
            uiObject= uicontrol(objHandler,'Style', style,...
            'Position', [20 5 50 20],'Units','normal');
            set(uiObject,'Position', position,'Enable','off');
        otherwise
            uiObject= uicontrol(objHandler,'Style', style, 'String', string,...
            'Position', [20 5 50 20],'Units','normal');
            set(uiObject,'Position', position);
    end

end

%======================================================================
%>@file nodeSelection.m
%>@brief Functions to create a GUI for boundary conditions and loads
%>definition
%>@see
%>MATLAB documentation<br>
%>[Callback definition](https://www.mathworks.com/help/matlab/creating_plots/callback-definition.html#buhztrr-8)
%>| [Nested functions](https://www.mathworks.com/help/matlab/matlab_prog/nested-functions.html)
%>| [brush](https://www.mathworks.com/help/matlab/ref/brush.html)
%>| [uiwait](https://www.mathworks.com/help/matlab/ref/uiwait.html)
%>@author Rúben Lourenço
%>@date 16-Oct-2018
%==========================