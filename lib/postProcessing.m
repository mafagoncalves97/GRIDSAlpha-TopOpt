%======================================================================
%>@brief Launches the GUI window for post-processing of results
%>
%>@param U (@b matrix) Displacement field
%>@param u1 (@b matrix) Horizontal displacements
%>@param u2 (@b matrix) Vertival displacements
%>@param numNodes (@b integer) Total number of nodes
%>@param connectInfo (@b matrix) Connectivity matrix
%>@param nodeCoord (@b matrix) Nodal coordinates matrix
%>
%>@retval none
%>
%>@details
%>Launches a figure window with a representation of the meshed structure,
%>allowing the user to choose differente field variables to be presented 
%>as contour plots. The corresponding nodal values can be inspected via
%>selection provided by the data brushin tool. The user may also toggle 
%>the deformed/undeformed state of the structure.
%>@todo
%>Reduce the number os parameters of the function. The horizontal and
%vertical displacements are readily available from <strong>`U`</strong>
%>@bug
%>There's a bug when the user chooses to represent the undeformed
%>structure, where the contour plot for the selected field variable remains
%>active. For this state, the contour plot should be blank.
%>
%======================================================================

function [] = postProcessing (U, u1, u2, S, numNodes, connectInfo, nodeCoord)

xData=nodeCoord(:,2);
yData=nodeCoord(:,3);

%Create figure
fPost=figure('Name','Visualization','NumberTitle','off',...
            'Position', [1 41 1440 700],'color','w');

%Create axes object
ax=axes(fPost,'Position',[0.05 0.07 0.75 0.86786]);

%scaleFact=1;
scaleFact=0.2*max(xData)/max((u1.^2+u2.^2).^0.5);

ax.XLim=[min(xData+u1*scaleFact)-max(xData+u1*scaleFact)*0.1 max(xData+u1*scaleFact)*1.1];
ax.YLim=[min(yData+u2*scaleFact)-max(yData+u2*scaleFact)*0.1 max(yData+u2*scaleFact)*1.1];

ratio=max(abs(xData))/max(abs(yData));
multiple=ratio/(max(xData+u1*scaleFact)/max(yData+u2*scaleFact));

if ratio<1
    ax.XLim=ax.XLim*multiple;
    ax.XLim=[min(ax.XLim)-(max(ax.XLim)-max(xData+u1*scaleFact))*0.3 max(ax.XLim)-(max(ax.XLim)-max(xData+u1*scaleFact))*0.3];
elseif ratio>1
    ax.YLim=ax.YLim*(1+multiple);
    ax.YLim=[min(ax.YLim)-(max(ax.YLim)-max(yData+u2*scaleFact))*0.3 max(ax.YLim)-(max(ax.YLim)-max(yData+u2*scaleFact))*0.3];
else
    
   axis equal
   
   
end

axis off
%Creating radiobutton group
bg = uibuttongroup(fPost,'Title','Select contour data...','BackgroundColor','w',...
                   'Position', [0.82 0.63 0.15 0.30],'Tag','radioGroup',...
                   'SelectionChangedFcn',@changeContourData);

cb=uicontrol(bg,'Style','checkbox','String','Plot undeformed','BackgroundColor','w', ...
                   'Value',0,'units','normalized','Position',[0.5 0.75 0.5 0.15],...
                    'Callback',@plotUndeformed);               
               
% Loop to create radio buttons inside button group.
radioLabels=["U1","U2","Sx","Sy","Sxy","Von Mises"];
radioBtns=cell(1,length(radioLabels));

for i=1:length(radioLabels)
    %Create object
    radioBtns{i}=createUIControl(bg,radioLabels(i),'radiobutton',[30 5 30 20]);
    set(radioBtns{i},'Units','normal','Tag',num2str(i));
    if i==1
        set(radioBtns{i},'Position',[0.1 0.75 0.2 0.15]);
    else
        set(radioBtns{i},'Position',[0.1 radioBtns{i-1}.Position(2)-0.15 0.5 0.15]);
    end
    
end

%Setting global variables
global nodes contour c h t

%Setting contour data according to default selection at first run
changeContourData (bg);

%Setting data brushing
b=brush(fPost);
set(b,'Color',[0 1 0],'Enable','on','ActionPostCallBack',@onBrushData);

%Wait until user closes figure
uiwait(fPost);

%----------------------------------------------------------------
%                  NESTED CALLBACK FUNCTIONS 
%----------------------------------------------------------------                                                   
    %======================================================================
    %>@brief Changes the contour plot 
    %>
    %>@param hObject (@b object) Handle to the object providing the interaction
    %>
    %>@retval none
    %>
    %>@details
    %>Provides the necessary logic in order for the user to be able to
    %>select different contour plots.
    %======================================================================
 
    function changeContourData (hObject,~)
        
        if ~isempty(findobj('Type','colorbar')) && ~isempty(findobj('Type','patch'))         
          delete([contour h]); 
        end
        
        if ~isempty(findobj('Type','uitable'))
            set(t,'Data',[],'ColumnName',{'Node',bg.SelectedObject.String});
        end
        
        selectedObjTag=hObject.SelectedObject.Tag;
        
        set(cb,'value',0);
        
        switch selectedObjTag
        
            case '1'
                contour=plotContour(ax,'U1');
            case '2'
                contour=plotContour(ax,'U2');
            case '3'
                contour=plotContour(ax,'Sx');
            case '4'
                contour=plotContour(ax,'Sy');
            case '5'
                contour=plotContour(ax,'Sxy');
            case '6'
                contour=plotContour(ax,'SVM');
        end
    
    end
%---------------------------------------------------
    
    %======================================================================
    %>@brief Plots the deformed/undeformed geometry.
    %>
    %>@param hObject (@b object) Object Handle to the object providing the interaction
    %>
    %>@retval none
    %>
    %>@details
    %>Allows the plotting of the deformed or undeformed geometry.
    %======================================================================
    function plotUndeformed(hObject, ~, ~)
    
        if (get(hObject,'Value') == get(hObject,'Max'))
              %Plot undeformed shape
              set(contour,'vertices',[xData(:) yData(:)]);              
              set(nodes,'XData',xData(:),'YData',yData(:));
        else
              %Plot deformed shape
              set(contour,'vertices',[xData(:)+u1*scaleFact yData(:)+u2*scaleFact]);
              set(nodes,'XData',xData(:)+u1*scaleFact,'YData',yData(:)+u2*scaleFact);
        end
        
    end
%-------------------------------------------------------

    %======================================================================
    %>@brief Plots the selected contour.
    %>
    %>@param axesHandler (@b object) Handler to the axes object where the
    %>the plot is to be drawn
    %>@param dataType (@b string) String identifying which field variable
    %is to be plotted
    %>
    %>@retval contour (@b contour) Handle to the patch object being created
    %>
    %>@details
    %>Allows for the user to choose which field variable is to be
    %>plotted. Creates a patch object whose faces are colored according to
    %>the data concerning the given field variable. This is done through
    %>the command <strong>`'FaceVertexCData'`</strong>, color interpolation
    %>is done automatically. A colormap with 12 colors is applied and a 
    %>corresponding legend is created.
    %======================================================================
 
    function [contour]=plotContour(axesHandler,dataType)
        
       switch dataType
       
           case 'U1'
                
               data=U(((1:numNodes)'-1)*2+1);
               %Set axes title
               ax.Title.String='Horizontal displacement (U1)';
               
           case 'U2'
               
               data=U(((1:numNodes)'-1)*2+2);
               %Set axes title
               ax.Title.String='Vertical displacement (U2)';
               
           case 'Sx'
               
               data=S(:,1);
               ax.Title.String='Normal stress (Sx)';
           
           case 'Sy'
               
               data=S(:,2);
               ax.Title.String='Normal stress (Sy)';
           
           case 'Sxy'
               
               data=S(:,3);
               ax.Title.String='Shear stress (Sxy)';
           case 'SVM'
               
               data=S(:,4);
               ax.Title.String='Von Mises stress (SVM)';
       
       end
      
      contour = patch(ax,'faces',connectInfo(:,2:end),'vertices',[xData(:)+u1*scaleFact yData(:)+u2*scaleFact],...
                      'facecolor','interp','FaceVertexCData',data,'FaceAlpha',1,'Marker','.',...
                      'MarkerFaceColor',[0 0 0],'EdgeAlpha',0.5);        
      
      hold on
           
      if ~isempty(findobj('Type','scatter'))
          
          delete(nodes)
          
      end
      
      nodes=scatter(ax,xData(:)+u1*scaleFact,yData(:)+u2*scaleFact,12,'o','filled','MarkerFaceAlpha',0);            
      numcolors=12; 
      c=colormap(ax,jet(numcolors));
      caxis(ax,[min(data) max(data)]);
      h=colorbar(ax,'Location','east');
      set(get(h,'Label'),'String',dataType,'FontWeight','bold','FontSize',12);
      set(h,'Ticks',linspace(min(data),max(data),numcolors+1));
      drawnow
    end
%-----------------------------------------
    %======================================================================
    %>@brief Handles the data brushing functionality
    %>
    %>@param hObject (@b object) Handle to the object providing the interaction
    %>
    %>@retval none
    %>
    %>@details
    %>Allows the user to probe the node values, for a given field variable,
    %>using the data brushing tool. A table is populated with the selected
    %>nodes' values.
    %>
    %======================================================================
 
    function onBrushData(hObject, ~)
       
        if ~isempty(nodes.BrushData)
            
            nodeLabels=(1:numNodes)';
            selectedNodes=find(nodes.BrushData);
                       
            switch bg.SelectedObject.Tag
        
                case '1'
                    var=u1(selectedNodes);
                case '2'
                    var=u2(selectedNodes);
                case '3'
                    var=S(selectedNodes,1);
                case '4'
                    var=S(selectedNodes,2);
                case '5'
                    var=S(selectedNodes,3);
                case '6'
                    var=S(selectedNodes,4);                    
            end                        
            
            data=[nodeLabels(selectedNodes),var];           
            
            if isempty(findobj('Type','uitable'))
                
                t = uitable(fPost,'Data',data,'Units','normal','Position',[0.82 0.07 0.15 0.55]);
                set(t,'ColumnName',{'Node',bg.SelectedObject.String},'RowName',{},'ColumnWidth',{55 100});
            
            else
                
                set(t,'Data',data,'ColumnName',{'Node',bg.SelectedObject.String});
                
            end    
            
        end
               
    end


end


function [uiObject]=createUIControl(objHandler,string,style,position)
    
    switch style
        case 'radiobutton'
            uiObject= uicontrol(objHandler,'Style', style, 'String', string,...
            'Position', position,'Units','normal','Enable','on','BackgroundColor','w');
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
%>@file postProcessing.m
%>@brief Functions to create a GUI for post-processing of results
%>@details
%>
%>@author Rúben Lourenço
%>@date 2-Nov-2018
%======================================================================