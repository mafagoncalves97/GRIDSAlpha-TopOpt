%======================================================================
%>@brief Plots the different eigenmodes for a given finite element
%>formulation
%>
%>@param vects (@b matrix) Matrix of eigenvectors
%>@param vals (@b matrix) Matrix of eigenvalues
%>@param nodeCoord (@b matrix) Matrix of nodal coordinates
%>@param connectInfo (@b matrix) Element connectivity matrix
%>
%>@retval none
%>
function plotEigen(V,lambdas,strainEnergy,nodeCoord,connectInfo,intPts)

nodeCoord(:,1)=[];
connectInfo(:,1)=[];

if size(intPts,1)>1
    intPts=intPts([3 4 2 1],:);
end

nModes=size(V,2);

xData=real(V(1:2:nModes,:));
yData=real(V(2:2:nModes,:));

x=xData+nodeCoord(:,1);
y=yData+nodeCoord(:,2);
nCols=3;
nRows=round(nModes/nCols);

fig=figure(1);
set(fig,'color','w');
hTabGroup=uitabgroup(fig,'TabLocation','left');

for i=1:nModes+1

    if i==1
        tab(i)=uitab(hTabGroup,'Title','Eigenmodes','BackgroundColor','w');
        createUIControl(tab(i),sprintf('rank(K)=%d',nnz(lambdas)),'text',[0 0.9 0.15 0.05]);
        axes('parent',tab(i));
        
        for k=1:nModes
            subplot(nRows,nCols,k);
            patch('faces',connectInfo,'vertices',[nodeCoord(:,1) nodeCoord(:,2)],'facecolor','none','edgecolor',[0 0 0],'linestyle','--');
            patch('faces',connectInfo,'vertices',[x(:,k) y(:,k)],'facecolor','none','edgecolor',[1 0 0]);
            title(sprintf('\\lambda_%d=%0.3f',k,lambdas(k)));
            pbaspect([1 1 1]); axis off;
            xlim([-2 2])
            ylim([-2 2])
        end
        
    else
        tab(i)=uitab(hTabGroup,'Title', sprintf('Mode %i', i-1),'BackgroundColor','w');
        axes('parent',tab(i));
        patch('faces',connectInfo,'vertices',[nodeCoord(:,1) nodeCoord(:,2)],'facecolor','none','edgecolor',[0 0 0],'linestyle','--');
        patch('faces',connectInfo,'vertices',[x(:,i-1) y(:,i-1)],'facecolor','none','edgecolor',[1 0 0]);
        tx=text(intPts(:,1), intPts(:,2), num2str(strainEnergy(:,i-1)),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize', 8);
        title(sprintf('\\lambda_%d=%0.3f',i-1,lambdas(i-1)));
        pbaspect([1 1 1]);
        xlim([-2 2])
        ylim([-2 2])
    end
    
end

end

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
            'Position', [20 5 50 20],'Units','normal','BackgroundColor','w','FontSize',11);
            set(uiObject,'Position', position);
    end

end

%======================================================================
%>@file plotEigen.m
%>@brief Algorithm to plot eigenmodes.
%>@details
%>
%>@author Rúben Lourenço
%>@date 25-Dez-2018
%======================================================================
