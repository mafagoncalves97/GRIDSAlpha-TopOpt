%======================================================================
%>@file
%>@brief Routine to perform a 2D FEM analysis under plane stress or plane strain
%>@todo
%>Introduce the necessary modifications to the routine in order to convert 
%>it into a funtion. This should make it easier to launch the analysis from 
%>a full fledged GUI, to be developed in the future.
%>@todo
%>Review code for multiple-part models
%>
%>@author Rúben Lourenço
%>@date 17-Oct-2018
%======================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
%                            PLANE STRESS/STRAIN                        %
%                                                                       %
%                                                                       %
%   Purpose: sub-routine for plane stress/strain calculation using FEM  %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Workspace cleaning
%(TODO - turn this into a function in the future)
% function planeStress(pathName,fileName)
clear all
drawnow; pause(0.05)
close all
drawnow; pause(0.05)
clc
drawnow; pause(0.05)
%% Path folders

%Adding functions directory
addpath('../lib')

%% Routine presentation
disp('***************************************************************')
disp('                     PLANE STRESS/STRAIN                       ')
disp('***************************************************************')
drawnow; pause(0.05)

%% Input file reading
disp('> Initiating Input File Reading Module...')
disp(' ')
%Gets model information from input file
model=readInput();
termoelastico=0;elastico=0;termico=0;TopOptElastic=0;TopOptTherm=0;TopOptMO=0;TopOpt=0;fem=0;TopOptElasticDelta=0;
%% Model variables

%(TODO - review code for models with more than one part)

part=model.Part;                 %Getting part from model 

nodeCoord=part.Node;             %Getting nodal coordinates from model
numNodes=size(nodeCoord,1);      %Total number of nodes

connectInfo=part.Connect.Elements;   %Getting connectivity data from model
numElements=size(connectInfo,1);     %Total number of elements

elType=part.Connect.Type;            %Getting element type
elNodes=model.Part.Connect.Nodes;    %Total number of nodes per element

disp('> Parsing element node coordinates...')
elNodeCoord=elementNodeCoord(nodeCoord,connectInfo,elNodes,numElements);
%% Tipo de Análise - FEM ou TopOpt
choiceInt = questdlg('What kind of analysis will be made?','Integration scheme',...
                     'Finite Element Analysis', 'Topology Optimization','Topology Optimization');
waitfor(choiceInt)
switch choiceInt       
    case 'Finite Element Analysis'
        fem=1;
        choiceProblemType = questdlg('What kind of analysis will be made?','Integration scheme',...
                            "Elástico", "Termoelástico","Térmico","Térmico");
        waitfor(choiceProblemType)
        switch choiceProblemType       
            case "Elástico"
                elastico=1;   
            case "Termoelástico"
                termoelastico=1;
            case "Térmico"
                termico=1; 
        end
    case 'Topology Optimization'
        TopOpt=1; 
        choiceInt2 = questdlg('What kind of analysis will be made?','Integration scheme',...
                     'Linear Elasticity', 'Thermal Analysis','Multi-objective','Elasticity+DeltaT');
        waitfor(choiceInt2)
        switch choiceInt2       
            case 'Linear Elasticity'
                choiceInt3 = questdlg('What kind of analysis will be made?','Integration scheme',...
                     'Linear Elasticity', 'Thermoelasticity','Topology Optimization');
                waitfor(choiceInt3)
                switch choiceInt3
                     case "Linear Elasticity"
                        TopOptElastic=1;     
                     case "Thermoelasticity"
                        TopOptElasticDelta=1;
                end
            case 'Thermal Analysis'
                TopOptTherm=1; 
            case 'Multi-objective'
                TopOptMO=1; 
        end
end

%%%%%%%%%%%%%%%%% Problema elástico ou termoelástico %%%%%%%%%%%%%%%%%%%%%%
if ((fem==1)&&(elastico==1||termoelastico==1))||((TopOpt==1)&&(TopOptElastic==1))||((TopOpt==1)&&(TopOptMO==1))||((TopOpt==1)&&(TopOptElasticDelta==1))
%% Defining material properties
    dOf=size(elNodeCoord,2);
    disp('> Waiting for user to input material properties...')
    try
        %Prompts the user for input
        materialProps = materialPrompt();  
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;  
    end
    E=str2double(materialProps{1});
    v=str2double(materialProps{2});
    t=str2double(materialProps{3});
    alpha=[11*10^(-6);11*10^(-6);0];

%% Input analysis and finite element formulation
    disp('> Waiting for user to select analysis type...')
    try
        [D, choiceElement, intMode] = inputAnalysis(E,v);
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
%% Nodal shape functions and integration points
    disp('> Calculating nodal shape functions...')
    nodeShapeFun=shapeFunctions(elType,elNodes);

    disp('> Calculating integration points...')
    %Calculation integration points
    %(TODO - include code logic to deal with reduced integration schemes.)
    intPts=intPoints(elType,elNodes);
%% Boundary conditions
    %Launching GUI for the user to define BC's and loads
    [fixedNodes,sSnodes,loads]=nodeSelection(nodeCoord, connectInfo);
    %Creates the global mechanical load vector
    F=setNodalForces(loads,numNodes,dOf);
    
    if (termoelastico==1)
        deltaT=7*ones(1,4);
        %Creates the global thermal load vector
        [ elF ] = assemGlobalF( nodeShapeFun, intPts, elNodeCoord, dOf,elNodes, numNodes, connectInfo, D,alpha,deltaT, t);
        %Creates the global vector
        f=F+elF;
    elseif (elastico==1||TopOptElastic==1||TopOptMO==1)
        f=F;
    end
    %Creates the global degrees of freedom vector
    activeDoF=setFreeNodes(fixedNodes,sSnodes,numNodes,dOf);
end
if ((fem==1)&&(elastico==1||termoelastico==1))
%% Stiffness Matrix
    disp('> Assembling global stiffness matrix...')
    if strcmp(choiceElement,'PS/FI - Q4')|| strcmp(choiceElement,'PE/FI - Q4')

        K = assemGlobalStiff(nodeShapeFun,intPts,elNodeCoord,dOf,...
            elNodes,numNodes,connectInfo,D,t);

    elseif strcmp(choiceElement,'PE/SRI - Q4')

        K = assemGlobalStiffSRI(nodeShapeFun,intPts,elNodeCoord,dOf,...
            elNodes,numNodes,connectInfo,D,t);

    else
        [K,kAlpha,kAlphaD] = assemGlobalStiffEnhanced(nodeShapeFun,intPts,elNodeCoord,dOf,...
            elNodes,numNodes,connectInfo,D,t,choiceElement);
    end

%% Global System of Equations
    U=K(activeDoF,activeDoF)\f(activeDoF);
    U1=zeros(dOf*numNodes,1);
    U1(activeDoF)=U;

%% Displacements
    U=U1;
    u1=U(1:dOf:dOf*numNodes);
    u2=U(2:dOf:dOf*numNodes);

%% Reaction forces
    F=K*U;
    f1=F(1:dOf:dOf*numNodes);
    f1(abs(f1)<eps)=0;
    f2=F(2:dOf:dOf*numNodes);
    f2(abs(f2)<eps)=0;

%% Results 
    if elastico==1
        if strcmp(choiceElement,'PS/FI - Q4')|| strcmp(choiceElement,'PE/FI - Q4')...
                || strcmp(choiceElement,'PE/SRI - Q4')

            S=computeStress(nodeShapeFun,intPts,elNodeCoord,dOf,...
                elNodes,numNodes,connectInfo,D,U);
        else
            alphas=computeAlphas(dOf,connectInfo,U,choiceElement,kAlpha,kAlphaD);

            S=computeStressEnhanced(nodeShapeFun, intPts, elNodeCoord, dOf,...
                elNodes, numNodes, connectInfo, D, U, alphas, choiceElement);

        end

    elseif termoelastico==1

        [S]=computeStressDelta( nodeShapeFun, intPts, elNodeCoord, dOf,...
            elNodes, numNodes, connectInfo, D, U,deltaT,alpha);

    end

%% Pós - Processamento
    postProcessing (U, u1, u2, S, numNodes, connectInfo, nodeCoord);
elseif ((fem==1)&&(termico==1))
%%%%%%%%%%%%%%%%%%%%%%% Problema Térmico %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('> Waiting for user to input material properties...')
    try
        %Prompts the user for input
        [materialPropsThermal]=materialPromptThermal();
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    E=0;v=0;
    k=str2double(materialPropsThermal{1});
    t=str2double(materialPropsThermal{3});
    
    try
        [D, choiceElement, intMode] = inputAnalysis(E,v);
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    disp('> Calculating nodal shape functions...')
    nodeShapeFun=shapeFunctions(elType,elNodes);
    
    disp('> Calculating integration points...');
    %Calculation integration points
    intPts=intPoints(elType,elNodes);
    [Q,T1]=ThermalAnalysis(elNodeCoord,elNodes,connectInfo,numNodes,nodeShapeFun,intPts,k,t,nodeCoord);
    
%%%%%%%%%%%%%%%%%%%%%%%% TopOpt Mecânico %%%%%%%%%%%%%%%%%%%%%%%%%    
elseif (TopOpt==1)&&(TopOptElastic==1)
    
    TopOptMechanical(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,D,activeDoF,nodeShapeFun,intPts,f);
%%%%%%%%%%%%%%%%%%%%%%%% TopOpt Mecânico + Delta T %%%%%%%%%%%%%%%%%%%%%%%%%    
elseif (TopOpt==1)&&(TopOptElasticDelta==1)

    TopOptMecDelta(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,D,activeDoF,nodeShapeFun,intPts,F);
%%%%%%%%%%%%%%%%%%%%%%%% TopOpt Térmico %%%%%%%%%%%%%%%%%%%%%%%%%      
elseif (TopOpt==1)&&(TopOptTherm==1)
    
    disp('> Waiting for user to input material properties...')
    try
        %Prompts the user for input
        [materialPropsThermal]=materialPromptThermal();
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    E=0;v=0;
    k=str2double(materialPropsThermal{1});
    t=str2double(materialPropsThermal{3});

    try
        [D, choiceElement, intMode] = inputAnalysis(E,v);
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    disp('> Calculating nodal shape functions...')
    nodeShapeFun=shapeFunctions(elType,elNodes);
    
    disp('> Calculating integration points...');
    %Calculation integration points
    intPts=intPoints(elType,elNodes);
    
    TopOptThermal(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,k,nodeShapeFun,intPts);
%%%%%%%%%%%%%%%%%%%%%%%% TopOpt Multi-objective %%%%%%%%%%%%%%%%%     
elseif (TopOpt==1)&&(TopOptMO==1)
    disp('> Waiting for user to input material properties...')
    try
        %Prompts the user for input
        [materialPropsThermal]=materialPromptThermal();
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    k=str2double(materialPropsThermal{1});
    t=str2double(materialPropsThermal{3});

    try
        [D, choiceElement, intMode] = inputAnalysis(E,v);
    catch e
        errorMsg=sprintf('> %s',e.message);
        disp(errorMsg)
        return;
    end
    disp('> Calculating nodal shape functions...')
    nodeShapeFun=shapeFunctions(elType,elNodes);
    
    disp('> Calculating integration points...');
    intPts=intPoints(elType,elNodes);
    
    TopOptMultiObjective(nodeCoord,numNodes,connectInfo,numElements,elNodes,elNodeCoord,t,D,k,activeDoF,nodeShapeFun,intPts,f);
end 