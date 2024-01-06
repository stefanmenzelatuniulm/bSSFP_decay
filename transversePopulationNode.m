classdef transversePopulationNode < populationNode

    properties (Access=public)
        
        transverseChild1 emptyNode; %not inverted phase
        transverseChild2 emptyNode; %inverted phase
        longitudinalChild1 emptyNode; %not inverted phase storage
        longitudinalChild2 emptyNode; %inverted phase storage
    
    end

    methods

        %Constructor
        function transversePopulationNode = transversePopulationNode(parent,transverseChild1,transverseChild2,longitudinalChild1,longitudinalChild2,label,xpos,ypos,Meq,dephasingDegree,amplitude)

            if nargin > 1

                transversePopulationNode.parent = parent;
                transversePopulationNode.label = label;

                if isa(parent,"populationNode")
                    transversePopulationNode.level = parent.level+1;
                else
                    transversePopulationNode.level = 0;
                end

                transversePopulationNode.xpos = xpos;
                transversePopulationNode.ypos = ypos;
                transversePopulationNode.Meq = Meq;
                transversePopulationNode.dephasingDegree = dephasingDegree;
                transversePopulationNode.amplitude = amplitude;
                transversePopulationNode.transverseChild1 = transverseChild1;
                transversePopulationNode.longitudinalChild1 = longitudinalChild1;  
                transversePopulationNode.transverseChild2 = transverseChild2;
                transversePopulationNode.longitudinalChild2 = longitudinalChild2;  

            else
                
                transversePopulationNode.parent = emptyNode();
                transversePopulationNode.label = "";
                transversePopulationNode.level = 0;
                transversePopulationNode.xpos = 0;
                transversePopulationNode.ypos = 0;
                transversePopulationNode.Meq = 1;
                transversePopulationNode.dephasingDegree = 0;
                transversePopulationNode.amplitude = sym(1);
                transversePopulationNode.transverseChild1 = emptyNode();
                transversePopulationNode.longitudinalChild1 = emptyNode();   
                transversePopulationNode.transverseChild2 = emptyNode();
                transversePopulationNode.longitudinalChild2 = emptyNode(); 
            
            end

        end

        function [transverseBottomNodes,longitudinalBottomNodes,transversePopulationNodeObject] = applyPulse(transversePopulationNodeObject,a,TR,f,height)

            transverseBottomNodes=[];
            longitudinalBottomNodes=[];
            
            if transversePopulationNodeObject.level==height %Pulse only changes nodes at the bottom of the tree

                syms T2 T2p;
                
                %Not inverted phase
                
                oldDephasingDegree=transversePopulationNodeObject.dephasingDegree;
                dephasingDegreeNotInverted=oldDephasingDegree+f*TR;

                E2=exp(-f*TR/T2);

                if oldDephasingDegree<=0 && dephasingDegreeNotInverted<=0

                    E2p=exp(f*TR/T2p);

                elseif oldDephasingDegree<=0 && dephasingDegreeNotInverted>=0

                    a=-TR*oldDephasingDegree/(dephasingDegreeNotInverted-oldDephasingDegree);
                    b=TR*dephasingDegreeNotInverted/(dephasingDegreeNotInverted-oldDephasingDegree);

                    E2p=exp(a*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeNotInverted>0

                    E2p=exp(-f*TR/T2p);

                end

                %Transverse child 1 %not inverted phase
                if ~isa(transversePopulationNodeObject.transverseChild1,"populationNode") %only change empty children
                    amplitude=transversePopulationNodeObject.amplitude*E2*E2p*cosd(a/2)^2;
                    %Care: distinguish between function
                    %transversePopulationNode and object
                    %transversePopulationNodeObject
                    transversePopulationNodeObject.transverseChild1=transversePopulationNode(transversePopulationNodeObject,emptyNode(),emptyNode(),emptyNode(),emptyNode(),transversePopulationNodeObject.label+"1_",transversePopulationNodeObject.xpos+f*TR,transversePopulationNodeObject.ypos+0.5*f*TR,transversePopulationNodeObject.Meq,dephasingDegreeNotInverted,amplitude);
                    transverseBottomNodes=[transverseBottomNodes,transversePopulationNodeObject.transverseChild1.label+"#"+string(transversePopulationNodeObject.transverseChild1.amplitude)+"#"+string(transversePopulationNodeObject.transverseChild1.dephasingDegree)];
                end

                %Longitudinal child 1 %not inverted phase storage
                if ~isa(transversePopulationNodeObject.longitudinalChild1,"populationNode")
                    amplitude=transversePopulationNodeObject.amplitude*E2*E2p*sind(a)*1i/2;
                    transversePopulationNodeObject.longitudinalChild1=longitudinalPopulationNode(transversePopulationNodeObject,emptyNode(),emptyNode(),transversePopulationNodeObject.label+"0_",transversePopulationNodeObject.xpos+f*TR,transversePopulationNodeObject.ypos+0.5*f*TR,transversePopulationNodeObject.Meq,dephasingDegreeNotInverted,amplitude);
                    longitudinalBottomNodes=[longitudinalBottomNodes,transversePopulationNodeObject.longitudinalChild1.label+"#"+string(transversePopulationNodeObject.longitudinalChild1.amplitude)+"#"+string(transversePopulationNodeObject.longitudinalChild1.dephasingDegree)];
                end

                %Inverted phase

                oldDephasingDegree=-oldDephasingDegree;
                dephasingDegreeInverted=oldDephasingDegree+f*TR;

                E2=exp(-f*TR/T2);

                if oldDephasingDegree<=0 && dephasingDegreeInverted<=0

                    E2p=exp(f*TR/T2p);

                elseif oldDephasingDegree<=0 && dephasingDegreeInverted>=0

                    a=-TR*oldDephasingDegree/(dephasingDegreeInverted-oldDephasingDegree);
                    b=TR*dephasingDegreeInverted/(dephasingDegreeInverted-oldDephasingDegree);

                    E2p=exp(a*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeInverted>0

                    E2p=exp(-f*TR/T2p);

                end

                %Transverse child 2 %inverted phase
                if ~isa(transversePopulationNodeObject.transverseChild2,"populationNode")
                    amplitude=transversePopulationNodeObject.amplitude*E2*E2p*sind(a/2)^2;
                    transversePopulationNodeObject.transverseChild2=transversePopulationNode(transversePopulationNodeObject,emptyNode(),emptyNode(),emptyNode(),emptyNode(),transversePopulationNodeObject.label+"-1_",transversePopulationNodeObject.xpos+f*TR,transversePopulationNodeObject.ypos+0.5*f*TR,transversePopulationNodeObject.Meq,dephasingDegreeInverted,amplitude);
                    transverseBottomNodes=[transverseBottomNodes,transversePopulationNodeObject.transverseChild2.label+"#"+string(transversePopulationNodeObject.transverseChild2.amplitude)+"#"+string(transversePopulationNodeObject.transverseChild2.dephasingDegree)];
                end

                %Longitudinal child 2 %inverted phase storage
                if ~isa(transversePopulationNodeObject.longitudinalChild2,"populationNode")
                    amplitude=-transversePopulationNodeObject.amplitude*E2*E2p*sind(a)*1i/2;
                    transversePopulationNodeObject.longitudinalChild2=longitudinalPopulationNode(transversePopulationNodeObject,emptyNode(),emptyNode(),transversePopulationNodeObject.label+"0*_",transversePopulationNodeObject.xpos+f*TR,transversePopulationNodeObject.ypos+0.5*f*TR,transversePopulationNodeObject.Meq,dephasingDegreeInverted,amplitude);
                    longitudinalBottomNodes=[longitudinalBottomNodes,transversePopulationNodeObject.longitudinalChild2.label+"#"+string(transversePopulationNodeObject.longitudinalChild2.amplitude)+"#"+string(transversePopulationNodeObject.longitudinalChild2.dephasingDegree)];
                end

            else

                [transverseBottomNodes1,longitudinalBottomNodes1,~]=applyPulse(transversePopulationNodeObject.transverseChild1,a,TR,f,height);
                [transverseBottomNodes2,longitudinalBottomNodes2,~]=applyPulse(transversePopulationNodeObject.longitudinalChild1,a,TR,f,height);
                [transverseBottomNodes3,longitudinalBottomNodes3,~]=applyPulse(transversePopulationNodeObject.transverseChild2,a,TR,f,height);
                [transverseBottomNodes4,longitudinalBottomNodes4,~]=applyPulse(transversePopulationNodeObject.longitudinalChild2,a,TR,f,height);

                transverseBottomNodes=cat(2,transverseBottomNodes1,transverseBottomNodes2,transverseBottomNodes3,transverseBottomNodes4);
                longitudinalBottomNodes=cat(2,longitudinalBottomNodes1,longitudinalBottomNodes2,longitudinalBottomNodes3,longitudinalBottomNodes4);

            end
        
        end

        %Prunes node with label
        function transversePopulationNode = prune(transversePopulationNode,label)

            if isa(transversePopulationNode.transverseChild1,"populationNode") && transversePopulationNode.transverseChild1.label==label

                transversePopulationNode.transverseChild1=emptyNode();

            elseif isa(transversePopulationNode.transverseChild2,"populationNode") && transversePopulationNode.transverseChild2.label==label

                transversePopulationNode.transverseChild2=emptyNode();

            elseif isa(transversePopulationNode.longitudinalChild1,"populationNode") && transversePopulationNode.longitudinalChild1.label==label

                transversePopulationNode.longitudinalChild1=emptyNode();

            elseif isa(transversePopulationNode.longitudinalChild2,"populationNode") && transversePopulationNode.longitudinalChild2.label==label

                transversePopulationNode.longitudinalChild2=emptyNode();

            else

                transversePopulationNode = prune(transversePopulationNode.transverseChild1,label);
                transversePopulationNode = prune(transversePopulationNode.transverseChild2,label);
                transversePopulationNode = prune(transversePopulationNode.longitudinalChild1,label);
                transversePopulationNode = prune(transversePopulationNode.longitudinalChild2,label);

            end

        end

        function transversePopulationNode = updateAmplitudeLabel(transversePopulationNode,updateLabel,summedAmplitudes,newLabel)

            if isa(transversePopulationNode.transverseChild1,"populationNode") && transversePopulationNode.transverseChild1.label==updateLabel

                transversePopulationNode.transverseChild1.label=newLabel;
                transversePopulationNode.transverseChild1.amplitude=summedAmplitudes;

            elseif isa(transversePopulationNode.longitudinalChild1,"populationNode") && transversePopulationNode.longitudinalChild1.label==updateLabel

                transversePopulationNode.longitudinalChild1.label=newLabel;
                transversePopulationNode.longitudinalChild1.amplitude=summedAmplitudes;

            elseif isa(transversePopulationNode.transverseChild2,"populationNode") && transversePopulationNode.transverseChild2.label==updateLabel

                transversePopulationNode.transverseChild2.label=newLabel;
                transversePopulationNode.transverseChild2.amplitude=summedAmplitudes;

            elseif isa(transversePopulationNode.longitudinalChild2,"populationNode") && transversePopulationNode.longitudinalChild2.label==updateLabel

                transversePopulationNode.longitudinalChild2.label=newLabel;
                transversePopulationNode.longitudinalChild2.amplitude=summedAmplitudes;
             
            else
                    
                transversePopulationNode = transversePopulationNode.transverseChild1.updateAmplitudeLabel(updateLabel,summedAmplitudes,newLabel);
                transversePopulationNode = transversePopulationNode.longitudinalChild1.updateAmplitudeLabel(updateLabel,summedAmplitudes,newLabel);
                transversePopulationNode = transversePopulationNode.transverseChild2.updateAmplitudeLabel(updateLabel,summedAmplitudes,newLabel);
                transversePopulationNode = transversePopulationNode.longitudinalChild2.updateAmplitudeLabel(updateLabel,summedAmplitudes,newLabel);

            end

        end

        function plotNode(transversePopulationNode,plotDigits)

            oldDigits = digits;
            digits(plotDigits);
            
            plot(transversePopulationNode.xpos,transversePopulationNode.ypos,'.','color',[0,0.65,0],'MarkerSize',2);
            text(transversePopulationNode.xpos,transversePopulationNode.ypos,transversePopulationNode.label+" "+string(vpa(transversePopulationNode.amplitude)),'FontSize',3);
            
            if isa(transversePopulationNode.longitudinalChild1,"populationNode")
                line([transversePopulationNode.xpos,transversePopulationNode.longitudinalChild1.xpos],[transversePopulationNode.ypos,transversePopulationNode.longitudinalChild1.ypos],'color',[0.8,0.8,0.8]);
            end

            if isa(transversePopulationNode.longitudinalChild2,"populationNode")
                line([transversePopulationNode.xpos,transversePopulationNode.longitudinalChild2.xpos],[transversePopulationNode.ypos,transversePopulationNode.longitudinalChild2.ypos],'color',[0.8,0.8,0.8]);
            end

            if isa(transversePopulationNode.transverseChild1,"populationNode")
                line([transversePopulationNode.xpos,transversePopulationNode.transverseChild1.xpos],[transversePopulationNode.ypos,transversePopulationNode.transverseChild1.ypos],'color',[0.8,0.8,0.8]);
            end

            if isa(transversePopulationNode.transverseChild2,"populationNode")
                line([transversePopulationNode.xpos,transversePopulationNode.transverseChild2.xpos],[transversePopulationNode.ypos,transversePopulationNode.transverseChild2.ypos],'color',[0.8,0.8,0.8]);
            end

            digits(oldDigits);

        end

    end
    
end