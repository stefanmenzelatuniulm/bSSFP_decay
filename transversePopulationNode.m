classdef transversePopulationNode < populationNode

    properties (Access=private)
        
        data transversePopulation;
        transverseChild1 transversePopulationNode; %not inverted phase
        transverseChild2 transversePopulationNode; %inverted phase
        longitudinalChild1 longitudinalPopulationNode; %not inverted phase storage
        longitudinalChild2 longitudinalPopulationNode; %inverted phase storage
        label string;
        level int64;
    
    end

    methods

        %Constructor
        function transversePopulationNode = transversePopulationNode(parent,transverseChild1,transverseChild2,longitudinalChild1,longitudinalChild2,data,label)

            if nargin > 0
                transversePopulationNode = populationNode(parent);
                transversePopulationNode.transverseChild1 = transverseChild1;
                transversePopulationNode.transverseChild2 = transverseChild2;
                transversePopulationNode.longitudinalChild1 = longitudinalChild1;
                transversePopulationNode.longitudinalChild2 = longitudinalChild2;
                transversePopulationNode.data = data;
                transversePopulationNode.label = label;
                transversePopulationNode.level = parent.level+1;
            else
                transversePopulationNode = populationNode();
                transversePopulationNode.transverseChild1 = emptyPopulationNode();
                transversePopulationNode.transverseChild2 = emptyPopulationNode();
                transversePopulationNode.longitudinalChild1 = emptyPopulationNode();
                transversePopulationNode.longitudinalChild2 = emptyPopulationNode();
                transversePopulationNode.data = population();
                transversePopulationNode.label = "";
                transversePopulationNode.level=-1;
            end
    
        end

        function [transverseBottomNodes,longitudinalBottomNodes,transversePopulationNode] = applyPulse(transversePopulationNode,a,TR,f,height)

            transverseBottomNodes=[];
            longitudinalBottomNodes=[];
            
            if transversePopulationNode.level==height %Pulse only changes nodes at the bottom of the tree

                syms T2 T2p;
                
                %Not inverted phase
                
                oldDephasingDegree=transversePopulationNode.data.dephasingDegree;
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
                if isa(transversePopulationNode.transverseChild1,"transversePopulationNode") %do nothing for empty nodes
                    amplitude=transversePopulationNode.data.amplitude*E2*E2p*cosd(a/2)^2;
                    transversePopulationNode.transverseChild1=transversePopulationNode(transversePopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegreeNotInverted,amplitude),transversePopulationNode.label+"_1");
                    transverseBottomNodes=[transverseBottomNodes,transversePopulationNode.transverseChild1.label+"#"+string(transversePopulationNode.transverseChild1.amplitude)+"#"+string(transversePopulationNode.transverseChild1.dephasingDegree)];
                end

                %Longitudinal child 1 %not inverted phase storage
                if isa(transversePopulationNode.longitudinalChild1,"longitudinalPopulationNode")
                    amplitude=transversePopulationNode.data.amplitude*E2*E2p*sind(a)*1i/2;
                    transversePopulationNode.longitudinalChild1=longitudinalPopulationNode(transversePopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegreeNotInverted,amplitude),transversePopulationNode.label+"_0");
                    longitudinalBottomNodes=[longitudinalBottomNodes,transversePopulationNode.longitudinalChild1.label+"#"+string(transversePopulationNode.longitudinalChild1.amplitude)+"#"+string(transversePopulationNode.longitudinalChild1.dephasingDegree)];
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
                if isa(transversePopulationNode.transverseChild2,"transversePopulationNode")
                    amplitude=transversePopulationNode.data.amplitude*E2*E2p*sind(a/2)^2;
                    transversePopulationNode.transverseChild2=transversePopulationNode(transversePopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegreeInverted,amplitude),transversePopulationNode.label+"_-1");
                    transverseBottomNodes=[transverseBottomNodes,transversePopulationNode.transverseChild2.label+"#"+string(transversePopulationNode.transverseChild2.amplitude)+"#"+string(transversePopulationNode.transverseChild2.dephasingDegree)];
                end

                %Longitudinal child 2 %inverted phase storage
                if isa(transversePopulationNode.longitudinalChild2,"longitudinalPopulationNode")
                    amplitude=-transversePopulationNode.data.amplitude*E2*E2p*sind(a)*1i/2;
                    transversePopulationNode.longitudinalChild2=longitudinalPopulationNode(transversePopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegreeInverted,amplitude),transversePopulationNode.label+"_0*");
                    longitudinalBottomNodes=[longitudinalBottomNodes,transversePopulationNode.longitudinalChild2.label+"#"+string(transversePopulationNode.longitudinalChild2.amplitude)+"#"+string(transversePopulationNode.longitudinalChild2.dephasingDegree)];
                end

            else

                [transverseBottomNodes1,longitudinalBottomNodes1,~]=applyPulse(transversePopulationNode.transverseChild1,a,TR,f,height);
                [transverseBottomNodes2,longitudinalBottomNodes2,~]=applyPulse(transversePopulationNode.longitudinalChild1,a,TR,f,height);
                [transverseBottomNodes3,longitudinalBottomNodes3,~]=applyPulse(transversePopulationNode.transverseChild2,a,TR,f,height);
                [transverseBottomNodes4,longitudinalBottomNodes4,~]=applyPulse(transversePopulationNode.longitudinalChild2,a,TR,f,height);

                transverseBottomNodes=cat(2,transverseBottomNodes1,transverseBottomNodes2,transverseBottomNodes3,transverseBottomNodes4);
                longitudinalBottomNodes=cat(2,longitudinalBottomNodes1,longitudinalBottomNodes2,longitudinalBottomNodes3,longitudinalBottomNodes4);

            end
        
        end

        %Prunes node with label
        function transversePopulationNode = prune(transversePopulationNode,label)

            if transversePopulationNode.transverseChild1.label==label

                transversePopulationNode.transverseChild1=emptyPopulationNode();

            elseif transversePopulationNode.transverseChild2.label==label

                transversePopulationNode.transverseChild2=emptyPopulationNode();

            elseif transversePopulationNode.longitudinalChild1.label==label

                transversePopulationNode.longitudinalChild1=emptyPopulationNode();

            elseif transversePopulationNode.longitudinalChild2.label==label

                transversePopulationNode.longitudinalChild2=emptyPopulationNode();

            else

                transversePopulationNode = prune(transversePopulationNode.transverseChild1,label);
                transversePopulationNode = prune(transversePopulationNode.transverseChild2,label);
                transversePopulationNode = prune(transversePopulationNode.longitudinalChild1,label);
                transversePopulationNode = prune(transversePopulationNode.longitudinalChild2,label);

            end

        end

    end
    
end