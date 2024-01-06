classdef longitudinalPopulationNode < populationNode

    properties (Access=private)
        
        data longitudinalPopulation;
        transverseChild transversePopulationNode; 
        longitudinalChild longitudinalPopulationNode;
        label string;
        level int64;
    
    end

    methods

        %Constructor
        function longitudinalPopulationNode = longitudinalPopulationNode(parent,transverseChild,longitudinalChild,data,label)

            if nargin > 0
                longitudinalPopulationNode = populationNode(parent);
                longitudinalPopulationNode.transverseChild = transverseChild;
                longitudinalPopulationNode.longitudinalChild = longitudinalChild;
                longitudinalPopulationNode.data = data;
                longitudinalPopulationNode.label = label;
                longitudinalPopulationNode.level = parent.level+1;
            else
                longitudinalPopulationNode = populationNode();
                longitudinalPopulationNode.transverseChild = emptyPopulationNode();
                longitudinalPopulationNode.longitudinalChild = emptyPopulationNode();
                longitudinalPopulationNode.data = population();
                longitudinalPopulationNode.label = "";
                longitudinalPopulationNode.level = -1;
            end
    
        end

        function [transverseBottomNodes,longitudinalBottomNodes,longitudinalPopulationNode] = applyPulse(longitudinalPopulationNode,a,TR,f,height)
            
            transverseBottomNodes=[];
            longitudinalBottomNodes=[];
            
            if longitudinalPopulationNode.level==height %Pulse only changes nodes at the bottom of the tree
 
                syms T1 Meq;
                E1=exp(-f*TR/T1);
                F=Meq*(1-E1);
                
                oldDephasingDegree=longitudinalPopulationNode.data.dephasingDegree;
                dephasingDegree=oldDephasingDegree;

                %Transverse child 
                if isa(longitudinalPopulationNode.transverseChild,"transversePopulationNode") %do nothing for empty nodes
                    amplitude=1i*sind(a)*(longitudinalPopulationNode.data.amplitude*E1+F);
                    longitudinalPopulationNode.transverseChild=transversePopulationNode(longitudinalPopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegree,amplitude),longitudinalPopulationNode.label+"_1");
                    transverseBottomNodes=[transverseBottomNodes,longitudinalPopulationNode.transverseChild.label+"#"+string(longitudinalPopulationNode.transverseChild.amplitude)+"#"+string(longitudinalPopulationNode.transverseChild.dephasingDegree)];
                end

                %Longitudinal child 
                if isa(longitudinalPopulationNode.longitudinalChild,"longitudinalPopulationNode") 
                    amplitude=cosd(a)*(longitudinalPopulationNode.data.amplitude*E1+F);
                    longitudinalPopulationNode.longitudinalChild=longitudinalPopulationNode(longitudinalPopulationNode,emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),emptyPopulationNode(),population(dephasingDegree,amplitude),longitudinalPopulationNode.label+"_0");
                    longitudinalBottomNodes=[longitudinalBottomNodes,longitudinalPopulationNode.longitudinalChild.label+"#"+string(longitudinalPopulationNode.longitudinalChild.amplitude)+"#"+string(longitudinalPopulationNode.longitudinalChild.dephasingDegree)];
                end
         
            else

                [transverseBottomNodes1,longitudinalBottomNodes1,~]=applyPulse(longitudinalPopulationNode.transverseChild,a,TR,f,height);
                [transverseBottomNodes2,longitudinalBottomNodes2,~]=applyPulse(longitudinalPopulationNode.longitudinalChild,a,TR,f,height);
                
                transverseBottomNodes=cat(2,transverseBottomNodes1,transverseBottomNodes2);
                longitudinalBottomNodes=cat(2,longitudinalBottomNodes1,longitudinalBottomNodes2);

            end
        
        end

        %Prunes node with label
        function longitudinalPopulationNode = prune(longitudinalPopulationNode,label)

            if longitudinalPopulationNode.transverseChild.label==label

                longitudinalPopulationNode.transverseChild=emptyPopulationNode();

            elseif longitudinalPopulationNode.longitudinalChild.label==label

                longitudinalPopulationNode.transverseChild=emptyPopulationNode();

            else

                longitudinalPopulationNode = prune(longitudinalPopulationNode.transverseChild,label);
                longitudinalPopulationNode = prune(longitudinalPopulationNode.longitudinalChild,label);

            end

        end

    end
    
end