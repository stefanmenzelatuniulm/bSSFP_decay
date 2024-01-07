classdef longitudinalPopulationNode < populationNode

    properties (Access = public)

        transverseChild emptyNode; 
        longitudinalChild emptyNode;

    end

    methods

        %Constructor
        function longitudinalPopulationNode = longitudinalPopulationNode(parent, transverseChild, longitudinalChild, label, xpos, ypos, dephasingDegree, amplitude, amplitudeLabel)

            if nargin > 1

                longitudinalPopulationNode.parent = parent;
                longitudinalPopulationNode.label = label;

                if isa(parent, "populationNode")
                    longitudinalPopulationNode.level = parent.level+1;
                else
                    longitudinalPopulationNode.level = 0;
                end

                longitudinalPopulationNode.xpos = xpos;
                longitudinalPopulationNode.ypos = ypos;
                longitudinalPopulationNode.dephasingDegree = dephasingDegree;
                longitudinalPopulationNode.amplitude = amplitude;
                longitudinalPopulationNode.amplitudeLabel = amplitudeLabel;
                longitudinalPopulationNode.transverseChild = transverseChild;
                longitudinalPopulationNode.longitudinalChild = longitudinalChild;  

            else
                
                longitudinalPopulationNode.parent = emptyNode();
                longitudinalPopulationNode.label = "";
                longitudinalPopulationNode.level = 0;
                longitudinalPopulationNode.xpos = 0;
                longitudinalPopulationNode.ypos = 0;
                longitudinalPopulationNode.dephasingDegree = 0;
                syms Meq;
                longitudinalPopulationNode.amplitude = Meq;
                longitudinalPopulationNode.amplitudeLabel = Meq;
                longitudinalPopulationNode.transverseChild = emptyNode();
                longitudinalPopulationNode.longitudinalChild = emptyNode();         
            
            end
    
        end

        function [transverseBottomNodes, longitudinalBottomNodes, longitudinalPopulationNodeObject] = applyPulse(longitudinalPopulationNodeObject, a_, TR_, f, height, yScale)
            
            transverseBottomNodes = [];
            longitudinalBottomNodes = [];
            
            if longitudinalPopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree
 
                syms T1 T2 T2p TR a Meq;
                E1 = exp(-f*TR/T1);
                E2 = exp(-f*TR/T2);
                
                %Not inverted phase
                oldDephasingDegree = longitudinalPopulationNodeObject.dephasingDegree;
                dephasingDegreeNotInverted = subs(oldDephasingDegree+f*TR, TR, TR_);

                if oldDephasingDegree <= 0 && dephasingDegreeNotInverted <= 0

                    E2p = exp(f*TR/T2p);

                elseif oldDephasingDegree <= 0 && dephasingDegreeNotInverted >= 0

                    as = -oldDephasingDegree/(dephasingDegreeNotInverted-oldDephasingDegree);
                    b = dephasingDegreeNotInverted/(dephasingDegreeNotInverted-oldDephasingDegree);

                    E2p = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeNotInverted>0

                    E2p = exp(-f*TR/T2p);

                end

                %Transverse child 
                if ~isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") %only change empty children
                    amplitude = subs(subs(longitudinalPopulationNodeObject.amplitude*1i*sind(a)*E2*E2p, TR, TR_), a, a_);
                    amplitudeLabel = longitudinalPopulationNodeObject.amplitudeLabel*1i*sind(a)*E2*E2p;
                    if height>0
                        newLabel = longitudinalPopulationNodeObject.label+"_1";
                    else
                        newLabel = "1";
                    end
                    longitudinalPopulationNodeObject.transverseChild = transversePopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*dephasingDegreeNotInverted, dephasingDegreeNotInverted, amplitude, amplitudeLabel);
                    transverseBottomNodes = [transverseBottomNodes, longitudinalPopulationNodeObject.transverseChild.label+"#"+string(longitudinalPopulationNodeObject.transverseChild.amplitude)+"#"+string(longitudinalPopulationNodeObject.transverseChild.amplitudeLabel)+"#"+string(longitudinalPopulationNodeObject.transverseChild.dephasingDegree)];
                end

                %Longitudinal child
                if ~isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") 
                    amplitude = subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitude-Meq)*E1+Meq, TR, TR_), a, a_);
                    amplitudeLabel = (cosd(a)*longitudinalPopulationNodeObject.amplitudeLabel-Meq)*E1+Meq;
                    %Care: distinguish between function
                    %longitudinalPopulationNode and object
                    %longitudinalPopulationNodeObject
                    if height>0
                        newLabel = longitudinalPopulationNodeObject.label+"_0";
                    else
                        newLabel = "0";
                    end
                    longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*oldDephasingDegree, oldDephasingDegree, amplitude, amplitudeLabel);
                    longitudinalBottomNodes = [longitudinalBottomNodes, longitudinalPopulationNodeObject.longitudinalChild.label+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.amplitude)+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel)+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.dephasingDegree)];
                end
         
            else

                [transverseBottomNodes1, longitudinalBottomNodes1, ~] = applyPulse(longitudinalPopulationNodeObject.transverseChild, a_, TR_, f, height, yScale);
                [transverseBottomNodes2, longitudinalBottomNodes2, ~] = applyPulse(longitudinalPopulationNodeObject.longitudinalChild, a_, TR_, f, height, yScale);
                
                transverseBottomNodes = cat(2, transverseBottomNodes1, transverseBottomNodes2);
                longitudinalBottomNodes = cat(2, longitudinalBottomNodes1, longitudinalBottomNodes2);

            end
        
        end

        %Prunes node with label
        function longitudinalPopulationNode = prune(longitudinalPopulationNode, label)

            if isa(longitudinalPopulationNode.transverseChild, "populationNode") && longitudinalPopulationNode.transverseChild.label == label

                longitudinalPopulationNode.transverseChild = emptyNode();

            elseif isa(longitudinalPopulationNode.longitudinalChild, "populationNode") && longitudinalPopulationNode.longitudinalChild.label == label

                longitudinalPopulationNode.transverseChild = emptyNode();

            else

                longitudinalPopulationNode = prune(longitudinalPopulationNode.transverseChild, label);
                longitudinalPopulationNode = prune(longitudinalPopulationNode.longitudinalChild, label);

            end

        end
        
        %Updates node with label
        function longitudinalPopulationNode = updateAmplitudeLabel(longitudinalPopulationNode, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel)

            if isa(longitudinalPopulationNode.transverseChild, "populationNode") && longitudinalPopulationNode.transverseChild.label == updateLabel

                longitudinalPopulationNode.transverseChild.label = newLabel;
                longitudinalPopulationNode.transverseChild.amplitude = summedAmplitudes;
                longitudinalPopulationNode.transverseChild.amplitude = summedAmplitudeLabels;

            elseif isa(longitudinalPopulationNode.longitudinalChild, "populationNode") && longitudinalPopulationNode.longitudinalChild.label == updateLabel

                longitudinalPopulationNode.longitudinalChild.label = newLabel;
                longitudinalPopulationNode.longitudinalChild.amplitude = summedAmplitudes;
                longitudinalPopulationNode.longitudinalChild.amplitude = summedAmplitudeLabels;

            else
                    
                longitudinalPopulationNode = longitudinalPopulationNode.transverseChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);
                longitudinalPopulationNode = longitudinalPopulationNode.longitudinalChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);

            end

        end

        function plotNode(longitudinalPopulationNode)
            
            plot(longitudinalPopulationNode.xpos, longitudinalPopulationNode.ypos, '.', 'color', [0, 0.65, 0], 'MarkerSize', 2);
            text(longitudinalPopulationNode.xpos, longitudinalPopulationNode.ypos, longitudinalPopulationNode.label+newline+string(longitudinalPopulationNode.amplitudeLabel), 'FontSize', 3);
            
            if isa(longitudinalPopulationNode.longitudinalChild, "populationNode")
                line([longitudinalPopulationNode.xpos, longitudinalPopulationNode.longitudinalChild.xpos], [longitudinalPopulationNode.ypos, longitudinalPopulationNode.longitudinalChild.ypos], 'color', [0.8, 0.8, 0.8]);
            end

            if isa(longitudinalPopulationNode.transverseChild, "populationNode")
                line([longitudinalPopulationNode.xpos, longitudinalPopulationNode.transverseChild.xpos], [longitudinalPopulationNode.ypos, longitudinalPopulationNode.transverseChild.ypos], 'color', [0.8, 0.8, 0.8]);
            end

        end

    end
    
end