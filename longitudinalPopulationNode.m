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
                syms M_eq;
                longitudinalPopulationNode.amplitude = M_eq;
                longitudinalPopulationNode.amplitudeLabel = M_eq;
                longitudinalPopulationNode.transverseChild = emptyNode();
                longitudinalPopulationNode.longitudinalChild = emptyNode();         
            
            end
    
        end

        function [transverseBottomNodes, longitudinalBottomNodes, longitudinalPopulationNodeObject] = applyPulse(longitudinalPopulationNodeObject, a_, TR_, f, height, yScale)
            
            transverseBottomNodes = [];
            longitudinalBottomNodes = [];
            
            if longitudinalPopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree
 
                syms T1 T2 T2p TR a M_eq;
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
                    amplitude = simplify(subs(subs(longitudinalPopulationNodeObject.amplitude*1i*sind(a)*E2*E2p, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify(longitudinalPopulationNodeObject.amplitudeLabel*1i*sind(a)*E2*E2p, "IgnoreAnalyticConstraints", true);
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
                    amplitude = simplify(subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitude-M_eq)*E1+M_eq, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify((cosd(a)*longitudinalPopulationNodeObject.amplitudeLabel-M_eq)*E1+M_eq, "IgnoreAnalyticConstraints", true);
                    if height>0
                        newLabel = longitudinalPopulationNodeObject.label+"_0";
                    else
                        newLabel = "0";
                    end
                    longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*oldDephasingDegree, oldDephasingDegree, amplitude, amplitudeLabel);
                    longitudinalBottomNodes = [longitudinalBottomNodes, longitudinalPopulationNodeObject.longitudinalChild.label+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.amplitude)+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel)+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.dephasingDegree)];
                end
         
            else

                [transverseBottomNodes1, longitudinalBottomNodes1, longitudinalPopulationNodeObject.transverseChild] = longitudinalPopulationNodeObject.transverseChild.applyPulse(a_, TR_, f, height, yScale);
                [transverseBottomNodes2, longitudinalBottomNodes2, longitudinalPopulationNodeObject.longitudinalChild] = longitudinalPopulationNodeObject.longitudinalChild.applyPulse(a_, TR_, f, height, yScale);
                
                transverseBottomNodes = cat(2, transverseBottomNodes1, transverseBottomNodes2);
                longitudinalBottomNodes = cat(2, longitudinalBottomNodes1, longitudinalBottomNodes2);

            end
        
        end

        %Prunes node with label
        function longitudinalPopulationNodeObject = prune(longitudinalPopulationNodeObject, label)

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") && longitudinalPopulationNodeObject.transverseChild.label == label

                longitudinalPopulationNodeObject.transverseChild = emptyNode();

            elseif isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") && longitudinalPopulationNodeObject.longitudinalChild.label == label

                longitudinalPopulationNodeObject.transverseChild = emptyNode();

            else

                longitudinalPopulationNodeObject.transverseChild = longitudinalPopulationNodeObject.transverseChild.prune(label);
                longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNodeObject.longitudinalChild.prune(label);

            end

        end
        
        %Updates node with label
        function longitudinalPopulationNodeObject = updateAmplitudeLabel(longitudinalPopulationNodeObject, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel)

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") && longitudinalPopulationNodeObject.transverseChild.label == updateLabel

                longitudinalPopulationNodeObject.transverseChild.label = newLabel;
                longitudinalPopulationNodeObject.transverseChild.amplitude = summedAmplitudes;
                longitudinalPopulationNodeObject.transverseChild.amplitudeLabel = summedAmplitudeLabels;

            elseif isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") && longitudinalPopulationNodeObject.longitudinalChild.label == updateLabel

                longitudinalPopulationNodeObject.longitudinalChild.label = newLabel;
                longitudinalPopulationNodeObject.longitudinalChild.amplitude = summedAmplitudes;
                longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel = summedAmplitudeLabels;

            else
                    
                longitudinalPopulationNodeObject.transverseChild = longitudinalPopulationNodeObject.transverseChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);
                longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNodeObject.longitudinalChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);

            end

        end

        function plotNode(longitudinalPopulationNodeObject, textOffsetX, textOffsetY)

            hold on;
            
            plot(longitudinalPopulationNodeObject.xpos, longitudinalPopulationNodeObject.ypos, '.', 'color', [0 0.4470 0.7410], 'MarkerSize', 20);

            hold on;

            if isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode")
                line([longitudinalPopulationNodeObject.xpos, longitudinalPopulationNodeObject.longitudinalChild.xpos], [longitudinalPopulationNodeObject.ypos, longitudinalPopulationNodeObject.longitudinalChild.ypos], 'color', [0.2 0.2 0.2], 'LineStyle', ':', 'LineWidth', 1);
                text(longitudinalPopulationNodeObject.longitudinalChild.xpos, textOffsetY+longitudinalPopulationNodeObject.longitudinalChild.ypos, strrep(strrep(longitudinalPopulationNodeObject.longitudinalChild.label, "_", "\_"), "&", "\&"), 'FontSize', 12, 'Color', [0.4940 0.1840 0.5560], 'Interpreter', 'latex'); %, 'BackgroundColor', [1 1 1]);
                text(-textOffsetX+(longitudinalPopulationNodeObject.xpos+longitudinalPopulationNodeObject.longitudinalChild.xpos)/2, (longitudinalPopulationNodeObject.ypos+longitudinalPopulationNodeObject.longitudinalChild.ypos)/2, string("$"+latex(simplify(longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980], 'Interpreter', 'latex'); %, 'BackgroundColor', [1 1 1]);               
                hold on;
                longitudinalPopulationNodeObject.longitudinalChild.plotNode(textOffsetX, textOffsetY);
                hold on;
            end

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode")
                line([longitudinalPopulationNodeObject.xpos, longitudinalPopulationNodeObject.transverseChild.xpos], [longitudinalPopulationNodeObject.ypos, longitudinalPopulationNodeObject.transverseChild.ypos], 'color', [0.2 0.2 0.2], 'LineStyle', '--', 'LineWidth', 1);
                text(longitudinalPopulationNodeObject.transverseChild.xpos, textOffsetY+longitudinalPopulationNodeObject.transverseChild.ypos, strrep(strrep(longitudinalPopulationNodeObject.transverseChild.label, "_", "\_"), "&", "\&"), 'FontSize', 12, 'Color', [0.4940 0.1840 0.5560], 'Interpreter', 'latex'); %, 'BackgroundColor', [1 1 1]);
                text(-textOffsetX+(longitudinalPopulationNodeObject.xpos+longitudinalPopulationNodeObject.transverseChild.xpos)/2, (longitudinalPopulationNodeObject.ypos+longitudinalPopulationNodeObject.transverseChild.ypos)/2, string("$"+latex(simplify(longitudinalPopulationNodeObject.transverseChild.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', 12, 'Color', [0.8500 0.3250 0.0980], 'Interpreter', 'latex'); %, 'BackgroundColor', [1 1 1]);               
                hold on;
                longitudinalPopulationNodeObject.transverseChild.plotNode(textOffsetX, textOffsetY);
                hold on;
            end

            hold on;

        end

    end
    
end