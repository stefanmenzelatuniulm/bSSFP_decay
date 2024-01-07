classdef longitudinalPopulationNode < populationNode

    properties (Access = public)

        transverseChild emptyNode; 
        longitudinalChild emptyNode;

    end

    methods

        %Constructor
        function longitudinalPopulationNode = longitudinalPopulationNode(parent, transverseChild, longitudinalChild, label, xpos, ypos, Meq, dephasingDegree, amplitude, labelAmplitude)

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
                longitudinalPopulationNode.Meq = Meq;
                longitudinalPopulationNode.dephasingDegree = dephasingDegree;
                longitudinalPopulationNode.amplitude = amplitude;
                longitudinalPopulationNode.labelAmplitude = labelAmplitude;
                longitudinalPopulationNode.transverseChild = transverseChild;
                longitudinalPopulationNode.longitudinalChild = longitudinalChild;  

            else
                
                longitudinalPopulationNode.parent = emptyNode();
                longitudinalPopulationNode.label = "";
                longitudinalPopulationNode.level = 0;
                longitudinalPopulationNode.xpos = 0;
                longitudinalPopulationNode.ypos = 0;
                longitudinalPopulationNode.Meq = 1;
                longitudinalPopulationNode.dephasingDegree = 0;
                longitudinalPopulationNode.amplitude = sym(1);
                longitudinalPopulationNode.labelAmplitude = sym(1);
                longitudinalPopulationNode.transverseChild = emptyNode();
                longitudinalPopulationNode.longitudinalChild = emptyNode();         
            
            end
    
        end

        function [transverseBottomNodes, longitudinalBottomNodes, longitudinalPopulationNodeObject] = applyPulse(longitudinalPopulationNodeObject, a_, TR_, f, height, Meq, yScale)
            
            transverseBottomNodes = [];
            longitudinalBottomNodes = [];
            
            if longitudinalPopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree
 
                syms T1 T2 T2p TR a;
                E1 = exp(-f*TR/T1);
                F = Meq*(1-E1);
                E2 = exp(-f*TR/T2);
                
                %Not inverted phase
                oldDephasingDegree = transversePopulationNodeObject.dephasingDegree;
                dephasingDegreeNotInverted = oldDephasingDegree+f*TR;

                if oldDephasingDegree <= 0 && dephasingDegreeNotInverted <= 0

                    E2p = exp(f*TR/T2p);

                elseif oldDephasingDegree <= 0 && dephasingDegreeNotInverted >= 0

                    a = -TR*oldDephasingDegree/(dephasingDegreeNotInverted-oldDephasingDegree);
                    b = TR*dephasingDegreeNotInverted/(dephasingDegreeNotInverted-oldDephasingDegree);

                    E2p = exp(a*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeNotInverted>0

                    E2p = exp(-f*TR/T2p);

                end

                %Transverse child 
                if ~isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") %only change empty children
                    amplitude = subs(subs(longitudinalPopulationNodeObject.amplitude*1i*sind(a)*E2*E2p, TR, TR_), a, a_);
                    labelAmplitude = longitudinalPopulationNodeObject.labelAmplitude*1i*sind(a)*E2*E2p;
                    if height>0
                        newLabel = longitudinalPopulationNodeObject.label+"_1";
                    else
                        newLabel = "1";
                    end
                    longitudinalPopulationNodeObject.transverseChild = transversePopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.xpos+f*TR, yScale*dephasingDegree, Meq, dephasingDegree, amplitude, labelAmplitude);
                    transverseBottomNodes = [transverseBottomNodes, longitudinalPopulationNodeObject.transverseChild.label+"#"+string(longitudinalPopulationNodeObject.transverseChild.amplitude)+"#"+string(longitudinalPopulationNodeObject.transverseChild.dephasingDegree)];
                end

                %Longitudinal child 
                if ~isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") 
                    amplitude = subs(subs(cosd(a)*(longitudinalPopulationNodeObject.amplitude*E1+F), TR, TR_), a, a_);
                    labelAmplitude = cosd(a)*(longitudinalPopulationNodeObject.labelAmplitude*E1+F);
                    %Care: distinguish between function
                    %longitudinalPopulationNode and object
                    %longitudinalPopulationNodeObject
                    if height>0
                        newLabel = longitudinalPopulationNodeObject.label+"_0";
                    else
                        newLabel = "0";
                    end
                    longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.xpos+f*TR, yScale*oldDephasingDegree, Meq, oldDephasingDegree, amplitude, labelAmplitude);
                    longitudinalBottomNodes = [longitudinalBottomNodes, longitudinalPopulationNodeObject.longitudinalChild.label+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.amplitude)+"#"+string(longitudinalPopulationNodeObject.longitudinalChild.dephasingDegree)];
                end
         
            else

                [transverseBottomNodes1, longitudinalBottomNodes1, ~] = applyPulse(longitudinalPopulationNodeObject.transverseChild, a_, TR_, f, height, Meq, yScale);
                [transverseBottomNodes2, longitudinalBottomNodes2, ~] = applyPulse(longitudinalPopulationNodeObject.longitudinalChild, a_, TR_, f, height, Meq, yScale);
                
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
        function longitudinalPopulationNode = updateAmplitudeLabel(longitudinalPopulationNode, updateLabel, summedAmplitudes, newLabel)

            if isa(longitudinalPopulationNode.transverseChild, "populationNode") && longitudinalPopulationNode.transverseChild.label == updateLabel

                longitudinalPopulationNode.transverseChild.label = newLabel;
                longitudinalPopulationNode.transverseChild.amplitude = summedAmplitudes;

            elseif isa(longitudinalPopulationNode.longitudinalChild, "populationNode") && longitudinalPopulationNode.longitudinalChild.label == updateLabel

                longitudinalPopulationNode.longitudinalChild.label = newLabel;
                longitudinalPopulationNode.longitudinalChild.amplitude = summedAmplitudes;

            else
                    
                longitudinalPopulationNode = longitudinalPopulationNode.transverseChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, newLabel);
                longitudinalPopulationNode = longitudinalPopulationNode.longitudinalChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, newLabel);

            end

        end

        function plotNode(longitudinalPopulationNode, plotDigits)

            oldDigits = digits;
            digits(plotDigits);
            
            plot(longitudinalPopulationNode.xpos, longitudinalPopulationNode.ypos, '.', 'color', [0, 0.65, 0], 'MarkerSize', 2);
            text(longitudinalPopulationNode.xpos, longitudinalPopulationNode.ypos, longitudinalPopulationNode.label+" "+string(vpa(longitudinalPopulationNode.amplitude)), 'FontSize', 3);
            
            if isa(longitudinalPopulationNode.longitudinalChild, "populationNode")
                line([longitudinalPopulationNode.xpos, longitudinalPopulationNode.longitudinalChild.xpos], [longitudinalPopulationNode.ypos, longitudinalPopulationNode.longitudinalChild.ypos], 'color', [0.8, 0.8, 0.8]);
            end

            if isa(longitudinalPopulationNode.transverseChild, "populationNode")
                line([longitudinalPopulationNode.xpos, longitudinalPopulationNode.transverseChild.xpos], [longitudinalPopulationNode.ypos, longitudinalPopulationNode.transverseChild.ypos], 'color', [0.8, 0.8, 0.8]);
            end

            digits(oldDigits);

        end

    end
    
end