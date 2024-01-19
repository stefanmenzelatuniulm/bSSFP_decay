classdef longitudinalPopulationNode < populationNode

    properties (Access = public)

        transverseChild emptyNode; 
        longitudinalChild emptyNode;

    end

    methods

        %Constructor
        function longitudinalPopulationNode = longitudinalPopulationNode(parent, transverseChild, longitudinalChild, label, totalTime, coherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulse, dephasingTimeDirectlyAfterPulse, dephasingTime)

            if nargin > 1

                longitudinalPopulationNode.parent = parent;
                longitudinalPopulationNode.label = label;

                if isa(parent, "populationNode")
                    longitudinalPopulationNode.level = parent.level+1;
                else
                    longitudinalPopulationNode.level = 0;
                end

                longitudinalPopulationNode.totalTime = totalTime;
                longitudinalPopulationNode.coherenceDegree = coherenceDegree;
                longitudinalPopulationNode.amplitude = amplitude;
                longitudinalPopulationNode.amplitudeDirectlyAfterPulse = amplitudeDirectlyAfterPulse;
                longitudinalPopulationNode.amplitudeWithoutT2p = amplitudeWithoutT2p;
                longitudinalPopulationNode.amplitudeDirectlyAfterPulseWithoutT2p = amplitudeDirectlyAfterPulseWithoutT2p;
                longitudinalPopulationNode.amplitudeLabel = amplitudeLabel;
                longitudinalPopulationNode.transverseChild = transverseChild;
                longitudinalPopulationNode.longitudinalChild = longitudinalChild;
                longitudinalPopulationNode.coherenceDegreeDirectlyAfterPulse = coherenceDegreeDirectlyAfterPulse;
                longitudinalPopulationNode.dephasingTimeDirectlyAfterPulse = dephasingTimeDirectlyAfterPulse;
                longitudinalPopulationNode.dephasingTime = dephasingTime;

            else
                
                longitudinalPopulationNode.parent = emptyNode();
                longitudinalPopulationNode.label = "";
                longitudinalPopulationNode.level = 0;
                longitudinalPopulationNode.totalTime = 0;
                longitudinalPopulationNode.coherenceDegree = 0;
                syms M_eq real;
                longitudinalPopulationNode.amplitude = M_eq;
                longitudinalPopulationNode.amplitudeDirectlyAfterPulse = sym(1);
                longitudinalPopulationNode.amplitudeLabel = M_eq;
                longitudinalPopulationNode.amplitudeWithoutT2p = sym(1);
                longitudinalPopulationNode.amplitudeDirectlyAfterPulseWithoutT2p = sym(1);
                longitudinalPopulationNode.transverseChild = emptyNode();
                longitudinalPopulationNode.longitudinalChild = emptyNode();    
                longitudinalPopulationNode.coherenceDegreeDirectlyAfterPulse = 0;
                longitudinalPopulationNode.dephasingTimeDirectlyAfterPulse = 0;
                longitudinalPopulationNode.dephasingTime = 0;
            
            end
    
        end

        function [transverseBottomNodes, longitudinalBottomNodes, longitudinalPopulationNodeObject] = applyPulse(longitudinalPopulationNodeObject, a_, TR_, f, height)

            transverseBottomNodes = [];
            longitudinalBottomNodes = [];

            if longitudinalPopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree

                if height == 0
                    aFactor = 0.5;
                else
                    aFactor = (-1)^height;
                end
 
                syms T1 T2 T2p TR a M_eq w real;
                E1 = exp(-f*TR/T1);
                E2 = exp(-f*TR/T2);
                dephasing = exp(1i*2*pi*w*TR*f);
                                
                %Not inverted coherenceDegree
                w0 = 1;
                oldCoherenceDegree = longitudinalPopulationNodeObject.coherenceDegree;
                coherenceDegreeNotInverted = subs(oldCoherenceDegree+f*TR*w0*360, TR, TR_);
                coherenceDegreeDirectlyAfterPulseNotInverted = longitudinalPopulationNodeObject.coherenceDegree;

                if oldCoherenceDegree <= 0 && coherenceDegreeNotInverted <= 0

                    E2p = exp(f*TR/T2p);

                elseif oldCoherenceDegree <= 0 && coherenceDegreeNotInverted >= 0

                    as = -oldCoherenceDegree/(coherenceDegreeNotInverted-oldCoherenceDegree);
                    b = coherenceDegreeNotInverted/(coherenceDegreeNotInverted-oldCoherenceDegree);

                    E2p = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldCoherenceDegree>0 -> coherenceDegreeNotInverted>0

                    E2p = exp(-f*TR/T2p);

                end

                %Transverse child 
                if ~isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") %only change empty children
                    if isa(longitudinalPopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+longitudinalPopulationNodeObject.label+" -> 1");
                    end

                    amplitudeLabel = dephasing*longitudinalPopulationNodeObject.amplitudeLabel*1i*sind(aFactor*a)*E2*E2p;
                    amplitude = subs(subs(dephasing*longitudinalPopulationNodeObject.amplitude*1i*sind(a)*E2*E2p, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulse = subs(subs(longitudinalPopulationNodeObject.amplitude*1i*sind(a), TR, TR_), a, a_);
                    amplitudeWithoutT2p = subs(subs(dephasing*longitudinalPopulationNodeObject.amplitudeWithoutT2p*1i*sind(a)*E2, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulseWithoutT2p = subs(subs(longitudinalPopulationNodeObject.amplitudeWithoutT2p*1i*sind(a), TR, TR_), a, a_);
                    dephasingTimeDirectlyAfterPulse = longitudinalPopulationNodeObject.dephasingTime;
                    dephasingTime = longitudinalPopulationNodeObject.dephasingTime+TR_*f;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(longitudinalPopulationNodeObject.label, "+"));
                        for k = 1:length(pathways)
                            pathways(k) = pathways(k)+", 1";
                            if k~=1
                                newLabel = newLabel+" + "+pathways(k);
                            else
                                newLabel = pathways(k);
                            end
                        end
                    else
                        newLabel = "1";
                    end
                    longitudinalPopulationNodeObject.transverseChild = transversePopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.totalTime+subs(f*TR, TR, TR_), coherenceDegreeNotInverted, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseNotInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if longitudinalPopulationNodeObject.transverseChild.level == height+1

                        transverseBottomNodes = [transverseBottomNodes, longitudinalPopulationNodeObject.transverseChild];

                    end

                end

                %Longitudinal child
                if ~isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") 
                    if isa(longitudinalPopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+longitudinalPopulationNodeObject.label+" -> 0");
                    end

                    %CARE: Longitudinal magnetization with phase memory decays to 0, not to Meq!
                    if contains(longitudinalPopulationNodeObject.label,"1") || contains(longitudinalPopulationNodeObject.label,"0*") || contains(longitudinalPopulationNodeObject.label,"-1")
                        amplitudeLabel = cosd(aFactor*a)*longitudinalPopulationNodeObject.amplitudeLabel*E1;
                        amplitude = subs(subs(cosd(a)*longitudinalPopulationNodeObject.amplitude*E1, TR, TR_), a, a_);
                        amplitudeDirectlyAfterPulse = subs(subs(cosd(a)*longitudinalPopulationNodeObject.amplitude, TR, TR_), a, a_);
                        amplitudeWithoutT2p = subs(subs(cosd(a)*longitudinalPopulationNodeObject.amplitudeWithoutT2p*E1, TR, TR_), a, a_);
                        amplitudeDirectlyAfterPulseWithoutT2p = subs(subs(cosd(a)*longitudinalPopulationNodeObject.amplitudeWithoutT2p, TR, TR_), a, a_);
                    else 
                        amplitudeLabel = (cosd(aFactor*a)*longitudinalPopulationNodeObject.amplitudeLabel-M_eq)*E1+M_eq;
                        amplitude = subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitude-M_eq)*E1+M_eq, TR, TR_), a, a_);
                        amplitudeDirectlyAfterPulse = subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitude-M_eq)+M_eq, TR, TR_), a, a_);
                        amplitudeWithoutT2p = subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitudeWithoutT2p-M_eq)*E1+M_eq, TR, TR_), a, a_);
                        amplitudeDirectlyAfterPulseWithoutT2p = subs(subs((cosd(a)*longitudinalPopulationNodeObject.amplitudeWithoutT2p-M_eq)+M_eq, TR, TR_), a, a_);
                    end
                     
                    dephasingTimeDirectlyAfterPulse = longitudinalPopulationNodeObject.dephasingTime;
                    dephasingTime = longitudinalPopulationNodeObject.dephasingTime;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(longitudinalPopulationNodeObject.label, "+"));
                        for k = 1:length(pathways)
                            pathways(k) = pathways(k)+", 0";
                            if k~=1
                                newLabel = newLabel+" + "+pathways(k);
                            else
                                newLabel = pathways(k);
                            end
                        end
                    else
                        newLabel = "0";
                    end
                    longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNode(longitudinalPopulationNodeObject, emptyNode(), emptyNode(), newLabel, longitudinalPopulationNodeObject.totalTime+subs(f*TR, TR, TR_), oldCoherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseNotInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if longitudinalPopulationNodeObject.longitudinalChild.level == height+1
                    
                        longitudinalBottomNodes = [longitudinalBottomNodes, longitudinalPopulationNodeObject.longitudinalChild];

                    end

                end
         
            else

                if isa(longitudinalPopulationNodeObject.transverseChild, "transversePopulationNode")
                    [transverseBottomNodes1, longitudinalBottomNodes1, longitudinalPopulationNodeObject.transverseChild] = longitudinalPopulationNodeObject.transverseChild.applyPulse(a_, TR_, f, height);
                else
                    transverseBottomNodes1 = [];
                    longitudinalBottomNodes1 = [];
                end

                 if isa(longitudinalPopulationNodeObject.longitudinalChild, "longitudinalPopulationNode")               
                    [transverseBottomNodes2, longitudinalBottomNodes2, longitudinalPopulationNodeObject.longitudinalChild] = longitudinalPopulationNodeObject.longitudinalChild.applyPulse(a_, TR_, f, height);
                 else
                    transverseBottomNodes2 = [];
                    longitudinalBottomNodes2 = [];
                 end

                transverseBottomNodes = cat(2, transverseBottomNodes1, transverseBottomNodes2);
                longitudinalBottomNodes = cat(2, longitudinalBottomNodes1, longitudinalBottomNodes2);

            end
        
        end

        %Prunes node with label
        function longitudinalPopulationNodeObject = prune(longitudinalPopulationNodeObject, label)

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") && longitudinalPopulationNodeObject.transverseChild.label == label
    
                disp("Pruning "+longitudinalPopulationNodeObject.transverseChild.label);
                longitudinalPopulationNodeObject.transverseChild = prunedNode();

            elseif isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") && longitudinalPopulationNodeObject.longitudinalChild.label == label

                disp("Pruning "+longitudinalPopulationNodeObject.longitudinalChild.label);
                longitudinalPopulationNodeObject.longitudinalChild = prunedNode();

            else

                longitudinalPopulationNodeObject.transverseChild = longitudinalPopulationNodeObject.transverseChild.prune(label);
                longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNodeObject.longitudinalChild.prune(label);

            end

        end
        
        %Updates node with label
        function longitudinalPopulationNodeObject = updateAmplitudeLabel(longitudinalPopulationNodeObject, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p)

            if updateLabel == newLabel
                return;
            end

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode") && longitudinalPopulationNodeObject.transverseChild.label == updateLabel
                
                disp("Updating "+updateLabel+" with "+newLabel);
                longitudinalPopulationNodeObject.transverseChild.label = newLabel;
                longitudinalPopulationNodeObject.transverseChild.amplitude = summedAmplitudes;
                longitudinalPopulationNodeObject.transverseChild.amplitudeLabel = summedAmplitudeLabels;
                longitudinalPopulationNodeObject.transverseChild.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                longitudinalPopulationNodeObject.transverseChild.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                longitudinalPopulationNodeObject.transverseChild.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;

            elseif isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode") && longitudinalPopulationNodeObject.longitudinalChild.label == updateLabel

                disp("Updating "+updateLabel+" with "+newLabel);
                longitudinalPopulationNodeObject.longitudinalChild.label = newLabel;
                longitudinalPopulationNodeObject.longitudinalChild.amplitude = summedAmplitudes;
                longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel = summedAmplitudeLabels;
                longitudinalPopulationNodeObject.longitudinalChild.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                longitudinalPopulationNodeObject.longitudinalChild.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                longitudinalPopulationNodeObject.longitudinalChild.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;

            else
                    
                longitudinalPopulationNodeObject.transverseChild = longitudinalPopulationNodeObject.transverseChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);
                longitudinalPopulationNodeObject.longitudinalChild = longitudinalPopulationNodeObject.longitudinalChild.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);

            end

        end

        function plotPathway(longitudinalPopulationNodeObject, TRnum, fnum, maxNodeLevel)

            hold on;

            if longitudinalPopulationNodeObject.level > maxNodeLevel
                return;
            end

            if isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode")
                line([longitudinalPopulationNodeObject.totalTime, longitudinalPopulationNodeObject.longitudinalChild.totalTime], [longitudinalPopulationNodeObject.coherenceDegree, longitudinalPopulationNodeObject.longitudinalChild.coherenceDegree], 'color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 0.5);           
                hold on;
                longitudinalPopulationNodeObject.longitudinalChild.plotPathway(TRnum, fnum, maxNodeLevel);
            end

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode")
                line([longitudinalPopulationNodeObject.totalTime, longitudinalPopulationNodeObject.transverseChild.totalTime], [longitudinalPopulationNodeObject.coherenceDegree, longitudinalPopulationNodeObject.transverseChild.coherenceDegree], 'color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 0.5);
                hold on;
                longitudinalPopulationNodeObject.transverseChild.plotPathway(TRnum, fnum, maxNodeLevel);
            end

            hold on;

        end

        function [plottedTransverseNodes, plottedLongitudinalNodes] = plotNode(longitudinalPopulationNodeObject, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel)
            
            if longitudinalPopulationNodeObject.level > maxNodeLevel+1
                return;
            end

            hold on;
            
            if coherenceDegreeIsInNodeList(plottedTransverseNodes, longitudinalPopulationNodeObject.coherenceDegree, longitudinalPopulationNodeObject.level)
                c = [0.8196 0 1.0000];
            else
                c = [0 0.4470 0.7410];
            end

            plot(longitudinalPopulationNodeObject.totalTime, longitudinalPopulationNodeObject.coherenceDegree, '.', 'color', c, 'MarkerSize', 15);
            plottedLongitudinalNodes = [plottedLongitudinalNodes, longitudinalPopulationNodeObject];

            hold on;

            if isa(longitudinalPopulationNodeObject.longitudinalChild, "populationNode")

                if coherenceDegreeIsInNodeList(plottedTransverseNodes, longitudinalPopulationNodeObject.longitudinalChild.coherenceDegree, longitudinalPopulationNodeObject.longitudinalChild.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(longitudinalPopulationNodeObject.longitudinalChild.amplitudeLabel);
                pathwayString = longitudinalPopulationNodeObject.longitudinalChild.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end

                if longitudinalPopulationNodeObject.longitudinalChild.level == 1 || abs(f-0.5)>labelOverlapThreshold
                    alignTextOnPathway(text((longitudinalPopulationNodeObject.totalTime+longitudinalPopulationNodeObject.longitudinalChild.totalTime)/2, longitudinalPopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, true, false); 
                    alignTextOnPathway(text(longitudinalPopulationNodeObject.longitudinalChild.totalTime, longitudinalPopulationNodeObject.longitudinalChild.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, false, noOverlap); 
                else
                    alignTextOnPathway(text((longitudinalPopulationNodeObject.totalTime+longitudinalPopulationNodeObject.longitudinalChild.totalTime)/2, longitudinalPopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, true, false); 
                    alignTextOnPathway(text(longitudinalPopulationNodeObject.longitudinalChild.totalTime, longitudinalPopulationNodeObject.longitudinalChild.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, false, noOverlap); 
                end  

                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = longitudinalPopulationNodeObject.longitudinalChild.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

            if isa(longitudinalPopulationNodeObject.transverseChild, "populationNode")

                if coherenceDegreeIsInNodeList(plottedLongitudinalNodes, longitudinalPopulationNodeObject.transverseChild.coherenceDegree, longitudinalPopulationNodeObject.transverseChild.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(longitudinalPopulationNodeObject.transverseChild.amplitudeLabel);
                pathwayString = longitudinalPopulationNodeObject.transverseChild.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end                

                if ~contains(longitudinalPopulationNodeObject.label, "0") && ~contains(longitudinalPopulationNodeObject.label, "-1")
                    alignTextOnPathway(text((longitudinalPopulationNodeObject.totalTime+longitudinalPopulationNodeObject.transverseChild.totalTime)/2, (longitudinalPopulationNodeObject.coherenceDegree+longitudinalPopulationNodeObject.transverseChild.coherenceDegree)/2, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, true, true);
                    alignTextOnPathway(text(longitudinalPopulationNodeObject.transverseChild.totalTime, longitudinalPopulationNodeObject.transverseChild.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, false, noOverlap);                     
                else
                    alignTextOnPathway(text((longitudinalPopulationNodeObject.totalTime+longitudinalPopulationNodeObject.transverseChild.totalTime)/2, (longitudinalPopulationNodeObject.coherenceDegree+longitudinalPopulationNodeObject.transverseChild.coherenceDegree)/2, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, true, false);
                    alignTextOnPathway(text(longitudinalPopulationNodeObject.transverseChild.totalTime, longitudinalPopulationNodeObject.transverseChild.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, false, noOverlap); 
                end          

                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = longitudinalPopulationNodeObject.transverseChild.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

        end

    end
    
end