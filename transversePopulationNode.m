classdef transversePopulationNode < populationNode

    properties (Access = public)
        
        transverseChild1 emptyNode; %not inverted phase
        transverseChild2 emptyNode; %inverted phase
        longitudinalChild1 emptyNode; %not inverted phase storage
        longitudinalChild2 emptyNode; %inverted phase storage
    
    end

    methods

        %Constructor
        function transversePopulationNode = transversePopulationNode(parent, transverseChild1, transverseChild2, longitudinalChild1, longitudinalChild2, label, xpos, ypos, dephasingDegree, amplitude, amplitudeLabel)

            if nargin > 1

                transversePopulationNode.parent = parent;
                transversePopulationNode.label = label;

                if isa(parent, "populationNode")
                    transversePopulationNode.level = parent.level+1;
                else
                    transversePopulationNode.level = 0;
                end

                transversePopulationNode.xpos = xpos;
                transversePopulationNode.ypos = ypos;
                transversePopulationNode.dephasingDegree = dephasingDegree;
                transversePopulationNode.amplitude = amplitude;
                transversePopulationNode.amplitudeLabel = amplitudeLabel;
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
                transversePopulationNode.dephasingDegree = 0;
                transversePopulationNode.amplitude = sym(1);
                transversePopulationNode.amplitudeLabel = sym(1);
                transversePopulationNode.transverseChild1 = emptyNode();
                transversePopulationNode.longitudinalChild1 = emptyNode();   
                transversePopulationNode.transverseChild2 = emptyNode();
                transversePopulationNode.longitudinalChild2 = emptyNode(); 
            
            end

        end

        function [transverseBottomNodes, longitudinalBottomNodes, transversePopulationNodeObject] = applyPulse(transversePopulationNodeObject, a_, TR_, f, height, yScale)

            transverseBottomNodes = [];
            longitudinalBottomNodes = [];
            
            if transversePopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree

                syms T1 T2 T2p TR a;
                E1 = exp(-f*TR/T1);
                E2 = exp(-f*TR/T2);
                
                %Not inverted phase
                oldDephasingDegree = transversePopulationNodeObject.dephasingDegree;
                dephasingDegreeNotInverted = subs(oldDephasingDegree+f*TR, TR, TR_);

                if oldDephasingDegree <= 0 && dephasingDegreeNotInverted <= 0

                    E2pNotInverted = exp(f*TR/T2p);

                elseif oldDephasingDegree <= 0 && dephasingDegreeNotInverted >= 0

                    as = -oldDephasingDegree/(dephasingDegreeNotInverted-oldDephasingDegree);
                    b = dephasingDegreeNotInverted/(dephasingDegreeNotInverted-oldDephasingDegree);

                    E2pNotInverted = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeNotInverted>0

                    E2pNotInverted = exp(-f*TR/T2p);

                end

                %Transverse child 1
                if ~isa(transversePopulationNodeObject.transverseChild1, "populationNode") %only change empty children
                    amplitude = simplify(subs(subs(transversePopulationNodeObject.amplitude*E2*E2pNotInverted*cosd(a/2)^2, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify(transversePopulationNodeObject.amplitudeLabel*E2*E2pNotInverted*cosd(a/2)^2, "IgnoreAnalyticConstraints", true);
                    if height>0
                        newLabel = transversePopulationNodeObject.label+"_1";
                    else
                        newLabel = "1";
                    end
                    transversePopulationNodeObject.transverseChild1 = transversePopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*dephasingDegreeNotInverted, dephasingDegreeNotInverted, amplitude, amplitudeLabel);
                    transverseBottomNodes = [transverseBottomNodes, transversePopulationNodeObject.transverseChild1.label+"#"+string(transversePopulationNodeObject.transverseChild1.amplitude)+"#"+string(transversePopulationNodeObject.transverseChild1.amplitudeLabel)+"#"+string(transversePopulationNodeObject.transverseChild1.dephasingDegree)];
                end                 

                %Longitudinal child 1
                if ~isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")
                    amplitude = simplify(subs(subs((1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitude, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify((1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitudeLabel, "IgnoreAnalyticConstraints", true);
                    if height>0
                        newLabel = transversePopulationNodeObject.label+"_0";
                    else
                        newLabel = "0";
                    end
                    transversePopulationNodeObject.longitudinalChild1 = longitudinalPopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*oldDephasingDegree, oldDephasingDegree, amplitude, amplitudeLabel);
                    longitudinalBottomNodes = [longitudinalBottomNodes, transversePopulationNodeObject.longitudinalChild1.label+"#"+string(transversePopulationNodeObject.longitudinalChild1.amplitude)+"#"+string(transversePopulationNodeObject.longitudinalChild1.amplitudeLabel)+"#"+string(transversePopulationNodeObject.longitudinalChild1.dephasingDegree)];
                end               

                %Inverted phase
                oldDephasingDegree = -oldDephasingDegree;
                dephasingDegreeInverted = subs(oldDephasingDegree+f*TR, TR, TR_);

                if oldDephasingDegree <= 0 && dephasingDegreeInverted <= 0

                    E2pInverted = exp(f*TR/T2p);

                elseif oldDephasingDegree <= 0 && dephasingDegreeInverted >= 0

                    as = -oldDephasingDegree/(dephasingDegreeInverted-oldDephasingDegree);
                    b = dephasingDegreeInverted/(dephasingDegreeInverted-oldDephasingDegree);

                    E2pInverted = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldDephasingDegree>0 -> dephasingDegreeInverted>0

                    E2pInverted = exp(-f*TR/T2p);

                end

                %Transverse child 2
                if ~isa(transversePopulationNodeObject.transverseChild2, "populationNode")
                    amplitude = simplify(subs(subs(transversePopulationNodeObject.amplitude*E2*E2pInverted*sind(a/2)^2, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify(transversePopulationNodeObject.amplitudeLabel*E2*E2pInverted*sind(a/2)^2, "IgnoreAnalyticConstraints", true);
                    if height>0
                        newLabel = transversePopulationNodeObject.label+"_-1";
                    else
                        newLabel = "-1";
                    end
                    transversePopulationNodeObject.transverseChild2 = transversePopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*dephasingDegreeInverted, dephasingDegreeInverted, amplitude, amplitudeLabel);
                    transverseBottomNodes = [transverseBottomNodes, transversePopulationNodeObject.transverseChild2.label+"#"+string(transversePopulationNodeObject.transverseChild2.amplitude)+"#"+string(transversePopulationNodeObject.transverseChild2.amplitudeLabel)+"#"+string(transversePopulationNodeObject.transverseChild2.dephasingDegree)];
                end

                %Longitudinal child 2
                if ~isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")
                    amplitude = simplify(subs(subs(-(1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitude, TR, TR_), a, a_), "IgnoreAnalyticConstraints", true);
                    amplitudeLabel = simplify(-(1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitudeLabel, "IgnoreAnalyticConstraints", true);
                    if height>0
                        newLabel = transversePopulationNodeObject.label+"_0*";
                    else
                        newLabel = "_0*";
                    end
                    transversePopulationNodeObject.longitudinalChild2 = longitudinalPopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.xpos+subs(f*TR, TR, TR_), yScale*oldDephasingDegree, oldDephasingDegree, amplitude, amplitudeLabel);
                    longitudinalBottomNodes = [longitudinalBottomNodes, transversePopulationNodeObject.longitudinalChild2.label+"#"+string(transversePopulationNodeObject.longitudinalChild2.amplitude)+"#"+string(transversePopulationNodeObject.longitudinalChild2.amplitudeLabel)+"#"+string(transversePopulationNodeObject.longitudinalChild2.dephasingDegree)];
                end

            else

                [transverseBottomNodes1, longitudinalBottomNodes1, transversePopulationNodeObject.transverseChild1] = transversePopulationNodeObject.transverseChild1.applyPulse(a_, TR_, f, height, yScale);
                [transverseBottomNodes2, longitudinalBottomNodes2, transversePopulationNodeObject.longitudinalChild1] = transversePopulationNodeObject.longitudinalChild1.applyPulse(a_, TR_, f, height, yScale);
                [transverseBottomNodes3, longitudinalBottomNodes3, transversePopulationNodeObject.transverseChild2] = transversePopulationNodeObject.transverseChild2.applyPulse(a_, TR_, f, height, yScale);
                [transverseBottomNodes4, longitudinalBottomNodes4, transversePopulationNodeObject.longitudinalChild2] = transversePopulationNodeObject.longitudinalChild2.applyPulse(a_, TR_, f, height, yScale);

                transverseBottomNodes = cat(2, transverseBottomNodes1, transverseBottomNodes2, transverseBottomNodes3, transverseBottomNodes4);
                longitudinalBottomNodes = cat(2, longitudinalBottomNodes1, longitudinalBottomNodes2, longitudinalBottomNodes3, longitudinalBottomNodes4);

            end
        
        end

        %Prunes node with label
        function transversePopulationNodeObject = prune(transversePopulationNodeObject, label)

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode") && transversePopulationNodeObject.transverseChild1.label == label

                transversePopulationNodeObject.transverseChild1 = prunedNode();

            elseif isa(transversePopulationNodeObject.transverseChild2, "populationNode") && transversePopulationNodeObject.transverseChild2.label == label

                transversePopulationNodeObject.transverseChild2 = prunedNode();

            elseif isa(transversePopulationNodeObject.longitudinalChild1, "populationNode") && transversePopulationNodeObject.longitudinalChild1.label == label

                transversePopulationNodeObject.longitudinalChild1 = prunedNode();

            elseif isa(transversePopulationNodeObject.longitudinalChild2, "populationNode") && transversePopulationNodeObject.longitudinalChild2.label == label

                transversePopulationNodeObject.longitudinalChild2 = prunedNode();

            else

                transversePopulationNodeObject.transverseChild1 = transversePopulationNodeObject.transverseChild1.prune(label);
                transversePopulationNodeObject.transverseChild2 = transversePopulationNodeObject.transverseChild2.prune(label);
                transversePopulationNodeObject.longitudinalChild1 = transversePopulationNodeObject.longitudinalChild1.prune(label);
                transversePopulationNodeObject.longitudinalChild2 = transversePopulationNodeObject.longitudinalChild2.prune(label);

            end

        end

        function transversePopulationNodeObject = updateAmplitudeLabel(transversePopulationNodeObject, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel)

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode") && transversePopulationNodeObject.transverseChild1.label == updateLabel

                transversePopulationNodeObject.transverseChild1.label = newLabel;
                transversePopulationNodeObject.transverseChild1.amplitude = summedAmplitudes;
                transversePopulationNodeObject.transverseChild1.amplitudeLabel = summedAmplitudeLabels;

            elseif isa(transversePopulationNodeObject.longitudinalChild1, "populationNode") && transversePopulationNodeObject.longitudinalChild1.label == updateLabel

                transversePopulationNodeObject.longitudinalChild1.label = newLabel;
                transversePopulationNodeObject.longitudinalChild1.amplitude = summedAmplitudes;
                transversePopulationNodeObject.longitudinalChild1.amplitudeLabel = summedAmplitudeLabels;

            elseif isa(transversePopulationNodeObject.transverseChild2, "populationNode") && transversePopulationNodeObject.transverseChild2.label == updateLabel

                transversePopulationNodeObject.transverseChild2.label = newLabel;
                transversePopulationNodeObject.transverseChild2.amplitude = summedAmplitudes;
                transversePopulationNodeObject.transverseChild2.amplitudeLabel = summedAmplitudeLabels;

            elseif isa(transversePopulationNodeObject.longitudinalChild2, "populationNode") && transversePopulationNodeObject.longitudinalChild2.label == updateLabel

                transversePopulationNodeObject.longitudinalChild2.label = newLabel;
                transversePopulationNodeObject.longitudinalChild2.amplitude = summedAmplitudes;
                transversePopulationNodeObject.longitudinalChild2.amplitudeLabel = summedAmplitudeLabels;
             
            else
                    
                transversePopulationNodeObject.transverseChild1 = transversePopulationNodeObject.transverseChild1.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);
                transversePopulationNodeObject.longitudinalChild1 = transversePopulationNodeObject.longitudinalChild1.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);
                transversePopulationNodeObject.transverseChild2 = transversePopulationNodeObject.transverseChild2.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);
                transversePopulationNodeObject.longitudinalChild2 = transversePopulationNodeObject.longitudinalChild2.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);

            end

        end

        function plotPathway(transversePopulationNodeObject, TRnum, fnum)

            hold on;

            if isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")
                line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.longitudinalChild1.xpos], [transversePopulationNodeObject.ypos, transversePopulationNodeObject.longitudinalChild1.ypos], 'color', [0 0 0], 'LineStyle', ':', 'LineWidth', 0.5);           
                hold on;
                transversePopulationNodeObject.longitudinalChild1.plotPathway(TRnum, fnum);
            end

            if isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")
                line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.longitudinalChild2.xpos], [-transversePopulationNodeObject.ypos, transversePopulationNodeObject.longitudinalChild2.ypos], 'color', [0 0 0], 'LineStyle', ':', 'LineWidth', 0.5);
                if transversePopulationNodeObject.dephasingDegree + 0.001 > TRnum*fnum+TRnum*(transversePopulationNodeObject.level-1)
                    hold on;
                    line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.xpos], [transversePopulationNodeObject.ypos, -transversePopulationNodeObject.ypos], 'color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 0.5);               
                end
                hold on;
                transversePopulationNodeObject.longitudinalChild2.plotPathway(TRnum, fnum)
            end

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode")
                line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.transverseChild1.xpos], [transversePopulationNodeObject.ypos, transversePopulationNodeObject.transverseChild1.ypos], 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 0.5);
                hold on;
                transversePopulationNodeObject.transverseChild1.plotPathway(TRnum, fnum);
            end

            if isa(transversePopulationNodeObject.transverseChild2, "populationNode")
                line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.transverseChild2.xpos], [-transversePopulationNodeObject.ypos, transversePopulationNodeObject.transverseChild2.ypos], 'color', [0 0 0], 'LineStyle', '--', 'LineWidth', 0.5);
                if transversePopulationNodeObject.dephasingDegree + 0.001 > TRnum*fnum+TRnum*(transversePopulationNodeObject.level-1)
                    hold on;
                    line([transversePopulationNodeObject.xpos, transversePopulationNodeObject.xpos], [transversePopulationNodeObject.ypos, -transversePopulationNodeObject.ypos], 'color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.5);
                end
                hold on;
                transversePopulationNodeObject.transverseChild2.plotPathway(TRnum, fnum);
            end

            hold on;

        end

        function [plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees] = plotNode(transversePopulationNodeObject, textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees, height)

            hold on;
            
            if ismembertol(transversePopulationNodeObject.dephasingDegree, plottedLongitudinalDephasingDegrees, 0.001)
                c = [0.7969 0.0469 0.7500];
            else
                c = [0.6350 0.0780 0.1840];
            end

            plot(transversePopulationNodeObject.xpos, transversePopulationNodeObject.ypos, '.', 'color', c, 'MarkerSize', 20);
            plottedTransverseDephasingDegrees = [plottedTransverseDephasingDegrees, transversePopulationNodeObject.dephasingDegree];

            hold on;

            if isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")
                text(textOffsetPathwayX+transversePopulationNodeObject.longitudinalChild1.xpos, -textOffsetPathwayY+transversePopulationNodeObject.longitudinalChild1.ypos, strrep(strrep(transversePopulationNodeObject.longitudinalChild1.label, "_", ","), "&", "\&"), 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');
                text(textOffsetLabelX+(transversePopulationNodeObject.xpos+transversePopulationNodeObject.longitudinalChild1.xpos)/2, -textOffsetLabelY+transversePopulationNodeObject.ypos, string("$"+latex(simplify(transversePopulationNodeObject.longitudinalChild1.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');               
                hold on;
                [plottedTransverseDephasingDegreesChild, plottedLongitudinalDephasingDegreesChild] = transversePopulationNodeObject.longitudinalChild1.plotNode(textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees, height);
                plottedTransverseDephasingDegrees = [plottedTransverseDephasingDegrees, plottedTransverseDephasingDegreesChild];
                plottedLongitudinalDephasingDegrees = [plottedLongitudinalDephasingDegrees, plottedLongitudinalDephasingDegreesChild];
                hold on;
            end

            if isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")
                if transversePopulationNodeObject.longitudinalChild2.level ~= height
                    text(textOffsetLabelX+(transversePopulationNodeObject.xpos+transversePopulationNodeObject.longitudinalChild2.xpos)/2, -textOffsetLabelY-transversePopulationNodeObject.ypos, string("$"+latex(simplify(transversePopulationNodeObject.longitudinalChild2.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');
                    text(textOffsetPathwayX+transversePopulationNodeObject.longitudinalChild2.xpos, -textOffsetPathwayY+transversePopulationNodeObject.longitudinalChild2.ypos, strrep(strrep(transversePopulationNodeObject.longitudinalChild2.label, "_", ","), "&", "\&"), 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');
                else
                    text(textOffsetLabelX+(transversePopulationNodeObject.xpos+transversePopulationNodeObject.longitudinalChild2.xpos)/2, textOffsetLabelY-transversePopulationNodeObject.ypos, string("$"+latex(simplify(transversePopulationNodeObject.longitudinalChild2.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');
                    text(textOffsetPathwayX+transversePopulationNodeObject.longitudinalChild2.xpos, textOffsetPathwayY+transversePopulationNodeObject.longitudinalChild2.ypos, strrep(strrep(transversePopulationNodeObject.longitudinalChild2.label, "_", ","), "&", "\&"), 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex');
                end               
                hold on;
                [plottedTransverseDephasingDegreesChild, plottedLongitudinalDephasingDegreesChild] = transversePopulationNodeObject.longitudinalChild2.plotNode(textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees, height);
                plottedTransverseDephasingDegrees = [plottedTransverseDephasingDegrees, plottedTransverseDephasingDegreesChild];
                plottedLongitudinalDephasingDegrees = [plottedLongitudinalDephasingDegrees, plottedLongitudinalDephasingDegreesChild];
                hold on;
            end

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode")
                text(textOffsetPathwayX+transversePopulationNodeObject.transverseChild1.xpos, textOffsetPathwayY+transversePopulationNodeObject.transverseChild1.ypos, strrep(strrep(transversePopulationNodeObject.transverseChild1.label, "_", ","), "&", "\&"), 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex');
                text(textOffsetLabelX+(transversePopulationNodeObject.xpos+transversePopulationNodeObject.transverseChild1.xpos)/2, textOffsetLabelY+(transversePopulationNodeObject.ypos+transversePopulationNodeObject.transverseChild1.ypos)/2, string("$"+latex(simplify(transversePopulationNodeObject.transverseChild1.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex');               
                hold on;
                [plottedTransverseDephasingDegreesChild, plottedLongitudinalDephasingDegreesChild] = transversePopulationNodeObject.transverseChild1.plotNode(textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees, height);
                plottedTransverseDephasingDegrees = [plottedTransverseDephasingDegrees, plottedTransverseDephasingDegreesChild];
                plottedLongitudinalDephasingDegrees = [plottedLongitudinalDephasingDegrees, plottedLongitudinalDephasingDegreesChild];
                hold on;
            end

            if isa(transversePopulationNodeObject.transverseChild2, "populationNode")
                text(textOffsetPathwayX+transversePopulationNodeObject.transverseChild2.xpos, textOffsetPathwayY+transversePopulationNodeObject.transverseChild2.ypos, strrep(strrep(transversePopulationNodeObject.transverseChild2.label, "_", ","), "&", "\&"), 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex');
                text(textOffsetLabelX+(transversePopulationNodeObject.xpos+transversePopulationNodeObject.transverseChild2.xpos)/2, textOffsetLabelY+(-transversePopulationNodeObject.ypos+transversePopulationNodeObject.transverseChild2.ypos)/2, string("$"+latex(simplify(transversePopulationNodeObject.transverseChild2.amplitudeLabel, "IgnoreAnalyticConstraints", true))+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex');               
                hold on;
                [plottedTransverseDephasingDegreesChild, plottedLongitudinalDephasingDegreesChild] = transversePopulationNodeObject.transverseChild2.plotNode(textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees, height);
                plottedTransverseDephasingDegrees = [plottedTransverseDephasingDegrees, plottedTransverseDephasingDegreesChild];
                plottedLongitudinalDephasingDegrees = [plottedLongitudinalDephasingDegrees, plottedLongitudinalDephasingDegreesChild];
                hold on;
            end

            hold on;

        end


    end
    
end