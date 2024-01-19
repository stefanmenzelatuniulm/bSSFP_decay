classdef transversePopulationNode < populationNode

    properties (Access = public)
        
        transverseChild1 emptyNode; %not inverted coherenceDegree
        transverseChild2 emptyNode; %inverted coherenceDegree
        longitudinalChild1 emptyNode; %not inverted coherenceDegree storage
        longitudinalChild2 emptyNode; %inverted coherenceDegree storage
    
    end

    methods

        %Constructor
        function transversePopulationNode = transversePopulationNode(parent, transverseChild1, transverseChild2, longitudinalChild1, longitudinalChild2, label, totalTime, coherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulse, dephasingTimeDirectlyAfterPulse, dephasingTime)

            if nargin > 1

                transversePopulationNode.parent = parent;
                transversePopulationNode.label = label;

                if isa(parent, "populationNode")
                    transversePopulationNode.level = parent.level+1;
                else
                    transversePopulationNode.level = 0;
                end

                transversePopulationNode.totalTime = totalTime;
                transversePopulationNode.coherenceDegree = coherenceDegree;
                transversePopulationNode.amplitude = amplitude;
                transversePopulationNode.amplitudeDirectlyAfterPulse = amplitudeDirectlyAfterPulse;
                transversePopulationNode.amplitudeDirectlyAfterPulseWithoutT2p = amplitudeDirectlyAfterPulseWithoutT2p;
                transversePopulationNode.amplitudeWithoutT2p = amplitudeWithoutT2p;
                transversePopulationNode.amplitudeLabel = amplitudeLabel;
                transversePopulationNode.transverseChild1 = transverseChild1;
                transversePopulationNode.longitudinalChild1 = longitudinalChild1;  
                transversePopulationNode.transverseChild2 = transverseChild2;
                transversePopulationNode.longitudinalChild2 = longitudinalChild2;
                transversePopulationNode.coherenceDegreeDirectlyAfterPulse = coherenceDegreeDirectlyAfterPulse;
                transversePopulationNode.dephasingTimeDirectlyAfterPulse = dephasingTimeDirectlyAfterPulse;
                transversePopulationNode.dephasingTime = dephasingTime;

            else
                
                transversePopulationNode.parent = emptyNode();
                transversePopulationNode.label = "";
                transversePopulationNode.level = 0;
                transversePopulationNode.totalTime = 0;
                transversePopulationNode.coherenceDegree = 0;
                transversePopulationNode.amplitude = sym(1);
                transversePopulationNode.amplitudeDirectlyAfterPulse = sym(1);
                transversePopulationNode.amplitudeDirectlyAfterPulseWithoutT2p = sym(1);
                transversePopulationNode.amplitudeWithoutT2p = sym(1);
                transversePopulationNode.amplitudeLabel = sym(1);
                transversePopulationNode.transverseChild1 = emptyNode();
                transversePopulationNode.longitudinalChild1 = emptyNode();   
                transversePopulationNode.transverseChild2 = emptyNode();
                transversePopulationNode.longitudinalChild2 = emptyNode(); 
                transversePopulationNode.coherenceDegreeDirectlyAfterPulse = 0;
                transversePopulationNode.dephasingTimeDirectlyAfterPulse = 0;
                transversePopulationNode.dephasingTime = 0;
            
            end

        end

        function [transverseBottomNodes, longitudinalBottomNodes, transversePopulationNodeObject] = applyPulse(transversePopulationNodeObject, a_, TR_, f, height)

            transverseBottomNodes = [];
            longitudinalBottomNodes = [];
            
            if transversePopulationNodeObject.level == height %Pulse only changes nodes at the bottom of the tree

                if height == 0
                    aFactor = 0.5;
                else
                    aFactor = (-1)^height;
                end



                syms T1 T2 T2p TR a w real;
                E1 = exp(-f*TR/T1);
                E2 = exp(-f*TR/T2);
                dephasing = exp(1i*2*pi*w*TR*f);
                
                %Not inverted coherenceDegree
                w0 = 1;
                oldCoherenceDegree = transversePopulationNodeObject.coherenceDegree;
                coherenceDegreeNotInverted = subs(oldCoherenceDegree+f*TR*w0*360, TR, TR_);
                coherenceDegreeDirectlyAfterPulseNotInverted = transversePopulationNodeObject.coherenceDegree;

                if oldCoherenceDegree <= 0 && coherenceDegreeNotInverted <= 0

                    E2pNotInverted = exp(f*TR/T2p);

                elseif oldCoherenceDegree <= 0 && coherenceDegreeNotInverted >= 0

                    as = -oldCoherenceDegree/(coherenceDegreeNotInverted-oldCoherenceDegree);
                    b = coherenceDegreeNotInverted/(coherenceDegreeNotInverted-oldCoherenceDegree);

                    E2pNotInverted = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldCoherenceDegree>0 -> coherenceDegreeNotInverted>0

                    E2pNotInverted = exp(-f*TR/T2p);

                end

                %Transverse child 1
                if ~isa(transversePopulationNodeObject.transverseChild1, "populationNode") %only change empty children
                    if isa(transversePopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+transversePopulationNodeObject.label+" -> 1");
                    end

                    amplitudeLabel = dephasing*transversePopulationNodeObject.amplitudeLabel*E2*E2pNotInverted*cosd(aFactor*a/2)^2;       
                    amplitude = subs(subs(dephasing*transversePopulationNodeObject.amplitude*E2*E2pNotInverted*cosd(a/2)^2, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulse = subs(subs(transversePopulationNodeObject.amplitude*cosd(a/2)^2, TR, TR_), a, a_);
                    amplitudeWithoutT2p = subs(subs(dephasing*transversePopulationNodeObject.amplitudeWithoutT2p*E2*cosd(a/2)^2, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulseWithoutT2p = subs(subs(transversePopulationNodeObject.amplitudeWithoutT2p*cosd(a/2)^2, TR, TR_), a, a_);
                    dephasingTimeDirectlyAfterPulse = transversePopulationNodeObject.dephasingTime;
                    dephasingTime = transversePopulationNodeObject.dephasingTime+TR_*f;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(transversePopulationNodeObject.label, "+"));
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
                    transversePopulationNodeObject.transverseChild1 = transversePopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.totalTime+subs(f*TR, TR, TR_), coherenceDegreeNotInverted, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseNotInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if transversePopulationNodeObject.transverseChild1.level == height+1

                        transverseBottomNodes = [transverseBottomNodes, transversePopulationNodeObject.transverseChild1];

                    end

                end                 

                %Longitudinal child 1
                if ~isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")
                    if isa(transversePopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+transversePopulationNodeObject.label+" -> 0");
                    end 

                    amplitudeLabel = (1i/2)*sind(aFactor*a)*E1*transversePopulationNodeObject.amplitudeLabel;
                    amplitude = subs(subs((1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitude, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulse = subs(subs((1i/2)*sind(a)*transversePopulationNodeObject.amplitude, TR, TR_), a, a_);
                    amplitudeWithoutT2p = subs(subs((1i/2)*sind(a)*E1*transversePopulationNodeObject.amplitudeWithoutT2p, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulseWithoutT2p = subs(subs((1i/2)*sind(a)*transversePopulationNodeObject.amplitudeWithoutT2p, TR, TR_), a, a_);
                    dephasingTimeDirectlyAfterPulse = transversePopulationNodeObject.dephasingTime;
                    dephasingTime = transversePopulationNodeObject.dephasingTime;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(transversePopulationNodeObject.label, "+"));
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
                    transversePopulationNodeObject.longitudinalChild1 = longitudinalPopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.totalTime+subs(f*TR, TR, TR_), oldCoherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseNotInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if transversePopulationNodeObject.longitudinalChild1.level == height+1

                        longitudinalBottomNodes = [longitudinalBottomNodes, transversePopulationNodeObject.longitudinalChild1];

                    end

                end               

                %Inverted coherenceDegree
                w0 = 1;
                coherenceDegreeDirectlyAfterPulseInverted = -oldCoherenceDegree;
                oldCoherenceDegree = -oldCoherenceDegree;
                coherenceDegreeInverted = subs(oldCoherenceDegree+f*TR*w0*360, TR, TR_);

                if oldCoherenceDegree <= 0 && coherenceDegreeInverted <= 0

                    E2pInverted = exp(f*TR/T2p);

                elseif oldCoherenceDegree <= 0 && coherenceDegreeInverted >= 0

                    as = -oldCoherenceDegree/(coherenceDegreeInverted-oldCoherenceDegree);
                    b = coherenceDegreeInverted/(coherenceDegreeInverted-oldCoherenceDegree);

                    E2pInverted = exp(as*f*TR/T2p)*exp(-b*f*TR/T2p);

                else %oldCoherenceDegree>0 -> coherenceDegreeInverted>0

                    E2pInverted = exp(-f*TR/T2p);

                end

                %Transverse child 2
                if ~isa(transversePopulationNodeObject.transverseChild2, "populationNode")
                    if isa(transversePopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+transversePopulationNodeObject.label+" -> -1");
                    end  

                    amplitudeLabel = dephasing*conj(transversePopulationNodeObject.amplitudeLabel)*E2*E2pInverted*sind(aFactor*a/2)^2;                   
                    amplitude = subs(subs(dephasing*conj(transversePopulationNodeObject.amplitude)*E2*E2pInverted*sind(a/2)^2, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulse = subs(subs(conj(transversePopulationNodeObject.amplitude)*sind(a/2)^2, TR, TR_), a, a_);
                    amplitudeWithoutT2p = subs(subs(dephasing*conj(transversePopulationNodeObject.amplitudeWithoutT2p)*E2*sind(a/2)^2, TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulseWithoutT2p = subs(subs(conj(transversePopulationNodeObject.amplitudeWithoutT2p)*sind(a/2)^2, TR, TR_), a, a_);
                    dephasingTimeDirectlyAfterPulse = -transversePopulationNodeObject.dephasingTime;
                    dephasingTime = -transversePopulationNodeObject.dephasingTime+TR_*f;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(transversePopulationNodeObject.label, "+"));
                        for k = 1:length(pathways)
                            pathways(k) = pathways(k)+", -1";
                            if k~=1
                                newLabel = newLabel+" + "+pathways(k);
                            else
                                newLabel = pathways(k);
                            end
                        end
                    else
                        newLabel = "-1";
                    end
                    transversePopulationNodeObject.transverseChild2 = transversePopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.totalTime+subs(f*TR, TR, TR_), coherenceDegreeInverted, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if transversePopulationNodeObject.transverseChild2.level == height+1

                        transverseBottomNodes = [transverseBottomNodes, transversePopulationNodeObject.transverseChild2];

                    end

                end

                %Longitudinal child 2
                if ~isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")
                    if isa(transversePopulationNodeObject, "populationNode")
                    	disp("Applying pulse "+transversePopulationNodeObject.label+" -> 0*");
                    end   

                    amplitudeLabel = -(1i/2)*sind(aFactor*a)*E1*conj(transversePopulationNodeObject.amplitudeLabel);                    
                    amplitude = subs(subs(-(1i/2)*sind(a)*E1*conj(transversePopulationNodeObject.amplitude), TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulse = subs(subs(-(1i/2)*sind(a)*conj(transversePopulationNodeObject.amplitude), TR, TR_), a, a_);  
                    amplitudeWithoutT2p = subs(subs(-(1i/2)*sind(a)*E1*conj(transversePopulationNodeObject.amplitudeWithoutT2p), TR, TR_), a, a_);
                    amplitudeDirectlyAfterPulseWithoutT2p = subs(subs(-(1i/2)*sind(a)*conj(transversePopulationNodeObject.amplitudeWithoutT2p), TR, TR_), a, a_);  
                    dephasingTimeDirectlyAfterPulse = -transversePopulationNodeObject.dephasingTime;
                    dephasingTime = -transversePopulationNodeObject.dephasingTime;

                    if height>0
                        newLabel = "";
                        pathways = strtrim(strsplit(transversePopulationNodeObject.label, "+"));
                        for k = 1:length(pathways)
                            pathways(k) = pathways(k)+", 0*";
                            if k~=1
                                newLabel = newLabel+" + "+pathways(k);
                            else
                                newLabel = pathways(k);
                            end
                        end
                    else
                        newLabel = "0*";
                    end
                    transversePopulationNodeObject.longitudinalChild2 = longitudinalPopulationNode(transversePopulationNodeObject, emptyNode(), emptyNode(), newLabel, transversePopulationNodeObject.totalTime+subs(f*TR, TR, TR_), oldCoherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulseInverted, dephasingTimeDirectlyAfterPulse, dephasingTime);
                    
                    if transversePopulationNodeObject.longitudinalChild2.level == height+1

                        longitudinalBottomNodes = [longitudinalBottomNodes, transversePopulationNodeObject.longitudinalChild2];

                    end

                end

            else
                
                if isa(transversePopulationNodeObject.transverseChild1, "transversePopulationNode")
                    [transverseBottomNodes1, longitudinalBottomNodes1, transversePopulationNodeObject.transverseChild1] = transversePopulationNodeObject.transverseChild1.applyPulse(a_, TR_, f, height);
                else
                    transverseBottomNodes1 = [];
                    longitudinalBottomNodes1 = [];
                end

                if isa(transversePopulationNodeObject.longitudinalChild1, "longitudinalPopulationNode")
                    [transverseBottomNodes2, longitudinalBottomNodes2, transversePopulationNodeObject.longitudinalChild1] = transversePopulationNodeObject.longitudinalChild1.applyPulse(a_, TR_, f, height);
                else
                    transverseBottomNodes2 = [];
                    longitudinalBottomNodes2 = [];
                end

                if isa(transversePopulationNodeObject.transverseChild2, "transversePopulationNode")
                    [transverseBottomNodes3, longitudinalBottomNodes3, transversePopulationNodeObject.transverseChild2] = transversePopulationNodeObject.transverseChild2.applyPulse(a_, TR_, f, height);
                else
                    transverseBottomNodes3 = [];
                    longitudinalBottomNodes3 = [];
                end                

                if isa(transversePopulationNodeObject.longitudinalChild2, "longitudinalPopulationNode")
                    [transverseBottomNodes4, longitudinalBottomNodes4, transversePopulationNodeObject.longitudinalChild2] = transversePopulationNodeObject.longitudinalChild2.applyPulse(a_, TR_, f, height);
                else
                    transverseBottomNodes4 = [];
                    longitudinalBottomNodes4 = [];
                end

                transverseBottomNodes = cat(2, transverseBottomNodes1, transverseBottomNodes2, transverseBottomNodes3, transverseBottomNodes4);
                longitudinalBottomNodes = cat(2, longitudinalBottomNodes1, longitudinalBottomNodes2, longitudinalBottomNodes3, longitudinalBottomNodes4);

            end
        
        end

        %Prunes node with label
        function transversePopulationNodeObject = prune(transversePopulationNodeObject, label)

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode") && transversePopulationNodeObject.transverseChild1.label == label

                disp("Pruning "+transversePopulationNodeObject.transverseChild1.label);
                transversePopulationNodeObject.transverseChild1 = prunedNode();

            elseif isa(transversePopulationNodeObject.transverseChild2, "populationNode") && transversePopulationNodeObject.transverseChild2.label == label
                
                disp("Pruning "+transversePopulationNodeObject.transverseChild2.label);
                transversePopulationNodeObject.transverseChild2 = prunedNode();

            elseif isa(transversePopulationNodeObject.longitudinalChild1, "populationNode") && transversePopulationNodeObject.longitudinalChild1.label == label

                disp("Pruning "+transversePopulationNodeObject.longitudinalChild1.label);
                transversePopulationNodeObject.longitudinalChild1 = prunedNode();

            elseif isa(transversePopulationNodeObject.longitudinalChild2, "populationNode") && transversePopulationNodeObject.longitudinalChild2.label == label

                disp("Pruning "+transversePopulationNodeObject.longitudinalChild2.label);                
                transversePopulationNodeObject.longitudinalChild2 = prunedNode();

            else

                transversePopulationNodeObject.transverseChild1 = transversePopulationNodeObject.transverseChild1.prune(label);
                transversePopulationNodeObject.transverseChild2 = transversePopulationNodeObject.transverseChild2.prune(label);
                transversePopulationNodeObject.longitudinalChild1 = transversePopulationNodeObject.longitudinalChild1.prune(label);
                transversePopulationNodeObject.longitudinalChild2 = transversePopulationNodeObject.longitudinalChild2.prune(label);

            end

        end

        %Updates node with label
        function transversePopulationNodeObject = updateAmplitudeLabel(transversePopulationNodeObject, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p)

            if updateLabel == newLabel
                return;
            end            

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode") && transversePopulationNodeObject.transverseChild1.label == updateLabel

                disp("Updating "+updateLabel+" with "+newLabel);
                transversePopulationNodeObject.transverseChild1.label = newLabel;
                transversePopulationNodeObject.transverseChild1.amplitude = summedAmplitudes;
                transversePopulationNodeObject.transverseChild1.amplitudeLabel = summedAmplitudeLabels;
                transversePopulationNodeObject.transverseChild1.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                transversePopulationNodeObject.transverseChild1.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                transversePopulationNodeObject.transverseChild1.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;

            elseif isa(transversePopulationNodeObject.longitudinalChild1, "populationNode") && transversePopulationNodeObject.longitudinalChild1.label == updateLabel

                disp("Updating "+updateLabel+" with "+newLabel);
                transversePopulationNodeObject.longitudinalChild1.label = newLabel;
                transversePopulationNodeObject.longitudinalChild1.amplitude = summedAmplitudes;
                transversePopulationNodeObject.longitudinalChild1.amplitudeLabel = summedAmplitudeLabels;
                transversePopulationNodeObject.longitudinalChild1.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                transversePopulationNodeObject.longitudinalChild1.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                transversePopulationNodeObject.longitudinalChild1.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;                

            elseif isa(transversePopulationNodeObject.transverseChild2, "populationNode") && transversePopulationNodeObject.transverseChild2.label == updateLabel

                disp("Updating "+updateLabel+" with "+newLabel);
                transversePopulationNodeObject.transverseChild2.label = newLabel;
                transversePopulationNodeObject.transverseChild2.amplitude = summedAmplitudes;
                transversePopulationNodeObject.transverseChild2.amplitudeLabel = summedAmplitudeLabels;
                transversePopulationNodeObject.transverseChild2.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                transversePopulationNodeObject.transverseChild2.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                transversePopulationNodeObject.transverseChild2.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;                

            elseif isa(transversePopulationNodeObject.longitudinalChild2, "populationNode") && transversePopulationNodeObject.longitudinalChild2.label == updateLabel

                disp("Updating "+updateLabel+" with "+newLabel);
                transversePopulationNodeObject.longitudinalChild2.label = newLabel;
                transversePopulationNodeObject.longitudinalChild2.amplitude = summedAmplitudes;
                transversePopulationNodeObject.longitudinalChild2.amplitudeLabel = summedAmplitudeLabels;
                transversePopulationNodeObject.longitudinalChild2.amplitudeDirectlyAfterPulse = summedAmplitudesDirectlyAfterPulse;
                transversePopulationNodeObject.longitudinalChild2.amplitudeWithoutT2p = summedAmplitudesWithoutT2p;
                transversePopulationNodeObject.longitudinalChild2.amplitudeDirectlyAfterPulseWithoutT2p = summedAmplitudesDirectlyAfterPulseWithoutT2p;                
             
            else
                    
                transversePopulationNodeObject.transverseChild1 = transversePopulationNodeObject.transverseChild1.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);
                transversePopulationNodeObject.longitudinalChild1 = transversePopulationNodeObject.longitudinalChild1.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);
                transversePopulationNodeObject.transverseChild2 = transversePopulationNodeObject.transverseChild2.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);
                transversePopulationNodeObject.longitudinalChild2 = transversePopulationNodeObject.longitudinalChild2.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel, summedAmplitudesDirectlyAfterPulse, summedAmplitudesWithoutT2p, summedAmplitudesDirectlyAfterPulseWithoutT2p);

            end

        end

        function plotPathway(transversePopulationNodeObject, TRnum, fnum, maxNodeLevel)

            if transversePopulationNodeObject.level > maxNodeLevel
                return;
            end

            hold on;

            w0 = 1;

            if isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")
                line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.longitudinalChild1.totalTime], [transversePopulationNodeObject.coherenceDegree, transversePopulationNodeObject.longitudinalChild1.coherenceDegree], 'color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 0.5);           
                hold on;
                transversePopulationNodeObject.longitudinalChild1.plotPathway(TRnum, fnum, maxNodeLevel);
            end

            if isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")
                line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.longitudinalChild2.totalTime], [-transversePopulationNodeObject.coherenceDegree, transversePopulationNodeObject.longitudinalChild2.coherenceDegree], 'color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 0.5);
                if transversePopulationNodeObject.coherenceDegree + 0.001 > TRnum*fnum*w0*360+TRnum*(transversePopulationNodeObject.level-1)*w0*360
                    hold on;
                    line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.totalTime], [transversePopulationNodeObject.coherenceDegree, -transversePopulationNodeObject.coherenceDegree], 'color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 0.5);               
                end
                hold on;
                transversePopulationNodeObject.longitudinalChild2.plotPathway(TRnum, fnum, maxNodeLevel)
            end

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode")
                line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.transverseChild1.totalTime], [transversePopulationNodeObject.coherenceDegree, transversePopulationNodeObject.transverseChild1.coherenceDegree], 'color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 0.5);
                hold on;
                transversePopulationNodeObject.transverseChild1.plotPathway(TRnum, fnum, maxNodeLevel);
            end

            if isa(transversePopulationNodeObject.transverseChild2, "populationNode")
                line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.transverseChild2.totalTime], [-transversePopulationNodeObject.coherenceDegree, transversePopulationNodeObject.transverseChild2.coherenceDegree], 'color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 0.5);
                if transversePopulationNodeObject.coherenceDegree + 0.001 > TRnum*fnum*w0*360+TRnum*(transversePopulationNodeObject.level-1)*w0*360
                    hold on;
                    line([transversePopulationNodeObject.totalTime, transversePopulationNodeObject.totalTime], [transversePopulationNodeObject.coherenceDegree, -transversePopulationNodeObject.coherenceDegree], 'color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 0.5);
                end
                hold on;
                transversePopulationNodeObject.transverseChild2.plotPathway(TRnum, fnum, maxNodeLevel);
            end

            hold on;

        end

        function [plottedTransverseNodes, plottedLongitudinalNodes] = plotNode(transversePopulationNodeObject, pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel)

            if transversePopulationNodeObject.level > maxNodeLevel+1
                return;
            end
            
            hold on;

            if coherenceDegreeIsInNodeList(plottedLongitudinalNodes, transversePopulationNodeObject.coherenceDegree, transversePopulationNodeObject.level)
                c = [0.8196 0 1.0000];
            else
                c = [0.6350 0.0780 0.1840];
            end

            plot(transversePopulationNodeObject.totalTime, transversePopulationNodeObject.coherenceDegree, '.', 'color', c, 'MarkerSize', 15);
            plottedTransverseNodes = [plottedTransverseNodes, transversePopulationNodeObject];

            hold on;

            if isa(transversePopulationNodeObject.longitudinalChild1, "populationNode")

                if coherenceDegreeIsInNodeList(plottedTransverseNodes, transversePopulationNodeObject.longitudinalChild1.coherenceDegree, transversePopulationNodeObject.longitudinalChild1.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(transversePopulationNodeObject.longitudinalChild1.amplitudeLabel);
                pathwayString = transversePopulationNodeObject.longitudinalChild1.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end                

                if abs(f-0.5)>labelOverlapThreshold
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.longitudinalChild1.totalTime)/2, transversePopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, true, false); 
                    alignTextOnPathway(text(transversePopulationNodeObject.longitudinalChild1.totalTime, transversePopulationNodeObject.longitudinalChild1.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, false, noOverlap); 
                else
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.longitudinalChild1.totalTime)/2, transversePopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, true, false); 
                    alignTextOnPathway(text(transversePopulationNodeObject.longitudinalChild1.totalTime, transversePopulationNodeObject.longitudinalChild1.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, false, noOverlap); 
                end  

                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = transversePopulationNodeObject.longitudinalChild1.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

            if isa(transversePopulationNodeObject.longitudinalChild2, "populationNode")

                if coherenceDegreeIsInNodeList(plottedTransverseNodes, transversePopulationNodeObject.longitudinalChild2.coherenceDegree, transversePopulationNodeObject.longitudinalChild2.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(transversePopulationNodeObject.longitudinalChild2.amplitudeLabel);
                pathwayString = transversePopulationNodeObject.longitudinalChild2.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end             

                if (~contains(transversePopulationNodeObject.label, "0") && ~contains(transversePopulationNodeObject.label, "-1")) || abs(f-0.5)>labelOverlapThreshold
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.longitudinalChild2.totalTime)/2, -transversePopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, true, false);
                    alignTextOnPathway(text(transversePopulationNodeObject.longitudinalChild2.totalTime, transversePopulationNodeObject.longitudinalChild2.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), 1, false, noOverlap);                    
                else
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.longitudinalChild2.totalTime)/2, -transversePopulationNodeObject.coherenceDegree, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, true, false);
                    alignTextOnPathway(text(transversePopulationNodeObject.longitudinalChild2.totalTime, transversePopulationNodeObject.longitudinalChild2.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0 0.4470 0.7410], 'Interpreter', 'latex'), -1, false, noOverlap);                    
                end  

                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = transversePopulationNodeObject.longitudinalChild2.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

            if isa(transversePopulationNodeObject.transverseChild1, "populationNode")

                if coherenceDegreeIsInNodeList(plottedLongitudinalNodes, transversePopulationNodeObject.transverseChild1.coherenceDegree, transversePopulationNodeObject.transverseChild1.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(transversePopulationNodeObject.transverseChild1.amplitudeLabel);
                pathwayString = transversePopulationNodeObject.transverseChild1.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end                  

                if ~contains(transversePopulationNodeObject.label, "0") && ~contains(transversePopulationNodeObject.label, "-1")
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.transverseChild1.totalTime)/2, (transversePopulationNodeObject.coherenceDegree+transversePopulationNodeObject.transverseChild1.coherenceDegree)/2, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, true, true); 
                    alignTextOnPathway(text(transversePopulationNodeObject.transverseChild1.totalTime, transversePopulationNodeObject.transverseChild1.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, false, noOverlap);                    
                else
                    alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.transverseChild1.totalTime)/2, (transversePopulationNodeObject.coherenceDegree+transversePopulationNodeObject.transverseChild1.coherenceDegree)/2, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, true, false); 
                    alignTextOnPathway(text(transversePopulationNodeObject.transverseChild1.totalTime, transversePopulationNodeObject.transverseChild1.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, false, noOverlap);
                end

                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = transversePopulationNodeObject.transverseChild1.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

            if isa(transversePopulationNodeObject.transverseChild2, "populationNode")

                if coherenceDegreeIsInNodeList(plottedLongitudinalNodes, transversePopulationNodeObject.transverseChild2.coherenceDegree, transversePopulationNodeObject.transverseChild2.level)
                    noOverlap = false;
                else
                    noOverlap = true;
                end

                amplitudeString = latex(transversePopulationNodeObject.transverseChild2.amplitudeLabel);
                pathwayString = transversePopulationNodeObject.transverseChild2.label;
                if strlength(amplitudeString)>1100
                    amplitudeString = "CharLimit";
                end
                if strlength(pathwayString)>1100
                    pathwayString = "CharLimit";
                end  

                alignTextOnPathway(text(transversePopulationNodeObject.transverseChild2.totalTime, transversePopulationNodeObject.transverseChild2.coherenceDegree, pathwayString, 'FontSize', pathwayLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, false, noOverlap);
                alignTextOnPathway(text((transversePopulationNodeObject.totalTime+transversePopulationNodeObject.transverseChild2.totalTime)/2, (-transversePopulationNodeObject.coherenceDegree+transversePopulationNodeObject.transverseChild2.coherenceDegree)/2, string("$"+amplitudeString+"$"), 'FontSize', amplitudeLabelFontsize, 'Color', [0.6350 0.0780 0.1840], 'Interpreter', 'latex'), 1, true, false);               
                
                hold on;
                [plottedTransverseNodes, plottedLongitudinalNodes] = transversePopulationNodeObject.transverseChild2.plotNode(pathwayLabelFontsize, amplitudeLabelFontsize, plottedTransverseNodes, plottedLongitudinalNodes, height, f, labelOverlapThreshold, maxNodeLevel);
                hold on;

            end

        end

    end
    
end