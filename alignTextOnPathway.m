function alignTextOnPathway(textObject, shiftSign, isAmplitudeLabel, noOverlap)
    
    drawnow; %important for correct value of extent

    extent = get(textObject).Extent;
    width = extent(3);
    height = extent(4); 

    if isAmplitudeLabel

        if noOverlap

            set(textObject, 'Position', get(textObject).Position+[-width/2 0 0]);

        else

            set(textObject, 'Position', get(textObject).Position+[-width/2 height*shiftSign 0]);

        end

    else %isPathwayLabel

        s = get(textObject).String;

        if ~isempty(s)
            width = width/length(s);
        end
        
        if noOverlap

            set(textObject, 'Position', get(textObject).Position+[width height/2 0]);

        else

            set(textObject, 'Position', get(textObject).Position+[width -height/2 0]);

        end
  
    end

end