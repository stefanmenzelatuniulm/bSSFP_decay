function [updateIndices, pruneIndices, summedAmplitudes, summedAmplitudeLabels, updateLabels] = populationTreeMergePruneHelper(phases, amplitudes, amplitudeLabels, labels)
    
    values = uniquetol(phases);
    counts = histcounts(phases, [transpose(values), 2*max(values)+1]);
    repeatedElements = values(counts >=  1);
    transverseDuplicateIndices = zeros(1, length(phases));

    iter = 0;
    for k = 1:length(repeatedElements)
        addElements = find(ismembertol(phases,repeatedElements(k)));
        for m = 1:length(addElements)
            iter = iter+1;
            transverseDuplicateIndices(iter) = addElements(m);
        end
    end

    if iter<length(transverseDuplicateIndices)
        transverseDuplicateIndices(iter+1:end) = [];
    end

    [transverseGC, transverseGR] = groupcounts(phases(transverseDuplicateIndices)); 
    transverseGC = transpose(transverseGC); %number of occurrences of duplicate elements
    transverseGR = transpose(transverseGR); %duplicate elements
    sums = sym(zeros(1, length(unique(phases))));
    sums2 = sym(zeros(1, length(unique(phases))));
    updateLabels = strrep(string(zeros(1, length(unique(phases)))), "0", "");

    iter1 = 0;
    iter2 = 0;
    for k = 1:length(transverseGR)

        for m = 1:transverseGC(k) 
            iter1 = iter1+1;
            index = transverseDuplicateIndices(iter1);
            sums(k) = sums(k)+amplitudes(index);
            sums2(k) = sums2(k)+amplitudeLabels(index);
            if m~= transverseGC(k)
                updateLabels(k) = updateLabels(k)+labels(index)+" + ";
            else
                updateLabels(k) = updateLabels(k)+labels(index);
            end
        end

        for m = 1:transverseGC(k) 
            iter2 = iter2+1;
            index = transverseDuplicateIndices(iter2);
            amplitudes(index) = sums(k); %A now contains the sum over all occurences of each element
            amplitudeLabels(index) = sums2(k);
            labels(index) = updateLabels(k);
        end

    
    end

    [uniqueAmplitudes, updateIndices, ~] = unique(amplitudes);
    updateIndices = transpose(updateIndices);
    pruneIndices = setdiff(linspace(1, length(amplitudes), length(amplitudes)), updateIndices);
    summedAmplitudes = uniqueAmplitudes;
    summedAmplitudeLabels = amplitudeLabels(updateIndices);
    updateLabels = labels(updateIndices);

end