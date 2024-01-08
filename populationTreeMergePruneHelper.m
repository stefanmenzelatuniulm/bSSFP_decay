function [updateIndices, pruneIndices, summedAmplitudes, summedAmplitudeLabels, updateLabels] = populationTreeMergePruneHelper(dephasingDegrees, amplitudes, amplitudeLabels, labels)
    
    edges = min(dephasingDegrees) : max(dephasingDegrees)+1;
    [counts, values] = histcounts(dephasingDegrees, edges);
    repeatedElements = values(counts >=  1);
    transverseDuplicateIndices = zeros(1, length(dephasingDegrees));

    iter = 0;
    for k = 1:length(repeatedElements)
        addElements = find(dephasingDegrees == repeatedElements(k));
        for m = 1:length(addElements)
            iter = iter+1;
            transverseDuplicateIndices(iter) = addElements(m);
        end
    end

    if iter<length(transverseDuplicateIndices)
        transverseDuplicateIndices(iter+1:end) = [];
    end

    [transverseGC, transverseGR] = groupcounts(dephasingDegrees(transverseDuplicateIndices)); 
    transverseGC = transpose(transverseGC); %number of occurrences of duplicate elements
    transverseGR = transpose(transverseGR); %duplicate elements
    sums = sym(zeros(1, length(unique(dephasingDegrees))));
    sums2 = sym(zeros(1, length(unique(dephasingDegrees))));
    updateLabels = strrep(string(zeros(1, length(unique(dephasingDegrees)))), "0", "");

    iter1 = 0;
    iter2 = 0;
    for k = 1:length(transverseGR)

        for m = 1:transverseGC(k) 
            iter1 = iter1+1;
            index = transverseDuplicateIndices(iter1);
            sums(k) = sums(k)+amplitudes(index);
            sums2(k) = sums2(k)+amplitudeLabels(index);
            if m~= transverseGC(k)
                updateLabels(k) = updateLabels(k)+labels(index)+"&";
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