function [updateIndices, pruneIndices, summedAmplitudes, summedAmplitudeLabels, updateLabels] = populationTreeMergePruneHelper(coherenceDegrees, amplitudes, amplitudeLabels, labels)
    
    values = uniquetol(coherenceDegrees);
    counts = histcounts(coherenceDegrees, [transpose(values), 2*max(values)+1]);
    repeatedElements = values(counts >=  1);
    transverseDuplicateIndices = zeros(1, length(coherenceDegrees));

    iter = 0;
    for k = 1:length(repeatedElements)
        addElements = find(ismembertol(coherenceDegrees,repeatedElements(k)));
        for m = 1:length(addElements)
            iter = iter+1;
            transverseDuplicateIndices(iter) = addElements(m);
        end
    end

    if iter<length(transverseDuplicateIndices)
        transverseDuplicateIndices(iter+1:end) = [];
    end

    [transverseGC, transverseGR] = groupcounts(coherenceDegrees(transverseDuplicateIndices)); 
    transverseGC = transpose(transverseGC); %number of occurrences of duplicate elements
    transverseGR = transpose(transverseGR); %duplicate elements
    sums = sym(zeros(1, length(unique(coherenceDegrees))));
    sums2 = sym(zeros(1, length(unique(coherenceDegrees))));
    updateLabels = strrep(string(zeros(1, length(unique(coherenceDegrees)))), "0", "");

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

    [~, updateIndices, ~] = unique(coherenceDegrees);
    updateIndices = transpose(updateIndices);
    pruneIndices = setdiff(linspace(1, length(amplitudes), length(amplitudes)), updateIndices);
    summedAmplitudes = amplitudes(updateIndices);
    summedAmplitudeLabels = amplitudeLabels(updateIndices);
    updateLabels = labels(updateIndices);

end