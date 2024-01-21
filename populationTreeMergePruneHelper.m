function [updateIndices, pruneIndices, summedAmplitudeLabels, updateLabels, summedAmplitudesWithoutT2s, summedAmplitudesDirectlyAfterPulseWithoutT2s] = populationTreeMergePruneHelper(coherenceDegrees, amplitudeLabels, labels, amplitudesWithoutT2s, amplitudesDirectlyAfterPulseWithoutT2s)
    
    values = uniquetol(coherenceDegrees);
    counts = histcounts(coherenceDegrees, [transpose(values), 2*max(values)+1]);
    repeatedElements = values(counts >=   1);
    duplicateIndices = zeros(1, length(coherenceDegrees));

    iter = 0;
    for k = 1:length(repeatedElements)
        addElements = find(ismembertol(coherenceDegrees, repeatedElements(k)));
        for m = 1:length(addElements)
            iter = iter+1;
            duplicateIndices(iter) = addElements(m);
        end
    end

    if iter<length(duplicateIndices)
        duplicateIndices(iter+1:end) = [];
    end

    [transverseGC, transverseGR] = groupcounts(coherenceDegrees(duplicateIndices)); 
    transverseGC = transpose(transverseGC); %number of occurrences of duplicate elements
    transverseGR = transpose(transverseGR); %duplicate elements
    sums = sym(zeros(1, length(unique(coherenceDegrees))));
    sums2 = sym(zeros(1, length(unique(coherenceDegrees))));
    sums3 = sym(zeros(1, length(unique(coherenceDegrees))));
    updateLabels = strrep(string(zeros(1, length(unique(coherenceDegrees)))), "0", "");

    iter1 = 0;
    iter2 = 0;
    for k = 1:length(transverseGR)

        for m = 1:transverseGC(k) 
            iter1 = iter1+1;
            index = duplicateIndices(iter1);
            sums(k) = sums(k)+amplitudeLabels(index);
            sums2(k) = sums2(k)+amplitudesWithoutT2s(index);
            sums3(k) = sums3(k)+amplitudesDirectlyAfterPulseWithoutT2s(index);
            if m~=  transverseGC(k)
                updateLabels(k) = updateLabels(k)+labels(index)+" + ";
            else
                updateLabels(k) = updateLabels(k)+labels(index);
            end
        end

        for m = 1:transverseGC(k) 
            iter2 = iter2+1;
            index = duplicateIndices(iter2);
            amplitudeLabels(index) = sums(k);
            amplitudesWithoutT2s(index) = sums2(k);
            amplitudesDirectlyAfterPulseWithoutT2s(index) = sums3(k);
            labels(index) = updateLabels(k);
        end

    end

    [~, updateIndices, ~] = unique(coherenceDegrees);
    updateIndices = transpose(updateIndices);
    pruneIndices = setdiff(linspace(1, length(amplitudeLabels), length(amplitudeLabels)), updateIndices);
    summedAmplitudeLabels = amplitudeLabels(updateIndices);
    summedAmplitudesWithoutT2s = amplitudesWithoutT2s(updateIndices);
    summedAmplitudesDirectlyAfterPulseWithoutT2s = amplitudesDirectlyAfterPulseWithoutT2s(updateIndices);
    updateLabels = labels(updateIndices);

end