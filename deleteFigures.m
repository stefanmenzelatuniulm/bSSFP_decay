%Deletes all non-directories in the Figure subdirectories

function deleteFigures(subfolder)

    % Specify the folder where the files live.
    myFolder = pwd+"\"+subfolder;
    % Check to make sure that folder actually exists.  Warn user if it
    % doesn't.
    if ~isfolder(myFolder)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
        uiwait(warndlg(errorMessage));
    return;
    end
    % Get a list of all files in the folder with the desired file name
    % pattern.
    filePattern = fullfile(myFolder, "**/*"); % Change to whatever pattern you need.
    theFiles = dir(filePattern);
    theFiles([theFiles.isdir]) = [];
    for k = 1 : length(theFiles)
        baseFileName = theFiles(k).name;
        baseFileFolder = theFiles(k).folder;
        fullFileName = string(baseFileFolder)+"\"+string(baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end

end