%% This script is to read the binary file from DT5730 digitiser's output.

tic;
clear all; clc;

    % variable
    
    % Identify the files for analysis
    fileName=dir('*001*'); % get all the text files end with . =dir('Run34 **/*.txt');
    if isempty(fileName); fprintf(' >!>!>! error in identifying the txt files in this directory!\n');end
    fileName={fileName(~[fileName.isdir]).name};

    % Read the binary file: 1.) convert the first six line into header. 2.)
    % Read the data into a four-column format. 
    % 1st col= timeStamp. 2nd col=Qlong. 3rd col= EXTRAS. 4th col=Qshort

        recordType = {'uint32' 'int16' 'uint32' 'int16'};
        recordLen = [4 2 4 2];
        R = cell(1,numel(recordType));

        %# read column-by-column
        fid = fopen(fileName{1},'rb'); %fseek(fid, 1*6, 'bof');
        Header=fread(fid, 6, 'uint32'); % Read the first six headerss in the ls bin. file
        for i=1:numel(recordType)
            %# seek to the first field of the first record
            fseek(fid, sum(recordLen(1:i-1)), 'bof');

            %# % read column with specified format, skipping required number of bytes
            R{i} = fread(fid, Inf, ['*' recordType{i}], sum(recordLen)-recordLen(i));
        end
        fclose(fid);
        fclose all;
        
        timeStamp=R{1};
        Qlong=R{2};
        Qshort=R{4};
        PSD=(minus(double(Qlong),double(Qshort))./double(Qlong));
     

        
        https://stackoverflow.com/questions/8096702/reading-multiple-precision-binary-files-through-fread-in-matlab
    