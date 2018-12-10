%% This script is to read the binary file from DT5730 digitiser's output.

tic;
clear all; clc;

    % variable
    
    % Identify the files for analysis
    fileName=dir('*_ls_4*'); % get all the text files end with . =dir('Run34 **/*.txt');
    if isempty(fileName); fprintf(' >!>!>! error in identifying the txt files in this directory!\n');end
    fileName={fileName(~[fileName.isdir]).name};

    % Read the binary file: 1.) convert the first six line into header. 2.)
    % Read the data into a four-column format. 
    % 1st col= timeStamp. 2nd col=Qlong. 3rd col= EXTRAS. 4th col=Qshort

iF=1; % loop through all files 
%for iF=1:1:length(fileName)
for iF=1:1:10
        recordType = {'uint32' 'int16' 'uint32' 'int16'};
        recordLen = [4 2 4 2];
        R = cell(1,numel(recordType));

        %# read column-by-column
        fid = fopen(fileName{1},'rb'); %fseek(fid, 1*6, 'bof');
        Header=fread(fid, 6, 'uint32'); % Read the first six headerss in the ls bin. file
        for i=1:numel(recordType) % Please find reference in https://stackoverflow.com/questions/8096702/reading-multiple-precision-binary-files-through-fread-in-matlab
            
            %# seek to the first field of the first record
            fseek(fid, sum(recordLen(1:i-1)), 'bof');

            %# % read column with specified format, skipping required number of bytes
            R{i} = fread(fid, Inf, ['*' recordType{i}], sum(recordLen)-recordLen(i));
        end
        fclose(fid);
        fclose all;
        
        timeStamp=R{1};
        Qlong=double(R{2});
        Qshort=double(R{4});
        PSD=(minus(double(Qlong),double(Qshort))./double(Qlong));
        
        nQlong=Qlong(PSD>0 & PSD <1 & Qlong>0);
        nQshort=Qshort(PSD>0 & PSD <1 & Qlong>0);
        nPSD=PSD(PSD>0 & PSD <1 & Qlong>0);      
        
        clear('Qlong','Qshort','PSD');
% Part1: *** ___ Calculate the acqTime ___ ***
         % Method 3: calculate the nCycle. tot_AcqT= nCycle*2^32
                nCycle=0;
                cnEvt=0;
                EvtspC=[];
                peaktS=[];
                troughtS=[];
                iS=2; % the 1st event
                for iS=2:1:length(timeStamp);

                        if timeStamp(iS-1)<timeStamp(iS);
                            cnEvt=cnEvt+1;  
                        else timeStamp(iS-1)>timeStamp(iS);
                            nCycle=nCycle+1;
                            EvtspC(nCycle, 1)=cnEvt; % number of events per cycle
                            troughtS(nCycle, 1)=timeStamp(iS); % trough values
                            peaktS(nCycle, 1)=timeStamp(iS-1); % peak values

                            cnEvt=0; % Count the noEvt per cycle again. 

                        end
                    iS=iS+1;
                end

                if nCycle==0;
                    AcqT3_s= (timeStamp(end)-timeStamp(1))/1E9;
                else
                    % nCycle * 2^32 + first cycle + last cycle
                AcqT3_s(iF)= (double(nCycle-1)*(2^32-1) + (2^32-1- double(timeStamp(1))) + double(timeStamp(end)))/1E9;

                end
                
                
                % plot timeStamp and save this figure
                f1=figure; plot(timeStamp,'LineWidth', 1 ); 
                str={strcat('AcqT3 (s)= ', num2str(round(AcqT3_s(iF),1)))};
%                 annotation(f1, 'textbox', [.3 .3 .3 .1], 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');
                annotation(f1, 'textbox','FitBoxToText', 'on', 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');             
                xlabel('th event (a.u.)');
                ylabel('timeStamp (a.u.)');
                set(gca,'FontWeight', 'bold', 'FontSize', 15, 'LineWidth', 2);
                pbaspect([1.5 1 1]);

                set(f1, 'PaperUnits', 'centimeter', 'PaperPosition', [0 0 15 10]);
                % extract the run number from the fileName
                runNo=string(extractBetween(fileName(iF), "Run__", "_ls"));
                tSfileName=sprintf('timeStamp_%s', runNo);
                saveas(gcf, tSfileName, 'png');
                
                close all;
% %             fileIDacqT=fopen('acqTime.txt', 'w');
% %             % write the acquisition time to a text file
% %             str=fprintf(fileIDacqT,'Acquisition time in second= %f', AcqTimeIns);
% %             fclose(fileIDacqT);

        % After acquiring the acquistion time, we need to free some
        % memory from Matlab.
        clear('f1','R','EvtspC','Header', 'i', 'iS', 'nCycle', 'peaktS', 'recordLen', 'recordType', 'troughtS', 'timeStamp', 'str', 'fid', 'cnEvt', 'tSfileName' )

% Part2: *** ___ Energy Histogram ___ ***
                eH_nbin=uint16(1000);
                max_chNo=2^15-1;
                eH_edge=linspace(0, max_chNo, eH_nbin);
                eH_xaxis=eH_edge(1:end-1)+ diff(eH_edge)./2;
                
                eH_hist=histcounts(nQlong, eH_edge);
                eH_err=sqrt(eH_hist);
                f_eH=figure;
                plot(eH_xaxis, eH_hist, 'bd:','LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerSize', 2);
%                 eB=errorbar(eH_xaxis, eH_hist, eH_err);
%                 eB.Color='black';
                str={strcat('1 bin = ', num2str(round(eH_edge(2),1)));
                     strcat('noEvts = ', num2str(length(nQlong)))};
%                 annotation(f1, 'textbox', [.3 .3 .3 .1], 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');
                annotation(f_eH, 'textbox', [0.5 0.7 0.3 0.1],'FitBoxToText', 'on', 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');             
                grid on;
                xlabel('Qlong (a.u.)'); 
                ylabel('Counts');
                axis tight;
                FontWeight= 'bold';
                FontSize= 50;
                set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2);
                pbaspect([1.5 1 1]);

                 set(f_eH, 'PaperUnits', 'centimeter', 'PaperPosition', [0 0 15 10]);
                % extract the run number from the fileName
                %runNo=string(extractBetween(fileName(iF), "Run__", "_ls"));
                eHfileName=sprintf('EnHist_%s', runNo);
                saveas(gcf, eHfileName, 'fig');
                
                close all;
                clear('f_eH','eHfileName','FontSize', 'str', 'eB', 'eHfileName');
                
% Part3: *** ___ FoM calculation ___ ***   
                foM_psd_nBin=1000;
                foM_psd_edge=linspace(0, 1, foM_psd_nBin);
                foM_psd_xaxis=foM_psd_edge(1:end-1)+diff(foM_psd_edge)./2;
                foM_hist=histcounts(nPSD, foM_psd_edge);
                %gaussEqn='a*exp(-((x-b)/c).^2) +d';
                foM_fit=fit(foM_psd_xaxis.', foM_hist.', 'gauss1');
                yfoM_fit=feval(foM_fit, foM_psd_xaxis);
                
                f_foM= figure; 
                histogram(nPSD, foM_psd_nBin); hold on
                plot(foM_psd_xaxis.', yfoM_fit, 'r--', 'LineWidth', 2);
                xlabel('PSD (a.u.)');
                ylabel('Frequency (a.u.)');
                xlim([0 0.5]);
                grid on;
                set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2);
                pbaspect([1.5 1 1]);
                

                
                close all;
                clear('f_foM','eH*');
                
% Part4: *** ___ 2D scatter plot ___ ***
                TD_nBin=500;   %% control the number of bins inbetween 0 and 1

                xb=linspace(0, max_chNo, TD_nBin);  
                yb=linspace(0,1,TD_nBin);  
                
                Data1= [nQlong, nPSD];
                %Data1= [En, nPSD1];
                h1=hist3(Data1, {xb yb});
                h1(h1==0)=NaN;  %% the color plot will not include 0 value
                hr1=h1';

                f_TDSP=figure;

            %     hplot=pcolor(xb, yb, hr1/acqT); hold on;
                hplot=pcolor(xb, yb, hr1); hold on;
                set(hplot, 'EdgeColor', 'none');
                str={strcat('Qlong TD nBin = ', num2str(TD_nBin));
                     strcat('PSD TD nBin = ', num2str(TD_nBin));
                     strcat('noEvts = ', num2str(length(nQlong)))};
%                 annotation(f1, 'textbox', [.3 .3 .3 .1], 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');
                annotation(f_TDSP, 'textbox', [0.5 0.7 0.3 0.1],'FitBoxToText', 'on', 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');             
                xlim([0 inf]);
                ylim([0 0.5]);
                %caxis([0 0.15])
                grid on;
                %grid minor;
                %title( 'RUN8 AL 500ns');
                colorbar;
                xlabel('Qlong (a.u.)');
                ylabel('PSD (a.u.)');
                set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2)
                pbaspect([1.5 1 1]);
                
                set(f_TDSP, 'PaperUnits', 'centimeter', 'PaperPosition', [0 0 15 10]);
                % extract the run number from the fileName
                %runNo=string(extractBetween(fileName(iF), "Run__", "_ls"));
                eHfileName=sprintf('TDSP_%s', runNo);
                saveas(gcf, eHfileName, 'fig');
                
                close all;
                clear('f_TDSP','Data1','eHfileName','h1', 'hplot', 'hr1', 'xb', 'yb', 'str', 'TD_nBin');
                
          
% Part5: *** ___ Output the Acquisition time in a text file ___ ***
                C=cat(2,fileName(1:length(AcqT3_s))', num2cell(round(AcqT3_s,1))' );
                acqT_table=cell2table(C, 'VariableNames', {'fileName', 'AcqT3_s'});
                writetable(acqT_table, 'Acq Time.txt', 'Delimiter', '\t');

                iF=iF+1;     
end