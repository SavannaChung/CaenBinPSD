%% This script is to read the binary file from DT5730 digitiser's output.

tic;
clear ; clc;


    % variable
    
    % Identify the files for analysis
    fileName=dir('*_ls_0*'); % get all the text files end with . =dir('Run34 **/*.txt');
    if isempty(fileName); fprintf(' >!>!>! error in identifying the txt files in this directory!\n');end
    fileName={fileName(~[fileName.isdir]).name};
    
    detNo=str2num(string(extractBetween(fileName(1), "ls_", ".dat")));
    if detNo==4;
        Det="NPL";
    elseif detNo==0;
        Det="UCL";
    end
        

    % Read the binary file: 1.) convert the first six line into header. 2.)
    % Read the data into a four-column format. 
    % 1st col= timeStamp. 2nd col=Qlong. 3rd col= EXTRAS. 4th col=Qshort

iF=1; % loop through all files 
%for iF=1:1:length(fileName)
for iF=1:1:length(fileName)
        recordType = {'uint32' 'int16' 'uint32' 'int16'};
        recordLen = [4 2 4 2];
        R = cell(1,numel(recordType));

        %# read column-by-column
        fid = fopen(fileName{iF},'rb'); %fseek(fid, 1*6, 'bof');
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
        
        clear('Qlong','Qshort','PSD', 'R');
        
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
                    % nCycle * 2^32 + first cycle + last cycle (V1751)
                % AcqT3_s(iF)= (double(nCycle-1)*(2^32-1) + (2^32-1- double(timeStamp(1))) + double(timeStamp(end)))/1E9;
                    % nCycle * 2^32 + first cycle + last cycle (DT5730)
                    % first cycle
                    if timeStamp(2)>timeStamp(1)
                        firstCycle= 2^31-1- double(timeStamp(1))*2;
                    else
                        firstCycle= 2;
                    end
                    
                    if timeStamp(end)>timeStamp(end-1)
                        lastCycle= double(timeStamp(end))*2;
                    else
                        lastCycle= 2;
                    end
                    
                AcqT3_s(iF)= (double(nCycle-1)*(2^31-1)*2 + firstCycle + lastCycle)/1E9;

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
                tSfileName=sprintf('%s_timeStamp_%s', Det, runNo);
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
                eH_nbin=uint16(2000);
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
                     strcat('noEvts = ', num2str(length(nQlong)));
                     strcat('noEvts = ', num2str(eH_nbin))};
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
                eHfileName=sprintf('%s_EnHist_%s', Det, runNo);
                saveas(gcf, eHfileName, 'fig');
                
                % plot the EnHist in log scale
                f_log_eH=figure;
                semilogy(eH_xaxis, eH_hist, 'bd:','LineWidth', 2, 'MarkerEdgeColor', 'k', 'MarkerSize', 2);
                annotation(f_log_eH, 'textbox', [0.5 0.7 0.3 0.1],'FitBoxToText', 'on', 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');             
                grid on;
                xlabel('Qlong (a.u.)'); 
                ylabel('Counts');
                axis tight;
                FontWeight= 'bold';
                FontSize= 50;
                set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2);
                pbaspect([1.5 1 1]);

                 set(f_log_eH, 'PaperUnits', 'centimeter', 'PaperPosition', [0 0 15 10]);
                % extract the run number from the fileName
                %runNo=string(extractBetween(fileName(iF), "Run__", "_ls"));
                eHfileName=sprintf('%s_LogEnHist_%s', Det, runNo);
                saveas(gcf, eHfileName, 'fig');
                close all;
                clear('f_*','eHfileName','FontSize', 'str', 'eB', 'eHfileName');
                
% Part3: *** ___ FoM calculation ___ ***   
                foM_psd_nBin=2000;
                foM_psd_edge=linspace(0, 1, foM_psd_nBin);
                foM_psd_xaxis=foM_psd_edge(1:end-1)+diff(foM_psd_edge)./2;
                foM_hist=histcounts(nPSD, foM_psd_edge);
                %gaussEqn='a*exp(-((x-b)/c).^2) +d';
                foM_fit=fit(foM_psd_xaxis.', foM_hist.', 'gauss2');
                yfoM_fit=feval(foM_fit, foM_psd_xaxis);
                
                gauFitCoeff=coeffvalues(foM_fit);
                y_eqt1=gauFitCoeff(1)*exp(-((foM_psd_xaxis-gauFitCoeff(2))/gauFitCoeff(3)).^2);
                y_eqt2=gauFitCoeff(4)*exp(-((foM_psd_xaxis-gauFitCoeff(5))/gauFitCoeff(6)).^2);
                
                f_foM= figure; 
                histogram(nPSD, foM_psd_nBin); hold on
                plot(foM_psd_xaxis.', yfoM_fit, 'r--', 'LineWidth', 2);
                plot(foM_psd_xaxis, y_eqt1, 'g:', 'LineWidth', 2);
                plot(foM_psd_xaxis, y_eqt2, 'M:', 'LineWidth', 2);
                xlabel('PSD (a.u.)');
                ylabel('Frequency (a.u.)');
                xlim([0 0.5]);
                grid on;
                set(gca, 'FontSize', 18, 'FontWeight', 'bold', 'LineWidth', 2);
                pbaspect([1.5 1 1]);
                str={strcat('Tot. evts= ', num2str(length(nQlong)), 'nBin = ', num2str(length(foM_psd_nBin))), ...
                     strcat('Mean1= ', num2str(round(gauFitCoeff(2),4)),'|',...
                     'Mean2= ', num2str(round(gauFitCoeff(5),4)),'|'),...
                     strcat('FWHM1= ', num2str(round(2.35*gauFitCoeff(3)/2, 4)),'|',...
                     'FWHM2= ', num2str(round(2.35*gauFitCoeff(6)/2, 4)), '|'), ...
                     strcat('Separation =',num2str(round(abs(minus(gauFitCoeff(2),gauFitCoeff(5))),4)), '|',...
                     'FoM =',num2str(round(abs(minus(gauFitCoeff(2),gauFitCoeff(5)))/(2.35*(gauFitCoeff(3)+gauFitCoeff(6))/2),4)))
                      };
                     %strcat('FoM =',num2str(round(abs(minus(gauFitCoeff(2),gauFitCoeff(5)))/(2.35*(gauFitCoeff(3)+gauFitCoeff(6))/2),4)))};
                %annotation(f_foM, 'textbox', [.5 .96 .1 .1], 'FitBoxToText', 'on', 'String', str , 'FontSize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w');
                title(str, 'FontSize', 10)
                %lgd=legend(str, 'Location','southoutside', 'FontSize', 10 );
                
                set(f_foM, 'PaperUnits', 'centimeter', 'PaperPosition', [0 0 15 10]);
                % extract the run number from the fileName
                %runNo=string(extractBetween(fileName(iF), "Run__", "_ls"));
                eHfileName=sprintf('%s_FoM_%s', Det, runNo);
                saveas(gcf, eHfileName, 'fig');
                
                close all;
                clear('f_foM','eH*', 'str', 'foM*', 'gauFitCoeff');
                
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
                eHfileName=sprintf('%s_TDSP_%s', Det, runNo);
                saveas(gcf, eHfileName, 'fig');
                
                close all;
                clear('f_TDSP','Data1','eHfileName','h1', 'hplot', 'hr1', 'xb', 'yb', 'str', 'TD_nBin', 'nQlong', 'nQshort', 'nPSD');
                
          
               
                runTime=round(toc);
                fprintf('%d file done. runtTime=%d \n', iF, runTime);
                iF=iF+1;     
end

% Part5: *** ___ Output the Acquisition time in a text file ___ ***
                C=cat(2,fileName(1:length(AcqT3_s))', num2cell(round(AcqT3_s,1))' );
                acqT_table=cell2table(C, 'VariableNames', {'fileName', 'AcqT3_s'});
                acq_fNx=sprintf('%s Acq Time.xlsx', Det); % output excel
                writetable(acqT_table, acq_fNx, 'sheet', 1);
                
                acq_fN=sprintf('%s Acq Time.txt', Det);
                writetable(acqT_table, acq_fN, 'Delimiter', '\t'); % output txt
                