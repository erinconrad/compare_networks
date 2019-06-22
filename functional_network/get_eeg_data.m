function data = get_eeg_data(dataName,indices,loginname,pwname)

n = 0;

%while n == 0
    %try
        session = IEEGSession(dataName, loginname, pwname);
        channelLabels = session.data.channelLabels;
        values = session.data.getvalues(indices,':');
        fs = session.data.sampleRate;
        n = 1; 
    %catch  
    %    error('what\n');
    %    fprintf('Failed to retrieve ieeg.org data, trying again...\n');        
    %end
    
%end


data.values = values;
data.ch_labels = channelLabels;
data.fs = fs;

session.delete;
clearvars -except data

end