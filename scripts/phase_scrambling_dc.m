% performs phase scrambling of concentration signals (dc) in the input data

function data_new = phase_scrambling_dc(data_all)
    % initialize output data
    data_new = data_all;
    % select dc data
    data_dc = data_all.procResult.dc;
    for type=1:size(data_dc,2) % loop over chromophores (ohb,hhb)
        for ch=1:size(data_dc,3) % loop over channels
            % extract data of a single chromophore
            data = data_dc(:,type, ch);
            % calculate Fast Forurier Transform (fft)
            data_fft = fft(data);
            % extract fft magnitude and phase
            mag = abs(data_fft);
            phase = angle(data_fft);
            % randomize fft phase between -pi and pi
            phase = -pi + (pi+pi)*rand(size(phase)); 
            % recompose fft
            scrambled_data_fft = mag .* exp(i * phase);
            % calculate the inverse fft (ifft) to reconstruct the signal
            scrambled_data = abs(ifft(scrambled_data_fft));
            % assign phase scrambled data to output file 
            data_new.procResult.dc(:,type,ch) = scrambled_data;
        end
    end
end