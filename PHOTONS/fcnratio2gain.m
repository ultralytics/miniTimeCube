function gain = fcnratio2gain(ratio,type)
%gain (dB)
%ratio (amplitudeOut/amplitudeIn)
%type = 'amplitude' or 'power' ratio

switch type
    case 'amplitude' %amplitudeRatio
        gain = 10*log10(ratio.^2);
    case 'power' %powerRatio
        gain = 10*log10(ratio);
end


