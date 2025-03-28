% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function ratio = fcngain2ratio(gain,type)
%gain (dB)
%ratio (amplitudeOut/amplitudeIn)
%type = 'amplitude' or 'power' ratio

switch type
    case 'amplitude' %amplitudeRatio
        ratio = exp(gain*log(10)/10).^.5;
    case 'power' %powerRatio
        ratio = exp(gain*log(10)/10);
end


