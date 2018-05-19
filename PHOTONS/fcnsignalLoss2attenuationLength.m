function al = fcnsignalLoss2attenuationLength(dBpm)
%This function converts signal loss (dB/meter) into attenuation length (meters)
%Example spectra to apply this function to on page 6 of PDF: http://www.crystals.saint-gobain.com/sites/imdf.crystals.com/files/documents/fiber-brochure.pdf
%https://en.wikipedia.org/wiki/Attenuation_length
%https://en.wikipedia.org/wiki/Attenuation

al = 1./(-log(1./exp(dBpm*log(10)/10))); %(m) attenuation length (distance at which 63% of particles are stopped)