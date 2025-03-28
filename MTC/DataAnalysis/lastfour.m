% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function i=lastfour(w,i)
i = i( w>(max(w)-4) );
end
