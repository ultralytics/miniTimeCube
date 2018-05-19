function [] = testOrderStatistic(input)
A=input.table.smearExp;




%randcdfc(Material(mi(1)).R.scintT,n)

x=A.x;
PDF = A.pdf;
CDF = cumsum(A.pdf)*(A.x(2)-A.x(1));

fig
for n = [1 2 3 5 10 50]
    k = 1; %sample
    fk = orderStatistic(k,n,PDF,CDF);
    [mu, sigma] = weightedMean(x,fk);
    plot(x,fk,'-','Display',sprintf('%2g PE \\mu = %.1f \\pm %.1f ns',n,mu,sigma)); xyzlabel('T (ns)','First PE Arrival Time (PDF)')
end
fcntight; set(gca,'xlim',[0 10]); fcnlinewidth(2); legend show; title('Order Statistic (Scintillator: 850 ps rise, 2000 ps fall)')
%fig; plot(x,cumsum(fk)/sum(fk) - .5)

fig; h=plotyy(x,PDF,x,CDF);
xyzlabel(h(1),'T (ns)','PDF'); ylabel(h(2),'CDF')

end

function fk = orderStatistic(k,n,PDF,CDF)
%k = order of sample. i.e. k=1 means first sample out of n total
%https://en.wikipedia.org/wiki/Order_statistic

%fk = (factorial(n)/factorial(k-1)/factorial(n-k))*CDF.^(k-1).*(1-CDF).^(n-k).*PDF;
fk = (prod((n-k+1):n)/factorial(k-1))*CDF.^(k-1).*(1-CDF).^(n-k).*PDF;  %fixes numerical issues with factorials>20 http://study.com/academy/lesson/division-of-factorials-definition-lesson-quiz.html
end