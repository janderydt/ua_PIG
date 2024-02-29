% test
clear all;

n=[100];
x = [1:21];
figure; hold on;

for jj=1:length(x)
    
    for ii=1:length(n)

        I = [1:n(ii)];

        terms = x(jj).^I.*(-1).^(I+1)./factorial(I);
        buildsum = cumsum(terms);

        TS(jj,ii) = buildsum(end);
        F(jj) = sinh(x(jj))-cosh(x(jj))+1;

    end
    
end

plot(x,TS(:,1),'ok');
plot(x,F,'r');