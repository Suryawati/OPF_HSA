
%function [totalcost]=gencost(Pgg, cost)
if exist('Pgg')~=1
Pgg=input('Enter the scheduled real power gen. in row matrix ');
else,end
if exist('cost')~=1
cost = input('Enter the cost function matrix ');
else, end
ngg = length(cost(:,1));
Pmt = [ones(1,ngg); Pgg; Pgg.^2];
for i = 1:ngg
costv(i) = cost(i,:)*Pmt(:,i);
end
totalcost=sum(costv);
%fprintf('\nTotal generation cost = % 10.2f $/h \n', totalcost)
