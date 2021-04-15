mkdir('pings')
clear; close all; clc
read = 0;
if read
    cd('..\PS4')
    [x,xt] = xlsread('cps09mar.xlsx','Sheet1');
    save 'data'
    cd('..\PS5')
else
    cd('..\PS4')
    load 'data'
    cd('..\PS5')
end

men = logical(1-x(:,2)); % first matlab question
Y = categorical(x(men,8));
X = x(men,[1 4 3]);
race = x(men,11);
black = race==2; % is this right?
X = [X black];
[B,~,stats] = mnrfit(X,Y);
tab = table(B,stats.se,stats.t,'VariableNames',{'Coeff','SE','T'}, ...
    'RowNames',{'intercept','age','education','hispanic','black'});
table2latex(tab,'tab1.tex')

women = logical(x(:,2)); % second matlab question
college = logical(x(:,4)>13);
woco = logical(women.*college);
Y = categorical(x(woco,12)<4);
age = x(woco,1);
ov30 = 0*age;
ov30(age>30,:) = age(age>30,:)-30;
ov40 = 0*age;
ov40(age>40,:) = age(age>40,:)-40;
ov50 = 0*age;
ov50(age>50,:) = age(age>50,:)-50;
ov60 = 0*age;
ov60(age>60,:) = age(age>60,:)-60;
ov70 = 0*age;
ov70(age>70,:) = age(age>70,:)-70;
X = [age ov30 ov40 ov50 ov60 ov70];
%X = [age ov30 ov30.^2 ov50 ov50.^2 ov70 ov70.^2];
B = mnrfit(X,Y);
aG = (min(age):max(age))';
aG30 = 0*aG;
aG30(aG>30,:) = aG(aG>30,:)-30;
aG40 = 0*aG;
aG40(aG>40,:) = aG(aG>40,:)-40;
aG50 = 0*aG;
aG50(aG>50,:) = aG(aG>50,:)-50;
aG60 = 0*aG;
aG60(aG>60,:) = aG(aG>60,:)-60;
aG70 = 0*aG;
aG70(aG>70,:) = aG(aG>70,:)-70;
G = [ones(size(aG)) aG aG30 aG40 aG50 aG60 aG70];
%G = [ones(size(aG)) aG aG30 aG30.^2  aG50 aG50.^2  aG70 aG70.^2];
fit = (1+exp(G*B)).^(-1);
emp = 0*fit;
for i=1:length(emp)
    emp(i) = mean(double(Y(age==aG(i)))-1);
end
figure
plot(aG,fit,'k')
hold on
plot(aG,emp,'r-.')
hold off
title('estimated marriage probabilities for women with degrees')
xlabel('age')
ylabel('marriage probability')
legend('fitted probabilities, spline','empirical mean conditioning on age')
set(gcf,'Color',[1 1 1])
cd('pings')
saveas(gcf,'fig1.png')
cd('..')