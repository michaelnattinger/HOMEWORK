clc; % Clears the output in the command window
clear; % Clears the memory/workspace

%% Introduction to MATLAB
% A simple introduction to MATLAB for ECON712 - UW Madison
% Written by Eirik Eylands Brandsaas, updated by Jason Choi for 2018
% Always write your program in the editor
% To run the entire program, F5
% To run just highlighted parts F9
% You can always find out what a code means by searching at https://www.mathworks.com/help/

%% MATLAB is a matrix language
s=1 % you created a variable called "s" that has value of 1! 
x=[1 2 3] % now you've created a vector x!
Y=[1 2 3; 4 5 6] % this one's a matrix! 

%% I recommend using lots of comments in your code
%% Using two %% creates a section, so you can run one and one section by pressing ctrl+enter, ctrl+shift+enter
%% Use ; to suppress output
Z=[1 2 3; 4 5 6];

%% Basic matrix operation
Y*Z'
size(Y)
W=inv(Y*Z')
A=[Y Z] % Stack matrices
eye(2) % "eyedentity matrix"
ones(1,4) 
zeros(4,1) 
B=1:0.01:2 % Creates evenly spaced vector
B=linspace(1,2,100)


%% Elementwise operations
Y.*Z
Y./Z
1./Z

%% Random numbers
X=rand(10,1)
Y=rand(10,1)

%% Loops
for iter=1:10
    XY(iter)=X(iter,1)*Y(iter,1);
end
XY'
X.*Y

tol=0.01;
iter=0;
error=1
while error>tol  % Loop stops whenever the error becomes less than the tolerance
    iter=iter+1
    X=X.*X;
    error=max(X);
end

%% What to do if program never stops? (ctrl+c in the command window)
% This loop will never stop, since the error is strictly increasing
tol=0.00001;
iter=0;
error=1
while error>tol
    iter=iter+1
    X=X+X;
    error=max(X);
end

%% Use maxiter and the AND statement
tol=0.00001;
iter=0;
error=1
maxiter=10^4
while error>tol && iter<=maxiter
    iter=iter+1
    X=X+X;
    error=max(X);
end

%% Using if, elseif, else.
% Note that ' indicates a string is coming, and when it ends
x=rand
if x<0.5
    disp('Macro is cooler than micro')
elseif x<0.99
    disp('Macro is WAY cooler than micro')
else
    disp('Micro is cooler than macro!')
end

%% Functions: probability that micro is cooler than macro
x=0; % Initialize
iter=0; % Reset counter
while x<0.995 && iter<=maxiter %(Remember to set a max iteration number)
    iter=iter+1;
    x=rand;
end
disp(['The simulated probability of micro being cooler than macro is ',num2str(1/iter)])
% Note: disp only takes one object (a matrix, a vector a scalar) as an
% input. So use the brackets to indicate a matrix, , to separate the two
% elements and finally the num2str function translates the number into a
% string
%% But maybe we often want to check this probability with different thresholds
% You should write functions! (I never do, but you should!) See the
% attached function for the syntax
[prob, iter]=prob_macro_micro(0.9999,100)
[prob, iter]=prob_macro_micro(0.9999,10000)
    
%% Figures
t=1997:2007
y=log(1:1:11)
figure(1);
plot(t,y);
title('Capital');
xlabel('Year'); ylabel('Capital');
legend('K_t','Location','Southeast');
saveas(1,'test','png'); %Saves the figure, usefull for homeworks!
%PS: To make many graphs in one figure, use command subplot.

%%
figure(2)
subplot(2,1,1)
plot(t,sqrt(1:1:11))
title('This is subplot 1 of a figure with 2 rows and 1 column')
subplot(2,1,2)
Y=1:1:11;
Y=Y'*y;
surf(Y)
title('This is a 3d subplot 2 of a figure with 2 rows and 1 column')

