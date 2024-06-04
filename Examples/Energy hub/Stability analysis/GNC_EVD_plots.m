%% GNC and EVD MATLAB plots from frequency-domain data
% Author: Francisco Javier Cifuentes Garcia
close all;clear all;clc
set(0,'defaulttextinterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
set(0,'defaultAxesFontSize',16)
format compact; format long

%% Import the GNC data and define standard colors
filename = 'stable_GNC.txt';
num_vars = 16; % Number of variables: matrix size is num_vars x num_vars
lambdas = importmatrix(filename, num_vars, "\t");
f = real(lambdas(:,1)); % Extract the frequency at the 1st column
lambdas = lambdas(:,2:end); % The rest of the data
N = size(lambdas,2);
c = []; % Standard colors as needed
for i = 1:ceil(N/12)
    c = [c;colororder("gem12")];
end

%% Nyquist plot
figure('Name','Nyquist')
hold on; box on
title('Nyquist Plot','interpreter','latex');
xlabel('$\Re$','interpreter','latex')
ylabel('$\Im$','interpreter','latex')
N = size(lambdas,2);
for L =1:N
plot(real(lambdas(:,L)),imag(lambdas(:,L)),"-",'Color',c(L,:),"LineWidth",2);hold on
end
for L =1:N
plot(real(lambdas(:,L)),-imag(lambdas(:,L)),"-",'Color',c(L,:),"LineWidth",2)
end
plot(-1,0,'kx','MarkerSize',10,'LineWidth',1);grid on
plot(cos(-pi:1e-3:pi),sin(-pi:1e-3:pi),'k--','Linewidth',1)

axes('position',[.62 .3 .25 .3]);box on;hold on
for L =1:N
plot(real(lambdas(:,L)),imag(lambdas(:,L)),"-",'Color',c(L,:),"LineWidth",2);hold on
end
for L =1:N
plot(real(lambdas(:,L)),-imag(lambdas(:,L)),"-",'Color',c(L,:),"LineWidth",2)
end
plot(-1,0,'kx')
plot(cos(-pi:1e-3:pi),sin(-pi:1e-3:pi),'k--','Linewidth',1)
xlim([-2 2]);ylim([-2 2]);xticks(-2:1:2);yticks(-2:1:2)
title('Critical point','interpreter','latex');


%% Import the EVD data and define standard colors
filename = 'stable_EVD.txt';
num_vars = 16; % Number of variables: matrix size is num_vars x num_vars
evd = importmatrix(filename, num_vars, "\t");
f = real(evd(:,1)); % Extract the frequency
evd = evd(:,2:end); % The rest of the data
N = size(evd,2);
c = []; % Standard colors
for i = 1:ceil(N/12)
    c = [c;colororder("gem12")];
end

%% Plot the EVD
figure('Name','EVD')
hold on; box on; grid on
title('Oscillation mode identification','interpreter','latex');
set(gca,'xScale', 'log');set(gca,'yScale', 'log');
xlabel("Frequency [Hz]");ylabel("$$|\lambda_i(\mathbf{Z}_{bus}^{cl})|$$")

set(gcf,'Units','centimeters','position',[15,5,25,10])
for L =1:N
plot(f,abs(evd(:,L)),"-",'Color',c(L,:),"LineWidth",2);hold on
end
grid off; grid


%% Function for data import
function y = importmatrix(filename, num_vars, delimiter)
%% Input handling
% If dataLines is not specified, define defaults
if nargin < 3
    fprintf('\n Missing inputs! \n')
    num_vars = 1;
end
dataLines = [2, Inf];

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", num_vars);
opts = setvartype(opts, 'string'); % Set variable type to string to handle text processing

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = delimiter;

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";

% Import the data as text
textData = readmatrix(filename, opts, 'OutputType', 'string');

% Remove the unwanted symbols ( and )
cleanedTextData = regexprep(textData, '[\(\)]', '');

% Convert the cleaned text data to numeric values
y = str2double(cleanedTextData);

end