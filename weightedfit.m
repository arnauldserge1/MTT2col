function Results = weightedfit(data, do_plot)

% This code fits makes a linear fit to a data set (using y =bx+a) where each data point
% has a different or constant standard deviation. Your data should have three or two columns.
% The first column should be the independent variable(x) and the second
% column should be the dependent variable(y). Column three should contain
% your standard deviations for each datapoint. In the situations where you
% do not specify a column three, the code assigns a weight of one to all
% data points and this corresponds to the regular linear fits.
%==========
% INPUTS
%==========
%data = 3 columns; column 1 = x, column2 = y and column 3 = standard dev.

%==========
%OUTPUTS
%==========
%Result.slope= b; Fitted slope
%Result.Intercept = a; Fitted intercept

%Coded by Ebo Ewusi-Annan (University of Florida, 2011)
%============
%REFERENCES
%===========
%1. Willam H. Press, Saul A. Teukolsky and Willan T. Vetterling (1997).
%Numerical Recipes in Fortran.
%2. Philip R. Bevington and D. Keith Robinson (2003). Data Reduction and
%Error Analysis for the Physical Sciences.

if nargin<2, do_plot = 1; end

x = data(:,1);
y = data(:,2);
[s, t]= size(data);
stdv = ones(s,1);
if t==3
    stdv = data(:,3);
end
w = 1./stdv.^2;
S = sum(w);
Sx = sum(w.*x);
Sy = sum(w.*y);
Sxx= sum(w.*x.^2);
Sxy= sum(w.*x.*y);
Delta = S*Sxx - (Sx)^2;
a = (Sxx*Sy - Sx*Sxy)./Delta;
b = (S*Sxy - Sx*Sy)./Delta;
Results.slope = b;
Results.Intercept = a;

if do_plot
    fprintf('\n slope=%f Int=%f \n',b,a)

    y_fit = a + b*x;
    set(gcf,'color','w');
    if t==2
        h=plot(x,y,'bs','MarkerFaceColor','b');
        title('Unweighted fit')
    else
        h=errorbar(x,y,stdv,'bs','MarkerFaceColor','b');
        title('Weighted fit')
    end
    hold on; q=plot(x,y_fit,'r.--','linewidth',2);
    xlabel('x (Column 1)')
    ylabel('y (Column 2)')
% % %     legend([h(1),q(1)],'Data',sprintf('\nSlope=%f\nIntercept=%f\n',b,a),'location','Southeast')
end

end