function [output] = myplot(x,y,myAxisNames,DataSeriesNames,varargin)
%MYPLOT plots data in a designated figure in standardised format
%   P = myplot(x,y,myAxisNames,DataSeriesNames) for vector x and matrix y
%   will plot each column in y as a function of x on axis labeled from
%   myAxisNames in the form {xlabel, ylabel}. Legend entries are derived
%   from DataSeriesNames with each value in the cell representing one name
%   
%   P = myplot(...,myLegendLocation) defines the legend position (default:
%   'eastoutside')
%
%   P = myplot(...,myLegendLocation,myLineStyles) gives the chance to
%   define custom linestyles. Set to 'colour' to cycle through colours
%
%   P = myplot(...,myLegendLocation,myLineStyles,myInterpreter) sets
%   Interpreter for text in plot (default: 'none')
%
%   P = myplot(...,myLegendLocation,myLineStyles,myInterpreter,myFontSize)
%   defines font sizes for the plot [Legend, AxisLabels, TickLabels, Title]
%   (default: Matlab default)
%
%   P = myplot(...,myLegendLocation,myLineStyles,myInterpreter,myFontSize,
%   Title) sets title
%   P = myplot(...,myLegendLocation,myLineStyles,myInterpreter,myFontSize,
%   Title,myax) sets current axes

if nargin >= 10 && ~isempty(varargin{6})
    ax  =   varargin{6};
else
    ax  =   axes;
end

if nargin >= 5
    if ~isempty(varargin{1})
        myLegendPosition        =   varargin{1};
    else
        myLegendPosition        =   'eastoutside';
    end
    if (nargin >= 6) && ~isempty(varargin{2})
        if ischar(varargin{2})
            ax.ColorOrder       = 	linspecer(size(y,2));
            ax.LineStyleOrder   =   '-x';
        elseif isa(varargin{2},'double')
            ax.ColorOrder       = 	varargin{2};
            ax.LineStyleOrder   =   '-x';
        else
            ax.ColorOrder       =   [0 0 0];
            ax.LineStyleOrder   =   varargin{2};
        end
    else
        ax.ColorOrder       =   [0 0 0];
        ax.LineStyleOrder       =   {'-','--','-.',':'};
    end
    if (nargin >= 7) && ~isempty(varargin{3})
        myInt                   =   varargin{3};
    else            
        myInt                   =   'none';
    end
    if (nargin >= 8) && ~isempty(varargin{4})
        FS                      =   varargin{4};
    else            
        FS                      =   [11 11 11 11];
    end
    if (nargin >= 9) && ~isempty(varargin{5})
        myTitle                 =   varargin{5};
    else            
        myTitle                 =   [];
    end
    if (nargin >= 12) && ~isempty(varargin{8})
        ax.ColorOrder           =   varargin{8};
    else            
        ax.ColorOrder           =   [0 0 0];
    end
else
    myLegendPosition        =   'eastoutside';
    ax.ColorOrder           =   [0 0 0];
    ax.LineStyleOrder       =   {'-','--','-.',':'};
    myInt                   =   'none';
    FS                      =   [11 11 11 11];
    myTitle                 =   [];
end

if size(x,3) > 1
    x_new   =   x(:,:,1);
    x_err   =   x(:,:,2);
    y_new   =   y(:,:,1);
    y_err   =   y(:,:,2);
    x       =   x_new;
    y       =   y_new;
    errors  =   1;
elseif size(x,3) == 1 && size(y,3) > 1
    x_err   =   zeros(size(x));
    y_new   =   y(:,:,1);
    y_err   =   y(:,:,2);
    y       =   y_new;
    errors  =   1;
else
    errors  =   0;
end

hold on

if errors
    for i = 1:size(y,2)
        errorbar(x,y(:,i),y_err(:,i),y_err(:,i),x_err,x_err,'DisplayName',DataSeriesNames{i},'LineWidth',1)
    end
else
    for i = 1:size(y,2)
        plot(x,y(:,i),'DisplayName',DataSeriesNames{i},'LineWidth',1)
    end
end
set(gca,'FontName','Arial')
grid on

if size(y,2) >= 2
    title(myTitle,'FontName','Arial','FontSize',FS(4),'Interpreter',myInt)
else
%     title(DataSeriesNames{i},'FontName','Arial','FontSize',FS(4),'Interpreter',myInt)
end
if ~isempty(myAxisNames)
    xlabel(myAxisNames{1},'FontName','Arial','FontSize',FS(2))
    ylabel(myAxisNames{2},'FontName','Arial','FontSize',FS(2))
end
if size(y,2) >= 2
    if nargin >= 11
        if varargin{7} && ~isempty(varargin{7})
            legend('Location',myLegendPosition,'Interpreter',myInt,'FontSize',FS(1))
        end
    else
        legend('Location',myLegendPosition,'Interpreter',myInt,'FontSize',FS(1))
    end
end
% ytickformat('%.2f')
output  = [];

end

