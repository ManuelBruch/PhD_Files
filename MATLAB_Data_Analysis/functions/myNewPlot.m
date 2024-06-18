function myNewPlot(x,y,NVArgs)
%   MYNEWPLOT this function has been designed by Manuel Bruch (M. Sc.) to
%   create plots in a standardised format and makes use of Name-Value-pair
%   arguments.
%   myNewPlot(x,y) will plot y vs x. X can be a vector, potentially
%   containing error values in its third dimension. Y can be a matrix that
%   has to agree with x in either the first or the second dimension and can
%   contain error values in the third dimension.
%
%   myNewPlot(x,y,Name=Value) can feed in additional formatting arguments
%   into the function. The default values of each pair can be found in the
%   "define input arguments" section
%   IMPORTANT NOTE: before ver R2021a use the myNewPlot(x,y,'Name','Value')
%   syntax instead (still works with R2021a as well)
%
%   myLegendPosition    -->     defines legend position according to
%                               MATLABs predefined cardinal locations. Set
%                               to 'none' if no legend is needed
%
%   myAxes              -->     defines axis to be used for the plot
%
%   myColorOrder        -->     defines colour order. Colour order will be
%                               permuted before the line styles. The order
%                               should be defined in an RGB matrix and can
%                               be obtained by using the linspecer function
%
%   myLSOrder           -->     defines a line style order
%
%   myInt               -->     sets the interpreter for legend, axes names
%                               and title
%
%   FS                  -->     define font sizes for x-axis, y-axis,
%                               legend and title. Define it as a
%                               four-element vector with a scalar number
%                               for each component
%
%   myTitle             -->     define a title for your plot
%
%   myAxisNames         -->     define a 1x2 cell with the axis names 
%                               {'x-axis name','y-axis name'}
%
%   DataSeriesNames     -->     define names for each data series in y as a
%                               cell array {'Name1','Name2',...}
%
%   myYTickFormat       -->     define a standardised tick format for the y
%                               axis (number of decimals)
%
%   myxLim              -->     sets lower and upper bound on x. define as
%                               a 1x2 vector of real numbers or +/- Inf.
%                               Must contain values for both boundaries
%
%   myyLim              -->     sets lower and upper bound on y. define as
%                               a 1x2 vector of real numbers or +/- Inf.
%                               Must contain values for both boundaries
%
%   myColor             -->     sets the colour for the line to be plotted
%
%   transpErr           -->     sets errorbars to transparent. set the alpha
%                               value of the transparency.

%% define input arguments
arguments
    x   double
    y   double
    NVArgs.myLegendPosition (1,:)   char    =   'eastoutside';
    NVArgs.myAxes                           =   gca;
    NVArgs.myColorOrder             double  =   [0 0 0];
    NVArgs.myLSOrder                        =   {'-x', '--x', ':x', '-.x'};
    NVArgs.myInt            (1,:)   char    =   'none';
    NVArgs.FS                       double  =   [11 11 11 11];
    NVArgs.myTitle                  char    =   [];
    NVArgs.myAxisNames              cell    =   {'',''};
    NVArgs.DataSeriesNames          cell    =   {};
    NVArgs.myYTickFormat            char    =   '%.2f';
    NVArgs.myxLim                           =   [-Inf,Inf];
    NVArgs.myyLim                           =   [-Inf,Inf];
    NVArgs.myColor                  double  =   [0 0 0];
    NVArgs.transpErr                double  =   [];
end

%% reshape x and y --> assume that when x is longer than high, values are to
% be plotted horizontally
if size(x,1) < size(x,2)
    x   =   permute(x,[2 1 3]);
    y   =   permute(y,[2 1 3]);
end

%% set up axes properties
ax                  =   NVArgs.myAxes;
ax.ColorOrder       =   NVArgs.myColorOrder;
ax.LineStyleOrder   =   NVArgs.myLSOrder;

%% decide upon error bars
if size(x,3) > 1
    x_new   =   x(:,:,1);
    x_err   =   x(:,:,2);
    y_new   =   y(:,:,1);
    y_err   =   y(:,:,2);
    x       =   x_new;
    y       =   y_new;
    errors  =   1;
    noXerrors   =   0;
elseif size(x,3) == 1 && size(y,3) > 1
%     x_err   =   zeros(size(x));
    y_new   =   y(:,:,1);
    y_err   =   y(:,:,2);
    y       =   y_new;
    errors  =   1;
    noXerrors   =   1;
else
    errors  =   0;
    noXerrors   =   0;
end

%% plot
hold on
if noXerrors
    for i = 1:size(y,2)
        if size(NVArgs.myColor,1) ~= 1
            myColor     =   NVArgs.myColor(i,:);
        else
            myColor     =   NVArgs.myColor;
        end
        if ~isempty(NVArgs.DataSeriesNames)
            h   =   errorbar(x,y(:,i),y_err(:,i),...
                    'DisplayName',NVArgs.DataSeriesNames{i},'LineWidth',1,...
                    'Color',myColor);
        else
            h   =   errorbar(x,y(:,i),y_err(:,i),...
                    'LineWidth',1,'Color',myColor);
        end
        if ~isempty(NVArgs.transpErr)
            % from https://uk.mathworks.com/matlabcentral/answers/473325-how-to-make-errorbar-transparent
            % Set transparency level (0:1)
            alpha = NVArgs.transpErr;
            % Set transparency (undocumented)
            set(h.Bar, 'ColorType', 'truecoloralpha', 'ColorData',...
                [h.Line.ColorData(1:3); 255*alpha])
            set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData',...
                [h.Cap.EdgeColorData(1:3); 255*alpha])
        end
    end
elseif errors && ~noXerrors
    for i = 1:size(y,2)
        if size(NVArgs.myColor,1) ~= 1
            myColor     =   NVArgs.myColor(i,:);
        else
            myColor     =   NVArgs.myColor;
        end
        if ~isempty(NVArgs.DataSeriesNames)
            h   =   errorbar(x,y(:,i),y_err(:,i),y_err(:,i),x_err,x_err,...
                    'DisplayName',NVArgs.DataSeriesNames{i},'LineWidth',1,...
                    'Color',myColor);
        else
            h   =   errorbar(x,y(:,i),y_err(:,i),y_err(:,i),x_err,x_err,...
                    'LineWidth',1,'Color',myColor);
        end
        if ~isempty(NVArgs.transpErr)
            % from https://uk.mathworks.com/matlabcentral/answers/473325-how-to-make-errorbar-transparent
            % Set transparency level (0:1)
            alpha = NVArgs.transpErr;
            % Set transparency (undocumented)
            set(h.Bar, 'ColorType', 'truecoloralpha', 'ColorData',...
                [h.Line.ColorData(1:3); 255*alpha])
            set(h.Cap, 'EdgeColorType', 'truecoloralpha', 'EdgeColorData',...
                [h.Cap.EdgeColorData(1:3); 255*alpha])
        end
    end
else
    for i = 1:size(y,2)
        if size(NVArgs.myColor,1) ~= 1
            myColor     =   NVArgs.myColor(i,:);
        else
            myColor     =   NVArgs.myColor;
        end
        if ~isempty(NVArgs.DataSeriesNames)
            plot(x,y(:,i),'DisplayName',NVArgs.DataSeriesNames{i},...
                'LineWidth',1,'Color',myColor)
        else
            plot(x,y(:,i),'LineWidth',1,'Color',myColor)
        end
    end
end

if ~strcmp(NVArgs.myxLim, 'none')
    xlim(NVArgs.myxLim)
end
if ~strcmp(NVArgs.myyLim, 'none')
    ylim(NVArgs.myyLim)
end

set(gca,'FontName','Arial')
grid on

if ~strcmp(NVArgs.myTitle, 'none')
    if size(y,2) >= 2
        title(NVArgs.myTitle,'FontName','Arial','FontSize',NVArgs.FS(4),...
            'Interpreter',NVArgs.myInt)
    elseif ~isempty(NVArgs.DataSeriesNames)
        title(NVArgs.DataSeriesNames{i},'FontName','Arial','FontSize',...
            NVArgs.FS(4),'Interpreter',NVArgs.myInt)
    end
end

xlabel(NVArgs.myAxisNames{1},'FontName','Arial','FontSize',NVArgs.FS(1))
ylabel(NVArgs.myAxisNames{2},'FontName','Arial','FontSize',NVArgs.FS(2))

if size(y,2) >= 2
    if ~strcmp(NVArgs.myLegendPosition,'none')
        legend('Location',NVArgs.myLegendPosition,'Interpreter',...
            NVArgs.myInt,'FontSize',NVArgs.FS(3))
    end
end

if ~isempty(NVArgs.myYTickFormat)
    ytickformat(NVArgs.myYTickFormat)
end

end

