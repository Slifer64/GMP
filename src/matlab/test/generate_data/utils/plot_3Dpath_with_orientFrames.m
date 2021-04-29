function ax = plot_3Dpath_with_orientFrames(Pos, Quat, varargin)
%% Plots 3D path with frames denoting the orientation at each point.
%  @param[in] Pos: 3xN matrix with position (x,y,z) at each column.
%  @param[in] Quat: 4xN matrix with orientation as a unit quaternion at each column.
%
%  Optional arguments:
%  @param[in] axes: The axis object where to plot the data (optional, default = gca, i.e. the current axes).
%
%  Variable argument Name-Value pairs:
%  @param[in] title: The title of the axis (optional, default = '').
%  @param[in] xlabel: The x-label of the axis (optional, default = '').
%  @param[in] ylabel: The y-label of the axis (optional, default = '').
%  @param[in] zlabel: The z-label of the axis (optional, default = '').
%  @param[in] LineWidth: The linewidth of the path line (optional, default = 1.0).
%  @param[in] LineStyle: The line style of the path line (optional, default = '-').
%  @param[in] LineColor: The color of the path line in rgb format (optional, default = [0.45 0.26 0.26]).
%  @param[in] Interpreter: The text interpreter (optional, default = 'latex').
%  @param[in] fontSize: The text fontsize (optional, default = 12).
%  @param[in] numberOfFrames: The number of orientation frames to plot. They will be equally spaced along the path (optional, default = 6).
%  @param[in] frameScale: The scaling of the orientation frames size (optional, default = 1.0).
%  @param[in] frameLineWidth: The linewidth of the orientation frame axes (optional, default = 1.0).
%  @param[in] frameLineStyle: The linestyle of the axes of the orientation frame.
%  @param[in] frameXAxisColor: The color of the x-axis of the orientation frame in rgb format (optional, default = [1 0 0]).
%  @param[in] frameYAxisColor: The color of the y-axis of the orientation frame in rgb format (optional, default = [0 1 0]).
%  @param[in] frameZAxisColor: The color of the z-axis of the orientation frame in rgb format (optional, default = [0 0 1]).
%  @param[in] animated: If true the path is plotted animated (optional, default = false).
%  @param[in] Time: 1xN matrix with timestamps (used to set the animation speed).


% inPars.addParameter('frameXAxisLegend', {}, @(x)assert_string(x));
% inPars.addParameter('frameYAxisLegend', {}, @(x)assert_string(x));
% inPars.addParameter('frameZAxisLegend', {}, @(x)assert_string(x));
% inPars.addParameter('LineLegend', {}, @(x)assert_string(x));

%
%  \note The Names of the Name-Value parameters are not case sensitive.
%  Partial matching is disabled.
%  Unmatched name parameters are ignored (no warnings or error messages are produced).
%     
    
    %% Parse the input arguments
    [inArgs, usingDefaults, unmatchedNames] = parseInputArguments(varargin{:});

    ax = inArgs.axes;
    axis_colors = cell({inArgs.FrameXAxisColor, inArgs.FrameYAxisColor, inArgs.FrameZAxisColor});
    Time = inArgs.Time;
    if (isscalar(Time))
        dt = Time;
        Time = (0:(size(Pos,2)-1))*dt;
    end
    
    
    %% Set view on the 'ax' object
    view(ax, inArgs.view);

    %% Extract position and orientation coordinates
    Axang = quat2axang(Quat')';
    X = Pos(1,:);   Y = Pos(2,:);    Z = Pos(3,:);
    U = Axang(1,:); V = Axang(2,:);  W = Axang(3,:);
plot_3Dpath_with_orientFrames
    
    %% Set axes limits
    w = max(max(X) - min(X), 1.0);
    a = w*0.1;
    % ax.XLim = [min(X)-a max(X)+a];
    ax.XLim = [min(ax.XLim(1),min(X)-a) max(ax.XLim(2),max(X)+a)];
    w = max(max(Y) - min(Y), 1.0);
    a = w*0.1;
    % ax.YLim = [min(Y)-a max(Y)+a];
    ax.YLim = [min(ax.YLim(1),min(Y)-a) max(ax.YLim(2),max(Y)+a)];
    w = max(max(Z) - min(Z), 1.0);
    a = w*0.1;
    % ax.ZLim = [min(Z)-a max(Z)+a];
    ax.ZLim = [min(ax.ZLim(1),min(Z)-a) max(ax.ZLim(2),max(Z)+a)];
    
    
    %% Find where to place the orientation frames so that they equally spaced
    m = inArgs.NumberOfFrames; % how many frames to put between the start and end frame
    dist = zeros(length(X),1);
    for i=2:length(dist)
        dist(i) = dist(i-1) + norm([X(i) Y(i) Z(i)]-[X(i-1) Y(i-1) Z(i-1)]);
    end
    d = dist(end)/(m+1);
    frames_ind = [1];
    
    j = 1;
    for i=1:length(dist)
        if (dist(i) >= d*j)
            frames_ind = [frames_ind i];
            j = j + 1;
        end
    end
    frames_ind = [frames_ind length(dist)];
    frames_ind = unique(frames_ind);
    
    
    %% Enable 'hold on' on the ax object
    NextPlot0 = ax.NextPlot;
    ax.NextPlot = 'add';
    
    
    %% set the text properties
    ax.XLabel.String = inArgs.xlabel;
    ax.XLabel.Interpreter = inArgs.Interpreter;
    ax.XLabel.FontSize = inArgs.FontSize;

    ax.YLabel.String = inArgs.ylabel;
    ax.YLabel.Interpreter = inArgs.Interpreter;
    ax.YLabel.FontSize = inArgs.FontSize;

    ax.ZLabel.String = inArgs.zlabel;
    ax.ZLabel.Interpreter = inArgs.Interpreter;
    ax.ZLabel.FontSize = inArgs.FontSize;
    
    ax.Title.String = inArgs.title;
    ax.Title.Interpreter = inArgs.Interpreter;
    ax.Title.FontSize = inArgs.FontSize;

    
    %% Initialize 3 quivers, one for each axis of the orientation frame
    quiv = cell(3,1); % three quivers, for x, y and z axis of orientation frame
    for j=1:3
        quiv{j} = quiver3(ax, 0,0,0,0,0,0, inArgs.FrameScale);
        quiv{j}.Color = axis_colors{j};
        quiv{j}.LineStyle = inArgs.FrameLineStyle;
        quiv{j}.LineWidth = inArgs.FrameLineWidth;    
        quiv{j}.AutoScale = 'on';
    end


    %% Initialize line object
    lineH = line(ax);
    lineH.XData = [];
    lineH.YData = [];
    lineH.ZData = [];
    lineH.Color = inArgs.LineColor;
    lineH.LineWidth = inArgs.LineWidth;
    lineH.LineStyle = inArgs.LineStyle;
    
    
    % Keep previous legends and corresponding graphic objects
    legends = {};
    gObgs = [];
    for i=1:length(ax.Children)
       if (~isempty(ax.Children(i).DisplayName))
           legends = [legends ax.Children(i).DisplayName];
           gObgs = [gObgs ax.Children(i)];
       end
    end
    
    
    if (~strcmpi(inArgs.LineLegend,''))
        legends = [legends inArgs.LineLegend];
        gObgs = [gObgs lineH];
    end
    
    FrameAxisLegend = {inArgs.FrameXAxisLegend, inArgs.FrameYAxisLegend, inArgs.FrameZAxisLegend};
    for i=1:length(FrameAxisLegend)
        if (~strcmpi(FrameAxisLegend{i},''))
            legends = [legends FrameAxisLegend{i}];
            gObgs = [gObgs quiv{i}];
        end
    end
    
    legH = legend(ax, gObgs, legends, 'Interpreter',inArgs.Interpreter, 'FontSize',inArgs.LegendFontSize, ...
        'TextColor',inArgs.LegendTextColor, 'FontSize',inArgs.LegendFontSize);
    
    legH.Location = inArgs.LegendLocation;
    legH.Orientation = inArgs.LegendOrientation;
    legH.Box = inArgs.LegendBoxDisplay;
    
    if (inArgs.VideoCapture)
        videoWriter = VideoWriter(inArgs.VideoFilename);
        videoWriter.open();
    end
    
    if (~isempty(inArgs.ScreenScale))
        unitsPrev = ax.Parent.Units;
        ax.Parent.Units = 'normalized';
        ax.Parent.OuterPosition = [0 0 inArgs.ScreenScale(1) inArgs.ScreenScale(2)];
        ax.Parent.Units = unitsPrev; % restore it to previous value
    end
    
    %% Plot (animated) path with frames
    Time = [Time Time(end)];
    for i=1:(length(Time)-1)
        lineH.XData = [lineH.XData X(i)];
        lineH.YData = [lineH.YData Y(i)];
        lineH.ZData = [lineH.ZData Z(i)];
        
        if (inArgs.animated)

            T = makehgtform('translate',[X(i) Y(i) Z(i)]) * makehgtform('axisrotate',[U(i) V(i) W(i)],Axang(4,i));
            for j=1:3
                quiv{j}.XData = [quiv{j}.XData X(i)];
                quiv{j}.YData = [quiv{j}.YData Y(i)];
                quiv{j}.ZData = [quiv{j}.ZData Z(i)];

                quiv{j}.UData = [quiv{j}.UData T(1,j)];
                quiv{j}.VData = [quiv{j}.VData T(2,j)];
                quiv{j}.WData = [quiv{j}.WData T(3,j)];
            end

            drawnow;
            pause(Time(i+1)-Time(i));

            if (inArgs.VideoCapture)
                frame = getframe(ax);
                videoWriter.writeVideo(frame);
            end
            
            for j=1:3
                quiv{j}.XData = quiv{j}.XData(1:end-1);
                quiv{j}.YData = quiv{j}.YData(1:end-1);
                quiv{j}.ZData = quiv{j}.ZData(1:end-1);
                
                quiv{j}.UData = quiv{j}.UData(1:end-1);
                quiv{j}.VData = quiv{j}.VData(1:end-1);
                quiv{j}.WData = quiv{j}.WData(1:end-1);
            end
        end
    end
    
    for k=1:length(frames_ind)
        i = frames_ind(k);

        T = makehgtform('translate',[X(i) Y(i) Z(i)]) * makehgtform('axisrotate',[U(i) V(i) W(i)],Axang(4,i));
        for j=1:3
            quiv{j}.XData = [quiv{j}.XData X(i)];
            quiv{j}.YData = [quiv{j}.YData Y(i)];
            quiv{j}.ZData = [quiv{j}.ZData Z(i)];

            quiv{j}.UData = [quiv{j}.UData T(1,j)];
            quiv{j}.VData = [quiv{j}.VData T(2,j)];
            quiv{j}.WData = [quiv{j}.WData T(3,j)];
        end
        
        if (inArgs.VideoCapture)
            frame = getframe(ax);
            videoWriter.writeVideo(frame);
        end
    end
    
    if (inArgs.VideoCapture)
        videoWriter.close();
    end
    
    for j=1:3  
        quiv{j}.AutoScale = 'on';
    end

    
    %% restore it to its previous state
    ax.NextPlot = NextPlot0; 

end


function [inArgs, usingDefaults, unmatchedNames] = parseInputArguments(varargin)

    % initialize parser with the names and default values of the input arguments
    inPars = inputParser;
    
    inPars.KeepUnmatched = true;
    inPars.PartialMatching = false;
    inPars.CaseSensitive = false;
    
    inPars.addOptional('axes', gca, @(x)assert_axes(x));
    
    inPars.addParameter('NumberOfFrames', 6, @(x)assert_numeric_scalar_positive(x));
    inPars.addParameter('FrameScale', 1.0, @(x)assert_numeric_scalar_nonnegative(x));
    inPars.addParameter('FrameLineWidth', 1.0, @(x)assert_numeric_scalar_positive(x));
    inPars.addParameter('FrameLineStyle', '-', @(x)assert_string(x));
    
    inPars.addParameter('FrameXAxisColor', [1.0 0.0 0.0], @(x)assert_rgbColor(x));
    inPars.addParameter('FrameYAxisColor', [0.0 1.0 0.0], @(x)assert_rgbColor(x));
    inPars.addParameter('FrameZAxisColor', [0.0 0.0 1.0], @(x)assert_rgbColor(x));
    
    inPars.addParameter('title', '', @(x)assert_string(x));
    inPars.addParameter('xlabel', '', @(x)assert_string(x));
    inPars.addParameter('ylabel', '', @(x)assert_string(x));
    inPars.addParameter('zlabel', '', @(x)assert_string(x));
    
    inPars.addParameter('FrameXAxisLegend', '', @(x)assert_string(x));
    inPars.addParameter('FrameYAxisLegend', '', @(x)assert_string(x));
    inPars.addParameter('FrameZAxisLegend', '', @(x)assert_string(x));
    inPars.addParameter('LineLegend', '', @(x)assert_string(x));
    inPars.addParameter('LegendLocation', 'northeastoutside', @(x)assert_string(x));
    inPars.addParameter('LegendOrientation', 'vertical', @(x)assert_string(x));
    inPars.addParameter('LegendBoxDisplay', 'on', @(x)assert_string(x));
    inPars.addParameter('LegendTextColor', [0 0 0], @(x)assert_rgbColor(x));
    inPars.addParameter('LegendFontSize', 10, @(x)assert_numeric_scalar_positive(x));
    
    inPars.addParameter('LineWidth', 1.0, @(x)assert_numeric_scalar_positive(x));
    inPars.addParameter('LineColor', [0.45 0.26 0.26], @(x)assert_rgbColor(x));
    inPars.addParameter('LineStyle', '-', @(x)assert_string(x));
    
    inPars.addParameter('Interpreter', 'latex', @(x)assert_string(x));
    inPars.addParameter('FontSize', 10, @(x)assert_numeric_scalar_positive(x));
    inPars.addParameter('view', 3, @(x) assert( length(x)<=3, 'view must have length at most 3.'));
    
    inPars.addParameter('animated', false, @(x)assert_boolean(x));
    inPars.addParameter('Time', 0.01, @(x)assert_numeric_nonnegative_increasing(x));
    
    inPars.addParameter('VideoCapture', false, @(x)assert_boolean(x));
    inPars.addParameter('VideoFilename', '3Dpath_with_orientFrames.avi', @(x)assert_string(x));
    
    inPars.addParameter('ScreenScale', [], @(x) assert(isempty(x) || (isnumeric(x) && length(x)==2 && isempty(find(x>1 | x<=0))), 'Scale must be a 2D vector with values in the range (0 1]')); 
    
    % Parse input arguments
    inPars.parse(varargin{:});
    
    unmatchedNames = fieldnames(inPars.Unmatched);
    usingDefaults = inPars.UsingDefaults;
    
    inArgs = inPars.Results;
    
    if (~isempty(cellfun(@(x) strcmpi(x,'LegendFontSize'), usingDefaults)))
        inArgs.LegendFontSize = inArgs.FontSize;
    end
    
    if (~isempty(unmatchedNames))
        str = sprintf('plot_3Dpath_with_orientFrames: Found unmatched argument names:\n');
        for i=1:length(unmatchedNames)
            str = [str sprintf('%s\n', unmatchedNames{i})];
        end
        warning('%s', str); 
    end
    
%     disp('plot_3Dpath_with_orientFrames: Using defaults for:\n%s', usingDefaults{:});
    
end

function assert_axes(x)

 assert(strcmp(get(x, 'type'), 'axes'), 'Input must be of type axes.');

end

function assert_boolean(x)

 assert( islogical(x), 'Input must be boolean.');

end

function assert_numeric_nonnegative_increasing(x)

 assert(isnumeric(x) && isempty(find(x<0)) && isempty(find(diff(x)<0)), 'Input must be positive scalar or a numeric vector with non-negative increasing values');

end

function assert_numeric_scalar_nonnegative(x)

 assert(isnumeric(x) && isscalar(x) && (x >= 0), 'Value must be nonnegative, scalar, and numeric.');

end

function assert_numeric_scalar_positive(x)

 assert(isnumeric(x) && isscalar(x) && (x > 0), 'Value must be positive, scalar, and numeric.');

end

function assert_rgbColor(x)

 valid_color_ranges = true;
 for i=1:length(x)
	valid_color_ranges = valid_color_ranges && (0.0<=x(i)<=1.0);
 end

 assert(isnumeric(x) && length(x)==3 && valid_color_ranges, 'Input must be an rgb triplet with values in the range [0 1].');

end

function assert_string(x)

 assert( ischar(x), 'Input must be a string.');

end

