%% Colors

%% Susanne's blue/pink color map
n = round(100/4);               % number of colors
R_1 = linspace(0,0,n);  
G_1 = linspace(0.403,0.667,n);
B_1 = linspace(0.603,1,n);
R_2 = linspace(0,1,n);  
G_2 = linspace(0.667,1,n);
B_2 = linspace(1,1,n);
R_3 = linspace(1,1,n);  
G_3 = linspace(1,0,n);
B_3 = linspace(1,0.49,n);
R_4 = linspace(1,0.596,n);  
G_4 = linspace(0,0,n);
B_4 = linspace(0.49,0.298,n);
cmap_SH = [R_1(:), G_1(:), B_1(:)
           R_2(:), G_2(:), B_2(:)
           R_3(:), G_3(:), B_3(:)
           R_4(:), G_4(:), B_4(:)];

% To use:       
% colormap(cmap_SH);  % 

%% Yufan's pasteljet
cmap_YX = customcolormap([0 .25 .5 .75 1], ...
    {'#9d0142','#f66e45','#f5f9f3','#65c0ae','#5e4f9f'});

% To use: 
% colormap(cmap_YX);  % 

%%
% ----------------------------------------------------------------------- %
% FUNCTION "customcolormap" defines a customized colobar given the        %
% positions and the colors that are going to generate the gradients.      %
%                                                                         %
%   Input parameters:                                                     %
%       - positions:    Vector of positions, from 0 to 1. Note that the   %
%                       first position must be 0, and the last one must   %
%                       be 1.                                             %
%       - colors:       Colors to place in each position. This parameter  %
%                       can be specified as a RGB matrix (n_colors x 3), or
%                       as a cell vector, containing HTML values.         %
%                       For instance: {'#ffffff','#ff0000','#000000'} is  %
%                       equivalent to [1 1 1; 1 0 0; 0 0 0].              %
%       - m:            (Optional) Number of points (recommended: m > 64).%
%                                                                         %
%   Output variables:                                                     %
%       - J:            Colormap in RGB values (dimensions [mx3]).        %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       J = customcolormap([0 0.5 1], {'#ffffff','#ff0000','#000000'});   %
%       colorbar; colormap(J);                                            %
% ----------------------------------------------------------------------- %
%   Versions:                                                             %
%       - v1.0.:    (19/11/2018) Original script.                         %
% ----------------------------------------------------------------------- %
%       - Author:   Víctor Martínez-Cagigal                               %
%       - Date:     19/11/2018                                            %
%       - Version:  1.0                                                   %
%       - E-mail:   vicmarcag (at) gmail (dot) com                        %
%                                                                         %
%       Biomedical Engineering Group (University of Valladolid), Spain    %
% ----------------------------------------------------------------------- %
function J = customcolormap(positions, colors, m)

    % Error detection and defaults
    if nargin < 3
       f = get(groot,'CurrentFigure');
       if isempty(f)
          m = size(get(groot,'DefaultFigureColormap'),1);
       else
          m = size(f.Colormap,1);
       end
    end
    if ~isnumeric(m), error('Parameter m must be numeric.'); end
    if nargin < 2, error('Not enough parameters.'); end
    
    if iscell(colors)
        colors = colors(:);
        n_colors = length(colors);
        for i = 1:n_colors
            temp = colors{i};
            if ~ischar(temp)
                error(['Colors must be specified in HEX format (e.g., #FFFFFF).' ...
                ' Type "help colorbar" for further information']);
            elseif ~strcmp(temp(1),'#')
                error(['Character # is missing if %s.' ...
                ' Type "help colorbar" for further information'], temp);
            elseif length(temp)~=7
                error(['Not a valid color format: %s (use this format: #FFFFFF).' ...
                ' Type "help colorbar" for further information'], temp);
            end
        end
    elseif ismatrix(colors)
        n_colors = size(colors);
        if length(n_colors) ~= 2
            error(['RGB colors must be a 2D matrix.' ...
                ' Type "help colorbar" for further information']);
        elseif n_colors(2) ~= 3
            error(['RGB colors matrix must have 3 columns.' ...
                ' Type "help colorbar" for further information']);
        elseif min(colors(:))<0 || max(colors(:))>255
            error(['RGB colors matrix values must range from 0 to 255.' ...
                ' Type "help colorbar" for further information']);
        end
    else
        error(['Colors must be a cell vector or a matrix of RGB values.' ...
            ' Type "help colorbar" for further information']);
    end
    if ~isvector(positions)
        error(['Positions must be specified as a vector.' ...
            ' Type "help colorbar" for further information']);
    elseif min(positions)<0 || max(positions)>1
        error(['Positions must range from 0 to 1 in an ascending order.' ...
            ' Type "help colorbar" for further information']);
    elseif length(positions) ~= length(unique(positions))
        error(['Check the positions vector, there are some duplicates.' ...
            ' Type "help colorbar" for further information']);
    else
        positions = sort(positions, 'ascend');
        if positions(1)~=0
            error(['The first positions must be 0.' ...
            ' Type "help colorbar" for further information']);
        elseif positions(length(positions))~=1
            error(['The last positions must be 1.' ...
            ' Type "help colorbar" for further information']);
        elseif length(positions) ~= n_colors
            error(['The number of positions does not match the number of colors.' ...
            ' Type "help colorbar" for further information']);
        end
    end

    % Convert HEX colors into RGB colors if required
    if iscell(colors)
        hex_colors = colors;
        colors = NaN(n_colors,3);
        for i = 1:n_colors
            colors(i,:) = hex2rgb(hex_colors{i});
        end
    end
    
    % Compute positions along the samples
    color_samples = round((m-1)*positions)+1;

    % Make the gradients among colors
    J = zeros(m,3);
    J(color_samples,:) = colors;
    diff_samples = diff(color_samples)-1;
    for d = 1:1:length(diff_samples)
        if diff_samples(d)~=0
            color1 = colors(d,:);
            color2 = colors(d+1,:);
            G = zeros(diff_samples(d),3);
            for idx_rgb = 1:3
                g = linspace(color1(idx_rgb), color2(idx_rgb), diff_samples(d)+2);
                g([1, length(g)]) = [];
                G(:,idx_rgb) = g';
            end
            J(color_samples(d)+1:color_samples(d+1)-1,:) = G;
        end
    end
    J = flipud(J);
end

% Function that converts an HEX string (format: #FFFFFF) to the
% corresponding RGB color
function rgb = hex2rgb(hexString)
	if size(hexString,2) ~= 7
		error('Not a color! %s', hexString);
	else
		r = double(hex2dec(hexString(2:3)))/255;
		g = double(hex2dec(hexString(4:5)))/255;
		b = double(hex2dec(hexString(6:7)))/255;
		rgb = [r, g, b];
	end
end

