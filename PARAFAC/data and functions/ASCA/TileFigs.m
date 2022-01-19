function tilefigs(figs,nrows,ncols,border)
%   function tilefigs(figs,nrows,ncols,border)
%
%   Purpose:
%   Tile figure windows in screen.
%
%   Parameters:
%   figs    -   (optional)
%               Array with numbers of figures to tile.
%               Default: all open windows.
%   nrows   -   (optional)
%               Number of rows in tiling.
%               Default: square root of number figures to tile.
%   ncols   -   (optional)
%               Number of columns in tiling.
%               Default: square root of number figures to tile or
%                         large enough to plot all figures.
%   border  -   (optional)
%               Size of border of tile in pixels.
%               Default:0
%
%   References:
%
%   Original version by:
%   Charles Plum                    Nichols Research Corp.
%   <cplum@nichols.com>             70 Westview Street
%   Tel: (781) 862-9400             Kilnbrook IV
%   Fax: (781) 862-9485             Lexington, MA 02173

if ( nargin < 1 )
    hands   = get (0,'Children');   % locate all open figure handles
else
    if ( isempty(figs))
        hands   = get (0,'Children');   % locate all open figure handles
    else
        hands   = figs(:);
    end    
end
hands   = sort(hands);          % sort the figure handles.
numfigs = size(hands,1);        % number of figures to tile
maxfigs = 100;                  % Maximum number of figure to tile.
maxpos  = get (0,'screensize'); % Determine size of screen in pixels
maxpos(4) = maxpos(4) - 25;

if (numfigs > maxfigs)            % figure limit check
        error(sprinf('More than %3d figures requested',maxfigs))
end
if nargin < 2
    maxfactor   = sqrt(maxfigs);                % max number of figures per row or column
    sq          = [2:maxfactor].^2;             % vector of integer squares
    sq          = sq(find(sq>=numfigs));        % determine square grid size
    gridsize    = sq(1);                        % best grid size
    nrows       = sqrt(gridsize);               % figure size screen scale factor
    ncols       = nrows;                        % figure size screen scale factor
elseif nargin < 3 
    ncols       = ceil(numfigs / nrows);
end

if ( nargin < 4 )
  border = 0;
else
  maxpos(3) = maxpos(3) - 2*border;
  maxpos(4) = maxpos(4) - 2*border;
end
xlen = fix(maxpos(3)/ncols) - 30; % new tiled figure width
ylen = fix(maxpos(4)/nrows) - 45; % new tiled figure height

% Ttile figures by position 
% Location (1,1) is at bottom left corner of screen
pnum    = 0 ;
for ( iy = 1:nrows )
    ypos  = maxpos(4) - fix((iy)*maxpos(4)/nrows) + border +25; % figure location (row)
    for ( ix = 1:ncols )
        xpos = fix((ix-1)*maxpos(3)/ncols + 1) + border+7;     % figure location (column)
        pnum = pnum+1;
        if (pnum>numfigs)
            break
        else
            figure(hands(pnum))
            set(hands(pnum),'Position',[ xpos ypos xlen ylen ]); % move figure
        end
  end
end