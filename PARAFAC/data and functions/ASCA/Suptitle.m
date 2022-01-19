function hout=suptitle(str,fontsize)
% [hout] = suptitle(str,fontsize)
%
% Purpose:
% Puts a title above all subplots in a figure window.
%
% Parameters:
% str       - String to plot
%
% Description:
% suptitle('text') adds text to the top of the figure
% above all subplots (a "super title"). This function
% after all subplot commands have been issued.
%
% Warning:
% If the figure or axis units are non-default, this
% rouine will break down.
%
% References:
% Original version by:
% Drea Thomas 6/15/95 drea@mathworks.com
%
if ( nargin < 2 )
    % Fontsize for supertitle
    fontsize = get(gcf,'defaultaxesfontsize')+4;
end
if ( isempty(fontsize) )
    fontsize = get(gcf,'defaultaxesfontsize')+4;
end

% Parameters used to position the supertitle.

    % Amount of the figure window devoted to subplots
    plotregion = .92;

    % Y position of title in normalized coordinates
    titleypos  = .95;


% Fudge factor to adjust y spacing between subplots
fudge       = 1;
haold       = gca;
figunits    = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

    if (~strcmp(figunits,'pixels')),
        set(gcf,'units','pixels');
        pos = get(gcf,'position');
        set(gcf,'units',figunits);
    else,
        pos = get(gcf,'position');
    end
    ff = (fontsize-4)*1.27*5/pos(4)*fudge;
        % The 5 here reflects about 3 characters of height below
        % an axis and 2 above. 1.27 is pixels per point.
        
% Determine the bounding rectange for all the plots
    ch = get(gcf,'children');
    h=[];
    for i=1:length(ch),
        if strcmp(get(ch(i),'type'),'axes'),
            h=[h,ch(i)];
        end
    end
    max_y=0;
    min_y=1;

    oldtitle =0;
    for i=1:length(h),
        if (~strcmp(get(h(i),'Tag'),'suptitle')),
            pos = get(h(i),'pos');
            if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
            if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
        else,
            oldtitle = h(i);
        end
    end

    if max_y > plotregion,
        scale = (plotregion-min_y)/(max_y-min_y);
        for i=1:length(h),
            pos = get(h(i),'position');
            pos(2) = (pos(2)-min_y)*scale+min_y;
            pos(4) = pos(4)*scale-(1-scale)*ff/5*3;
            set(h(i),'position',pos);
        end
    end

    np  = get(gcf,'nextplot');
    set(gcf,'nextplot','add');
    if (oldtitle),
        delete(oldtitle);
    end
    ha  = axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
    ht  = text(.5,titleypos-1,str);set(ht,'horizontalalignment','center','fontsize',fontsize);
    set(gcf,'nextplot',np);
    axes(haold);
    if nargout,
        hout=ht;
    end
