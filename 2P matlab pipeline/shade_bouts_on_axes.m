function shade_bouts_on_axes(ax, bouts, colorRGB, alphaVal)
% Shade [on_s, off_s] time windows on existing axes.
    if isempty(bouts) || ~ishandle(ax), return; end
    yL = ylim(ax);
    for i = 1:size(bouts,1)
        x1 = bouts(i,1); x2 = bouts(i,2);
        if x2 <= x1, continue; end
        patch('Parent', ax, ...
              'XData',[x1 x2 x2 x1], ...
              'YData',[yL(1) yL(1) yL(2) yL(2)], ...
              'FaceColor', colorRGB, 'FaceAlpha', alphaVal, ...
              'EdgeColor','none', 'HitTest','off', 'PickableParts','none');
    end
end