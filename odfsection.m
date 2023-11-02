function [odfsection] = odfsection(ori, odf, slice, oriColors, mrd_max)

    fprintf('Plotting ODF section...')

    odf.SS=specimenSymmetry('orthorhombic');
    plotSection(odf, 'contourf', 'phi2', slice*degree, 'minmax');

%     [rgb] = add_white(viridis, 10);
%     colormap(rgb);
    CLim(gcm, [0 mrd_max])
    mtexColorbar ('location', 'southoutside', 'title', 'mrd', 'FontSize', 38);

    f = gcm;
    for i=1:length(f.children)
     sp = getappdata(f.children(i),'sphericalPlot');
     sp.TR.Position = [1.2, 1, 1];               % slice label position
     sp.TL.Position=[1.05, 0.3, 1];               % max label posn
     sp.BL.Visible='off';                         % get rid of min
    end
end