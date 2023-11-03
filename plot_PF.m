function [pf] = plot_PF(ori, odf, phase, oriColors, mrd_max)
    fprintf('Plotting %s pole figure...\n', phase)
        
    % plotting convention for rolling pole figures
    setMTEXpref('xAxisDirection','north');
    setMTEXpref('zAxisDirection','intoPlane');

    if contains(phase, 'alpha')
        %PF = figure();
        hkil = [Miller(0,0,0,2,odf.CS), Miller(1,0,-1,0,odf.CS), Miller(1,1,-2,0,odf.CS)]; % include hkil figures here
        plotPDF(odf, hkil,'antipodal', 'contourf', 0:0.1:mrd_max, 'minmax') % plot with contouring
%         plotPDF(ori, hkil, 'antipodal', 'property', oriColors, 'minmax'); % plot IPF color points
    elseif contains(phase, 'beta')
        %PF = figure();
        hkil = [Miller(0,0,1,odf.CS), Miller(1,1,0,odf.CS), Miller(1,1,1,odf.CS)]; % include hkil figures here
        plotPDF(odf, hkil,'antipodal', 'contourf', 0:0.1:mrd_max, 'minmax') % plot with contouring
%         plotPDF(ori, hkil, 'antipodal', 'property', oriColors, 'minmax'); % plot IPF color points
    end
    
    text(vector3d.X,'X','VerticalAlignment','bottom'); % moving the vector3d axis labels outside of the hemisphere boundary
    text(vector3d.Y,'Y','HorizontalAlignment','left');
    f = gcm; % moving up the hkil labels to make room for the rolling direction labels
    
    for i = 1:length(f.children)
        f.children(i).Title.Position=[1,1.25,1]; % use [-1,1.25,1] to make layout similar to Channel 5 and turn off minmax if you want to
        ax = getappdata(f.children(i),'sphericalPlot');
        maxval(i) = ax.maxData;
        ax.BL.Visible='off'; % omitt min
        
    end
    
    CLim(gcm, [0, mrd_max]) % define a colour giving the range (gcm, [min, max])
%     [rgb] = add_white(viridis,10); % Nick's colourbar
%     colormap(rgb);
    mtexColorbar ('location', 'southoutside', 'title', 'mrd'); % move colorbar to horizontal to avoid overlap
    set(gcf, 'PaperPositionMode', 'auto');
    %saveas (PF, output_filename, 'png');
    %close(PF);
end