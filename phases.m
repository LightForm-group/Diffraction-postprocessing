function [cs] = phases(phase)
    
    if contains(phase, "alpha") % HCP phases %        
        % IMPORTANT NOTE: DAMASK uses 10-10||Y and 11-20||X %
        cs = crystalSymmetry('6/mmm', [3 3 4.7], 'X||b', 'Y||a*', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0.53 0.81 0.98]);
    end
        
%     if isequal(phase, 'Mg') 
%         w = miller_indicies(4);
%         m_str = strcat(m_str, string(w));
%         %cs = 'INSERT Mg cs HERE';
%     end
        
    if contains(phase, 'beta')
        cs = crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0.56 0.74 0.56]);
    end

    if isequal(phase, 'Ni-gamma')
        cs = crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni-superalloy', 'color', [0.53 0.81 0.98]); % FCC
    end
    
    if isequal(phase, 'Al')
        cs = crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]); % FCC
    end
    
    fprintf("Got parameters for phase: %s\n", phase)
end

% Notes: %
% crystalSymmetry('laue group', [lattice parameters(unit cell size)], axes convention, 'mineral', 'ebsddatacolor', [R G B]);