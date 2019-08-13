function chs = get_soz_chs(pt,whichPt)
    chs = [];
    for j = 1:length(pt(whichPt).sz)
        elecs = pt(whichPt).sz(j).electrodes;
        for k = 1:length(elecs)
            for l = 1:length(pt(whichPt).new_elecs.names)
                if strcmp(elecs{k},pt(whichPt).new_elecs.names{l}) == 1
                    chs = [chs;l];
                end
            end
        end
    end
    
    chs = unique(chs);

end