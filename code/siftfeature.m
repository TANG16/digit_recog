function [f, d] = siftfeature(I, sizeImage)
    peakthresh = 5;
    edgethresh = 5;
    o = -1;
    I = single(reshape(I, sizeImage)) * 500;
    [f, d] = vl_sift(I, 'peakthresh', peakthresh, 'edgethresh', edgethresh, ...
                    'FirstOctave', o);
end