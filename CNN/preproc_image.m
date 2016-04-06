function out = preproc_image(Inorm)
    %Preprocess single image
    Inorm(~isfinite(Inorm)) = 1;
    Inorm = abs(Inorm-1)';
    out = zeros(32);
    out(3:30,3:30) = Inorm;
    if sum(out(:))>0,
        gain = 1./ std(out(:));
        out = (out - mean(out(:))).*gain;        
    end
end
