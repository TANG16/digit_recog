% HOG for digit
function Features = Hog_digit(I,c,r);
% I:  pXn matrix
% c,r: column and row of the image.
I = reshape(I, c, r,[]);
[~, ~, n ] = size(I);
for i = 1:n
    hog = vl_hog(single(I(:,:,i)), 4) ;
    Features(i,:) = hog(:);
end
    


