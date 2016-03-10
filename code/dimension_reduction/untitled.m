for i =1:10
    for j = 1:10
        t = 10*(i-1)+j;
        temp_im = reshape(fea(t,:),16,16);
        M(((i-1)*16+1):i*16,((j-1)*16+1):j*16) = temp_im;
    end
end