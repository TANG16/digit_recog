TestRG = V(:,1:40)*Test_low'+ repmat(mean(Test', 1)',1,1000);
for i = 0:9
    p = TestRG(:,i*100+1);
    start = mod(i,5);
    start2 = floor(i/5);
    toTest2((start2*28+1):(start2*28+28),start*28+1:(start+1)*28) = reshape(p,28,28);
end
I = uint8((1-toTest2)*255);
h = imagesc(I);
axis equal tight 
axis off

% eorigin = 1-size(find(predictlabel_origin~=GroundTruth))/1000
% e_lle = 1-size(find(predictlabel_lle~=GroundTruth))/1000
% epca = 1-size(find(predictlabel_pca~=GroundTruth))/1000
% epoly = 1-size(find(predictlabel_poly~=GroundTruth))/1000
% egau = 1-size(find(predictlabel_gaussian~=GroundTruth))/1000
% for i = 0:9
%     eorigin_on(i+1,:) = 1-size(find(predictlabel_origin((i*100+1):(i+1)*100)~=GroundTruth((i*100+1):(i+1)*100)))/100
%     e_lle_on(i+1,:) = 1-size(find(predictlabel_lle((i*100+1):(i+1)*100)~=GroundTruth((i*100+1):(i+1)*100)))/100
%     epca_on(i+1,:) = 1-size(find(predictlabel_pca((i*100+1):(i+1)*100)~=GroundTruth((i*100+1):(i+1)*100)))/100
%     epoly_on(i+1,:) = 1-size(find(predictlabel_poly((i*100+1):(i+1)*100)~=GroundTruth((i*100+1):(i+1)*100)))/100
%     egau_on(i+1,:) = 1-size(find(predictlabel_gaussian((i*100+1):(i+1)*100)~=GroundTruth((i*100+1):(i+1)*100)))/100
% end
% eorigin_ev = mean(eorigin_on);
% e_lle_on_ev = mean(e_lle_on);
% epca_on_ev = mean(epca_on);
% epoly_on_ev = mean(epoly_on);
% egau_on_ev = mean(egau_on);