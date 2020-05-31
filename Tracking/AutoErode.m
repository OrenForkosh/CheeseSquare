function map = AutoErode(map, ErodeRadios)
% Erode each blob with a map of itself enlarged by 'ErodeRadios'
%
%       Created by OREN FORKOSH
%
reg = bwconncomp(map);
conn = regionprops(reg, 'BoundingBox');
for j=1:length(conn)
    img = map(ceil(conn(j).BoundingBox(2):conn(j).BoundingBox(2)+conn(j).BoundingBox(4)-1), ceil(conn(j).BoundingBox(1):conn(j).BoundingBox(1)+conn(j).BoundingBox(3)-1));
    img = imerode(padarray(img, [ErodeRadios ErodeRadios]), imresize(img, ErodeRadios * [2 2] + 1));
    img = img(ErodeRadios + 1:end - ErodeRadios, ErodeRadios + 1:end - ErodeRadios);
    map(ceil(conn(j).BoundingBox(2):conn(j).BoundingBox(2)+conn(j).BoundingBox(4)-1), ceil(conn(j).BoundingBox(1):conn(j).BoundingBox(1)+conn(j).BoundingBox(3)-1)) = img;
end
