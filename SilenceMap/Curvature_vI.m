function [curv] = Curvature_vI(x,y)

dx  = gradient((x));
ddx = gradient(dx);
dy  = gradient((y));
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom .* denom .* denom;
curv = abs(num) ./ denom;
curv(denom == 0) = 0;


end