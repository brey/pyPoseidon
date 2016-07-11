function B = ft(A,b,c)
 %This function calculates the vertice of a triangle with the spherical tangent.
 %It uses as reference a triangle with vertices (A,B,C) and opposite sides (a,b,c).

 B=atan2(sin(A),sin(c)./tan(b)-cos(c).*cos(A));

 end