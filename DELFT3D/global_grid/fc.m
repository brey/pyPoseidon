function a = fc(A,b,c)
 %This function determines the length of a side of a spherical triangle. It
 %It is in reference to a triangle with vertices (A,B,C) and opposite sides
 %(a,b,c).

 a=acos(cos(A).*sin(b).*sin(c)+cos(b).*cos(c));

 end