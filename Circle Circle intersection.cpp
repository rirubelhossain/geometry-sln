/*
x^2 + y^2 = r1^2
(x - x2 )^2 + (y - y2 )^2 = r2^2
x (-2*x2 ) + y (-2*y2 ) + (x2^2 + y2^2 + r1^2 - r2^2 ) = 0

so, Ax + By + C = 0, A =-2x2 , B =-2y2 , C = x2^2 + y2^2 + r1^2 - r2^2
*/
double r, a, b, c; // ??????? ??????

double x0 = -a*c/(a*a+b*b),  y0 = -b*c/(a*a+b*b);
if (c*c > r*r*(a*a+b*b)+EPS)
	puts ("no points");
else if (abs (c*c - r*r*(a*a+b*b)) < EPS) {
	puts ("1 point");
	cout << x0 << ' ' << y0 << '\n';
}
else {
	double d = r*r - c*c/(a*a+b*b);
	double mult = sqrt (d / (a*a+b*b));
	double ax,ay,bx,by;
	ax = x0 + b * mult;
	bx = x0 - b * mult;
	ay = y0 - a * mult;
	by = y0 + a * mult;
	puts ("2 points");
	cout << ax << ' ' << ay << '\n' << bx << ' ' << by << '\n';
}