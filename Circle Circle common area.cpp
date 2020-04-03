double get_area(double x1,double y1,double x2,double y2,double r1,double r2) {
	double	a1,a2;
	double	ta1,ta2;
	double	th1,th2;
	double	x,d,h;

	d = hypot(y2 - y1,x2 - x1);
	x = (r1 * r1 - r2 * r2 + d * d) / (2 * d);
	h = sqrt(r1 * r1 - x * x);

	th1 = atan2(h,x);
	th2 = atan2(h,d - x);

	ta1 = 0.5 * h * x;
	ta2 = 0.5 * h * (d - x);

	a1 = 0.5 * th1 * r1 * r1;
	a2 = 0.5 * th2 * r2 * r2;

	a1 -= ta1;
	a2 -= ta2;

	return 2 * (a1 + a2);
}
