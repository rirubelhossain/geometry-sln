struct point { double x, y; };
point bcenter(point pnt[], int n){
	point p, s;
	double tp, area = 0, tpx = 0, tpy = 0;
	p.x = pnt[0].x; p.y = pnt[0].y;
	for (int i = 1; i <= n; ++i) { // point: 0 ~ n-1
		s.x = pnt[(i == n) ? 0 : i].x;
		s.y = pnt[(i == n) ? 0 : i].y;
		tp = (p.x * s.y - s.x * p.y); area += tp / 2;
		tpx += (p.x + s.x) * tp; tpy += (p.y + s.y) * tp;
		p.x = s.x; p.y = s.y;
	}
	s.x = tpx / (6 * area); s.y = tpy / (6 * area);
	return s;
}

///Filled polygon
#include <stdio.h>
#include <math.h>
typedef struct TPoint
{
double x;
double y;
}TPoint;
double triangleArea(TPoint p0, TPoint p1, TPoint p2)
{
//已知三角形三个顶点的坐标，求三角形的面积
double k = p0.x * p1.y + p1.x * p2.y
+ p2.x * p0.y - p1.x * p0.y
- p2.x * p1.y - p0.x * p2.y;
//if(k >= 0) return k / 2;
// else return -k / 2;
return k / 2;
}
int main()
{
int i, n, test;
TPoint p0, p1, p2, center;
double area, sumarea, sumx, sumy;
scanf("%d", &test);
while(test--){
scanf("%d", &n);
scanf("%lf%lf", &p0.x, &p0.y);
scanf("%lf%lf", &p1.x, &p1.y);
sumx = 0;
sumy = 0;
sumarea = 0;
for(i = 2;i < n;i++){
scanf("%lf%lf", &p2.x, &p2.y);
center.x = p0.x + p1.x + p2.x;
center.y = p0.y + p1.y + p2.y;
area = triangleArea(p0, p1, p2);
sumarea += area;
sumx += center.x * area;
sumy += center.y * area;
p1 = p2;
}
printf("%.2lf %.2lf\n", sumx / sumarea / 3, sumy / sumarea / 3);
}
return 0;
}