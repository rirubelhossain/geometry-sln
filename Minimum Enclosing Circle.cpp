//O(n) algorithm not sure
circle MINCIRCLE (Set I, Set B)
	If | I | == 0 or | B |> = 3
		returns the circle formed by the points of B
	
	random point belonging to the set P = I
	C = MINCIRCLE (I - {P}, B)
	If P is not included in C
		returns MINCIRCLE (I - {P}, B + {P})
	returns C


//O(nlgn)
by voronoi diagram

//O(n^2)
Select a random center,
Shrink the circle towards center so that it touches a point.
Shrink the circle towards the point if there is not 2 points on side already
Now shrink circle towards the edge
if all the points on the side are less than pi gap we are done
otherwise select such pair and again shrink toward it.
if the center is on the edge then break it.

//Implementation
#include <stdio.h>
#include <math.h>
const int maxn = 1005;
//const double eps = 1e-6;
struct TPoint
{
double x, y;
TPoint operator-(TPoint &a)
{
TPoint p1;
p1.x = x - a.x;
p1.y = y - a.y;
return p1;
}
};
struct TCircle
{
double r;
TPoint centre;
};
struct TTriangle
{
TPoint t[3];
};
TCircle c;
TPoint a[maxn];
double distance(TPoint p1, TPoint p2)
{
TPoint p3;
p3.x = p2.x - p1.x;
p3.y = p2.y - p1.y;
return sqrt(p3.x * p3.x + p3.y * p3.y);
}
double triangleArea(TTriangle t)
{
TPoint p1, p2;
p1 = t.t[1] - t.t[0];
p2 = t.t[2] - t.t[0];
return fabs(p1.x * p2.y - p1.y * p2.x) / 2;
}
TCircle circumcircleOfTriangle(TTriangle t)
{
//�����ε����Բ
TCircle tmp;
double a, b, c, c1, c2;
double xA, yA, xB, yB, xC, yC;
a = distance(t.t[0], t.t[1]);
b = distance(t.t[1], t.t[2]);
c = distance(t.t[2], t.t[0]);
//����S = a * b * c / R / 4;��뾶R
tmp.r = a * b * c / triangleArea(t) / 4;
xA = t.t[0].x; yA = t.t[0].y;
xB = t.t[1].x; yB = t.t[1].y;
xC = t.t[2].x; yC = t.t[2].y;
c1 = (xA * xA + yA * yA - xB * xB - yB * yB) / 2;
c2 = (xA * xA + yA * yA - xC * xC - yC * yC) / 2;
tmp.centre.x = (c1 * (yA - yC) - c2 * (yA - yB)) /
((xA - xB) * (yA - yC) - (xA - xC) * (yA - yB));
tmp.centre.y = (c1 * (xA - xC) - c2 * (xA - xB)) /
((yA - yB) * (xA - xC) - (yA - yC) * (xA - xB));
return tmp;
}
TCircle MinCircle2(int tce, TTriangle ce)
{
TCircle tmp;
if(tce == 0) tmp.r = -2;
else if(tce == 1)
{
tmp.centre = ce.t[0];
tmp.r = 0;
}
else if(tce == 2)
{
tmp.r = distance(ce.t[0], ce.t[1]) / 2;
tmp.centre.x = (ce.t[0].x + ce.t[1].x) / 2;
tmp.centre.y = (ce.t[0].y + ce.t[1].y) / 2;
}
else if(tce == 3) tmp = circumcircleOfTriangle(ce);
return tmp;
}
void MinCircle(int t, int tce, TTriangle ce)
{
int i, j;
TPoint tmp;
c = MinCircle2(tce, ce);
if(tce == 3) return;
for(i = 1;i <= t;i++)
{
if(distance(a[i], c.centre) > c.r)
{
ce.t[tce] = a[i];
MinCircle(i - 1, tce + 1, ce);
tmp = a[i];
for(j = i;j >= 2;j--)
{
a[j] = a[j - 1];
}
a[1] = tmp;
}
}
}
void run(int n)
{
TTriangle ce;
int i;
MinCircle(n, 0, ce);
printf("%.2lf %.2lf %.2lf\n", c.centre.x, c.centre.y, c.r);
}
int main()
{
freopen("circle.in", "r", stdin);
freopen("out.txt", "w", stdout);
int n;
while(scanf("%d", &n) != EOF && n)
{
for(int i = 1;i <= n;i++)
scanf("%lf%lf", &a[i].x, &a[i].y);
run(n);
}
return 0;
}
