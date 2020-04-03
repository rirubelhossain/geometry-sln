//Polygon
#include <stdlib.h>
#include <math.h>
#define MAXN 1000
#define offset 10000
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
#define _sign(x) ((x)>eps?1:((x)<-eps?2:0))
struct point{double x,y;};
struct line{point a,b;};
double xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
//Determine the convex polygon, vertex in a clockwise or counterclockwise to allow adjacent edges sharing the same line
int is_convex(int n,point* p){
int i,s[3]={1,1,1};
for (i=0;i<n&&s[1]|s[2];i++)
s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
return s[1]|s[2];
}
//Determine the convex polygon, vertex in a clockwise or counterclockwise to not allow adjacent edges sharing the same line
int is_convex_v2(int n,point* p){
int i,s[3]={1,1,1};
for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
s[_sign(xmult(p[(i+1)%n],p[(i+2)%n],p[i]))]=0;
return s[0]&&s[1]|s[2];
}
//Award points in a convex polygon or polygon edge, vertex in a clockwise or counterclockwise to
int inside_convex(point q,int n,point* p){
int i,s[3]={1,1,1};
for (i=0;i<n&&s[1]|s[2];i++)
s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
return s[1]|s[2];
}
//Liable to a point within the convex polygon, vertex in a clockwise or counterclockwise, on the edge of the polygon returns 0
int inside_convex_v2(point q,int n,point* p){
int i,s[3]={1,1,1};
for (i=0;i<n&&s[0]&&s[1]|s[2];i++)
s[_sign(xmult(p[(i+1)%n],q,p[i]))]=0;
return s[0]&&s[1]|s[2];
}
//Found within the freeform, vertices in a clockwise or counterclockwise to
//On_edge value that represents the point in polygon edge returns, offset maximum polygon coordinatesint inside_polygon(point q,int n,point* p,int on_edge=1){
point q2;
int i=0,count;
while (i<n)
for (count=i=0,q2.x=rand()+offset,q2.y=rand()+offset;i<n;i++)
if
(zero(xmult(q,p[i],p[(i+1)%n]))&&(p[i].x-q.x)*(p[(i+1)%n].x-q.x)<eps&&(p[i].y-q.y)*(p[(i+1)%
n].y-q.y)<eps)
return on_edge;
else if (zero(xmult(q,q2,p[i])))
break;
else if
(xmult(q,p[i],q2)*xmult(q,p[(i+1)%n],q2)<-eps&&xmult(p[i],q,p[(i+1)%n])*xmult(p[i],q2,p[(i+1)
%n])<-eps)
count++;
return count&1;
}
inline int opposite_side(point p1,point p2,point l1,point l2){
return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}
inline int dot_online_in(point p,point l1,point l2){
return zero(xmult(p,l1,l2))&&(l1.x-p.x)*(l2.x-p.x)<eps&&(l1.y-p.y)*(l2.y-p.y)<eps;
}
//Sentenced segments in a freeform, vertices in a clockwise or counterclockwise, cross the border to return 1
int inside_polygon(point l1,point l2,int n,point* p){
point t[MAXN],tt;
int i,j,k=0;
if (!inside_polygon(l1,n,p)||!inside_polygon(l2,n,p))
return 0;
for (i=0;i<n;i++)
if (opposite_side(l1,l2,p[i],p[(i+1)%n])&&opposite_side(p[i],p[(i+1)%n],l1,l2))
return 0;
else if (dot_online_in(l1,p[i],p[(i+1)%n]))
t[k++]=l1;
else if (dot_online_in(l2,p[i],p[(i+1)%n]))
t[k++]=l2;
else if (dot_online_in(p[i],l1,l2))
t[k++]=p[i];
for (i=0;i<k;i++)
for (j=i+1;j<k;j++){
tt.x=(t[i].x+t[j].x)/2;
tt.y=(t[i].y+t[j].y)/2;
if (!inside_polygon(tt,n,p))
return 0;
}
return 1;
}
point intersection(line u,line v){
	point ret=u.a;
double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
ret.x+=(u.b.x-u.a.x)*t;
ret.y+=(u.b.y-u.a.y)*t;
return ret;
}
point barycenter(point a,point b,point c){
line u,v;
u.a.x=(a.x+b.x)/2;
u.a.y=(a.y+b.y)/2;
u.b=c;
v.a.x=(a.x+c.x)/2;
v.a.y=(a.y+c.y)/2;
v.b=b;
return intersection(u,v);
}
//Polygon Center of gravity
point barycenter(int n,point* p){
point ret,t;
double t1=0,t2;
int i;
ret.x=ret.y=0;
for (i=1;i<n-1;i++)
if (fabs(t2=xmult(p[0],p[i],p[i+1]))>eps){
t=barycenter(p[0],p[i],p[i+1]);
ret.x+=t.x*t2;
ret.y+=t.y*t2;
t1+=t2;
}
if (fabs(t1)>eps)
ret.x/=t1,ret.y/=t1;
return ret;
}




//Polygon cutting
//For half-plane
#define MAXN 100
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
struct point{double x,y;};
double xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
int same_side(point p1,point p2,point l1,point l2){
return xmult(l1,p1,l2)*xmult(l1,p2,l2)>eps;
}
point intersection(point u1,point u2,point v1,point v2){
point ret=u1;
double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
ret.x+=(u2.x-u1.x)*t;
ret.y+=(u2.y-u1.y)*t;
return ret;
}
//Polygon side cutting of along straight lines cut in the side of the L1,L2 to determine, ensure that L1,L2,side are not collinear
void polygon_cut(int& n,point* p,point l1,point l2,point side){
point pp[100];
int m=0,i;
for (i=0;i<n;i++){
if (same_side(p[i],side,l1,l2))
pp[m++]=p[i];
if
(!same_side(p[i],p[(i+1)%n],l1,l2)&&!(zero(xmult(p[i],l1,l2))&&zero(xmult(p[(i+1)%n],l1,l2))))
pp[m++]=intersection(p[i],p[(i+1)%n],l1,l2);
}
for (n=i=0;i<m;i++)
if (!i||!zero(pp[i].x-pp[i-1].x)||!zero(pp[i].y-pp[i-1].y))
p[n++]=pp[i];
if (zero(p[n-1].x-p[0].x)&&zero(p[n-1].y-p[0].y))
n--;
if (n<3)
n=0;
}


//Floating-point functions
#include <math.h>
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
struct point{double x,y;};
struct line{point a,b;};
//Calculation cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
double xmult(double x1,double y1,double x2,double y2,double x0,double y0){
return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}
//Calculation dot product (P1-P0).(P2-P0)
double dmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.x-p0.x)+(p1.y-p0.y)*(p2.y-p0.y);
}
double dmult(double x1,double y1,double x2,double y2,double x0,double y0){
return (x1-x0)*(x2-x0)+(y1-y0)*(y2-y0);
}
//Two point distances
double distance(point p1,point p2){
return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
double distance(double x1,double y1,double x2,double y2){
return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
//Three points collinear
int dots_inline(point p1,point p2,point p3){
return zero(xmult(p1,p2,p3));
}
int dots_inline(double x1,double y1,double x2,double y2,double x3,double y3){
return zero(xmult(x1,y1,x2,y2,x3,y3));
}
//Liable to a point is on the line, including the endpoint
int dot_online_in(point p,line l){
return zero(xmult(p,l.a,l.b))&&(l.a.x-p.x)*(l.b.x-p.x)<eps&&(l.a.y-p.y)*(l.b.y-p.y)<eps;
}
int dot_online_in(point p,point l1,point l2){
return zero(xmult(p,l1,l2))&&(l1.x-p.x)*(l2.x-p.x)<eps&&(l1.y-p.y)*(l2.y-p.y)<eps;
}
int dot_online_in(double x,double y,double x1,double y1,double x2,double y2){
	return zero(xmult(x,y,x1,y1,x2,y2))&&(x1-x)*(x2-x)<eps&&(y1-y)*(y2-y)<eps;
}
//Liable to a point is on the segment, excluding endpoints
int dot_online_ex(point p,line l){
return
dot_online_in(p,l)&&(!zero(p.x-l.a.x)||!zero(p.y-l.a.y))&&(!zero(p.x-l.b.x)||!zero(p.y-l.b.y));
}
int dot_online_ex(point p,point l1,point l2){
return
dot_online_in(p,l1,l2)&&(!zero(p.x-l1.x)||!zero(p.y-l1.y))&&(!zero(p.x-l2.x)||!zero(p.y-l2.y));
}
int dot_online_ex(double x,double y,double x1,double y1,double x2,double y2){
return
dot_online_in(x,y,x1,y1,x2,y2)&&(!zero(x-x1)||!zero(y-y1))&&(!zero(x-x2)||!zero(y-y2));
}
//Two points online with the side, points on the line return 0
int same_side(point p1,point p2,line l){
return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)>eps;
}
int same_side(point p1,point p2,point l1,point l2){
return xmult(l1,p1,l2)*xmult(l1,p2,l2)>eps;
}
//ISO side two point in the segment, segment returns 0
int opposite_side(point p1,point p2,line l){
return xmult(l.a,p1,l.b)*xmult(l.a,p2,l.b)<-eps;
}
int opposite_side(point p1,point p2,point l1,point l2){
return xmult(l1,p1,l2)*xmult(l1,p2,l2)<-eps;
}
//Two lines run parallel
int parallel(line u,line v){
return zero((u.a.x-u.b.x)*(v.a.y-v.b.y)-(v.a.x-v.b.x)*(u.a.y-u.b.y));
}
int parallel(point u1,point u2,point v1,point v2){
return zero((u1.x-u2.x)*(v1.y-v2.y)-(v1.x-v2.x)*(u1.y-u2.y));
}
//Two linear vertical
int perpendicular(line u,line v){
return zero((u.a.x-u.b.x)*(v.a.x-v.b.x)+(u.a.y-u.b.y)*(v.a.y-v.b.y));
}
int perpendicular(point u1,point u2,point v1,point v2){
return zero((u1.x-u2.x)*(v1.x-v2.x)+(u1.y-u2.y)*(v1.y-v2.y));
}
//Two line segments intersect, including endpoint and partial coincidence
int intersect_in(line u,line v){
if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b))
return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
return dot_online_in(u.a,v)||dot_online_in(u.b,v)||dot_online_in(v.a,u)||dot_online_in(v.b,u);
}
int intersect_in(point u1,point u2,point v1,point v2){
if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
return
dot_online_in(u1,v1,v2)||dot_online_in(u2,v1,v2)||dot_online_in(v1,u1,u2)||dot_online_in(v2,u1,u
2);
}
//Two line segments intersect, not including endpoint and partial coincidence
int intersect_ex(line u,line v){
return opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}
int intersect_ex(point u1,point u2,point v1,point v2){
return opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
//Calculate the intersection of two lines, attention to determine whether lines are parallel in advance!
//Intersects a line segment intersection sentenced segment (also was to determine whether parallel!)
point intersection(line u,line v){
point ret=u.a;
double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
ret.x+=(u.b.x-u.a.x)*t;
ret.y+=(u.b.y-u.a.y)*t;
return ret;
}
point intersection(point u1,point u2,point v1,point v2){
point ret=u1;
double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
ret.x+=(u2.x-u1.x)*t;
ret.y+=(u2.y-u1.y)*t;
return ret;
}
//Point to the nearest point in a straight line
point ptoline(point p,line l){
point t=p;
t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
return intersection(p,t,l.a,l.b);
}
point ptoline(point p,point l1,point l2){
point t=p;
t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
return intersection(p,t,l1,l2);
}
//Point to straight line distance
double disptoline(point p,line l){
return fabs(xmult(p,l.a,l.b))/distance(l.a,l.b);
}
double disptoline(point p,point l1,point l2){
return fabs(xmult(p,l1,l2))/distance(l1,l2);
}
double disptoline(double x,double y,double x1,double y1,double x2,double y2){
return fabs(xmult(x,y,x1,y1,x2,y2))/distance(x1,y1,x2,y2);
}
//Point to the nearest point in the segment
point ptoseg(point p,line l){
point t=p;
t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
if (xmult(l.a,t,p)*xmult(l.b,t,p)>eps)
return distance(p,l.a)<distance(p,l.b)?l.a:l.b;
return intersection(p,t,l.a,l.b);
}
point ptoseg(point p,point l1,point l2){
point t=p;
t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
return distance(p,l1)<distance(p,l2)?l1:l2;
return intersection(p,t,l1,l2);
}
//Distance from points to lines
double disptoseg(point p,line l){
point t=p;
t.x+=l.a.y-l.b.y,t.y+=l.b.x-l.a.x;
if (xmult(l.a,t,p)*xmult(l.b,t,p)>eps)
return distance(p,l.a)<distance(p,l.b)?distance(p,l.a):distance(p,l.b);
return fabs(xmult(p,l.a,l.b))/distance(l.a,l.b);
}
double disptoseg(point p,point l1,point l2){
point t=p;
t.x+=l1.y-l2.y,t.y+=l2.x-l1.x;
if (xmult(l1,t,p)*xmult(l2,t,p)>eps)
return distance(p,l1)<distance(p,l2)?distance(p,l1):distance(p,l2);
return fabs(xmult(p,l1,l2))/distance(l1,l2);
}
//Vector v is a vertex with p counterclockwise angle and zoom scale times
point rotate(point v,point p,double angle,double scale){
point ret=p;
v.x-=p.x,v.y-=p.y;
p.x=scale*cos(angle);
p.y=scale*sin(angle);
ret.x+=v.x*p.x-v.y*p.y;
ret.y+=v.x*p.y+v.y*p.x;
return ret;
}




//Area
#include <math.h>
struct point{double x,y;};
//cross product (P1-P0)x(P2-P0)
double xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
double xmult(double x1,double y1,double x2,double y2,double x0,double y0){
return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}
double area_triangle(point p1,point p2,point p3){
return fabs(xmult(p1,p2,p3))/2;
}
double area_triangle(double x1,double y1,double x2,double y2,double x3,double y3){
return fabs(xmult(x1,y1,x2,y2,x3,y3))/2;
}
double area_triangle(double a,double b,double c){
double s=(a+b+c)/2;
return sqrt(s*(s-a)*(s-b)*(s-c));
}
double area_polygon(int n,point* p){
double s1=0,s2=0;
int i;
for (i=0;i<n;i++)
s1+=p[(i+1)%n].y*p[i].x,s2+=p[(i+1)%n].y*p[(i+2)%n].x;
return fabs(s1-s2)/2;
}




//Spherical
#include <math.h>
const double pi=acos(-1);

//Computing Central angle LAT latitude, -90<=w<=90,LNG longitude
//Returned corresponds to the two points where the big minor arc center point, 0<=angle<=PI
double angle(double lng1,double lat1,double lng2,double lat2){
double dlng=fabs(lng1-lng2)*pi/180;
while (dlng>=pi+pi)
dlng-=pi+pi;
if (dlng>pi)
dlng=pi+pi-dlng;
lat1*=pi/180,lat2*=pi/180;
return acos(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2));
}
//Calculate the distance, r-ball RADIUS
double line_dist(double r,double lng1,double lat1,double lng2,double lat2){
double dlng=fabs(lng1-lng2)*pi/180;
while (dlng>=pi+pi)
dlng-=pi+pi;
if (dlng>pi)
dlng=pi+pi-dlng;
lat1*=pi/180,lat2*=pi/180;
return r*sqrt(2-2*(cos(lat1)*cos(lat2)*cos(dlng)+sin(lat1)*sin(lat2)));
}
//Calculating spherical distance, r-ball RADIUS
inline double sphere_dist(double r,double lng1,double lat1,double lng2,double lat2){
return r*angle(lng1,lat1,lng2,lat2);
}



//Triangle
#include <math.h>
struct point{double x,y;};
struct line{point a,b;};
double distance(point p1,point p2){
return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
point intersection(line u,line v){
point ret=u.a;
double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
ret.x+=(u.b.x-u.a.x)*t;
ret.y+=(u.b.y-u.a.y)*t;
return ret;
}
//Circumcenter
point circumcenter(point a,point b,point c){
line u,v;
u.a.x=(a.x+b.x)/2;
u.a.y=(a.y+b.y)/2;
u.b.x=u.a.x-a.y+b.y;
u.b.y=u.a.y+a.x-b.x;
v.a.x=(a.x+c.x)/2;
v.a.y=(a.y+c.y)/2;
v.b.x=v.a.x-a.y+c.y;
v.b.y=v.a.y+a.x-c.x;
return intersection(u,v);
}
//Incenter
point incenter(point a,point b,point c){
line u,v;
double m,n;
u.a=a;
m=atan2(b.y-a.y,b.x-a.x);
n=atan2(c.y-a.y,c.x-a.x);
u.b.x=u.a.x+cos((m+n)/2);
u.b.y=u.a.y+sin((m+n)/2);
v.a=b;
m=atan2(a.y-b.y,a.x-b.x);
n=atan2(c.y-b.y,c.x-b.x);
v.b.x=v.a.x+cos((m+n)/2);
v.b.y=v.a.y+sin((m+n)/2);
return intersection(u,v);
}

point perpencenter(point a,point b,point c){
line u,v;
u.a=c;
u.b.x=u.a.x-a.y+b.y;
u.b.y=u.a.y+a.x-b.x;
v.a=b;
v.b.x=v.a.x-a.y+c.y;
v.b.y=v.a.y+a.x-c.x;
return intersection(u,v);
}

//Focus
//Three vertices of a triangle to the square of the distance, and the smallest of points
//Triangle triangular plot of maximum points
point barycenter(point a,point b,point c){
line u,v;
u.a.x=(a.x+b.x)/2;
u.a.y=(a.y+b.y)/2;
u.b=c;
v.a.x=(a.x+c.x)/2;
v.a.y=(a.y+c.y)/2;
v.b=b;
return intersection(u,v);
}
//Horse
//And smallest of the three vertices of a triangle from the point
point fermentpoint(point a,point b,point c){
point u,v;
double step=fabs(a.x)+fabs(a.y)+fabs(b.x)+fabs(b.y)+fabs(c.x)+fabs(c.y);
int i,j,k;
u.x=(a.x+b.x+c.x)/3;
u.y=(a.y+b.y+c.y)/3;
while (step>1e-10)
for (k=0;k<10;step/=2,k++)
for (i=-1;i<=1;i++)
for (j=-1;j<=1;j++){
v.x=u.x+step*i;
v.y=u.y+step*j;
if
(distance(u,a)+distance(u,b)+distance(u,c)>distance(v,a)+distance(v,b)+distance(v,c))
u=v;
}
return u;
}





//Three-dimensional geometry
#include <math.h>
#define eps 1e-8
#define zero(x) (((x)>0?(x):-(x))<eps)
struct point3{double x,y,z;};
struct line3{point3 a,b;};
struct plane3{point3 a,b,c;};

//Cross
point3 xmult(point3 u,point3 v){
point3 ret;
ret.x=u.y*v.z-v.y*u.z;
ret.y=u.z*v.x-u.x*v.z;
ret.z=u.x*v.y-u.y*v.x;
return ret;
}
//dot
double dmult(point3 u,point3 v){
return u.x*v.x+u.y*v.y+u.z*v.z;
}
//U-V
point3 subt(point3 u,point3 v){
point3 ret;
ret.x=u.x-v.x;
ret.y=u.y-v.y;
ret.z=u.z-v.z;
return ret;
}

//Get plane normal vector
point3 pvec(plane3 s){
return xmult(subt(s.a,s.b),subt(s.b,s.c));
}
point3 pvec(point3 s1,point3 s2,point3 s3){
return xmult(subt(s1,s2),subt(s2,s3));
}

//Two-point distance, orientation of single parameter size
double distance(point3 p1,point3 p2){
return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)+(p1.z-p2.z)*(p1.z-p2.z));
}

//Vector size
double vlen(point3 p){
return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}

//Sentenced to three points collinear
int dots_inline(point3 p1,point3 p2,point3 p3){
return vlen(xmult(subt(p1,p2),subt(p2,p3)))<eps;
}
//Four coplanar points
int dots_onplane(point3 a,point3 b,point3 c,point3 d){
return zero(dmult(pvec(a,b,c),subt(d,a)));
}
//Given point is on the line, including endpoint and sharing the same line
int dot_online_in(point3 p,line3 l){
return zero(vlen(xmult(subt(p,l.a),subt(p,l.b))))&&(l.a.x-p.x)*(l.b.x-p.x)<eps&&
(l.a.y-p.y)*(l.b.y-p.y)<eps&&(l.a.z-p.z)*(l.b.z-p.z)<eps;
}
int dot_online_in(point3 p,point3 l1,point3 l2){
return zero(vlen(xmult(subt(p,l1),subt(p,l2))))&&(l1.x-p.x)*(l2.x-p.x)<eps&&
(l1.y-p.y)*(l2.y-p.y)<eps&&(l1.z-p.z)*(l2.z-p.z)<eps;
}
//Liable to a point is on the segment, excluding endpoints
int dot_online_ex(point3 p,line3 l){
return dot_online_in(p,l)&&(!zero(p.x-l.a.x)||!zero(p.y-l.a.y)||!zero(p.z-l.a.z))&&
(!zero(p.x-l.b.x)||!zero(p.y-l.b.y)||!zero(p.z-l.b.z));
}
int dot_online_ex(point3 p,point3 l1,point3 l2){
return dot_online_in(p,l1,l2)&&(!zero(p.x-l1.x)||!zero(p.y-l1.y)||!zero(p.z-l1.z))&&
(!zero(p.x-l2.x)||!zero(p.y-l2.y)||!zero(p.z-l2.z));
}
//Found on a triangle point in space, including border, three points collinear is useless
int dot_inplane_in(point3 p,plane3 s){
return zero(vlen(xmult(subt(s.a,s.b),subt(s.a,s.c)))-vlen(xmult(subt(p,s.a),subt(p,s.b)))-
vlen(xmult(subt(p,s.b),subt(p,s.c)))-vlen(xmult(subt(p,s.c),subt(p,s.a))));
}
int dot_inplane_in(point3 p,point3 s1,point3 s2,point3 s3){
return zero(vlen(xmult(subt(s1,s2),subt(s1,s3)))-vlen(xmult(subt(p,s1),subt(p,s2)))-
vlen(xmult(subt(p,s2),subt(p,s3)))-vlen(xmult(subt(p,s3),subt(p,s1))));
}
//Found on a triangle point in space, not including borders, three points collinear is useless
int dot_inplane_ex(point3 p,plane3 s){
return dot_inplane_in(p,s)&&vlen(xmult(subt(p,s.a),subt(p,s.b)))>eps&&
vlen(xmult(subt(p,s.b),subt(p,s.c)))>eps&&vlen(xmult(subt(p,s.c),subt(p,s.a)))>eps;
}
int dot_inplane_ex(point3 p,point3 s1,point3 s2,point3 s3){
return dot_inplane_in(p,s1,s2,s3)&&vlen(xmult(subt(p,s1),subt(p,s2)))>eps&&
vlen(xmult(subt(p,s2),subt(p,s3)))>eps&&vlen(xmult(subt(p,s3),subt(p,s1)))>eps;
}
//Online sentenced to two points with the side, points on the line return 0, not total nonsense
int same_side(point3 p1,point3 p2,line3 l){
return dmult(xmult(subt(l.a,l.b),subt(p1,l.b)),xmult(subt(l.a,l.b),subt(p2,l.b)))>eps;
}
int same_side(point3 p1,point3 p2,point3 l1,point3 l2){
return dmult(xmult(subt(l1,l2),subt(p1,l2)),xmult(subt(l1,l2),subt(p2,l2)))>eps;
}
//ISO side two point in the segment, points on the line return 0, not total nonsense
int opposite_side(point3 p1,point3 p2,line3 l){
return dmult(xmult(subt(l.a,l.b),subt(p1,l.b)),xmult(subt(l.a,l.b),subt(p2,l.b)))<-eps;
}
int opposite_side(point3 p1,point3 p2,point3 l1,point3 l2){
return dmult(xmult(subt(l1,l2),subt(p1,l2)),xmult(subt(l1,l2),subt(p2,l2)))<-eps;
}
//Sentenced to two points in the plane with the side, points on the plane returned 0
int same_side(point3 p1,point3 p2,plane3 s){
return dmult(pvec(s),subt(p1,s.a))*dmult(pvec(s),subt(p2,s.a))>eps;
}
int same_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){
return dmult(pvec(s1,s2,s3),subt(p1,s1))*dmult(pvec(s1,s2,s3),subt(p2,s1))>eps;
}
//Different side sentenced to two points in the plane, on a point in the plane returns 0
int opposite_side(point3 p1,point3 p2,plane3 s){
return dmult(pvec(s),subt(p1,s.a))*dmult(pvec(s),subt(p2,s.a))<-eps;
}
int opposite_side(point3 p1,point3 p2,point3 s1,point3 s2,point3 s3){
return dmult(pvec(s1,s2,s3),subt(p1,s1))*dmult(pvec(s1,s2,s3),subt(p2,s1))<-eps;
}
//Two lines run parallel
int parallel(line3 u,line3 v){
return vlen(xmult(subt(u.a,u.b),subt(v.a,v.b)))<eps;
}
int parallel(point3 u1,point3 u2,point3 v1,point3 v2){
return vlen(xmult(subt(u1,u2),subt(v1,v2)))<eps;
}
//Sentenced two plane parallel
int parallel(plane3 u,plane3 v){
return vlen(xmult(pvec(u),pvec(v)))<eps;
}
int parallel(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
return vlen(xmult(pvec(u1,u2,u3),pvec(v1,v2,v3)))<eps;
}
//Sentenced to straight line and plane-parallel
int parallel(line3 l,plane3 s){
return zero(dmult(subt(l.a,l.b),pvec(s)));
}
int parallel(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
return zero(dmult(subt(l1,l2),pvec(s1,s2,s3)));
}
//Sentenced two lines vertical
int perpendicular(line3 u,line3 v){
return zero(dmult(subt(u.a,u.b),subt(v.a,v.b)));
}
int perpendicular(point3 u1,point3 u2,point3 v1,point3 v2){
	return zero(dmult(subt(u1,u2),subt(v1,v2)));
}

//Two planes vertical
int perpendicular(plane3 u,plane3 v){
return zero(dmult(pvec(u),pvec(v)));
}
int perpendicular(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
return zero(dmult(pvec(u1,u2,u3),pvec(v1,v2,v3)));
}
//Sentenced to straight line and plane-parallel
int perpendicular(line3 l,plane3 s){
return vlen(xmult(subt(l.a,l.b),pvec(s)))<eps;
}
int perpendicular(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
return vlen(xmult(subt(l1,l2),pvec(s1,s2,s3)))<eps;
}
//Two line segments intersect, including endpoint and partial coincidence
int intersect_in(line3 u,line3 v){
if (!dots_onplane(u.a,u.b,v.a,v.b))
return 0;
if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b))
return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
return dot_online_in(u.a,v)||dot_online_in(u.b,v)||dot_online_in(v.a,u)||dot_online_in(v.b,u);
}
int intersect_in(point3 u1,point3 u2,point3 v1,point3 v2){
if (!dots_onplane(u1,u2,v1,v2))
return 0;
if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
return
dot_online_in(u1,v1,v2)||dot_online_in(u2,v1,v2)||dot_online_in(v1,u1,u2)||dot_online_in(v2,u1,u
2);
}
//Two line segments intersect, not including endpoint and partial coincidence
int intersect_ex(line3 u,line3 v){
return dots_onplane(u.a,u.b,v.a,v.b)&&opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}
int intersect_ex(point3 u1,point3 u2,point3 v1,point3 v2){
return
dots_onplane(u1,u2,v1,v2)&&opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
//Sentenced segment and space intersect the triangle, including cross-border and (some) contain
int intersect_in(line3 l,plane3 s){
return !same_side(l.a,l.b,s)&&!same_side(s.a,s.b,l.a,l.b,s.c)&&
!same_side(s.b,s.c,l.a,l.b,s.a)&&!same_side(s.c,s.a,l.a,l.b,s.b);
}
int intersect_in(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
return !same_side(l1,l2,s1,s2,s3)&&!same_side(s1,s2,l1,l2,s3)&&
!same_side(s2,s3,l1,l2,s1)&&!same_side(s3,s1,l1,l2,s2);
}
//Sentenced segment and space intersect the triangle, not including cross-border and (some) contain
int intersect_ex(line3 l,plane3 s){
return opposite_side(l.a,l.b,s)&&opposite_side(s.a,s.b,l.a,l.b,s.c)&&
opposite_side(s.b,s.c,l.a,l.b,s.a)&&opposite_side(s.c,s.a,l.a,l.b,s.b);
}
int intersect_ex(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
return opposite_side(l1,l2,s1,s2,s3)&&opposite_side(s1,s2,l1,l2,s3)&&
opposite_side(s2,s3,l1,l2,s1)&&opposite_side(s3,s1,l1,l2,s2);
}
//Calculate the intersection of two lines, attention to determine whether lines are parallel and coplanar in advance!
//Intersects a line segment intersection sentenced segment (also was to determine whether parallel!)
point3 intersection(line3 u,line3 v){
point3 ret=u.a;
double t=((u.a.x-v.a.x)*(v.a.y-v.b.y)-(u.a.y-v.a.y)*(v.a.x-v.b.x))
/((u.a.x-u.b.x)*(v.a.y-v.b.y)-(u.a.y-u.b.y)*(v.a.x-v.b.x));
ret.x+=(u.b.x-u.a.x)*t;
ret.y+=(u.b.y-u.a.y)*t;
ret.z+=(u.b.z-u.a.z)*t;
return ret;
}
point3 intersection(point3 u1,point3 u2,point3 v1,point3 v2){
point3 ret=u1;
double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
ret.x+=(u2.x-u1.x)*t;
ret.y+=(u2.y-u1.y)*t;
ret.z+=(u2.z-u1.z)*t;
return ret;
}
//Calculated intersection of straight lines and planes, attention in advance to determine whether parallel and guaranteed three points are not collinear!
//Segment intersection space triangle and judging
point3 intersection(line3 l,plane3 s){
point3 ret=pvec(s);
double t=(ret.x*(s.a.x-l.a.x)+ret.y*(s.a.y-l.a.y)+ret.z*(s.a.z-l.a.z))/
(ret.x*(l.b.x-l.a.x)+ret.y*(l.b.y-l.a.y)+ret.z*(l.b.z-l.a.z));
ret.x=l.a.x+(l.b.x-l.a.x)*t;
ret.y=l.a.y+(l.b.y-l.a.y)*t;
ret.z=l.a.z+(l.b.z-l.a.z)*t;
return ret;
}
point3 intersection(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
point3 ret=pvec(s1,s2,s3);
double t=(ret.x*(s1.x-l1.x)+ret.y*(s1.y-l1.y)+ret.z*(s1.z-l1.z))/
(ret.x*(l2.x-l1.x)+ret.y*(l2.y-l1.y)+ret.z*(l2.z-l1.z));
ret.x=l1.x+(l2.x-l1.x)*t;
ret.y=l1.y+(l2.y-l1.y)*t;
ret.z=l1.z+(l2.z-l1.z)*t;
return ret;
}
//Calculating intersection, attention in advance to determine whether parallel and guaranteed three points are not collinear!
line3 intersection(plane3 u,plane3 v){
line3 ret;
ret.a=parallel(v.a,v.b,u.a,u.b,u.c)?intersection(v.b,v.c,u.a,u.b,u.c):intersection(v.a,v.b,u.a,u.b,u.
c);
ret.b=parallel(v.c,v.a,u.a,u.b,u.c)?intersection(v.b,v.c,u.a,u.b,u.c):intersection(v.c,v.a,u.a,u.b,u.
c);
return ret;
}
line3 intersection(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
line3 ret;
ret.a=parallel(v1,v2,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v1,v2,u1,u2,u3);
ret.b=parallel(v3,v1,u1,u2,u3)?intersection(v2,v3,u1,u2,u3):intersection(v3,v1,u1,u2,u3);
return ret;
}
//Point to straight line distance
double ptoline(point3 p,line3 l){
return vlen(xmult(subt(p,l.a),subt(l.b,l.a)))/distance(l.a,l.b);
}
double ptoline(point3 p,point3 l1,point3 l2){
return vlen(xmult(subt(p,l1),subt(l2,l1)))/distance(l1,l2);
}
//Points to the horizontal distance between
double ptoplane(point3 p,plane3 s){
return fabs(dmult(pvec(s),subt(p,s.a)))/vlen(pvec(s));
}
double ptoplane(point3 p,point3 s1,point3 s2,point3 s3){
return fabs(dmult(pvec(s1,s2,s3),subt(p,s1)))/vlen(pvec(s1,s2,s3));
}
//Linear to linear distance
double linetoline(line3 u,line3 v){
point3 n=xmult(subt(u.a,u.b),subt(v.a,v.b));
return fabs(dmult(subt(u.a,v.a),n))/vlen(n);
}
double linetoline(point3 u1,point3 u2,point3 v1,point3 v2){
point3 n=xmult(subt(u1,u2),subt(v1,v2));
return fabs(dmult(subt(u1,v1),n))/vlen(n);
}
//Angle between two linear COS values
double angle_cos(line3 u,line3 v){
return dmult(subt(u.a,u.b),subt(v.a,v.b))/vlen(subt(u.a,u.b))/vlen(subt(v.a,v.b));
}
double angle_cos(point3 u1,point3 u2,point3 v1,point3 v2){
return dmult(subt(u1,u2),subt(v1,v2))/vlen(subt(u1,u2))/vlen(subt(v1,v2));
}
//Angle between two planes COS values
double angle_cos(plane3 u,plane3 v){
return dmult(pvec(u),pvec(v))/vlen(pvec(u))/vlen(pvec(v));
}
double angle_cos(point3 u1,point3 u2,point3 u3,point3 v1,point3 v2,point3 v3){
return dmult(pvec(u1,u2,u3),pvec(v1,v2,v3))/vlen(pvec(u1,u2,u3))/vlen(pvec(v1,v2,v3));
}
//Linear plane angle Sin values
double angle_sin(line3 l,plane3 s){
return dmult(subt(l.a,l.b),pvec(s))/vlen(subt(l.a,l.b))/vlen(pvec(s));
}
double angle_sin(point3 l1,point3 l2,point3 s1,point3 s2,point3 s3){
return dmult(subt(l1,l2),pvec(s1,s2,s3))/vlen(subt(l1,l2))/vlen(pvec(s1,s2,s3));
}




//Grid
#define abs(x) ((x)>0?(x):-(x))
struct point{int x,y;};
int gcd(int a,int b){
return b?gcd(b,a%b):a;
}
//Number of grid points on the polygon
int grid_onedge(int n,point* p){
int i,ret=0;
for (i=0;i<n;i++)
ret+=gcd(abs(p[i].x-p[(i+1)%n].x),abs(p[i].y-p[(i+1)%n].y));
return ret;
}
//Grid points within the polygon number
int grid_inside(int n,point* p){
int i,ret=0;
for (i=0;i<n;i++)
ret+=p[(i+1)%n].y*(p[i].x-p[(i+2)%n].x);
return (abs(ret)-grid_onedge(n,p))/2+1;
}





//Circle
#include <math.h>
#define eps 1e-8
struct point{double x,y;};
double xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
double distance(point p1,point p2){
return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y));
}
double disptoline(point p,point l1,point l2){
return fabs(xmult(p,l1,l2))/distance(l1,l2);
}
point intersection(point u1,point u2,point v1,point v2){
point ret=u1;
double t=((u1.x-v1.x)*(v1.y-v2.y)-(u1.y-v1.y)*(v1.x-v2.x))
/((u1.x-u2.x)*(v1.y-v2.y)-(u1.y-u2.y)*(v1.x-v2.x));
ret.x+=(u2.x-u1.x)*t;
ret.y+=(u2.y-u1.y)*t;
return ret;
}
//Sentenced to intersect the lines and circles, including the tangent
int intersect_line_circle(point c,double r,point l1,point l2){
return disptoline(c,l1,l2)<r+eps;
}
//Sentenced to intersect the lines and circles, including the endpoints and the tangent
int intersect_seg_circle(point c,double r,point l1,point l2){
double t1=distance(c,l1)-r,t2=distance(c,l2)-r;
point t=c;
if (t1<eps||t2<eps)
return t1>-eps||t2>-eps;
t.x+=l1.y-l2.y;
t.y+=l2.x-l1.x;
return xmult(l1,c,t)*xmult(l2,c,t)<eps&&disptoline(c,l1,l2)-r<eps;
}
//Sentenced to intersect the circle and circle, including the tangent
int intersect_circle_circle(point c1,double r1,point c2,double r2){
return distance(c1,c2)<r1+r2+eps&&distance(c1,c2)>fabs(r1-r2)-eps;
}
//Calculated point p on the circle recently, such as p coincides with the Center, return p itself
point dot_to_circle(point c,double r,point p){
point u,v;
if (distance(p,c)<eps)
return p;
u.x=c.x+r*fabs(c.x-p.x)/distance(c,p);
u.y=c.y+r*fabs(c.y-p.y)/distance(c,p)*((c.x-p.x)*(c.y-p.y)<0?-1:1);
v.x=c.x-r*fabs(c.x-p.x)/distance(c,p);
v.y=c.y-r*fabs(c.y-p.y)/distance(c,p)*((c.x-p.x)*(c.y-p.y)<0?-1:1);
return distance(u,p)<distance(v,p)?u:v;
}
//Calculated intersection of lines and circles to ensure straight lines and circular intersection
//Calculate the intersection of line segments and circular sentenced after the functions available on the segment
void intersection_line_circle(point c,double r,point l1,point l2,point& p1,point& p2){
point p=c;
double t;
p.x+=l1.y-l2.y;
p.y+=l2.x-l1.x;
p=intersection(p,c,l1,l2);
t=sqrt(r*r-distance(p,c)*distance(p,c))/distance(l1,l2);
p1.x=p.x+(l2.x-l1.x)*t;
p1.y=p.y+(l2.y-l1.y)*t;
p2.x=p.x-(l2.x-l1.x)*t;
p2.y=p.y-(l2.y-l1.y)*t;
}
//Calculate the intersection of the circle and circle, ensure the intersection of circle and circle, circle do not coincide
void intersection_circle_circle(point c1,double r1,point c2,double r2,point& p1,point& p2){
point u,v;
double t;
t=(1+(r1*r1-r2*r2)/distance(c1,c2)/distance(c1,c2))/2;
u.x=c1.x+(c2.x-c1.x)*t;
u.y=c1.y+(c2.y-c1.y)*t;
v.x=u.x+c1.y-c2.y;
v.y=u.y-c1.x+c2.x;
intersection_line_circle(c1,r1,u,v,p1,p2);
}





//Integer function
//Integer geometry functions
//Note some integer arithmetic is out!
#define sign(a) ((a)>0?1:(((a)<0?-1:0)))
struct point{int x,y;};
struct line{point a,b;};
//cross
int xmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.y-p0.y)-(p2.x-p0.x)*(p1.y-p0.y);
}
int xmult(int x1,int y1,int x2,int y2,int x0,int y0){
return (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0);
}
//dot
int dmult(point p1,point p2,point p0){
return (p1.x-p0.x)*(p2.x-p0.x)+(p1.y-p0.y)*(p2.y-p0.y);
}
int dmult(int x1,int y1,int x2,int y2,int x0,int y0){
return (x1-x0)*(x2-x0)+(y1-y0)*(y2-y0);
}
//Sentenced to three points collinear
int dots_inline(point p1,point p2,point p3){
return !xmult(p1,p2,p3);
}
int dots_inline(int x1,int y1,int x2,int y2,int x3,int y3){
return !xmult(x1,y1,x2,y2,x3,y3);
}
//Given point is on the line, including endpoint and partial coincidence
int dot_online_in(point p,line l){
return !xmult(p,l.a,l.b)&&(l.a.x-p.x)*(l.b.x-p.x)<=0&&(l.a.y-p.y)*(l.b.y-p.y)<=0;
}
int dot_online_in(point p,point l1,point l2){
return !xmult(p,l1,l2)&&(l1.x-p.x)*(l2.x-p.x)<=0&&(l1.y-p.y)*(l2.y-p.y)<=0;
}
int dot_online_in(int x,int y,int x1,int y1,int x2,int y2){
return !xmult(x,y,x1,y1,x2,y2)&&(x1-x)*(x2-x)<=0&&(y1-y)*(y2-y)<=0;
}
//Liable to a point is on the segment, excluding endpoints
int dot_online_ex(point p,line l){
return dot_online_in(p,l)&&(p.x!=l.a.x||p.y!=l.a.y)&&(p.x!=l.b.x||p.y!=l.b.y);
}
int dot_online_ex(point p,point l1,point l2){
return dot_online_in(p,l1,l2)&&(p.x!=l1.x||p.y!=l1.y)&&(p.x!=l2.x||p.y!=l2.y);
}
int dot_online_ex(int x,int y,int x1,int y1,int x2,int y2){
return dot_online_in(x,y,x1,y1,x2,y2)&&(x!=x1||y!=y1)&&(x!=x2||y!=y2);
}
//Sentenced to two points in a straight line with the side, points on the line return 0
int same_side(point p1,point p2,line l){
return sign(xmult(l.a,p1,l.b))*xmult(l.a,p2,l.b)>0;
}
int same_side(point p1,point p2,point l1,point l2){
return sign(xmult(l1,p1,l2))*xmult(l1,p2,l2)>0;
}
//Abnormal side of two points in a straight line, on the line return 0
int opposite_side(point p1,point p2,line l){
return sign(xmult(l.a,p1,l.b))*xmult(l.a,p2,l.b)<0;
}
int opposite_side(point p1,point p2,point l1,point l2){
return sign(xmult(l1,p1,l2))*xmult(l1,p2,l2)<0;
}
//Two lines run parallel
int parallel(line u,line v){
return (u.a.x-u.b.x)*(v.a.y-v.b.y)==(v.a.x-v.b.x)*(u.a.y-u.b.y);
}
int parallel(point u1,point u2,point v1,point v2){
return (u1.x-u2.x)*(v1.y-v2.y)==(v1.x-v2.x)*(u1.y-u2.y);
}
//Sentenced two lines vertical
int perpendicular(line u,line v){
return (u.a.x-u.b.x)*(v.a.x-v.b.x)==-(u.a.y-u.b.y)*(v.a.y-v.b.y);
}
int perpendicular(point u1,point u2,point v1,point v2){
return (u1.x-u2.x)*(v1.x-v2.x)==-(u1.y-u2.y)*(v1.y-v2.y);
}
//Two line segments intersect, including endpoint and partial coincidence
int intersect_in(line u,line v){
if (!dots_inline(u.a,u.b,v.a)||!dots_inline(u.a,u.b,v.b))
return !same_side(u.a,u.b,v)&&!same_side(v.a,v.b,u);
return dot_online_in(u.a,v)||dot_online_in(u.b,v)||dot_online_in(v.a,u)||dot_online_in(v.b,u);
}
int intersect_in(point u1,point u2,point v1,point v2){
if (!dots_inline(u1,u2,v1)||!dots_inline(u1,u2,v2))
return !same_side(u1,u2,v1,v2)&&!same_side(v1,v2,u1,u2);
return
dot_online_in(u1,v1,v2)||dot_online_in(u2,v1,v2)||dot_online_in(v1,u1,u2)||dot_online_in(v2,u1,u
2);
}
//Two line segments intersect, not including endpoint and partial coincidence
int intersect_ex(line u,line v){
return opposite_side(u.a,u.b,v)&&opposite_side(v.a,v.b,u);
}
int intersect_ex(point u1,point u2,point v1,point v2){
return opposite_side(u1,u2,v1,v2)&&opposite_side(v1,v2,u1,u2);
}
