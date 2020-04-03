/*
Liuctic computational geometry libraries
p	-Lpoint ln,
l	- Lline 
ls	- Llineseg 
lr	- Lrad
Distance between two points on the plane					p2pdis
Return (P1-P0) * (P2-P0) the cross product.					xmulti
Determines whether two line segments intersecting			lsinterls
Tell if a point p is on the line segment l					ponls
Determine whether two points are equal to					Euqal_Point
Line segment intersects non-endpoint						lsinterls_A
Determines whether the point q in the polygon in Polygon	pinplg
Polygon area												area_of_polygon
Solve quadratic equation Ax^2+Bx+C=0						equa
Distance from points to lines								p2lndis
Intersection of lines and circles, known line will intersect the circle			lncrossc
Points directly to the ray														samedir
Ray intersection and rounded the first											lrcrossc
Seeking symmetry point of the point P1 on the line LN P2						mirror
Angle between two lines (ARC)													angle_LL
*/

#define infinity 1e20
#define EP 1e-10
const int MAXV = 300 ;
const double PI = 2.0*asin(1.0); //�߾�����PI
struct Lpoint {double x,y;}; //��
struct Llineseg{ Lpoint a,b;}; //�߶�
struct Ldir{double dx,dy;}; //��������
struct Lline{Lpoint p; Ldir dir;}; //ֱ��
struct Lrad{Lpoint Sp; Ldir dir;}; //����
struct Lround{Lpoint co; double r;};//Բ
| ��ƽ��������֮��ľ���
double p2pdis(Lpoint p1,Lpoint p2) {
return (sqrt((p1.x-p2.x) * (p1.x-p2.x) +
(p1.y-p2.y) * (p1.y-p2.y)));
}
double xmulti(Lpoint p1,Lpoint p2,Lpoint p0) {
return((p1.x-p0.x) * (p2.y-p0.y) -
(p2.x-p0.x) * (p1.y-p0.y));
}
(p2.x-p0.x) * (p1.y-p0.y));
double mx(double t1,double t2)
{
if(t1>t2) return t1;
return t2;
}
double mn(double t1,double t2)
{
if(t1<t2) return t1;
return t2;
}
int lsinterls(Llineseg u,Llineseg v)
{
return( (mx(u.a.x,u.b.x)>=mn(v.a.x,v.b.x))&&
(mx(v.a.x,v.b.x)>=mn(u.a.x,u.b.x))&&
(mx(u.a.y,u.b.y)>=mn(v.a.y,v.b.y))&&
(mx(v.a.y,v.b.y)>=mn(u.a.y,u.b.y))&&
(xmulti(v.a,u.b,u.a)*xmulti(u.b,v.b,u.a)>=0)&&
(xmulti(u.a,v.b,v.a)*xmulti(v.b,u.b,v.a)>=0));
}
| �жϵ�p�Ƿ����߶�l��
int ponls(Llineseg l,Lpoint p) {
return( (xmulti(l.b,p,l.a)==0) &&
( ((p.x-l.a.x)*(p.x-l.b.x)<0 ) ||
((p.y-l.a.y)*(p.y-l.b.y)<0 )) );
}
| �ж��������Ƿ����
int Euqal_Point(Lpoint p1,Lpoint p2) {
return((fabs(p1.x-p2.x)<EP)&&(fabs(p1.y-p2.y)<EP));
}
| �߶��ཻ�жϺ���
���ҽ���u,v�ཻ���ҽ��㲻��u,v�Ķ˵�ʱ����Ϊtrue;
int lsinterls_A(Llineseg u,Llineseg v) {
return((lsinterls(u,v)) && (!Euqal_Point(u.a,v.a))&&
(!Euqal_Point(u.a,v.b)) && (!Euqal_Point(u.b,v.a))&&
(!Euqal_Point(u.b,v.b)));
}
/*===============================================
| �жϵ�q�Ƿ��ڶ������
���ж�����������͹�򰼶���Σ�
Polygon�д�Ŷ���ε���ʱ�붥������
================================================*/
int pinplg(int vcount,Lpoint Polygon[],Lpoint q)
{
int c=0,i,n;
Llineseg l1,l2;
l1.a=q; l1.b=q; l1.b.x=infinity; n=vcount;
for (i=0;i<vcount;i++) {
l2.a=Polygon[i];
l2.b=Polygon[(i+1)%n];
if ((lsinterls_A(l1,l2))||
(
(ponls(l1,Polygon[(i+1)%n]))&&
(
(!ponls(l1,Polygon[(i+2)%n]))&&
(xmulti(Polygon[i],Polygon[(i+1)%n],l1.a) *
xmulti(Polygon[(i+1)%n],Polygon[(i+2)%n],l1.a)>0)
||
(ponls(l1,Polygon[(i+2)%n]))&&
(xmulti(Polygon[i],Polygon[(i+2)%n],l1.a) *
xmulti(Polygon[(i+2)%n],Polygon[(i+3)%n],l1.a)>0)
) ) ) c++;
}
return(c%2!=0);
}
/*==================================================*\
| �������ε����
| Ҫ������ʱ�뷽���������ζ���
| ������͹����λ򰼶����
\*==================================================*/
double areaofp(int vcount,double x[],double y[],Lpoint plg[])
{
int i;
double s;
if (vcount<3) return 0;
s=plg[0].y*(plg[vcount-1].x-plg[1].x);
for (i=1;i<vcount;i++)
s+=plg[i].y*(plg[(i-1)].x-plg[(i+1)%vcount].x);
return s/2;
}
/*********************\
| ����η��� Ax^2+Bx+C=0
����-1��ʾ�޽� ����1 ��ʾ�н�
\*********************/
int equa(double A,double B,double C,double& x1,double& x2)
{
double f=B*B-4*A*C;
if(f<0) return -1;
x1=(-B+sqrt(f))/(2*A);
x2=(-B-sqrt(f))/(2*A);
return 1;
}
| ����ֱ�ߵ�һ��ʽ Ax+By+C=0
void format(Lline ln,double& A,double& B,double& C)
{
A=ln.dir.dy;
B=-ln.dir.dx;
C=ln.p.y*ln.dir.dx-ln.p.x*ln.dir.dy;
}
| �㵽ֱ�߾���
double p2ldis(Lpoint a,Lline ln)
{
double A,B,C;
format(ln,A,B,C);
return(fabs(A*a.x+B*a.y+C)/sqrt(A*A+B*B));
}
| ֱ����Բ�Ľ��㣬��ֱ֪����Բ�ཻ
int lncrossc(Lline ln,Lround Y,Lpoint& p1,Lpoint& p2)
{
double A,B,C,t1,t2;
int zz=-1;
format(ln,A,B,C);
if(fabs(B)<1e-8)
{
p1.x=p2.x=-1.0*C/A;
zz=equa(1.0,-2.0*Y.co.y,Y.co.y*Y.co.y
+(p1.x-Y.co.x)*(p1.x-Y.co.x)-Y.r*Y.r,t1,t2);
p1.y=t1;p2.y=t2;
}
else if(fabs(A)<1e-8)
{
p1.y=p2.y=-1.0*C/B;
zz=equa(1.0,-2.0*Y.co.x,Y.co.x*Y.co.x
+(p1.y-Y.co.y)*(p1.y-Y.co.y)-Y.r*Y.r,t1,t2);
p1.x=t1;p2.x=t2;
}
else
{
zz=equa(A*A+B*B,2.0*A*C+2.0*A*B*Y.co.y
-2.0*B*B*Y.co.x,B*B*Y.co.x*Y.co.x+C*C+2*B*C*Y.co.y
+B*B*Y.co.y*Y.co.y-B*B*Y.r*Y.r,t1,t2);
p1.x=t1,p1.y=-1*(A/B*t1+C/B);
p2.x=t2,p2.y=-1*(A/B*t2+C/B);
}
return 0;
}
| ���Ƿ������ߵ�����
bool samedir(Lrad ln,Lpoint P)
{
double ddx,ddy;
ddx=P.x-ln.Sp.x;ddy=P.y-ln.Sp.y;
if((ddx*ln.dir.dx>0||fabs(ddx*ln.dir.dx)<1e-7)
&&(ddy*ln.dir.dy>0||(fabs(ddy*ln.dir.dy)<1e-7)))
return true;
else return false;
}
| ������Բ�ĵ�һ������
�Ѿ�ȷ����������ֱ����Բ�ཻ����-1��ʾ�������򽻵� �����򷵻�1
int lrcrossc(Lrad ln, Lround Y, Lpoint& P)
{
Lline ln2;
Lpoint p1,p2;
int res=-1;
double dis=1e20;
ln2.p=ln.Sp,ln2.dir=ln.dir;
lncrossc(ln2,Y,p1,p2);
if(samedir(ln,p1))
{
res=1;
if(p2pdis(p1,ln.Sp)<dis)
{
dis=p2pdis(p1,ln.Sp);
P=p1;
}
}
if(samedir(ln,p2))
{
res=1;
if(p2pdis(p2,ln.Sp)<dis)
{
dis=p2pdis(p2,ln.Sp);
P=p2;
}
}
return res;
}
| ���p1����ֱ��ln�ĶԳƵ�p2
Lpoint mirror(Lpoint P,Lline ln)
{
Lpoint Q;
double A,B,C;
format(ln,A,B,C);
Q.x=((B*B-A*A)*P.x-2*A*B*P.y-2*A*C)/(A*A+B*B);
Q.y=((A*A-B*B)*P.y-2*A*B*P.x-2*B*C)/(A*A+B*B);
return Q;
}
| ��ֱ�߼нǣ����ȣ�
double angle_LL(Lline line1, Lline line2)
{
double A1, B1, C1;
format(line1, A1, B1, C1);
double A2, B2, C2;
format(line2, A2, B2, C2);
if( A1*A2+B1*B2 == 0 ) return PI/2.0; // ��ֱ
else{
double t = fabs((A1*B2-A2*B1)/(A1*A2+B1*B2));
return atan(t);
}
}