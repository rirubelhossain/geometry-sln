typedef long long real;
struct point { real x,y; struct edge *e; };
struct edge {
point *o, *d;
struct edge *on, *op, *dn, *dp;
};
typedef struct point point;
typedef struct edge edge;
#define Op(e,p) ((e)->o==p?(e)->d:(e)->o)
#define Next(e,p) ((e)->o==p?(e)->on:(e)->dn)
#define Prev(e,p) ((e)->o==p?(e)->op:(e)->dp)
edge *make_edge(point *u, point *v){
edge *e = new edge();
e->on=e->op=e->dn=e->dp=e;
e->o=u;e->d=v;
if (u->e==NULL) u->e=e;
if (v->e==NULL) v->e=e;
return e;
}
void delete_edge(edge *e){
point *u=e->o, *v=e->d;
if (u->e==e) u->e=e->on;
if (v->e==e) v->e=e->dn;
if (e->on->o==u) e->on->op=e->op; else e->on->dp=e->op;
if (e->op->o==u) e->op->on=e->on; else e->op->dn=e->on;
if (e->dn->o==v) e->dn->op=e->dp; else e->dn->dp=e->dp;
if (e->dp->o==v) e->dp->on=e->dn; else e->dp->dn=e->dn;
delete e;
}
void splice(edge *a, edge *b, point *v){
edge *n;
if (a->o==v) { n=a->on; a->on=b; } else { n=a->dn; a->dn=b; }
if (n->o==v) n->op=b; else n->dp=b;
if (b->o==v) { b->on=n; b->op=a; } else { b->dn=n; b->dp=a; }
}
edge *join(edge *a, point *u, edge *b, point *v, int s){
edge *e = make_edge(u,v);
if (s == 0){ if (a->o==u) splice(a->op,e,u); else splice(a->dp,e,u); splice(b,e,v); }
else { splice(a,e,u); if (b->o==v) splice(b->op,e,v); else splice(b->dp,e,v); }
return e;
}
#define Vector(p1,p2,u,v) (u=p2->x-p1->x,v=p2->y-p1->y)
#define cross_v(u1,v1,u2,v2) (u1*v2-v1*u2)
#define cross_p(p1,p2,p3) ((p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x))
#define dot_v(u1,v1,u2,v2) (u1*u2+v1*v2)
const double eps = 1e-10;
void lower_tangent(edge *l, point *s, edge *r, point *u, edge **ll, point **llo, edge **rl, point **rlo){
point *ol=s, *dl=Op(l,s);
point *oor=u, *dr=Op(r,u);
while(1)
if (cross_p((*ol),(*dl),(*oor))>0) { l=Prev(l,dl); ol=dl; dl=Op(l,ol); }
else if (cross_p((*oor),(*dr),(*ol))<0) { r=Next(r, dr); oor=dr; dr=Op(r,oor); }
else break;
*ll=l; *rl=r; *llo=ol; *rlo=oor;
}
static void merge(edge *r_cw_l, point *s, edge *l_ccw_r, point *u, edge **l_tangent){
edge *b, *lc, *rc, *ll, *rl, *next, *prev;
point *ob, *db, *dlc, *drc, *orl, *oll, *dest_next, *dest_prev;
real ulcob, vlcob, ulcdb, vlcdb, urcob, vrcob, urcdb, vrcdb, cplc, cprc, dplc, dprc;
real unob, vnob, undb, vndb, cpn, dpn, upob, vpob, updb, vpdb, cpp, dpp;
real alc, arc, an, ap;
double crc, clc;
lower_tangent(r_cw_l, s, l_ccw_r, u, &ll, &oll, &rl, &orl);
b = join(ll, oll, rl, orl, 1);
ob = oll; db = orl;
*l_tangent = b;
do{
lc=Next(b,ob); rc=Prev(b,db); dlc=Op(lc,ob); drc=Op(rc,db);
Vector(dlc,ob,ulcob,vlcob); Vector(dlc,db,ulcdb,vlcdb); Vector(drc,ob,urcob,vrcob); Vector(drc,db,urcdb,vrcdb);
cplc=cross_v(ulcob, vlcob, ulcdb, vlcdb); cprc=cross_v(urcob, vrcob, urcdb, vrcdb);
alc=(cplc > 0); arc=(cprc > 0);
if (!alc && !arc) break;
if (alc){
dplc = dot_v(ulcob, vlcob, ulcdb, vlcdb);
clc = (double)dplc / cplc;
do{
next = Next(lc, ob);
dest_next = Op(next, ob);
Vector(dest_next, ob, unob, vnob);
Vector(dest_next, db, undb, vndb);
cpn = cross_v(unob, vnob, undb, vndb);
an = (cpn > 0);
if (!an) break;
dpn = dot_v(unob, vnob, undb, vndb);
double cn = (double)dpn / cpn;
if (cn > clc) break;
delete_edge(lc);
lc = next;
clc = cn;
}while(1);
}
if (arc){
dprc = dot_v(urcob, vrcob, urcdb, vrcdb);
crc = (double)dprc / cprc;
do{
prev = Prev(rc, db);
dest_prev = Op(prev, db);
Vector(dest_prev, ob, upob, vpob);
Vector(dest_prev, db, updb, vpdb);
cpp = cross_v(upob, vpob, updb, vpdb);
ap = (cpp > 0);
if (!ap) break;
dpp = dot_v(upob, vpob, updb, vpdb);
double cp = (double)dpp / cpp;
if (cp > crc) break;
delete_edge(rc);
rc = prev;
crc = cp;
} while (1);
}
dlc=Op(lc, ob); drc=Op(rc, db);
if (!alc || (alc && arc && crc < clc)){ b = join(b, ob, rc, drc, 1); db = drc; }
else { b = join(lc, dlc, b, db, 1); ob = dlc; }
}while(1);
}
void divide(point *p, int l, int r, edge **l_ccw, edge **r_cw){
int n=r-l+1;
edge *l_ccw_l, *r_cw_l, *l_ccw_r, *r_cw_r, *l_tangent, *a, *b, *c;
if (n == 2) {
*l_ccw = *r_cw = make_edge(&p[l], &p[r]);
}else if (n == 3) {
a = make_edge(&p[l], &p[l+1]);
b = make_edge(&p[l+1], &p[r]);
splice(a,b,&p[l+1]);
real c_p = cross_p(p[l], p[l+1], p[r]);
if (c_p>0){ c=join(a,&p[l],b,&p[r],1); *l_ccw=a; *r_cw=b; }
else if (c_p<0) { c=join(a,&p[l],b,&p[r],0); *l_ccw=c; *r_cw=c; }
else { *l_ccw=a; *r_cw=b; }
}else if (n > 3){
int split=(l+r)/2;
divide(p,l,split,&l_ccw_l,&r_cw_l);
divide(p,split+1,r,&l_ccw_r,&r_cw_r);
merge(r_cw_l,&p[split],l_ccw_r,&p[split+1],&l_tangent);
if(l_tangent->o == &p[l]) l_ccw_l=l_tangent;
if(l_tangent->d == &p[r]) r_cw_r=l_tangent;
*l_ccw=l_ccw_l; *r_cw=r_cw_r;
}
}
int n, m;
point *p;
int*eu, *ev, *pr, *r;
double *ew;
double dis(point p1, point p2){ return sqrt((p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y)); }
void enum_edges(int n){
edge *e_start, *e;
point *u, *v;
eu = new int[n*3]; ev = new int[n*3]; ew = new double[n*3];
m = 0;
for (int i = 0; i < n; i++) {
u = &p[i];
e_start = e = u->e;
do{
v = Op(e, u);
if (u < v){ eu[m] = u-p; ev[m] = v-p; ew[m]=dis(*u, *v); m++; }
e = Next(e, u);
} while (e!=e_start);
}
}
int cmp(const point& a, const point& b){ return (a.x < b.x || a.x == b.x && a.y < b.y); }
int cmp_e(const point& a, const point& b){ return a.x == b.x && a.y == b.y; }
int cmp2(const int& i, const int& j){ return ew[i]<ew[j]; }
int find(int x){ return x==pr[x]?x:pr[x]=find(pr[x]); }
double kruskal(){
double ans=0.0;
int i, j=0, e, u, v;
pr = new int[n]; r = new int[m];
for(i=0;i<n;i++) pr[i]=i;
for(i=0;i<m;i++) r[i]=i;
sort(r,r+m,cmp2);
for(i=0; i<n-1; i++){
	do{ e=r[j++]; u=find(eu[e]); v=find(ev[e]); }while(u==v);
ans+=ew[e]; pr[u]=v;
}
return ans;
}
int main(){
edge* l_cw, *r_ccw;
int i;
scanf("%d", &n);
p = new point[n];
for (i = 0; i < n; i++){
scanf("%I64d%I64d", &p[i].x, &p[i].y);
p[i].e = NULL;
}
sort(p, p+n, cmp);
n=unique(p, p+n, cmp_e)-p;
if(n == 1)
printf("%.4lf\n", 0.0);
else{
divide(p, 0, n-1, &l_cw, &r_ccw);
enum_edges(n);
printf("%.4lf\n", kruskal());
}
return 0;
}