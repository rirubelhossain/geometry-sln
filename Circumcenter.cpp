void Circumcenter (CPoint p0 , CPoint p1 , CPoint p2 , CPoint &cp)
{
double a1=p1.x-p0.x,b1 = p1.y - p0.y,c1 = (sqr(a1) + sqr(b1)) / 2;
double a2=p2.x-p0.x,b2 = p2.y - p0.y,c2 = (sqr(a2) + sqr(b2)) / 2;
double d = a1 * b2 - a2 * b1;
cp.x = p0.x + (c1 * b2 - c2 * b1) / d;
cp.y = p0.y + (a1 * c2 - a2 * c1) / d;
}

