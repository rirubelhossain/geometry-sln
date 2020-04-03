/////////////////////tested in TC
int isinside(vector<PDD> &V, PDD P)
{
	double area1 = 0, area2=0;
	int i,j;

	for(i=0;i<SZ(V);i++)
	{
		j = (i+1)%SZ(V);
		area1 += fabs(P.first*V[j].second + V[j].first*V[i].second + V[i].first*P.second
			  - P.second*V[j].first - V[j].second*V[i].first - V[i].second*P.first);
	}

	area1 = fabs(area1);

	for(i=0;i<SZ(V);i++)
	{
		j = (i+1)%SZ(V);

		area2 += V[i].first*V[j].second - V[j].first*V[i].second;
	}

	area2 = fabs(area2);

	if( fabs(area1-area2) < 1e-7 ) return 1;
	return 0;
}
////////////////////////////////////

struct pt {
	int x, y;
};

struct ang {
	int a, b;
};

bool operator < (const ang & p, const ang & q) {
	if (p.b == 0 && q.b == 0)
		return p.a < q.a;
	return p.a * 1ll * q.b < p.b * 1ll * q.a;
}

long long sq (pt & a, pt & b, pt & c) {
	return a.x*1ll*(b.y-c.y) + b.x*1ll*(c.y-a.y) + c.x*1ll*(a.y-b.y);
}

int main() {

	int n;
	cin >> n;
	vector<pt> p (n);
	int zero_id = 0;
	for (int i=0; i<n; ++i) {
		scanf ("%d%d", &p[i].x, &p[i].y);
		if (p[i].x < p[zero_id].x || p[i].x == p[zero_id].x && p[i].y < p[zero_id].y)
			zero_id = i;
	}
	pt zero = p[zero_id];
	rotate (p.begin(), p.begin()+zero_id, p.end());
	p.erase (p.begin());
	--n;

	vector<ang> a (n);
	for (int i=0; i<n; ++i) {
		a[i].a = p[i].y - zero.y;
		a[i].b = p[i].x - zero.x;
		if (a[i].a == 0)
			a[i].b = a[i].b < 0 ? -1 : 1;
	}

	for (;;) {
		pt q; // ????????? ??????
		bool in = false;
		if (q.x >= zero.x)
			if (q.x == zero.x && q.y == zero.y)
				in = true;
			else {
				ang my = { q.y-zero.y, q.x-zero.x };
				if (my.a == 0)
					my.b = my.b < 0 ? -1 : 1;
				vector<ang>::iterator it = upper_bound (a.begin(), a.end(), my);
				if (it == a.end() && my.a == a[n-1].a && my.b == a[n-1].b)
					it = a.end()-1;
				if (it != a.end() && it != a.begin()) {
					int p1 = int (it - a.begin());
					if (sq (p[p1], p[p1-1], q) <= 0)
						in = true;
				}
			}
		puts (in ? "INSIDE" : "OUTSIDE");
	}

}