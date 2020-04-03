typedef pair<double, double> PDD;

PDD operator+(PDD A, PDD B) { return PDD(A.first+B.first, A.second+B.second); }
PDD operator-(PDD A, PDD B) { return PDD(A.first-B.first, A.second-B.second); }
PDD operator*(double A, PDD B) { return PDD(A*B.first, A*B.second); }
PDD operator*(PDD B, double A) { return PDD(A*B.first, A*B.second); }
PDD operator/(PDD B, double A) { return PDD(B.first/A, B.second/A); }

//I = A + Bt, I = C + Ds
int line_intersect_vector(PDD A, PDD B, PDD C, PDD D, PDD &I)
{
	double hor = B.first*-D.second + D.first*B.second;
	double lob = -D.first*(A.second-C.second) + (A.first-C.first)*D.second;
	if(fabs(hor)<1e-7) return 0;
//	assert(!EQ(hor,0));
	double t = lob/hor;

	I = PDD(A.first+t*B.first,A.second+t*B.second);
	return 1;
}

//L1(A,B) X L2(C,D) = I
int line_intersect(PDD A, PDD B, PDD C, PDD D, PDD &I)
{
	return line_intersect_vector(A, PDD(B.first-A.first,B.second-A.second), C, PDD(D.first-C.first,D.second-C.second), I);
}

double dist(PDD A, PDD B)
{
	return sqrt( (A.first - B.first)*(A.first - B.first) + (A.second - B.second)*(A.second - B.second));
}

double sdist(PDD A, PDD B)
{
	return ( (A.first - B.first)*(A.first - B.first) + (A.second - B.second)*(A.second - B.second));
}
