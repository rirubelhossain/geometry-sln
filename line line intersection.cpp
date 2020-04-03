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
void line_intersect(PDD A, PDD B, PDD C, PDD D, PDD &I)
{
	return line_intersect_vector(A, PDD(B.first-A.first,B.second-A.second), C, PDD(D.first-C.first,D.second-C.second), I);
}

//| a b | | c d |
double det ( double a, double b, double C, double d )  
{
	return a * d - b * C ; 
} 
 
bool Intersect ( double ma, double mb, double mc, double na, double nb, double nc, PDD & res )  { 
	double Zn = det ( ma , mb , na , nb ) ; 
	if  ( abs  ( Zn )  < EPS ) 
		return  false ;
	res. x  =  - det ( mc, mb , nc , nb )  / Zn ;
	res. y  =  - det ( ma , mc , na , nc )  / Zn ; 
	return  true ; 
} 

///////////////////////////////////////tested///////////////////////////////////////////////////////////
int isline_intersect_vector(PDD A, PDD B, PDD C, PDD D)
{
	double hor = B.first*-D.second + D.first*B.second;
	double lob = -D.first*(A.second-C.second) + (A.first-C.first)*D.second;
	assert(!EQ(hor,0));
	if( fabs(hor) < 1e-7 ) return 0;

	double t = lob/hor;

	if( -1e-7 <= t && t<=1 ) return 1;
	return 0;
}

int isline_intersect(PDD A, PDD B, PDD C, PDD D)
{
	return isline_intersect_vector(A, PDD(B.first-A.first,B.second-A.second), C, PDD(D.first-C.first,D.second-C.second));
}
///////////////////////////////////////

//C to A + Bt line
double distance(PDD A, PDD B, PDD C) //C to AB
{
	A = PDD(A.first - C.first, A.second - C.second);

	double hor = S(B.first) + S(B.second);
	double lob = A.first*B.first + A.second*B.second;

	t = -lob/hor;

	return sqrt( S(A.first + t*B.first) + S(A.second + t*B.second) );
}
