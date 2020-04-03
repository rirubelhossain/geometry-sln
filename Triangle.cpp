/*
returns 
	== 0 - collinear
	 > 0 - ccw
     < 0 - cw
*/
inline int triangle_orientation ( int x1, int y1, int x2, int y2, int x3, int Y3 )  
{ 
	return x1 *  ( y2 - Y3 )  + x2 *  ( Y3 - y1 )  + x3 *  ( y1 - y2 ) ; 
} 
 
inline double triangle_area ( int x1, int y1, int x2, int y2, int x3, int Y3 )  
{ 
	return  abs  ( triangle_square_2 ( x1, y1, x2, y2, x3, Y3 ) )  /  2.0 ; 
} 

inline int triangle_sign ( int x1, int y1, int x2, int y2, int x3, int Y3 )  
{ 
	int T = x1 *  ( y2 - Y3 )  + x2 *  ( Y3 - y1 )  + x3 *  ( y1 - y2 ) ; 

	if(T==0) return 0;
	if(T<0) return -1;
	return 1;
} 
