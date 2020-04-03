double area_polygon()
{
	area = 0;
	for(i=0;i<sz;i++)
	{
		area += x[i]*y[(i+1)%sz] - y[i]*x[(i+1)%sz];
	}

	return area/2.;
}