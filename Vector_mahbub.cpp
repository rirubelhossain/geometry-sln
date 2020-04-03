struct Vector
{
	double x, y, z;

	Vector(){}
	Vector(double x, double y, double z) { this->x = x, this->y = y, this->z = z; }
	void input() { scanf("%lf%lf%lf",&x,&y,&z); }
	void output() { printf("%lf %lf %lf\n",x,y,z); }
	void set(double x, double y, double z) { this->x = x, this->y = y, this->z = z; }
	void unitify() { double d = length(); x/=d; y/=d; z/=d; }
	double length() { return sqrt( x*x + y*y + z*z ); }
	Vector operator+(Vector &a) { return Vector(a.x + x, a.y + y, a.z + z); }
	Vector operator-(Vector &a) { return Vector(x - a.x, y - a.y, z - a.z); }
	double operator*(Vector &a) { return (a.x*x + a.y*y + a.z*z); }
	Vector operator*(double a) { return Vector(a*x, a*y, a*z); }
	Vector operator^(Vector &a) { return Vector(y*a.z-z*a.y, z*a.x-x*a.z, x*a.y-y*a.x);	}
	void rotate_deg(Vector &a, double angle) { rotate_rad(a, angle/180*PI); }
	void rotate_rad(Vector w, double angle)
	{
		Vector a = *this;

		w.unitify();
		double d = a.length();
		a.unitify();
		
		Vector a_perp = w ^ a;
		Vector b = a*cos(angle) + a_perp*sin(angle);
		b = b*d;

		*this = b;
	}
	double dist(Vector &a) { return sqrt( S(a.x-x) + S(a.y-y) + S(a.z-z) ); }
};




struct Vector
{
	int x, y;

	Vector(){}
	Vector(int x, int y) { this->x = x, this->y = y;}
	void input() { scanf("%d%d",&x,&y); }
	void output() { printf("%d %d\n",x,y); }
	void set(int x, int y) { this->x = x, this->y = y; }
	void unitify() { double d = length(); x/=d; y/=d; }
	double length() { return sqrt( x*x + y*y ); }
	Vector operator+(Vector &a) { return Vector(a.x + x, a.y + y); }
	Vector operator-(Vector &a) { return Vector(x - a.x, y - a.y); }
	int operator*(Vector &a) { return (a.x*x + a.y*y); }
	Vector operator*(int a) { return Vector(a*x, a*y); }
	double dist(Vector &a) { return sqrt( S(a.x-x) + S(a.y-y) ); }
};
