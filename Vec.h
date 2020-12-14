#ifndef VEC_H
#define VEC_H
#include<vector>
#include<cmath>
using namespace std;

struct vec2{
	double x, y;
	vec2(double x, double y):x(x),y(y){}
	vec2(){}
	vec2 operator * (const double & d) const{
		return vec2(x*d, y*d);
	}
	vec2 operator / (const double & d) const{
		return vec2(x/d, y/d);
	}
	vec2 operator - (const vec2 & ano) const{
		return vec2(x-ano.x, y-ano.y);
	}
	vec2 operator + (const vec2 & ano) const{
		return vec2(x+ano.x, y+ano.y);
	}
	bool operator == (const vec2 & ano){
		return x == ano.x && y==ano.y;
	}
};


typedef vector<vec2> vec2l;

inline double dis(vec2 p1, vec2 p2){
	return hypot(p1.x-p2.x, p1.y-p2.y);
}

inline double norm(vec2 p){
	return p.x*p.x + p.y*p.y;
}

struct vec3{
	double x, y, z;
	vec3(double x, double y, double z):x(x), y(y), z(z){}
	vec3(){}
	vec3 operator * (const double & d) const{
		return vec3(x*d, y*d, z*d);
	}
	vec3 operator / (const double & d) const{
		return vec3(x/d, y/d, z/d);
	}
	vec3 operator + (const vec3 &ano) const{
		return vec3(x+ano.x, y+ano.y, z+ano.z);
	}
	vec3 operator - (const vec3 &ano) const{
		return vec3(x-ano.x, y-ano.y, z-ano.z);
	}
};

inline double dis(vec3 p1, vec3 p2){
	return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

inline double norm(vec3 p){
	return p.x * p.x + p.y * p.y + p.z * p.z;
}
#endif