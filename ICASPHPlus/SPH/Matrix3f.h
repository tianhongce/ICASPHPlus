#pragma once
#include <vector>
#include <math.h>
#include <algorithm>

//////////////////////////////////////////////////////////////////////////
//Matrix3f类
//作用：用于刚体模拟当中的旋转等数据，在SPH计算当中，没有使用该类。
//////////////////////////////////////////////////////////////////////////
class Matrix3f
{
public:
	double value[3][3];

public:
	Matrix3f(){}
	Matrix3f(Matrix3f& v)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				value[i][j] = v[i][j];
			}
		}
	}
	double* operator[](int x);
	const double* operator[](const int x)const;
	//重载赋值操作符
	Matrix3f& operator=(Matrix3f& a);
	Matrix3f& operator=(const Matrix3f& a);
};