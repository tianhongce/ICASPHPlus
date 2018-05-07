#include "Matrix3f.h"

double * Matrix3f::operator[](int x)
{
	return value[x];
}

const double * Matrix3f::operator[](const int x) const
{
	return value[x];
}

//对数据进行拷贝
Matrix3f & Matrix3f::operator=(Matrix3f & a)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			value[i][j] = a[i][j];
		}
	}
	return *this;
}

Matrix3f & Matrix3f::operator=(const Matrix3f & a)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			value[i][j] = a[i][j];
		}
	}
	return *this;
}
