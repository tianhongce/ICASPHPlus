#pragma once
#include <vector>
#include <math.h>
#include <algorithm>

//////////////////////////////////////////////////////////////////////////
//Matrix3f��
//���ã����ڸ���ģ�⵱�е���ת�����ݣ���SPH���㵱�У�û��ʹ�ø��ࡣ
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
	//���ظ�ֵ������
	Matrix3f& operator=(Matrix3f& a);
	Matrix3f& operator=(const Matrix3f& a);
};