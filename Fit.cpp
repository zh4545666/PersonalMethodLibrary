#include "stdafx.h"
#include "Fit.h"

CFit::CFit():ssr(0), sse(0), rmse(0)
{
	factorlen = 2;
	pfactor = new double[factorlen];
	fitedYslen = 0;
	pfitedYs = NULL;
}

CFit::~CFit()
{
	ClearData();
}

void CFit::ClearData()
{
	if (pfactor)
		delete[] pfactor;
	if (pfitedYs)
		delete[] pfitedYs;
	factorlen = 0;
	fitedYslen = 0;
}

double* CFit::getFactor(int& length)
{
	length = factorlen;
	if (factorlen <= 0)
		return NULL;

	double* pfac = new double[factorlen];
	memcpy(pfac, this->pfactor, sizeof(double)*factorlen);
	return pfac;
}

double* CFit::getFitedYs(int& length)
{
	length = fitedYslen;
	if (fitedYslen <= 0)
		return NULL;

	double* pYs = new double[fitedYslen];
	memcpy(pYs, this->pfitedYs, sizeof(double)*fitedYslen);
	return pYs;
}

string CFit::ToString()
{
	if (factorlen <= 0)
		return "";

	string s;
	ostringstream ostr;
	ostr << "y = ";
	for (int i = factorlen - 1; i >= 0; i--)
	{
		if (pfactor[i] == 0)
			continue;
		else if (pfactor[i] < 0)
			ostr << " - " << fabs(pfactor[i]);
		else
		{
			if (i == factorlen - 1)
				ostr << pfactor[i];
			else
				ostr << " + " << pfactor[i];
		}
		
		if (i > 1)
			ostr << "x^" << i;
		else if (i == 1)
			ostr << "x";
	}

	s = ostr.str();
	return s;
}

double CFit::getFactor(size_t i)
{
	if (0 == factorlen || i>=factorlen)
		return 0;
	return pfactor[i];
}

double CFit::getSlope()
{
	if (factorlen<2)
		return 0;
	return pfactor[1];
}

double CFit::getIntercept()
{
	if (0 == factorlen)
		return 0;
	return pfactor[0];
}

double CFit::getSSE()
{
	return sse;
}

double CFit::getSSR()
{
	return ssr;
}

double CFit::getRMSE()
{
	return rmse;
}

double CFit::getR_square()
{
	return 1 - (sse / (ssr + sse));
}

void CFit::FreeData(double **dat, double *d, int count)
{
	int i;
	free(d);
	for (i = 0; i < count; i++)
		free(dat[i]);
	free(dat);
}
// 解线性方程。data[count*(count+1)]矩阵数组；count：方程元数；
// Answer[count]：求解数组 。返回：0求解成功，-1无解或者无穷解
int CFit::LinearEquations(double *data, int count, double *Answer)
{
	int j, m, n;
	double tmp, **dat, *d = data;
	dat = (double**)malloc(count * sizeof(double*));
	for (m = 0; m < count; m++, d += (count + 1))
	{
		dat[m] = (double*)malloc((count + 1) * sizeof(double));
		memcpy(dat[m], d, (count + 1) * sizeof(double));
	}
	d = (double*)malloc((count + 1) * sizeof(double));
	for (m = 0; m < count - 1; m++)
	{
		// 如果主对角线元素为0，行交换
		for (n = m + 1; n < count && dat[m][m] == 0.0; n++)
		{
			if (dat[n][m] != 0.0)
			{
				memcpy(d, dat[m], (count + 1) * sizeof(double));
				memcpy(dat[m], dat[n], (count + 1) * sizeof(double));
				memcpy(dat[n], d, (count + 1) * sizeof(double));
			}
		}
		// 行交换后，主对角线元素仍然为0，无解，返回-1
		if (dat[m][m] == 0.0)
		{
			FreeData(dat, d, count);
			return -1;
		}
		// 消元
		for (n = m + 1; n < count; n++)
		{
			tmp = dat[n][m] / dat[m][m];
			for (j = m; j <= count; j++)
				dat[n][j] -= tmp * dat[m][j];
		}
	}
	for (j = 0; j < count; j++)
		d[j] = 0.0;
	// 求得count - 1的元
	Answer[count - 1] = dat[count - 1][count] / dat[count - 1][count - 1];
	// 逐行代入求各元
	for (m = count - 2; m >= 0; m--)
	{
		for (j = count - 1; j > m; j--)
			d[m] += Answer[j] * dat[m][j];
		Answer[m] = (dat[m][count] - d[m]) / dat[m][m];
	}
	FreeData(dat, d, count);
	return 0;
}
// 求多元回归方程：Y = B0 + B1X1 + B2X2 + ...BnXn
// data[rows*cols]二维数组；X1i,X2i,...Xni,Yi (i=0 to rows-1)
// rows：数据行数；cols数据列数；
// Answer[cols]：返回回归系数数组(B0,B1...Bn)
// SquarePoor[4]：返回方差分析指标: 回归平方和，剩余平方和，回归平方差，剩余平方差
// 返回值：0求解成功，-1错误
int CFit::MultipleRegression(double *data, int rows, int cols, double *Answer, double *SquarePoor)
{
	int m, n, i, count = cols - 1;
	double *dat, *p, a, b;
	if (data == 0 || Answer == 0 || rows < 2 || cols < 2)
		return -1;
	dat = (double*)malloc(cols * (cols + 1) * sizeof(double));
	dat[0] = (double)rows;
	for (n = 0; n < count; n++)                     // n = 0 to cols - 2
	{
		a = b = 0.0;
		for (p = data + n, m = 0; m < rows; m++, p += cols)
		{
			a += *p;
			b += (*p * *p);
		}
		dat[n + 1] = a;                              // dat[0, n+1] = Sum(Xn)
		dat[(n + 1) * (cols + 1)] = a;               // dat[n+1, 0] = Sum(Xn)
		dat[(n + 1) * (cols + 1) + n + 1] = b;       // dat[n+1,n+1] = Sum(Xn * Xn)
		for (i = n + 1; i < count; i++)             // i = n+1 to cols - 2
		{
			for (a = 0.0, p = data, m = 0; m < rows; m++, p += cols)
				a += (p[n] * p[i]);
			dat[(n + 1) * (cols + 1) + i + 1] = a;   // dat[n+1, i+1] = Sum(Xn * Xi)
			dat[(i + 1) * (cols + 1) + n + 1] = a;   // dat[i+1, n+1] = Sum(Xn * Xi)
		}
	}
	for (b = 0.0, m = 0, p = data + n; m < rows; m++, p += cols)
		b += *p;
	dat[cols] = b;                                   // dat[0, cols] = Sum(Y)
	for (n = 0; n < count; n++)
	{
		for (a = 0.0, p = data, m = 0; m < rows; m++, p += cols)
			a += (p[n] * p[count]);
		dat[(n + 1) * (cols + 1) + cols] = a;        // dat[n+1, cols] = Sum(Xn * Y)
	}
	n = CFit::LinearEquations(dat, cols, Answer);          // 计算方程式
	// 方差分析
	if (n == 0 && SquarePoor)
	{
		b = b / rows;                                // b = Y的平均值
		SquarePoor[0] = SquarePoor[1] = 0.0;
		p = data;
		for (m = 0; m < rows; m++, p++)
		{
			for (i = 1, a = Answer[0]; i < cols; i++, p++)
				a += (*p * Answer[i]);               // a = Ym的估计值
			SquarePoor[0] += ((a - b) * (a - b));    // U(回归平方和)
			SquarePoor[1] += ((*p - a) * (*p - a));  // Q(剩余平方和)(*p = Ym)
		}
		SquarePoor[2] = SquarePoor[0] / count;       // 回归方差
		if (rows - cols > 0.0)
			SquarePoor[3] = SquarePoor[1] / (rows - cols); // 剩余方差
		else
			SquarePoor[3] = 0.0;
	}
	free(dat);
	return n;
}

void CFit::Display(double *dat, double *Answer, double *SquarePoor, int rows, int cols)
{
	double v, *p;
	int i, j;
	printf("回归方程式:    Y = %.5lf", Answer[0]);
	for (i = 1; i < cols; i++)
		printf(" + %.5lf*X%d", Answer[i], i);
	printf(" \r\n");
	printf("回归显著性检验: \r\n");
	printf("回归平方和：%12.4lf  回归方差：%12.4lf \r\n", SquarePoor[0], SquarePoor[2]);
	printf("剩余平方和：%12.4lf  剩余方差：%12.4lf \r\n", SquarePoor[1], SquarePoor[3]);
	printf("离差平方和：%12.4lf  标准误差：%12.4lf \r\n", SquarePoor[0] + SquarePoor[1], sqrt(SquarePoor[3]));
	printf("F   检  验：%12.4lf  相关系数：%12.4lf \r\n", SquarePoor[2] / SquarePoor[3],
		sqrt(SquarePoor[0] / (SquarePoor[0] + SquarePoor[1])));
	printf("剩余分析: \r\n");
	printf("      观察值      估计值      剩余值    剩余平方 \r\n");
	for (i = 0, p = dat; i < rows; i++, p++)
	{
		v = Answer[0];
		for (j = 1; j < cols; j++, p++)
			v += *p * Answer[j];
		printf("%12.2lf%12.2lf%12.2lf%12.2lf \r\n", *p, v, *p - v, (*p - v) * (*p - v));
	}
	system("pause");
}

//pY = A + B*pM + C*pN
void CFit::GetParam(double* pM, double* pN, double* pY, int nCount, double& A, double& B, double& C)
{
	if (!pM || !pN || !pY || nCount<1)
		return;
	double a11 = 0, a12 = 0, a13 = 0, a21 = 0, a22 = 0, a23 = 0, a31 = 0, a32 = 0, a33 = 0, c1 = 0, c2 = 0, c3 = 0;
	for (int i = 0; i<nCount; i++)
	{
		double X1, X2, Y;
		X1 = pM[i];
		X2 = pN[i];
		Y = pY[i];

		a12 += X1;
		a13 += X2;
		a22 += pow(X1, 2);
		a23 += X1*X2;
		a33 += pow(X2, 2);
		c1 += Y;
		c2 += X1*Y;
		c3 += X2*Y;
	}

	a11 = nCount;
	a21 = a12;
	a31 = a13;
	a32 = a23;                     //求解 
	A = c1*(a22*a33 - a23*a32) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) - c2*(a12*a33 - a13*a32) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) + c3*(a12*a23 - a13*a22) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22);
	B = -c1*(a21*a33 - a23*a31) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) + c2*(a11*a33 - a13*a31) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) - c3*(a11*a23 - a13*a21) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22);
	C = c1*(a21*a32 - a22*a31) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) - c2*(a11*a32 - a12*a31) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22) + c3*(a11*a22 - a12*a21) / (a11*a22*a33 - a11*a23*a32 - a21*a12*a33 + a21*a13*a32 + a31*a12*a23 - a31*a13*a22);
}