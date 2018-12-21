#pragma once

#ifdef AFX_CLASS  
#define AFX_EX_CLASS _declspec(dllexport)  
#else  
#define AFX_EX_CLASS _declspec(dllimport)  
#endif 

class AFX_EX_CLASS CFit
{
public:
	CFit();
	~CFit();

private:
	size_t factorlen;
	double* pfactor;
	size_t fitedYslen;
	double* pfitedYs;

	double ssr;
	double sse;
	double rmse;
public:
	//
	/// \拟合y=a0+a1*x
	/// \brief 直线拟合-一元回归,拟合的结果可以使用getFactor获取，或者使用getSlope获取斜率，getIntercept获取截距
	/// \param x 观察值的x
	/// \param y 观察值的y
	/// \param length x,y数组的长度
	/// \param isSaveFitYs 拟合后的数据是否保存，默认否
	///
	template<typename T>
	bool linearFit(const T* x, const T* y, size_t length, bool isSaveFitYs = false)
	{
		ClearData();
		factorlen = 2;
		pfactor = new double[factorlen];

		if (0 == length || !x || !y)
			return false;

		typename T t1 = 0, t2 = 0, t3 = 0, t4 = 0;
		for (size_t i = 0; i<length; ++i)
		{
			t1 += x[i] * x[i];
			t2 += x[i];
			t3 += x[i] * y[i];
			t4 += y[i];
		}
		pfactor[1] = (t3*length - t2*t4) / (t1*length - t2*t2);
		pfactor[0] = (t1*t4 - t2*t3) / (t1*length - t2*t2);
		//////////////////////////////////////////////////////////////////////////
		//计算误差
		calcError(x, y, length, this->ssr, this->sse, this->rmse, isSaveFitYs);
		return true;
	}
	///
	/// \brief 多项式拟合，拟合y=a0+a1*x+a2*x^2+……+apoly_n*x^poly_n
	/// \param x 观察值的x
	/// \param y 观察值的y
	/// \param poly_n 期望拟合的阶数，若poly_n=2，则y=a0+a1*x+a2*x^2
	/// \param isSaveFitYs 拟合后的数据是否保存，默认是
	/// 
	template<typename T>
	void polyfit(const T* x, const T* y, size_t length, size_t poly_n, bool isSaveFitYs = false)
	{
		ClearData();
		factorlen = poly_n + 1;
		pfactor = new double[factorlen];

		if (0 == length || !x || !y)
			return;

		size_t i, j;
		std::vector<double> tempx(length, 1.0);
		std::vector<double> tempy(y, y + length);
		std::vector<double> sumxx(poly_n * 2 + 1);
		std::vector<double> ata((poly_n + 1)*(poly_n + 1));
		std::vector<double> sumxy(poly_n + 1);

		for (i = 0; i<2 * poly_n + 1; i++){
			for (sumxx[i] = 0, j = 0; j<length; j++)
			{
				sumxx[i] += tempx[j];
				tempx[j] *= x[j];
			}
		}
		for (i = 0; i<poly_n + 1; i++){
			for (sumxy[i] = 0, j = 0; j<length; j++)
			{
				sumxy[i] += tempy[j];
				tempy[j] *= x[j];
			}
		}
		for (i = 0; i<poly_n + 1; i++)
		for (j = 0; j<poly_n + 1; j++)
			ata[i*(poly_n + 1) + j] = sumxx[i + j];
		gauss_solve(poly_n + 1, &ata[0], &pfactor[0], &sumxy[0]);
		//计算拟合后的数据并计算误差
		calcError(&x[0], &y[0], length, this->ssr, this->sse, this->rmse, isSaveFitYs);
	}

	//拟合方程的系数
	double getFactor(size_t i);
	//返回x对应的y值
	template<typename T>
	double getY(const T x) const
	{
		double ans(0);
		for (size_t i = 0; i<factorlen; ++i)
			ans += pfactor[i] * pow((double)x, (int)i);
		return ans;
	}
	//斜率值
	double getSlope();
	//截距值
	double getIntercept();
	//剩余平方和
	double getSSE();
	//回归平方和
	double getSSR();
	//均方根误差
	double getRMSE();
	//确定系数，系数是0~1之间的数，是数理上判定拟合优度的一个量
	double getR_square();
	//均值
	template <typename T>
	static T Mean(const T* v, size_t length)
	{
		T total(0);
		for (size_t i = 0; i<length; ++i)
		{
			total += v[i];
		}
		return (total / length);
	}

	double* getFactor(int& length);
	double* getFitedYs(int& length);
	string ToString();
private:
	void ClearData();
	template<typename T>
	void calcError(const T* x
		, const T* y
		, size_t length
		, double& r_ssr
		, double& r_sse
		, double& r_rmse
		, bool isSaveFitYs = true
		)
	{
		fitedYslen = 0;
		if (pfitedYs)
			delete[] pfitedYs;
		if (0 == length || !x || !y)
			return;
		if (isSaveFitYs)
			pfitedYs = new double[length];

		T mean_y = Mean/*<T>*/(y, length);
		T yi(0);
		for (size_t i = 0; i<length; ++i)
		{
			yi = (T)getY(x[i]);
			r_ssr += ((yi - mean_y)*(yi - mean_y));//计算回归平方和
			r_sse += ((yi - y[i])*(yi - y[i]));//残差平方和
			if (isSaveFitYs)
				pfitedYs[i] = double(yi);
		}
		r_rmse = sqrt(r_sse / (double(length)));
	}
	template<typename T>
	void gauss_solve(int n, T* A, T* x, T* b)
	{
		int i, j, k, r;
		double max;
		for (k = 0; k<n - 1; k++)
		{
			max = fabs(A[k*n + k]); /*find maxmum*/
			r = k;
			for (i = k + 1; i<n - 1; i++){
				if (max<fabs(A[i*n + i]))
				{
					max = fabs(A[i*n + i]);
					r = i;
				}
			}
			if (r != k){
				for (i = 0; i<n; i++)         /*change array:A[k]&A[r] */
				{
					max = A[k*n + i];
					A[k*n + i] = A[r*n + i];
					A[r*n + i] = max;
				}
			}
			max = b[k];                    /*change array:b[k]&b[r]     */
			b[k] = b[r];
			b[r] = max;
			for (i = k + 1; i<n; i++)
			{
				for (j = k + 1; j<n; j++)
					A[i*n + j] -= A[i*n + k] * A[k*n + j] / A[k*n + k];
				b[i] -= A[i*n + k] * b[k] / A[k*n + k];
			}
		}

		for (i = n - 1; i >= 0; x[i] /= A[i*n + i], i--)
		for (j = i + 1, x[i] = b[i]; j<n; j++)
			x[i] -= A[i*n + j] * x[j];
	}

public:
	//二元一次拟合  pY = A + B*pM + C*pN
	static void GetParam(double* pM, double* pN, double* pY, int nCount, double& A, double& B, double& C);
	//
	//多元回归方程：Y = B0 + B1X1 + B2X2 + ...BnXn
	// data[rows*cols]二维数组；X1i,X2i,...Xni,Yi (i=0 to rows-1)
	// rows：数据行数；cols数据列数；
	// Answer[cols]：返回回归系数数组(B0,B1...Bn)
	// SquarePoor[4]：返回方差分析指标: 回归平方和，剩余平方和，回归平方差，剩余平方差
	// 返回值：0求解成功，-1错误
	//
	static int MultipleRegression(double *data, int rows, int cols, double *Answer, double *SquarePoor);
	static void Display(double *dat, double *Answer, double *SquarePoor, int rows, int cols);
private:
	static void FreeData(double **dat, double *d, int count);
	static int LinearEquations(double *data, int count, double *Answer);
}; 
