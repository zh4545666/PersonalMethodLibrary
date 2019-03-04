#include "stdafx.h"
#include "IntegralTransform.h"

#ifndef PI
#define PI 3.14159265358
#endif

CIntegralTransform::CIntegralTransform()
{
}


CIntegralTransform::~CIntegralTransform()
{
}


void CIntegralTransform::FFT(CComplex* TD, CComplex* FD, int r)
{
	if (r <= 0 || !TD || !FD)
		return;
	int count = 1 << r;
	CComplex* W = new (std::nothrow) CComplex[count / 2];
	CComplex* X1 = new (std::nothrow) CComplex[count];
	CComplex* X2 = new (std::nothrow) CComplex[count];
	CComplex* X = NULL;
	if (!W || !X1 || !X2)
		return;

	int  i, j, k;
	int  dist, p;
	double f = 2*PI / count;
	double a = 0;
	for (i = 0; i < count / 2; ++i)
	{
		W[i] = CComplex(cos(a), -sin(a));
		a += f;
	}

	for (i = 0; i < count; ++i)
	{
		X1[i] = TD[i];
	}

	for (k = 0; k < r; ++k)
	{
		for (j = 0; j < (1 << k); ++j)
		{
			dist = 1 << (r - k);
			for (i = 0; i < dist / 2; ++i)
			{
				p = j * dist;
				X2[i + p] = X1[i + p] + X1[i + p + dist / 2];
				X2[i + p + dist / 2] = (X1[i + p] - X1[i + p + dist / 2])* W[i * (1 << k)];
			}
		}
		X = X1;
		X1 = X2;
		X2 = X;
	}

	for (j = 0; j < count; ++j)
	{
		p = 0;
		for (i = 0; i < r; ++i)
		{
			if (j&(1 << i))
			{
				p += 1 << (r - i - 1);
			}
		}
		FD[j] = X1[p];
	}

	delete[] W;
	delete[] X1;
	delete[] X2;
}

void CIntegralTransform::FFT(double* TD, CComplex* FD, int r)
{
	if (r <= 0 || !TD || !FD)
		return;
	int count = 1 << r;
	CComplex* W = new (std::nothrow) CComplex[count / 2];
	CComplex* X1 = new (std::nothrow) CComplex[count];
	CComplex* X2 = new (std::nothrow) CComplex[count];
	CComplex* X = NULL;
	if (!W || !X1 || !X2)
		return;

	int  i, j, k;
	int  dist, p;
	double f = 2 * PI / count;
	double a = 0;
	for (i = 0; i < count / 2; ++i)
	{
		W[i] = CComplex(cos(a), -sin(a));
		a += f;
	}

	for (i = 0; i < count; ++i)
	{
		X1[i].SetReal(TD[i]);
		X1[i].SetImag(0.0);
	}

	for (k = 0; k < r; ++k)
	{
		for (j = 0; j < (1 << k); ++j)
		{
			dist = 1 << (r - k);
			for (i = 0; i < dist / 2; ++i)
			{
				p = j * dist;
				X2[i + p] = X1[i + p] + X1[i + p + dist / 2];
				X2[i + p + dist / 2] = (X1[i + p] - X1[i + p + dist / 2])* W[i * (1 << k)];
			}
		}
		X = X1;
		X1 = X2;
		X2 = X;
	}

	for (j = 0; j < count; ++j)
	{
		p = 0;
		for (i = 0; i < r; ++i)
		{
			if (j&(1 << i))
			{
				p += 1 << (r - i - 1);
			}
		}
		FD[j] = X1[p];
	}

	delete[] W;
	delete[] X1;
	delete[] X2;
}

void CIntegralTransform::IFFT(CComplex* FD, CComplex* TD, int r)
{
	if (r <= 0 || !TD || !FD)
		return;

	int count = 1 << r;
	CComplex X;
	
	if (FD != TD)
		memcpy(TD, FD, sizeof(CComplex)*count);

	int i, j, a, b, p;
	// 位反转置换 Bit-reversal Permutation
	for (i = 1, p = 0; i < count; i *= 2)
		p++;
	for (i = 0; i < count; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			X = TD[i];
			TD[i] = TD[b];
			TD[b] = X;
		}
	}

	double treal, timag, ureal, uimag, arg;
	CComplex* W = new (std::nothrow) CComplex[count / 2];
	if (!W)
		return;
	int m, k, t, index1, index2;
	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / count;
	treal = cos(arg);
	timag = sin(arg);
	W[0] = CComplex(1.0, 0);
	for (j = 1; j < count / 2; j++)
	{
		W[j].SetReal(W[j - 1].GetReal()*treal - W[j - 1].GetImag()*timag);
		W[j].SetImag(W[j - 1].GetReal()*timag + W[j - 1].GetImag()*treal);
	}

	for (m = 2; m <= count; m *= 2)
	{
		for (k = 0; k < count; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = (count / m) * j;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = W[t].GetReal()*TD[index2].GetReal() - W[t].GetImag()*TD[index2].GetImag();
				timag = W[t].GetReal()*TD[index2].GetImag() + W[t].GetImag()*TD[index2].GetReal();
				ureal = TD[index1].GetReal();
				uimag = TD[index1].GetImag();
				TD[index1].SetReal(ureal + treal);
				TD[index1].SetImag(uimag + timag);
				TD[index2].SetReal(ureal - treal);
				TD[index2].SetImag(uimag - timag);
			}
		}
	}

	for (j = 0; j < count; j++)
		TD[j] = TD[j] * CComplex(1 / count, 0);
}

void CIntegralTransform::IFFT(CComplex* FD, double* TD, int r)
{
	if (r <= 0 || !TD || !FD)
		return;

	int count = 1 << r;
	CComplex X;

	CComplex* Y = new CComplex[count];
	memcpy(Y, FD, sizeof(CComplex)*count);

	int i, j, a, b, p;
	// 位反转置换 Bit-reversal Permutation
	for (i = 1, p = 0; i < count; i *= 2)
		p++;
	for (i = 0; i < count; i++)
	{
		a = i;
		b = 0;
		for (j = 0; j < p; j++)
		{
			b = (b << 1) + (a & 1);    // b = b * 2 + a % 2;
			a >>= 1;        // a = a / 2;
		}
		if (b > i)
		{
			X = Y[i];
			Y[i] = Y[b];
			Y[b] = X;
		}
	}

	double treal, timag, ureal, uimag, arg;
	CComplex* W = new (std::nothrow) CComplex[count / 2];
	if (!W)
		return;
	int m, k, t, index1, index2;
	// 计算 1 的前 n / 2 个 n 次方根 Wj = wreal [j] + i * wimag [j] , j = 0, 1, ... , n / 2 - 1
	arg = 2 * PI / count;
	treal = cos(arg);
	timag = sin(arg);
	W[0] = CComplex(1.0, 0);
	for (j = 1; j < count / 2; j++)
	{
		W[j].SetReal(W[j - 1].GetReal()*treal - W[j - 1].GetImag()*timag);
		W[j].SetImag(W[j - 1].GetReal()*timag + W[j - 1].GetImag()*treal);
	}

	for (m = 2; m <= count; m *= 2)
	{
		for (k = 0; k < count; k += m)
		{
			for (j = 0; j < m / 2; j++)
			{
				index1 = k + j;
				index2 = index1 + m / 2;
				t = (count / m) * j;    // 旋转因子 w 的实部在 wreal [] 中的下标为 t
				treal = W[t].GetReal()*Y[index2].GetReal() - W[t].GetImag()*Y[index2].GetImag();
				timag = W[t].GetReal()*Y[index2].GetImag() + W[t].GetImag()*Y[index2].GetReal();
				ureal = Y[index1].GetReal();
				uimag = Y[index1].GetImag();
				Y[index1].SetReal(ureal + treal);
				Y[index1].SetImag(uimag + timag);
				Y[index2].SetReal(ureal - treal);
				Y[index2].SetImag(uimag - timag);
			}
		}
	}

	for (j = 0; j < count; j++)
	{
		double dValueReal = Y[j].GetReal() / count;
		double dValueImag = Y[j].GetImag() / count;

		TD[j] = dValueReal;
	}
}

void CIntegralTransform::DFT(CComplex* TD, CComplex* FD, size_t N)
{
	size_t i, j, k;
	if (!TD || !FD || N == 0)
		return;
	double* _cosQik = new double[N*N];
	double* _sinQik = new double[N*N];

	double Qk, Qik;
	double Q = 2 * PI / N;

	Qk = 0;
	for (j = k = 0; k<N; k++)
	{
		Qik = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			_cosQik[j] = cos(Qik);
			_sinQik[j] = sin(Qik);
			Qik += Qk;
		}
		Qk += Q;
	}

	for (j = k = 0; k<N; k++)
	{
		double dReal = 0.0;
		double dImag = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			dReal += TD[i].GetReal() * _cosQik[j];
			dImag -= TD[i].GetReal() * _sinQik[j];
		}

		FD[k].SetReal(dReal);
		FD[k].SetImag(dImag);
	}
	delete[] _cosQik;
	delete[] _sinQik;
}

void CIntegralTransform::DFT(double* TD, CComplex* FD, size_t N)
{
	size_t i, j, k;
	if (!TD || !FD || N == 0)
		return;
	double* _cosQik = new double[N*N];
	double* _sinQik = new double[N*N];

	double Qk, Qik;
	double Q = 2 * PI / N;

	Qk = 0;
	for (j = k = 0; k<N; k++)
	{
		Qik = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			_cosQik[j] = cos(Qik);
			_sinQik[j] = sin(Qik);
			Qik += Qk;
		}
		Qk += Q;
	}

	for (j = k = 0; k<N; k++)
	{
		double dReal = 0.0;
		double dImag = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			dReal += TD[i] * _cosQik[j];
			dImag -= TD[i] * _sinQik[j];
		}

		FD[k].SetReal(dReal);
		FD[k].SetImag(dImag);
	}
	delete[] _cosQik;
	delete[] _sinQik;
}

void CIntegralTransform::IDFT(CComplex* FD, CComplex* TD, size_t N)
{
	size_t i, j, k;
	if (!TD || !FD || N == 0)
		return;
	double* _cosQik = new double[N*N];
	double* _sinQik = new double[N*N];

	double Qk, Qik;
	double Q = 2 * PI / N;

	Qk = 0;
	for (j = k = 0; k<N; k++)
	{
		Qik = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			_cosQik[j] = cos(Qik);
			_sinQik[j] = sin(Qik);
			Qik += Qk;
		}
		Qk += Q;
	}

	for (j = k = 0; k<N; k++)
	{
		double sumX = 0.0;
		double sumY = 0.0;
		for (i = 0; i < N; i++, j++)
		{
			sumX += FD[i].GetReal() * _cosQik[j] - FD[i].GetImag() * _sinQik[j];
			sumY += FD[i].GetReal() * _sinQik[j] + FD[i].GetImag() * _cosQik[j];
		}

		TD[k].SetReal(sumX / N);
		TD[k].SetImag(sumY / N);
	}
	delete[] _cosQik;
	delete[] _sinQik;
}

void CIntegralTransform::IDFT(CComplex* FD, double* TD, size_t N)
{
	size_t i, j, k;
	if (!TD || !FD || N == 0)
		return;
	double* _cosQik = new double[N*N];
	double* _sinQik = new double[N*N];

	double Qk, Qik;
	double Q = 2 * PI / N;

	Qk = 0;
	for (j = k = 0; k<N; k++)
	{
		Qik = 0.0;
		for (i = 0; i<N; i++, j++)
		{
			_cosQik[j] = cos(Qik);
			_sinQik[j] = sin(Qik);
			Qik += Qk;
		}
		Qk += Q;
	}

	for (j = k = 0; k<N; k++)
	{
		double sumX = 0.0;
		for (i = 0; i<N; i++, j++)
			sumX += FD[i].GetReal() * _cosQik[j] - FD[i].GetImag() * _sinQik[j];

		TD[k] = sumX / N;
		if (fabs(TD[k]) < 0.0000000001)
			TD[k] = 0.0;

	}
	delete[] _cosQik;
	delete[] _sinQik;
}

void CIntegralTransform::HilBert(double* TD, CComplex* FD, size_t N)
{
	if (!TD || !FD || N <= 0)
		return;

	size_t i;
	CComplex* H = new CComplex[N];
	DFT(TD, H, N);

	for (i = 1; i < (N + 1) / 2; i++)
		H[i] = H[i] * CComplex(2.0, 0.0);
	for (i = N / 2 + 1; i < N; i++)
		H[i] = CComplex(0, 0);

	IDFT(H, FD, N);

	delete[] H;
}
void CIntegralTransform::HilBert(CComplex* TD, CComplex* FD, size_t N)
{
	if (!TD || !FD || N <= 0)
		return;

	size_t i;
	CComplex* H = new CComplex[N];
	DFT(TD, H, N);

	for (i = 1; i < (N + 1) / 2; i++)
		H[i] = H[i] * CComplex(2.0, 0.0);
	for (i = N / 2 + 1; i < N; i++)
		H[i] = CComplex(0, 0);

	IDFT(H, FD, N);

	delete[] H;
}

inline double CIntegralTransform::gauss(int n, int m)
{
	return exp(-2. * PI * PI * m * m / (n * n));
}

void CIntegralTransform::STransform(double* TD, CComplex* FD, size_t N)
{
	if (!TD || !FD || N <= 0)
		return;

	size_t i, k;
	double s = 0.0;
	for (i = 0; i < N; i++)
		s += TD[i];
	s /= N;

	CComplex* H = new CComplex[N];
	CComplex* G = new CComplex[N];
	double* g = new double[N];

	/* FFT. */
	DFT(TD, H, N);

	/* Hilbert transform. The upper half-circle gets multiplied by
	two, and the lower half-circle gets set to zero.  The real axis
	is left alone. */
	for (i = 1; i < (N + 1) / 2; i++)
		H[i] = H[i] * CComplex(2.0, 0);
	for (i = N / 2 + 1; i < N; i++)
		H[i] = CComplex(0, 0);

	size_t n = 0;
	for (i = 0; i < N; i++) 
		FD[i + n*N] = CComplex(s, 0.0);
	n++;

	while (n <= N/2) {

		/* Scale the FFT of the gaussian. Negative frequencies
		wrap around. */
		g[0] = gauss(n, 0);
		for (i = 1; i < N / 2 + 1; i++)
			g[i] = g[N - i] = gauss(n, i);


		for (i = 0; i < N; i++) {
			s = g[i];
			k = n + i;
			if (k >= N) k -= N;
			G[i] = H[k] * CComplex(s, 0);
		}

		/* Inverse FFT the result to get the next row. */
		CComplex* h = new CComplex[N];
		IDFT(G, h, N);
		for (i = 0; i < N; i++)
			FD[i + n*N] = CComplex(h[i].GetReal() / 2, h[i].GetImag() / 2);
		delete[] h;

		/* Go to the next row. */
		n++;
	}

	delete[] H;
	delete[] G;
	delete[] g;
}

void CIntegralTransform::ISTransform(CComplex* FD, double* TD, size_t N)
{
	if (!TD || !FD || N <= 0)
		return;

	size_t i,j;
	CComplex* H = new CComplex[N];
	CComplex* h = new CComplex[N];

	/* Sum the complex array across time. */
	for (i = 0; i <= N/2; i++)
	{
		double sumX = 0.0;
		double sumY = 0.0;
		for (j = 0; j < N; j++)
		{
			sumX += FD[i*N + j].GetReal();
			sumY += FD[i*N + j].GetImag();
		}
		H[i] = CComplex(sumX, sumY);
	}

	/* Invert the Hilbert transform. */
	for (i = 1; i < (N + 1) / 2; i++)
		H[i] = H[i] * CComplex(2, 0);
	for (i = N / 2 + 1; i < N; i++)
		H[i] = CComplex(H[N - i].GetReal(), -H[N - i].GetImag());

	/* Inverse FFT. */
	IDFT(H, h, N);
	for (i = 0; i < N; i++)
		TD[i] = h[i].GetReal() / 2;

	delete[] H;
	delete[] h;
}

void CIntegralTransform::DCT(double* TD, double* FD, size_t N)
{
	if (!TD || !FD || N == 0)
		return;

	double c1 = sqrt(1.0f / N);
	double c2 = sqrt(2.0f / N);

	for (size_t i = 0; i < N; i++)
	{
		FD[i] = 0;
		for (size_t j = 0; j < N; j++)
		{
			FD[i] += TD[j] * cos((2.0f * j + 1.0f) * PI * i / (2.0f * N));
		}

		if (i == 0)
			FD[i] *= c1;
		else
			FD[i] *= c2;
	}

}

void CIntegralTransform::IDCT(double* FD, double* TD, size_t N)
{
	if (!TD || !FD || N == 0)
		return;

	double c1 = sqrt(1.0 / N);
	double c2 = sqrt(2.0 / N);

	for (size_t i = 0; i < N; i++)
	{
		TD[i] = 0;
		for (size_t j = 0; j < N; j++)
		{
			if (j == 0)
				TD[i] += c1 * FD[j] * cos((2.0 * i + 1.0) * PI * j / (2.0 * N));
			else
				TD[i] += c2 * FD[j] * cos((2.0 * i + 1.0) * PI * j / (2.0 * N));
		}
	}
}

void CIntegralTransform::DCT2(CMatrix* TD, CMatrix* FD)
{
	if (!TD || !FD || TD->GetNumRows() != TD->GetNumColumns() || FD->GetNumRows() != FD->GetNumColumns() ||
		TD->GetNumRows() != FD->GetNumRows())
		return;

	size_t N = TD->GetNumRows();
	CMatrix* A = new CMatrix(N, N);


	for (size_t i = 0; i < N; i++)
	{
		double c;
		if (i == 0)
			c = sqrt(1.0 / N);
		else
			c = sqrt(2.0 / N);


		for (size_t j = 0; j < N; j++)
		{
			A->SetElement(i, j, c*cos((2.0 *j + 1) * PI * i / (2.0 * N)));
		}
	}

	*FD = (*A) * (*TD) * (A->Transpose());
}

void CIntegralTransform::IDCT2(CMatrix* FD, CMatrix* TD)
{
	if (!TD || !FD || TD->GetNumRows() != TD->GetNumColumns() || FD->GetNumRows() != FD->GetNumColumns() ||
		TD->GetNumRows() != FD->GetNumRows())
		return;

	size_t N = TD->GetNumRows();
	CMatrix* A = new CMatrix(N, N);


	for (size_t i = 0; i < N; i++)
	{
		double c;
		if (i == 0)
			c = sqrt(1.0 / N);
		else
			c = sqrt(2.0 / N);

		for (size_t j = 0; j < N; j++)
		{
			A->SetElement(i, j, c*cos((2.0 *j + 1) * PI * i / (2.0 * N)));
		}
	}

	*TD = (A->Transpose()) * (*FD) * (*A);
}