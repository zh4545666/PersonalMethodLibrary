#include "stdafx.h"
#include "BigNum.h"

#define MAXN 9999
#define MAXSIZE 10
#define DLEN 4
#include <map>

CBigNum::CBigNum(const int b)     //��һ��int���͵ı���ת��Ϊ����
{
	int c, d = b;
	len = 0;
	memset(a, 0, sizeof(a));
	while (d > MAXN)
	{
		c = d - (d / (MAXN + 1)) * (MAXN + 1);
		d = d / (MAXN + 1);
		a[len++] = c;
	}
	a[len++] = d;
}
CBigNum::CBigNum(const char*s, bool bIgnoreNonnumeric)     //��һ���ַ������͵ı���ת��Ϊ����
{

	if (bIgnoreNonnumeric)
	{
		memset(a, 0, sizeof(a));
		map<char, int> m; int i, j, cout = 0;
		for (i = 0; i < 10; i++)
			m['0' + i] = i;

		int t, index = 0;
		for (i = strlen(s) - 1; i >= 0;)
		{
			int cout = 0;
			t = 0;
			for (j = i; j >= 0; j--)
			{
				if (0 == m.count(s[j]))
					continue;

				int value = m[s[j]];
				int k = cout;
				while (k >= 1)
				{
					value *= 10;
					k--;
				}

				t += value;
				cout++;
				if (cout == 4)
				{
					a[index++] = t;
					t = 0; i = j - 1;
					break;
				}
			}
			if (t != 0)
			{
				a[index++] = t;
				break;
			}
		}
		len = index;

		return;
	}

	int t, k, index, l, i;
	memset(a, 0, sizeof(a));
	l = strlen(s);
	len = l / DLEN;
	if (l%DLEN)
		len++;
	index = 0;
	for (i = l - 1; i >= 0; i -= DLEN)
	{
		t = 0;
		k = i - DLEN + 1;
		if (k<0)
			k = 0;
		for (int j = k; j <= i; j++)
			t = t * 10 + s[j] - '0';
		a[index++] = t;
	}
}
CBigNum::CBigNum(const CBigNum & T) : len(T.len)  //�������캯��
{
	int i;
	memset(a, 0, sizeof(a));
	for (i = 0; i < len; i++)
		a[i] = T.a[i];
}
CBigNum & CBigNum::operator=(const CBigNum & n)   //���ظ�ֵ�����������֮����и�ֵ����
{
	int i;
	len = n.len;
	memset(a, 0, sizeof(a));
	for (i = 0; i < len; i++)
		a[i] = n.a[i];
	return *this;
}
istream& operator>>(istream & in, CBigNum & b)   //�������������
{
	char ch[MAXSIZE * 4];
	int i = -1;
	in >> ch;
	int l = strlen(ch);
	int count = 0, sum = 0;
	for (i = l - 1; i >= 0;)
	{
		sum = 0;
		int t = 1;
		for (int j = 0; j<4 && i >= 0; j++, i--, t *= 10)
		{
			sum += (ch[i] - '0')*t;
		}
		b.a[count] = sum;
		count++;
	}
	b.len = count++;
	return in;
}
ostream& operator<<(ostream& out, CBigNum& b)   //������������
{
	int i;
	cout << b.a[b.len - 1];
	for (i = b.len - 2; i >= 0; i--)
	{
		cout.width(DLEN);
		cout.fill('0');
		cout << b.a[i];
	}
	return out;
}
CBigNum CBigNum::operator+(const CBigNum & T) const   //��������֮����������
{
	CBigNum t(*this);
	int i, big;      //λ��   
	big = T.len > len ? T.len : len;
	for (i = 0; i < big; i++)
	{
		t.a[i] += T.a[i];
		if (t.a[i] > MAXN)
		{
			t.a[i + 1]++;
			t.a[i] -= MAXN + 1;
		}
	}
	if (t.a[big] != 0)
		t.len = big + 1;
	else
		t.len = big;
	return t;
}
CBigNum CBigNum::operator-(const CBigNum & T) const   //��������֮���������� 
{
	int i, j, big;
	bool flag;
	CBigNum t1, t2;
	if (*this>T)
	{
		t1 = *this;
		t2 = T;
		flag = 0;
	}
	else
	{
		t1 = T;
		t2 = *this;
		flag = 1;
	}
	big = t1.len;
	for (i = 0; i < big; i++)
	{
		if (t1.a[i] < t2.a[i])
		{
			j = i + 1;
			while (t1.a[j] == 0)
				j++;
			t1.a[j--]--;
			while (j > i)
				t1.a[j--] += MAXN;
			t1.a[i] += MAXN + 1 - t2.a[i];
		}
		else
			t1.a[i] -= t2.a[i];
	}
	t1.len = big;
	while (t1.a[len - 1] == 0 && t1.len > 1)
	{
		t1.len--;
		big--;
	}
	if (flag)
		t1.a[big - 1] = 0 - t1.a[big - 1];
	return t1;
}
CBigNum CBigNum::operator*(const CBigNum & T) const   //��������֮���������� 
{
	CBigNum ret;
	int i, j, up;
	int temp, temp1;
	for (i = 0; i < len; i++)
	{
		up = 0;
		for (j = 0; j < T.len; j++)
		{
			temp = a[i] * T.a[j] + ret.a[i + j] + up;
			if (temp > MAXN)
			{
				temp1 = temp - temp / (MAXN + 1) * (MAXN + 1);
				up = temp / (MAXN + 1);
				ret.a[i + j] = temp1;
			}
			else
			{
				up = 0;
				ret.a[i + j] = temp;
			}
		}
		if (up != 0)
			ret.a[i + j] = up;
	}
	ret.len = i + j;
	while (ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--;
	return ret;
}
CBigNum CBigNum::operator/(const int & b) const   //������һ�����������������
{
	CBigNum ret;
	int i, down = 0;
	for (i = len - 1; i >= 0; i--)
	{
		ret.a[i] = (a[i] + down * (MAXN + 1)) / b;
		down = a[i] + down * (MAXN + 1) - ret.a[i] * b;
	}
	ret.len = len;
	while (ret.a[ret.len - 1] == 0 && ret.len > 1)
		ret.len--;
	return ret;
}
int CBigNum::operator %(const int & b) const    //������һ��int���͵ı�������ȡģ����    
{
	int i, d = 0;
	for (i = len - 1; i >= 0; i--)
	{
		d = ((d * (MAXN + 1)) % b + a[i]) % b;
	}
	return d;
}
CBigNum CBigNum::operator^(const int & n) const    //������n�η�����
{
	CBigNum t, ret(1);
	int i;
	if (n<0)
		exit(-1);
	if (n == 0)
		return 1;
	if (n == 1)
		return *this;
	int m = n;
	while (m>1)
	{
		t = *this;
		for (i = 1; i << 1 <= m; i <<= 1)
		{
			t = t*t;
		}
		m -= i;
		ret = ret*t;
		if (m == 1)
			ret = ret*(*this);
	}
	return ret;
}
bool CBigNum::operator>(const CBigNum & T) const   //��������һ�������Ĵ�С�Ƚ�
{
	int ln;
	if (len > T.len)
		return true;
	else if (len == T.len)
	{
		ln = len - 1;
		while (a[ln] == T.a[ln] && ln >= 0)
			ln--;
		if (ln >= 0 && a[ln] > T.a[ln])
			return true;
		else
			return false;
	}
	else
		return false;
}
bool CBigNum::operator >(const int & t) const    //������һ��int���͵ı����Ĵ�С�Ƚ�
{
	CBigNum b(t);
	return *this>b;
}
string CBigNum::ToString()
{
	ostringstream ostr;
	int i;
	ostr << a[len - 1];
	for (i = len - 2; i >= 0; i--)
	{
		ostr.width(DLEN);
		ostr.fill('0');
		ostr << a[i];
	}
	return ostr.str();
}