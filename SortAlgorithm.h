#pragma once

namespace SortAlgorithm
{
	//����(Ĭ������)
	//Ԫ�ؽ���
	template<class T>
	inline void SwapT(T& t1, T& t2)
	{
		T temp = t1;
		t1 = t2;
		t2 = temp;
	}

	//ð���㷨���ȶ�����
	template<typename T>
	void BubbleSort(T* pData, size_t  N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		if (bflag)
		{
			for (size_t i = 0; i < N; i++)
			{
				for (size_t j = 0; j < N - 1 - i; j++)
				{
					if (pData[j] > pData[j + 1])
						SwapT(pData[j], pData[j + 1]);
				}
			}
		}
		else
		{
			for (size_t i = 0; i < N; i++)
			{
				for (size_t j = 0; j < N - 1 - i; j++)
				{
					if (pData[j] < pData[j + 1])
						SwapT(pData[j], pData[j + 1]);
				}
			}
		}
	}

	//�������򣨲��ȶ�����
	template<typename T>
	void QuickSort(T* pData, size_t N, bool bflag, size_t begin, size_t end)
	{
		size_t i, j;

		if (begin < end)
		{
			i = begin + 1;  // ��array[begin]��Ϊ��׼������˴�array[begin+1]��ʼ���׼���Ƚϣ�  
			j = end;        // array[end]����������һλ  

			if (bflag)
			{
				while (i < j)
				{
					if (pData[i] > pData[begin])  // ����Ƚϵ�����Ԫ�ش��ڻ�׼�����򽻻�λ�á�  
					{
						SwapT(pData[i], pData[j]);  // ����������  
						j--;
					}
					else
					{
						i++;  // �����������һλ���������׼���Ƚϡ�  
					}
				}

				/* ����whileѭ����i = j��
				* ��ʱ���鱻�ָ����������  -->  array[begin+1] ~ array[i-1] < array[begin]
				*                           -->  array[i+1] ~ array[end] > array[begin]
				* ���ʱ������array�ֳ��������֣��ٽ�array[i]��array[begin]���бȽϣ�����array[i]��λ�á�
				* ���array[i]��array[begin]���������������ָ�ֵ������Դ����ƣ�ֱ�����i = j�������������˳���
				*/

				if (pData[i] >= pData[begin])  // �������Ҫȡ�ȡ�>=������������Ԫ������ͬ��ֵʱ������ִ���  
				{
					i--;
				}
			}
			else
			{
				while (i < j)
				{
					if (pData[i] < pData[begin])  // ����Ƚϵ�����Ԫ��С�ڻ�׼�����򽻻�λ�á�  
					{
						SwapT(pData[i], pData[j]);  // ����������  
						j--;
					}
					else
					{
						i++;  // �����������һλ���������׼���Ƚϡ�  
					}
				}

				/* ����whileѭ����i = j��
				* ��ʱ���鱻�ָ����������  -->  array[begin+1] ~ array[i-1] < array[begin]
				*                           -->  array[i+1] ~ array[end] > array[begin]
				* ���ʱ������array�ֳ��������֣��ٽ�array[i]��array[begin]���бȽϣ�����array[i]��λ�á�
				* ���array[i]��array[begin]���������������ָ�ֵ������Դ����ƣ�ֱ�����i = j�������������˳���
				*/

				if (pData[i] <= pData[begin])  // �������Ҫȡ�ȡ�>=������������Ԫ������ͬ��ֵʱ������ִ���  
				{
					i--;
				}
			}

			SwapT(pData[begin], pData[i]);  // ����array[i]��array[begin]  

			QuickSort(pData, N, bflag, begin, i);
			QuickSort(pData, N, bflag, j, end);
		}
	}
	template<typename T>
	void QuickSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		QuickSort(pData, N, bflag, 0, N - 1);
	}

	//ѡ�����򣨷��ȶ�����
	template<typename T>
	void SelectionSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		size_t i, j;
		T fValue;
		if (bflag)
		{
			for (i = 0; i < N; i++)
			{
				size_t item = i;
				fValue = pData[item];
				for (j = i + 1; j < N; j++)
				{
					if (fValue>pData[j])
					{
						item = j;
						fValue = pData[item];
					}
				}
				if (item != i)
					SwapT(pData[i], pData[item]);
			}
		}
		else
		{
			for (i = 0; i < N; i++)
			{
				size_t item = i;
				fValue = pData[item];
				for (j = i + 1; j < N; j++)
				{
					if (fValue<pData[j])
					{
						item = j;
						fValue = pData[item];
					}
				}
				if (item != i)
					SwapT(pData[i], pData[item]);
			}
		}
	}

	//�����򣨷��ȶ�����
	template<typename T>
	void MinHeapify(T* pData, size_t N, int element, bool bflag)
	{
		size_t lchild = element * 2 + 1, rchild = lchild + 1;//��������
		if (!bflag)
		{
			while (rchild<N)//�������ڷ�Χ��
			{
				if (pData[element] <= pData[lchild] && pData[element] <= pData[rchild])//���������������С���������
				{
					return;
				}
				if (pData[lchild] <= pData[rchild])//��������С
				{
					SwapT(pData[element], pData[lchild]);//��������ᵽ����
					element = lchild;//ѭ��ʱ��������
				}
				else//����������С
				{
					SwapT(pData[element], pData[rchild]);//ͬ��
					element = rchild;
				}
				lchild = element * 2 + 1;
				rchild = lchild + 1;//���¼�������λ��
			}
			if (lchild<N&&pData[lchild]<pData[element])//ֻ��������������С���Լ�
			{
				SwapT(pData[lchild], pData[element]);
			}
		}
		else
		{
			while (rchild<N)//�������ڷ�Χ��
			{
				if (pData[element] >= pData[lchild] && pData[element] >= pData[rchild])//������������������������
				{
					return;
				}
				if (pData[lchild] >= pData[rchild])//���������
				{
					SwapT(pData[element], pData[lchild]);//��������ᵽ����
					element = lchild;//ѭ��ʱ��������
				}
				else//�����������
				{
					SwapT(pData[element], pData[rchild]);//ͬ��
					element = rchild;
				}
				lchild = element * 2 + 1;
				rchild = lchild + 1;//���¼�������λ��
			}
			if (lchild<N&&pData[lchild]>pData[element])//ֻ�������������������Լ�
			{
				SwapT(pData[lchild], pData[element]);
			}
		}
		return;
	}
	template<typename T>
	void HeapSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		int i;
		for (i = N - 1; i >= 0; i--)//��������ʼ������
		{
			MinHeapify(pData, N, i, bflag);
		}
		while (size>0)//�����
		{
			SwapT(pData[N - 1], pData[0]);//��������С����������ĩ����
			N--;//����С��С
			MinHeapify(pData, N, 0, bflag);//������
		}
		return;
	}

	//���������ȶ�����
	template<typename T>
	void InsertionSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		size_t i, j;
#ifdef _InsertionSort_Method_1
		if (bflag)
		{
			for (i = 1; i < N; i++)
			{
				for (j = 0; j < i; j++)
				{
					if (pData[i] < pData[j])
						break;
				}
				for (size_t m = i; m > j; m--)
					SwapT(pData[m], pData[m - 1]);
			}
		}
		else
		{
			for (i = 1; i < N; i++)
			{
				for (j = 0; j < i; j++)
				{
					if (pData[i] > pData[j])
						break;
				}
				for (size_t m = i; m > j; m--)
					SwapT(pData[m], pData[m - 1]);
			}
		}
#else
		if (bflag)
		{
			for (i = 1; i < N; i++)
			{
				for (j = 0; j < i; j++)
				{
					if (pData[i] < pData[j])
						SwapT(pData[i], pData[j]);
				}
			}
		}
		else
		{
			for (i = 1; i < N; i++)
			{
				for (j = 0; j < i; j++)
				{
					if (pData[i] > pData[j])
						SwapT(pData[i], pData[j]);
				}
			}
		}
#endif
	}

	//ϣ�����򣨲�������Ľ��棬���ȶ�����
	template<typename T>
	void ShellSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		if (bflag)
		{
			for (size_t div = N / 2; div >= 1; div = div / 2)//������div�������ϼ�С
			{
				for (size_t i = 0; i <= div; ++i)//�����div��
				{
					for (size_t j = i; j<N - div; j += div)//��ÿ����в�������
					for (size_t k = j; k<N; k += div)
					if (pData[j]>pData[k])
						SwapT(pData[j], pData[k]);//������������ֵ
				}
			}
		}
		else
		{
			for (size_t div = N / 2; div >= 1; div = div / 2)//������div�������ϼ�С
			{
				for (size_t i = 0; i <= div; ++i)//�����div��
				{
					for (size_t j = i; j<N - div; j += div)//��ÿ����в�������
					for (size_t k = j; k<N; k += div)
					if (pData[j]<pData[k])
						SwapT(pData[j], pData[k]);//������������ֵ
				}
			}
		}

	}

	//�鲢�����ȶ�����
	template<typename T>
	void merge(T* pData, size_t start, size_t end, T* result, bool bflag = true)
	{
		size_t left_length = (end - start + 1) / 2 + 1;//�󲿷����������Ԫ�صĸ���
		size_t left_index = start;
		size_t right_index = start + left_length;
		size_t result_index = start;
		if (bflag)
		{
			while (left_index < start + left_length && right_index < end + 1)
			{
				//�Էֱ��Ѿ��ź�������������������кϲ�
				if (pData[left_index] <= pData[right_index])
					result[result_index++] = pData[left_index++];
				else
					result[result_index++] = pData[right_index++];
			}
			while (left_index < start + left_length)
				result[result_index++] = pData[left_index++];
			while (right_index < end + 1)
				result[result_index++] = pData[right_index++];
		}
		else
		{
			while (left_index < start + left_length && right_index < end + 1)
			{
				//�Էֱ��Ѿ��ź�������������������кϲ�
				if (pData[left_index] >= pData[right_index])
					result[result_index++] = pData[left_index++];
				else
					result[result_index++] = pData[right_index++];
			}
			while (left_index < start + left_length)
				result[result_index++] = pData[left_index++];
			while (right_index < end + 1)
				result[result_index++] = pData[right_index++];
		}

	}
	template<typename T>
	void merge_sort(T* pData, size_t start, size_t end, T* result, bool bflag = true)
	{
		if (1 == end - start)//���������ֻ������Ԫ�أ����������Ԫ�ؽ�������
		{
			if (bflag)
			{
				if (pData[start] > pData[end])
					SwapT(pData[start], pData[end]);
			}
			else
			{
				if (pData[start] < pData[end])
					SwapT(pData[start], pData[end]);
			}
			return;
		}
		else if (0 == end - start)//���ֻ��һ��Ԫ�أ���������
			return;
		else
		{
			//�������������䣬�ֱ�������������������
			merge_sort(pData, start, (end - start + 1) / 2 + start, result, bflag);
			merge_sort(pData, (end - start + 1) / 2 + start + 1, end, result, bflag);
			//��ʼ�鲢�Ѿ��ź����start��end֮�������
			merge(pData, start, end, result, bflag);
			//���������������ݸ��Ƶ�ԭʼ������ȥ
			for (size_t i = start; i <= end; ++i)
				pData[i] = result[i];
		}
	}
	template<typename T>
	void MergeSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		T* pResult = new T[N];
		merge_sort(pData, 0, N - 1, pResult, bflag);
		delete[] pResult;
	}
}