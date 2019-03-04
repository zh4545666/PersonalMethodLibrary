#pragma once

namespace SortAlgorithm
{
	//排序(默认升序)
	//元素交换
	template<class T>
	inline void SwapT(T& t1, T& t2)
	{
		T temp = t1;
		t1 = t2;
		t2 = temp;
	}

	//冒泡算法（稳定排序）
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

	//快速排序（不稳定排序）
	template<typename T>
	void QuickSort(T* pData, size_t N, bool bflag, size_t begin, size_t end)
	{
		size_t i, j;

		if (begin < end)
		{
			i = begin + 1;  // 将array[begin]作为基准数，因此从array[begin+1]开始与基准数比较！  
			j = end;        // array[end]是数组的最后一位  

			if (bflag)
			{
				while (i < j)
				{
					if (pData[i] > pData[begin])  // 如果比较的数组元素大于基准数，则交换位置。  
					{
						SwapT(pData[i], pData[j]);  // 交换两个数  
						j--;
					}
					else
					{
						i++;  // 将数组向后移一位，继续与基准数比较。  
					}
				}

				/* 跳出while循环后，i = j。
				* 此时数组被分割成两个部分  -->  array[begin+1] ~ array[i-1] < array[begin]
				*                           -->  array[i+1] ~ array[end] > array[begin]
				* 这个时候将数组array分成两个部分，再将array[i]与array[begin]进行比较，决定array[i]的位置。
				* 最后将array[i]与array[begin]交换，进行两个分割部分的排序！以此类推，直到最后i = j不满足条件就退出！
				*/

				if (pData[i] >= pData[begin])  // 这里必须要取等“>=”，否则数组元素由相同的值时，会出现错误！  
				{
					i--;
				}
			}
			else
			{
				while (i < j)
				{
					if (pData[i] < pData[begin])  // 如果比较的数组元素小于基准数，则交换位置。  
					{
						SwapT(pData[i], pData[j]);  // 交换两个数  
						j--;
					}
					else
					{
						i++;  // 将数组向后移一位，继续与基准数比较。  
					}
				}

				/* 跳出while循环后，i = j。
				* 此时数组被分割成两个部分  -->  array[begin+1] ~ array[i-1] < array[begin]
				*                           -->  array[i+1] ~ array[end] > array[begin]
				* 这个时候将数组array分成两个部分，再将array[i]与array[begin]进行比较，决定array[i]的位置。
				* 最后将array[i]与array[begin]交换，进行两个分割部分的排序！以此类推，直到最后i = j不满足条件就退出！
				*/

				if (pData[i] <= pData[begin])  // 这里必须要取等“>=”，否则数组元素由相同的值时，会出现错误！  
				{
					i--;
				}
			}

			SwapT(pData[begin], pData[i]);  // 交换array[i]与array[begin]  

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

	//选择排序（非稳定排序）
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

	//堆排序（非稳定排序）
	template<typename T>
	void MinHeapify(T* pData, size_t N, int element, bool bflag)
	{
		size_t lchild = element * 2 + 1, rchild = lchild + 1;//左右子树
		if (!bflag)
		{
			while (rchild<N)//子树均在范围内
			{
				if (pData[element] <= pData[lchild] && pData[element] <= pData[rchild])//如果比左右子树都小，完成整理
				{
					return;
				}
				if (pData[lchild] <= pData[rchild])//如果左边最小
				{
					SwapT(pData[element], pData[lchild]);//把左面的提到上面
					element = lchild;//循环时整理子树
				}
				else//否则右面最小
				{
					SwapT(pData[element], pData[rchild]);//同理
					element = rchild;
				}
				lchild = element * 2 + 1;
				rchild = lchild + 1;//重新计算子树位置
			}
			if (lchild<N&&pData[lchild]<pData[element])//只有左子树且子树小于自己
			{
				SwapT(pData[lchild], pData[element]);
			}
		}
		else
		{
			while (rchild<N)//子树均在范围内
			{
				if (pData[element] >= pData[lchild] && pData[element] >= pData[rchild])//如果比左右子树都大，完成整理
				{
					return;
				}
				if (pData[lchild] >= pData[rchild])//如果左边最大
				{
					SwapT(pData[element], pData[lchild]);//把左面的提到上面
					element = lchild;//循环时整理子树
				}
				else//否则右面最大
				{
					SwapT(pData[element], pData[rchild]);//同理
					element = rchild;
				}
				lchild = element * 2 + 1;
				rchild = lchild + 1;//重新计算子树位置
			}
			if (lchild<N&&pData[lchild]>pData[element])//只有左子树且子树大于自己
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
		for (i = N - 1; i >= 0; i--)//从子树开始整理树
		{
			MinHeapify(pData, N, i, bflag);
		}
		while (size>0)//拆除树
		{
			SwapT(pData[N - 1], pData[0]);//将根（最小）与数组最末交换
			N--;//树大小减小
			MinHeapify(pData, N, 0, bflag);//整理树
		}
		return;
	}

	//插入排序（稳定排序）
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

	//希尔排序（插入排序改进版，非稳定排序）
	template<typename T>
	void ShellSort(T* pData, size_t N, bool bflag = true)
	{
		if (N <= 1 || !pData)
			return;
		if (bflag)
		{
			for (size_t div = N / 2; div >= 1; div = div / 2)//定增量div，并不断减小
			{
				for (size_t i = 0; i <= div; ++i)//分组成div组
				{
					for (size_t j = i; j<N - div; j += div)//对每组进行插入排序
					for (size_t k = j; k<N; k += div)
					if (pData[j]>pData[k])
						SwapT(pData[j], pData[k]);//交换两个数的值
				}
			}
		}
		else
		{
			for (size_t div = N / 2; div >= 1; div = div / 2)//定增量div，并不断减小
			{
				for (size_t i = 0; i <= div; ++i)//分组成div组
				{
					for (size_t j = i; j<N - div; j += div)//对每组进行插入排序
					for (size_t k = j; k<N; k += div)
					if (pData[j]<pData[k])
						SwapT(pData[j], pData[k]);//交换两个数的值
				}
			}
		}

	}

	//归并排序（稳定排序）
	template<typename T>
	void merge(T* pData, size_t start, size_t end, T* result, bool bflag = true)
	{
		size_t left_length = (end - start + 1) / 2 + 1;//左部分区间的数据元素的个数
		size_t left_index = start;
		size_t right_index = start + left_length;
		size_t result_index = start;
		if (bflag)
		{
			while (left_index < start + left_length && right_index < end + 1)
			{
				//对分别已经排好序的左区间和右区间进行合并
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
				//对分别已经排好序的左区间和右区间进行合并
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
		if (1 == end - start)//如果区间中只有两个元素，则对这两个元素进行排序
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
		else if (0 == end - start)//如果只有一个元素，则不用排序
			return;
		else
		{
			//继续划分子区间，分别对左右子区间进行排序
			merge_sort(pData, start, (end - start + 1) / 2 + start, result, bflag);
			merge_sort(pData, (end - start + 1) / 2 + start + 1, end, result, bflag);
			//开始归并已经排好序的start到end之间的数据
			merge(pData, start, end, result, bflag);
			//把排序后的区间数据复制到原始数据中去
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