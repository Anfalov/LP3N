#pragma once
#include<iostream>
#include<math.h>
#include<vector>
#include<algorithm>
#define Simplex
using namespace std;

int n, m;
vector<double> f, b, numCb, delta;
vector<vector<double>> a;

//����������
void Read()
{
	cout << "������� ����� ���������� � �������:" << endl;
	cin >> n >> m;

	f.resize(n + m + 1);
	b.resize(m);
	numCb.resize(m);
	delta.resize(n + m);
	a.resize(m, vector<double>(n + m));

	cout << "������� ������������ ������� ������� (������� ��������� ����):" << endl;
	for (int i = 0; i < n; i++)
		cin >> f[i];
	//��������� ������������ ��� �������������� ���������� ������
	for (int i = 0; i < m; i++)
		f[n + i] = 0;
	//��������� ��������� ����
	cin >> f[n + m];

	cout << "������� ������������ �������:" << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cin >> a[i][j];

		//������ ������� ���������� ���������� � ����������
		cin >> b[i];
		//��������� ������������ ��� �������������� ���������� ������
		for (int j = n; j < n + m; j++)
			a[i][j] = 0;
		//����������� ��� �������������� ���������� ������� ������ = 1
		a[i][n + i]++;
		numCb[i] = n + i;
	}
}
//�������������� ��������-�������
void GaussStep(int row, int column)
{
	//��������� � ������ ������
	numCb[row] = column;
	//������ �������������� ��������-�������
	//����� ������ �� ������� �������
	double tmp = a[row][column];
	b[row] /= tmp;
	for (int i = 0; i < n + m; i++)
		a[row][i] /= tmp;

	//�������� ��������� ������ �� ������ �����, ��� ��������� ��������� ��������� � ������� �������
	for (int i = 0; i < m; i++)
	{
		//���� ��� ���� ���������� ������, �� � �� �������
		if (i == row)
			continue;

		tmp = a[i][column];
		for (int j = 0; j < n + m; j++)
			a[i][j] -= a[row][j] * tmp;
		b[i] -= b[row] * tmp;
	}
}

//�������� ����, ���������� true, ���� ������� ������� � false, � �������� ������
bool DualSimplex()
{
	//�������� ����
	do
	{
		int row = -1, deltaNeg = -1, column = -1;
		//���� ����������� ������������� b
		for (int i = 0; i < m; i++)
			if (b[i] < 0 && (row == -1 || b[i] < b[row]))
				row = i;

		//������� ������
		for (int i = 0; i < n + m; i++)
		{
			//���� ��������� � ������, �� ������ ����� 0
			if (find(numCb.begin(), numCb.end(), i) != numCb.end())
			{
				delta[i] = 0;
				continue;
			}

			//������� ����������� -Ci � ������
			delta[i] = -f[i];
			//���������� � ������ Ai*Cb
			for (int j = 0; j < m; j++)
				delta[i] += a[j][i] * f[numCb[j]];

			//���� ����������� ������������� ������
			if (delta[i] < 0 && (deltaNeg == -1 || delta[i] < delta[deltaNeg]))
				deltaNeg = i;
		}

		//���� ������������� b ���
		if (row == -1)
		{
			//���� �� ������������ �������� �����, �� �� ����� �������
#ifndef Simplex
			return true;
#endif
			//����� ���������� ������ �������� �������
#ifdef Simplex
			//���� ��� ������������� ������, �� �� ����� �������
			if (deltaNeg == -1)
				return true;

			//�����, ������� ������� ���������� �� ����������� �� ������
			column = deltaNeg;

			//� ��������� �������, ���� ������ � ������������� ��������� Ai ��� �������� bi/Ai ����������
			for (int i = 0; i < m; i++)
				if (a[i][column]>0 && (row == -1 || b[i] / a[i][column] < b[row] / a[row][column]))
					row = i;

			//���� ��� �� ������ �������������� ��������, �� ������� ���
			if (row == -1)
				return false;
#endif
		}
		else
		{
			//��� ��������� ������, ���� ������� � ������������� ��������� Ai ��� �������� -delta/Ai ����������
			for (int i = 0; i < n + m; i++)
				if (a[row][i] < 0 && (column == -1 || -delta[i] / a[row][i] < -delta[column] / a[row][column]))
					column = i;

			//���� ��� �� ������ �������������� ��������, �� ������� ���
			if (column == -1)
				return false;
		}
		//�������� � ������ ������
		GaussStep(row, column);
	} while (true);
}

//����� ���������� �� �������
void Write(bool isSolution)
{
	//���� ������� �������, ������� ���
	if (isSolution)
	{
		//������� � ������� fmax, �������� �������� ���������� �����
		double fmax = f[n + m];
		for (int i = 0; i < m; i++)
			fmax += f[numCb[i]] * b[i];
		cout << "������������ �������� �������: " << fmax << endl;

		//������� ������ ������������ �������
		cout << "������ X���: (";
		for (int i = 0; i < n + m; i++)
		{
			//���������, ����� �� ������� � ����������� ����, ���� ��,
			//�� ���������� ����� ��� ������ �� �������-�������, � ������� ���
			int bas = -1;
			for (int j = 0; j < m; j++)
				if (numCb[j] == i)
				{
					bas = j;
					cout << b[bas];
				}
			//���� ���, �� ������� 0
			if (bas == -1)
				cout << 0;
			if (i != n + m - 1)
				cout << "; ";
			else
				cout << ")" << endl;
		}
	}
	//���� ������� �� �������, ������� ��������� �� ����
	else
		cout << "������� �� �������." << endl;
}

//���������� �������
void DualSimlexMethod()
{
	Read();
	Write(DualSimplex());
}

int main()
{
	setlocale(LC_ALL, "RUSSIAN");
	DualSimlexMethod();
	system("pause");
	return 0;
}

