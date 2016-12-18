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

//Считывание
void Read()
{
	cout << "Введите число переменных и условий:" << endl;
	cin >> n >> m;

	f.resize(n + m + 1);
	b.resize(m);
	numCb.resize(m);
	delta.resize(n + m);
	a.resize(m, vector<double>(n + m));

	cout << "Введите коэффициенты целевой функции (включая свободный член):" << endl;
	for (int i = 0; i < n; i++)
		cin >> f[i];
	//Заполняем коэффициенты при дополнительных переменных нулями
	for (int i = 0; i < m; i++)
		f[n + i] = 0;
	//Считываем свободный член
	cin >> f[n + m];

	cout << "Введите коэффициенты условий:" << endl;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
			cin >> a[i][j];

		//Правую сторону неравенств записываем в псевдоплан
		cin >> b[i];
		//Заполняем коэффициенты при дополнительных переменных нулями
		for (int j = n; j < n + m; j++)
			a[i][j] = 0;
		//Коэффициент при дополнительной переменной текущей строки = 1
		a[i][n + i]++;
		numCb[i] = n + i;
	}
}
//Преобразование симплекс-таблицы
void GaussStep(int row, int column)
{
	//Переходим к новому базису
	numCb[row] = column;
	//Делаем преобразование симплекс-таблицы
	//Делим строку на ведущий элемент
	double tmp = a[row][column];
	b[row] /= tmp;
	for (int i = 0; i < n + m; i++)
		a[row][i] /= tmp;

	//Вычитаем изменённую строку из других строк, для обнуления остальных элементов в ведущем столбце
	for (int i = 0; i < m; i++)
	{
		//Если это сама вычитаемая строка, то её не трогаем
		if (i == row)
			continue;

		tmp = a[i][column];
		for (int j = 0; j < n + m; j++)
			a[i][j] -= a[row][j] * tmp;
		b[i] -= b[row] * tmp;
	}
}

//Основной цикл, возвращает true, если решение найдено и false, в обратном случае
bool DualSimplex()
{
	//Основной цикл
	do
	{
		int row = -1, deltaNeg = -1, column = -1;
		//Ищем минимальный отрицательный b
		for (int i = 0; i < m; i++)
			if (b[i] < 0 && (row == -1 || b[i] < b[row]))
				row = i;

		//Считаем оценки
		for (int i = 0; i < n + m; i++)
		{
			//Если относится к базису, то оценка равна 0
			if (find(numCb.begin(), numCb.end(), i) != numCb.end())
			{
				delta[i] = 0;
				continue;
			}

			//Заранее присваиваем -Ci в оценку
			delta[i] = -f[i];
			//Прибавляем к оценке Ai*Cb
			for (int j = 0; j < m; j++)
				delta[i] += a[j][i] * f[numCb[j]];

			//Ищем минимальную отрицательную оценку
			if (delta[i] < 0 && (deltaNeg == -1 || delta[i] < delta[deltaNeg]))
				deltaNeg = i;
		}

		//Если отрицательных b нет
		if (row == -1)
		{
			//Если не использовать симплекс метод, то мы нашли решение
#ifndef Simplex
			return true;
#endif
			//Иначе продолжать искать симплекс методом
#ifdef Simplex
			//Если нет отрицательных оценок, то мы нашли решение
			if (deltaNeg == -1)
				return true;

			//Иначе, текущий столбец выбирается по минимальной из оценок
			column = deltaNeg;

			//В выбранном столбце, ищем строку с положительным элементом Ai для которого bi/Ai минимально
			for (int i = 0; i < m; i++)
				if (a[i][column]>0 && (row == -1 || b[i] / a[i][column] < b[row] / a[row][column]))
					row = i;

			//Если нет ни одного положительного элемента, то решения нет
			if (row == -1)
				return false;
#endif
		}
		else
		{
			//Для выбранной строки, ищем столбец с отрицательным элементом Ai для которого -delta/Ai минимально
			for (int i = 0; i < n + m; i++)
				if (a[row][i] < 0 && (column == -1 || -delta[i] / a[row][i] < -delta[column] / a[row][column]))
					column = i;

			//Если нет ни одного отрицательного элемента, то решения нет
			if (column == -1)
				return false;
		}
		//Переодим к новому базису
		GaussStep(row, column);
	} while (true);
}

//Вывод результата на консоль
void Write(bool isSolution)
{
	//Если решение найдено, выводим его
	if (isSolution)
	{
		//Считаем и выводим fmax, учитывая значение свободного члена
		double fmax = f[n + m];
		for (int i = 0; i < m; i++)
			fmax += f[numCb[i]] * b[i];
		cout << "Максимальное значение функции: " << fmax << endl;

		//Выводим вектор оптимального решения
		cout << "Вектор Xопт: (";
		for (int i = 0; i < n + m; i++)
		{
			//Проверяем, попал ли элемент в оптимальный план, если да,
			//то запоминаем номер его строки из симпекс-таблицы, и выводим его
			int bas = -1;
			for (int j = 0; j < m; j++)
				if (numCb[j] == i)
				{
					bas = j;
					cout << b[bas];
				}
			//Если нет, то выводим 0
			if (bas == -1)
				cout << 0;
			if (i != n + m - 1)
				cout << "; ";
			else
				cout << ")" << endl;
		}
	}
	//Если решение не найдено, выводим сообщение об этом
	else
		cout << "Решение не найдено." << endl;
}

//Вызываемая функция
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

