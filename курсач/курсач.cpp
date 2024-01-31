#include <stdio.h>
#include <math.h>
#include <conio.h>
#include <locale.h>
#include <Windows.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

int maxiter = 100; // максимальное количество итераций
double eps = 1e-15; // величина требуемой относительной невзяки
int n = 0; // количество конечных элементов
int nv = 0; // количество узлов в сетке
int nr = 0; // количество ребер по r
int nz = 0; // количество ребер по z
int ne = 0; // количество краевых условий
int nt = 0; // количество временных слоев
double r_min = 0; // наименьшая координата по r
double r_max = 0; // наибольшая координата по r
double z_min = 0; // наименьшая координата по z
double z_max = 0; // наибольшая координата по z
double times_max = 0; // наибольшее значение по времени
double hr = 0; // шаг по r
double hz = 0; // шаг по z
double ht = 0; // шаг по t
int n1 = 0; // количество 1 краевых
int n2 = 0; // количество 2 краевых
int n3 = 0; // количество 3 краевых

vector<vector<double>> G(4); // матрица жесткости
vector<vector<double>> M(4); // матрица массы
vector<vector<double>> localA(4); // локальная матрица
vector<double> localb(4); // локальный вектор
vector<vector<double>> nodes; // матрица координат узлов
vector<vector<int>> num_elem; // матрица глобальных номеров узлов
vector<vector<int>> edge1; // 1 краевые условия (узел1, узел2, формула)
vector<vector<int>> edge2; // 2 краевые условия (узел1, узел2, формула)
vector<vector<int>> edge3; // 3 краевые условия (узел1, узел2, формула)
vector<double> localb2(2); // для учета 2 и 3 краевых
vector<vector<double>> locala2(2); // для учета 3 краевых

vector<int> ig; // строки
vector<int> jg; // столбцы
vector<double> al; // нижний треугольник 
vector<double> au; // верхний треугольник
vector<double> di; // диагональ
vector<double> b; // вектор правой части
// для ЛОС
vector<double> L; // нижний треугольник 
vector<double> U; // верхний треугольник
vector<double> D; // диагональ
vector<double> x0;
vector<double> z;
vector<double> r;
vector<double> y;
vector<double> p;
vector<double> t;

vector<double> times; // временные слои
vector<double> u0; // решение на j-3 временном слое
vector<double> u1; // решение на j-2 временном слое
vector<double> u2; // решение на j-1 временном слое

// функция f (правая часть уравнения)
double Function(double r, double z, double t)
{
	return 1;
}

double Sigma(int num, double r, double z)
{
	switch (num)
	{
	case 0:
		return 1;
	default:
		break;
	}
}

//коэффициент диффузии
double Lambda(double r, double z)
{
	return 1;
}

//первые краевые условия
double ug(int num, float r, float z, double t)
{
	switch (num)
	{
	case 0:
		return t;
	default:
		break;
	}

}

// 2 краевые условия
double Teta(int num, float r, float z, double t)
{
	switch (num)
	{
	case 0:
		return 0;
	case 1:
		return 0;
	case 2:
		return 0;
	default:
		break;
	}
}

// константа для 3 краевых условий (коэффициент теплоотдачи)
double Betta(int num)
{
	switch (num)
	{
	case 0:
		return 1;
	case 1:
		return 1;
	default:
		break;
	}
}

// функция для 3 краевых условий 
double u_betta(int num, float r, float z, double t)
{
	switch (num)
	{
	case 0:
		return 1;
	case 1:
		return 1;
	default:
		break;
	}
}

// подсчет матрицы жесткости
void GetLocalG(double rp, double zs, double hr, double hz)
{
	double lambda1 = Lambda(rp, zs);
	double lambda2 = Lambda(rp + hr, zs);
	double lambda3 = Lambda(rp, zs + hz);
	double lambda4 = Lambda(rp + hr, zs + hz);

	double a1 = (lambda1 * hz * rp) / (6 * hr);
	double a2 = (lambda2 * hz) / 12;
	double a3 = (lambda3 * hr * rp) / (6 * hz);
	double a4 = (lambda4 * hr * hr) / (12 * hz);

	G[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + a4;
	G[0][1] = -2 * a1 - 2 * a2 + a3 + a4;
	G[0][2] = a1 + a2 - 2 * a3 - a4;
	G[0][3] = - a1 - a2 - a3 - a4;

	G[1][0] = -2 * a1 - 2 * a2 + a3 + a4;
	G[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
	G[1][2] = -a1 - a2 - a3 - a4;
	G[1][3] = a1 + a2 - 2 * a3 - 3 * a4;

	G[2][0] = a1 + a2 - 2 * a3 - a4;
	G[2][1] = -a1 - a2 - a3 - a4;
	G[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + a4;
	G[2][3] = -2 * a1 - 2 * a2 + a3 + a4;

	G[3][0] = -a1 - a2 - a3 - a4;
	G[3][1] = a1 + a2 - 2 * a3 - 3 * a4;
	G[3][2] = -2 * a1 - 2 * a2 + a3 + a4;
	G[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
}

// подсчет матрицы массы
void GetLocalM(double rp, double zs, double hr, double hz)
{

	M[0][0] = hz * ((((rp * hr) / 4 + (hr * hr) / 20)) / 4 + 
		(((rp * hr) / 12 + (hr * hr) / 30)) / 4 +
		(((rp * hr) / 4 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12);

	M[0][1] = hz * ((((rp * hr) / 12 + (hr * hr) / 30)) / 4 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 4 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12);

	M[0][2] = hz * ((((rp * hr) / 4 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12);

	M[0][3] = hz * ((((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12);

	M[1][0] = M[0][1];

	M[1][1] = hz * ((((rp * hr) / 12 + (hr * hr) / 20)) / 4 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 4 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 12);

	M[1][2] = hz * ((((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12);

	M[1][3] = hz * ((((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 12);

	M[2][0] = M[0][2];

	M[2][1] = M[1][2];

	M[2][2] = hz * ((((rp * hr) / 4 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 20)) / 4 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 4);

	M[2][3] = hz * ((((rp * hr) / 12 + (hr * hr) / 30)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 30)) / 4 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 4);

	M[3][0] = M[0][3];

	M[3][1] = M[1][3];

	M[3][2] = M[2][3];

	M[3][3] = hz * ((((rp * hr) / 12 + (hr * hr) / 20)) / 12 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 12 +
		(((rp * hr) / 12 + (hr * hr) / 20)) / 4 +
		(((rp * hr) / 4 + (hr * hr) / 5)) / 4);
}

// подсчет локального вектора правой части
void GetLocalB(int num, double rp, double zs, double hr, double hz, int numt)
{
	double t = times[numt];
	double f1 = Function(rp, zs, t);
	double f2 = Function(rp + hr, zs, t);
	double f3 = Function(rp, zs + hz, t);
	double f4 = Function(rp + hr, zs + hz, t);
	// добавки от изначальной задачи
	localb[0] = f1 * ((hz * rp * hr) / 9 + (hz * hr * hr) / 36) +
		f2 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 36) +
		f3 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 72) +
		f4 * ((hz * rp * hr) / 36 + (hz * hr * hr) / 72);

	localb[1] = f1 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 36) +
		f2 * ((hz * rp * hr) / 9 + (hz * hr * hr) / 12) +
		f3 * ((hz * rp * hr) / 36 + (hz * hr * hr) / 72) +
		f4 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 24);

	localb[2] = f1 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 72) +
		f2 * ((hz * rp * hr) / 36 + (hz * hr * hr) / 72) +
		f3 * ((hz * rp * hr) / 9 + (hz * hr * hr) / 36) +
		f4 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 36);

	localb[3] = f1 * ((hz * rp * hr) / 36 + (hz * hr * hr) / 72) +
		f2 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 24) +
		f3 * ((hz * rp * hr) / 18 + (hz * hr * hr) / 36) +
		f4 * ((hz * rp * hr) / 9 + (hz * hr * hr) / 12);
	// коэффициенты по времени
	double t0, t1, t2;
	double delta_t0 = (times[numt - 3] - times[numt - 2]) * (times[numt - 3] - times[numt - 1]) * (times[numt - 3] - times[numt]);
	double delta_t1 = (times[numt - 2] - times[numt - 3]) * (times[numt - 2] - times[numt - 1]) * (times[numt - 2] - times[numt]);
	double delta_t2 = (times[numt - 1] - times[numt - 3]) * (times[numt - 1] - times[numt - 2]) * (times[numt - 1] - times[numt]);

	t0 = (t * t - times[numt - 1] * t - times[numt - 2] * t + times[numt - 2] * times[numt - 1]) / delta_t0;
	t1 = (t * t - times[numt - 1] * t - times[numt - 3] * t + times[numt - 3] * times[numt - 1]) / delta_t1;
	t2 = (t * t - times[numt - 3] * t - times[numt - 2] * t + times[numt - 2] * times[numt - 3]) / delta_t2;
	// используемые узлы
	int a = num_elem[num][0];
	int b = num_elem[num][1];
	int c = num_elem[num][2];
	int d = num_elem[num][3];
	double s = Sigma(0, rp, zs);
	// добавки параболической задачи
	localb[0] -= s * (t0 * (M[0][0] * u0[a] + M[0][1] * u0[b] + M[0][2] * u0[c] + M[0][3] * u0[d])
		+ t1 * (M[0][0] * u1[a] + M[0][1] * u1[b] + M[0][2] * u1[c] + M[0][3] * u1[d])
		+ t2 * (M[0][0] * u2[a] + M[0][1] * u2[b] + M[0][2] * u2[c] + M[0][3] * u2[d]));

	localb[1] -= t0 * (M[1][0] * u0[a] + M[1][1] * u0[b] + M[1][2] * u0[c] + M[1][3] * u0[d])
		+ t1 * (M[1][0] * u1[a] + M[1][1] * u1[b] + M[1][2] * u1[c] + M[1][3] * u1[d])
		+ t2 * (M[1][0] * u2[a] + M[1][1] * u2[b] + M[1][2] * u2[c] + M[1][3] * u2[d]);

	localb[2] -= s * (t0 * (M[2][0] * u0[a] + M[2][1] * u0[b] + M[2][2] * u0[c] + M[2][3] * u0[d])
		+ t1 * (M[2][0] * u1[a] + M[2][1] * u1[b] + M[2][2] * u1[c] + M[2][3] * u1[d])
		+ t2 * (M[2][0] * u2[a] + M[2][1] * u2[b] + M[2][2] * u2[c] + M[2][3] * u2[d]));

	localb[3] -= s * (t0 * (M[3][0] * u0[a] + M[3][1] * u0[b] + M[3][2] * u0[c] + M[3][3] * u0[d])
		+ t1 * (M[3][0] * u1[a] + M[3][1] * u1[b] + M[3][2] * u1[c] + M[3][3] * u1[d])
		+ t2 * (M[3][0] * u2[a] + M[3][1] * u2[b] + M[3][2] * u2[c] + M[3][3] * u2[d]));
}
// cетка
void GridGeneration()
{
	ifstream input("matrix.txt");
	ofstream node("node.txt");
	ofstream elem("elem.txt");
	input >> r_min >> z_min >> r_max >> z_max >> nr >> nz >> times_max >> nt;
	hr = (r_max - r_min) / nr;
	hz = (z_max - z_min) / nz;
	ht = times_max / nt;

	int index = 0;
	nv = (nr + 1) * (nz + 1); // количество узлов
	node << nv << endl;
	nodes.resize(nv);
	for(double z = z_min; z <= z_max; z += hz)
		for (double r = r_min; r <= r_max; r += hr)
		{
			nodes[index].resize(2);
			nodes[index][0] = r;
			nodes[index][1] = z;
			node << index << " " << nodes[index][0] << " " << nodes[index][1] << endl;
			index++;
		}
	n = nr * nz;
	elem << n << endl; // количество областей
	num_elem.resize(n);
	for (int j = 0; j < nz; j++)
		for (int i = 0; i < nr; i++)
		{
			num_elem[nr * j + i].resize(4);
			num_elem[nr * j + i][0] = i + j * (nr + 1); // 1 вершина
			num_elem[nr * j + i][1] = i + j * (nr + 1) + 1; // 2 вершина
			num_elem[nr * j + i][2] = i + (j + 1) * (nr + 1); // 3 вершина
			num_elem[nr * j + i][3] = i + (j + 1) * (nr + 1) + 1; // 4 вершина
			elem << nr * j + i << " " << num_elem[nr * j + i][0] << " " << num_elem[nr * j + i][1] << " " << num_elem[nr * j + i][2] << " " << num_elem[nr * j + i][3] << endl;
		}

	times.resize(nt + 1);
	times[0] = 0;
	for (int i = 1; i < nt + 1; i++)
		times[i] = times[i - 1] + ht;
}

// портрет матрицы
void MatrixPortrait()
{
	vector<vector<int>> list(nv); // вспомогательный массив по количеству узлов
	vector<int> tmp(4); // вспомогательный вектор вершин одного элемента
	for (int i = 0; i < nv; i++)
		list[i];
	int number = 0; // количество элементов для массива ig
	for (int i = 0; i < n; i++) //идем по конечным элементам
	{
		for (int j = 0; j < 4; j++)
			tmp[j] = num_elem[i][j];
		reverse(tmp.begin(), tmp.end()); //сортируем по убывaнию, т.к. портрет симметричный (в нашем случае переворачиваем т.к. все изначально отсортировано)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				int flag = 1;
				for (int p = 0; p < list[tmp[j]].size() && flag; p++)
					if (list[tmp[j]][p] == tmp[k]) flag = 0; // если элемент уже есть в списке смежности
				if (flag)
				{
					list[tmp[j]].push_back(tmp[k]); // если в списке еще нет такого элемента, то добавляем
					number++; // увеличиваем количество для jg
				}
			}
	}
	for (int i = 0; i < nv; i++)
		sort(list[i].begin(), list[i].end()); //сортируем по возрастанию 

	ig.resize(nv + 1);
	
	ig[0] = 0;
	ig[1] = 0;
	for (int i = 2; i < nv + 1; i++)
		ig[i] = ig[i - 1] + list[i - 1].size();

	for (int i = 0; i < nv; i++)
		for (int j = 0; j < list[i].size(); j++)
		jg.push_back(list[i][j]);

	//выделение памяти для работы с векторами
	al.resize(number);
	au.resize(number);
	di.resize(nv);
	b.resize(nv);
	//ЛОС
	L.resize(number);
	U.resize(number);
	D.resize(nv);
	x0.resize(nv);
	y.resize(nv);
	z.resize(nv);
	r.resize(nv);
	p.resize(nv);
	t.resize(nv);
	u0.resize(nv);
	u1.resize(nv);
	u2.resize(nv);

	locala2[0].resize(2);
	locala2[1].resize(2);

	for (int i = 0; i < 4; i++)
	{
		G[i].resize(4);
		M[i].resize(4);
		localA[i].resize(4);
	}
}

// подсчет локальной матрицы
void LocalMatrix(int num, int numt)
{
	//вершины прямоугольника
	int a = num_elem[num][0];
	int b = num_elem[num][1];
	int c = num_elem[num][2];
	int d = num_elem[num][3];
	double rp, zs, hr, hz;
	rp = nodes[a][0];
	zs = nodes[a][1];
	hr = nodes[b][0] - nodes[a][0];
	hz = nodes[c][1] - nodes[b][1];

	GetLocalG(rp, zs, hr, hz);
	GetLocalM(rp, zs, hr, hz);
	GetLocalB(num, rp, zs, hr, hz, numt);
	// коэффициент от параболической задачи
	double delta_t = (times[numt] - times[numt - 3]) * (times[numt] - times[numt - 2]) * (times[numt] - times[numt - 1]);
	double t3 = (3 * times[numt] * times[numt] - 2 * times[numt] * (times[numt - 1] + times[numt - 2] + times[numt - 3]) + times[numt - 1] * times[numt - 2] + times[numt - 1] * times[numt - 3] + times[numt - 2] * times[numt - 3]) / delta_t;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			localA[i][j] = G[i][j] + t3 * Sigma(0, rp, zs) * M[i][j];
		}
}

// добавление локальной матрицы в глобальную
void AddInGlobal(int num)
{
	// диагональ
	for (int i = 0; i < 4; i++)
		di[num_elem[num][i]] += localA[i][i];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < i; j++)
		{
			int igg = num_elem[num][i];
			int jgg = num_elem[num][j];
			if (igg < jgg)
			{
				igg = jgg;
				jgg = num_elem[num][i];
			}
			int index = ig[igg];
			int flag = 1;
			for (; index < ig[igg + 1] && flag; index++)
				if (jg[index] == jgg) flag = 0;
			index--;
			al[index] += localA[i][j];
			au[index] += localA[j][i];
		}
	}
	// вектор правой части
	for (int i = 0; i < 4; i++)
		b[num_elem[num][i]] += localb[i];
}

// чтение краевых из файлов
void EdgeInput()
{
	ifstream e1("e1.txt");
	e1 >> n1;
	edge1.resize(n1);
	for (int i = 0; i < n1; i++)
	{
		edge1[i].resize(3);
		e1 >> edge1[i][0] >> edge1[i][1] >> edge1[i][2];
	}

	ifstream e2("e2.txt");
	e2 >> n2;
	edge2.resize(n2);
	for (int i = 0; i < n2; i++)
	{
		edge2[i].resize(3);
		e2 >> edge2[i][0] >> edge2[i][1] >> edge2[i][2];
	}

	ifstream e3("e3.txt");
	e3 >> n3;
	edge3.resize(n3);
	for (int i = 0; i < n3; i++)
	{
		edge3[i].resize(3);
		e3 >> edge3[i][0] >> edge3[i][1] >> edge3[i][2];
	}

}

// учет вторых краевых условий
void Second(int num, int numt)
{
	int i = edge2[num][0]; // 1 узел
	int j = edge2[num][1]; // 2 узел
	int formula = edge2[num][2]; // формула
	double t = times[numt];
	// параллельно z
	if (nodes[i][0] == nodes[j][0])
	{
		double rp = nodes[i][0];
		double h = nodes[j][1] - nodes[i][1];
		double t1 = Teta(formula, nodes[i][0], nodes[i][1], t);
		double t2 = Teta(formula, nodes[j][0], nodes[j][1], t);
		localb2[0] = rp * h * (2 * t1 + t2) / 6;
		localb2[1] = rp * h * (t1 + 2 * t2) / 6;
	}
	// параллельно r
	else
	{
		double rp = nodes[i][0];
		double h = nodes[j][0] - nodes[i][0];
		double t1 = Teta(formula, nodes[i][0], nodes[i][1], t);
		double t2 = Teta(formula, nodes[j][0], nodes[j][1], t);

		localb2[0] = rp * h * (2 * t1 + t2) / 6 + h * h * (t1 + t2) / 12;
		localb2[1] = rp * h * (t1 + 2 * t2) / 6 + h * h * (t1 + 3 * t2) / 12;
	}
	b[i] += localb2[0];
	b[j] += localb2[1];
}

// учет третьих краевых условий
void Third(int num, int numt)
{
	int i = edge3[num][0]; // 1 узел
	int j = edge3[num][1]; // 2 узел
	int formula = edge3[num][2]; // формула
	double t = times[numt];
	// параллельно z
	if (nodes[i][0] == nodes[j][0])
	{
		double rp = nodes[i][0];
		double h = nodes[j][1] - nodes[i][1];
		double u1 = u_betta(formula, nodes[i][0], nodes[i][1], t);
		double u2 = u_betta(formula, nodes[j][0], nodes[j][1], t);

		locala2[0][0] = Betta(formula) * rp * h / 3;
		locala2[1][0] = Betta(formula) * rp * h / 6;
		locala2[0][1] = Betta(formula) * rp * h / 6;
		locala2[1][1] = Betta(formula) * rp * h / 3;

		localb2[0] = u1 * locala2[0][0] + u2 * locala2[0][1];
		localb2[0] = u1 * locala2[0][0] + u2 * locala2[1][1];

	}
	// параллельно r
	else
	{
		double rp = nodes[i][0];
		double h = nodes[j][0] - nodes[i][0];
		double u1 = u_betta(formula, nodes[i][0], nodes[i][1], t);
		double u2 = u_betta(formula, nodes[j][0], nodes[j][1], t);

		locala2[0][0] = Betta(formula) * rp * h / 3 + Betta(formula) * h * h / 12;
		locala2[0][1] = Betta(formula) * rp * h / 6 + Betta(formula) * h * h / 12;
		locala2[1][0] = Betta(formula) * rp * h / 6 + Betta(formula) * h * h / 12;
		locala2[1][1] = Betta(formula) * rp * h / 3 + Betta(formula) * h * h / 4;

		localb2[0] = u1 * locala2[0][0] + u2 * locala2[0][1];
		localb2[0] = u1 * locala2[0][0] + u2 * locala2[1][1];
	}
	b[i] += localb2[0];
	b[j] += localb2[1];

	di[i] += locala2[0][0];
	di[j] += locala2[1][1];

	int ia, ja;
	ia = j;
	ja = i;
	int index = ig[ia];
	int flag = 1;
	for (; index < ig[ia + 1] && flag; index++)
		if (jg[index] == ja) flag = 0;
	index--;
	al[index] += locala2[0][1];
	au[index] += locala2[1][0];
}

// учет первых краевых условий
void First(int num, int numt)
{
	int i = edge1[num][0]; // 1 узел
	int j = edge1[num][1]; // 2 узел
	int formula = edge1[num][2]; // формула
	double t = times[numt];
	b[i] = ug(formula, nodes[i][0], nodes[i][1], t);
	b[j] = ug(formula, nodes[j][0], nodes[j][1], t);

	di[i] = 1;
	di[j] = 1;

	for (int k = ig[i]; k < ig[i + 1]; k++)
		al[k] = 0;

	for (int k = ig[i + 1]; k < ig[nv]; k++)
		if (jg[k] == i) au[k] = 0;

	for (int k = ig[j]; k < ig[j + 1]; k++)
		al[k] = 0;

	for (int k = ig[j + 1]; k < ig[nv]; k++)
		if (jg[k] == j) au[k] = 0;
}

void LUsq()
{
	for (int i = 0; i < ig[nv]; i++)
	{
		L[i] = al[i];
		U[i] = au[i];
	}
	for (int i = 0; i < nv; i++)
		D[i] = di[i];

	for (int i = 0; i < nv; i++)
	{
		int ia0 = ig[i];
		int ia1 = ig[i + 1];
		double sd = 0;
		for (int k = ia0; k < ia1; k++)
		{
			double s1 = 0;
			double s2 = 0;
			int j = jg[k]; // номер столбца(строки) элемента al(au)
			int ja0 = ig[j]; // номер, с которого начинается элементы столбца(строки)
			int ja1 = ig[j + 1]; // номер, с которого начинаются элементы нового столбца(строки)
			int ki = ia0; 
			int kj = ja0;
			while (ki < k && kj < ja1)
				if (jg[ki] == jg[kj])
				{
					s1 += L[ki] * U[kj];
					s2 += L[kj] * U[ki];
					ki++;
					kj++;
				}
				else if (jg[ki] > jg[kj])
					kj++;
				else
					ki++;
			L[k] = L[k] - s1;
			U[k] = U[k] - s2;
			L[k] = L[k] / D[jg[k]];
			U[k] = U[k] / D[jg[k]];
			sd += U[k] * L[k];
		}
		D[i] -= sd;
		D[i] = sqrt(D[i]);
	}
}

void LSolve(vector<double>& x1, vector<double> x2)
{
	x1[0] = x2[0] / D[0];
	for (int i = 1; i < nv; i++)
	{
		double s = 0;
		for (int j = ig[i]; j < ig[i + 1]; j++)
			s += L[j] * x1[jg[j]];
		x1[i] = x2[i] - s;
		x1[i] = x1[i] / D[i];
	}
}

void USolve(vector<double>& x1, vector<double> x2)
{
	for (int i = 0; i < nv; i++)
		x1[i] = x2[i];
	x1[nv - 1] = x1[nv - 1] / D[nv - 1];
	for (int i = nv - 1; i >= 1; i--)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
			x1[jg[j]] -= U[j] * x1[i];
		x1[i - 1] = x1[i - 1] / D[i - 1];
	}
}

// умножения вектора на матрицу A
void multiplication_matrix_on_vector(vector<double> a, vector<double>& b)
{
	for (int i = 0; i < nv; i++)
		b[i] = 0;
	for (int i = 0; i < nv; i++)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			b[i] += al[j] * a[jg[j]];
			b[jg[j]] += au[j] * a[i];
		}
		b[i] += di[i] * a[i];
	}
}

// умножение двух векторов
double vectors_multiplication(vector<double> v1, vector<double> v2)
{
	double s = 0;
	for (int i = 0; i < nv; i++)
		s += v1[i] * v2[i];
	return s;
}

// норма вектора
double norma(vector<double> vector)
{
	double s = 0;
	for (int i = 0; i < nv; i++)
		s += vector[i] * vector[i];
	return(sqrt(s));
}
// d = a + b * c
void summ(vector<double> a, double b, vector<double> c, vector<double>& d)
{
	for (int i = 0; i < nv; i++)
		d[i] = a[i] + b * c[i];
}

void LOS()
{
	double alpha, betta;
	double pk_1_rk_1, pk_1_pk_1;
	// нулевое начальное приближение
	for (int i = 0; i < nv; i++)
		x0[i] = 0;
	multiplication_matrix_on_vector(x0, y); // Ax0
	summ(b, -1, y, r); // r0 = f - Ax0
	LSolve(r, r); // r0 = (f - Ax0) / L
	USolve(z, r); // z0 = r0 / U
	multiplication_matrix_on_vector(z, y); // Az0
	LSolve(p, y); // p0 = Az0 / L
	for (int k = 1; k < maxiter; k++)
	{
		pk_1_rk_1 = vectors_multiplication(p, r); // (p_(k-1),r_(k-1))
		pk_1_pk_1 = vectors_multiplication(p, p); // (p_(k-1),p_(k-1))
		alpha = pk_1_rk_1 / pk_1_pk_1; // alpha_k = (p_(k-1),r_(k-1)) / (p_(k-1),p_(k-1))
		summ(x0, alpha, z, x0); // x_k = x_(k-1) + alpha_k * z_(k-1)
		summ(r, -alpha, p, r); // r_k = r_(k-1) - alpha_k * p_(k-1)

		USolve(t, r); // t = r_k/ U
		multiplication_matrix_on_vector(t, y); // y = A * r_k / U
		LSolve(t, y); // t = A * r_k / UL
		betta = -vectors_multiplication(p, t) / pk_1_pk_1; // betta_k = (p_k-1,L_1*A*U_1*r_k) / (p_(k-1),p_(k-1))
		USolve(y, r); // y = r_k / U
		summ(y, betta, z, z);//	z_k = r_k / U + betta_k * z_(k-1)
		summ(t, betta, p, p);//p_k = L_1*A*U_1* r_k + betta_k * p_(k-1)

		if (vectors_multiplication(r, r) < eps) // (r_k, r_k) < e
		{
			return;
		}
	}
}

// объединение всех функций 
void Method()
{
	GridGeneration();
	MatrixPortrait();
	// первые три слоя 
	for (int i = 0; i < nv; i++)
	{
		u0[i] = ug(0, nodes[i][0], nodes[i][1], 0);
		u1[i] = ug(0, nodes[i][0], nodes[i][1], ht);
		u2[i] = ug(0, nodes[i][0], nodes[i][1], 2 * ht);
	}
	cout << "t = 0" << endl;
	for (int i = 0; i < nv; i++)
		cout << u0[i] << endl;
	cout << endl;

	cout << "t = " << ht << endl;
	for (int i = 0; i < nv; i++)
		cout << u1[i] << endl;
	cout << endl;

	cout << "t = " << 2 * ht << endl;
	for (int i = 0; i < nv; i++)
		cout << u2[i] << endl;
	cout << endl;
	// остальные слои
	for (int j = 3; j < nt + 1; j++)
	{
		for (int i = 0; i < n; i++)
		{
			LocalMatrix(i, j);
			AddInGlobal(i);
		}
		EdgeInput();
		for (int i = 0; i < n2; i++)
			Second(i, j);
		for (int i = 0; i < n3; i++)
			Third(i, j);
		for (int i = 0; i < n1; i++)
			First(i, j);
		LUsq();
		LOS();
		cout << "t = " << times[j] << endl;
		//cout.precision(15);
		for (int i = 0; i < nv; i++)
			cout << x0[i] << endl;
		cout << endl;
		// меняем значения для нового слоя 
		for (int i = 0; i < nv; i++)
		{
			u0[i] = u1[i];
			u1[i] = u2[i];
			u2[i] = x0[i];
			di[i] = 0;
			b[i] = 0;
		}
		for (int i = 0; i < al.size(); i++)
		{
			al[i] = 0;
			au[i] = 0;
		}
	}

}


int main()
{
	Method();
}
