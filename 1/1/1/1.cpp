#include <iostream>
#include <iomanip>               
#include <math.h>
#include <conio.h>

using namespace std;

void inputM(double** WorkMatrix, double** SaveMatrix, int n)
{
	cout << "Заполните матрицу: " << endl;

	for (int i = 0; i < n; i++)
	{
		for (int k = 0; k < n; k++)
		{
			cin >> WorkMatrix[i][k];
			//a[i][k] = rand()%100;
			SaveMatrix[i][k] = WorkMatrix[i][k];
		}
	}
	cout << endl;
}

void inputV(double** WorkMatrix, double** SaveMatrix, int n)
{
	cout << "Введите вектор" << endl;
	for (int i = 0; i < n; i++) //заполнение числового вектора происходит за пределами определения числовой матрицы, поэтому 
	{
		cin >> WorkMatrix[i][n];			  //записываем матричный вектор 
		//a[i][n] = rand() % 100;
		SaveMatrix[i][n] = WorkMatrix[i][n];
	}
}

void outputM(double** WorkMatrix, int n)		//Вывод исходной матрицы
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			cout << setw(12) << WorkMatrix[i][j];
		}
		cout << endl;
	}
}

void Gauss(double** WorkMatrix, double* x, double* x2, double** SaveMatrix, int n, int& checkpoint)
{
	for (int i = 0; i < n; i++)			 //Прямой ход метода Гаусса
	{										 //поиск максимального элемента в первом столбце
		double maxElem = fabs(WorkMatrix[i][i]);		 //определим главный элемент i - ого столбца
		int lineNumb = i;

		for (int j = i; j < n; j++)
		{
			if (fabs(WorkMatrix[j][i]) > maxElem)
			{
				maxElem = fabs(WorkMatrix[j][i]);		//нашли макс элемент в 'i' строке и закунилуи его на верх
				lineNumb = j;
			}
		}

		if (lineNumb != i)					//перестановка строк
		{
			double* buf = WorkMatrix[i];
			WorkMatrix[i] = WorkMatrix[lineNumb];
			WorkMatrix[lineNumb] = buf;
		}

		double leadElem = WorkMatrix[i][i];             //определение ведущего элемента

		for (int z = i; z < n + 1; z++)			//делю строку на первый свобоный элемент
		{
			WorkMatrix[i][z] = WorkMatrix[i][z] / leadElem;
		}

		for (int j = i + 1; j < n; j++)			//зануляю нижнюю строку
		{
			double down_number = WorkMatrix[j][i];
			for (int m = i; m < n + 1; m++)
			{
				WorkMatrix[j][m] = WorkMatrix[j][m] - WorkMatrix[i][m] * down_number;
			}
		}

		cout << "Преобразования матрицы: " << endl;
		outputM(WorkMatrix, n);

	}

	for (int i = n - 1; i > 0; i--)			//обратный ход Гаусса
	{
		for (int j = i - 1; j >= 0; j--)
		{
			WorkMatrix[j][n] = WorkMatrix[j][n] - WorkMatrix[j][i] * WorkMatrix[i][n];
		}
	}

	//невязка и ~b высчитываются один раз
	if (checkpoint == 0)
	{
		for (int i = 0; i < n; i++)		//согласно F=a*x-b
		{									//вектор невязки - разность между точным и приближенным значениями
			x[i] = -SaveMatrix[i][n];
			for (int j = 0; j < n; j++)
			{
				x[i] = (SaveMatrix[i][j] * WorkMatrix[j][n] + x[i]);
			}
			cout << "Значение вектора невязки " << i + 1 << " = " << x[i] << endl;	//вектор невязки
		}

		//A*x=b~
		for (int i = 0; i < n; i++)		//согласно ~b=a*x       x2 = b~
		{
			x2[i] = 0;
			for (int j = 0; j < n; j++)
			{
				x2[i] = x2[i] + SaveMatrix[i][j] * WorkMatrix[j][n];
			}
		}

		//b~~ save in SaveMatrix
		for (int i = 0; i < n; i++) {
			SaveMatrix[i][n] = x2[i];
		}


		int i = 0;

		double norm = fabs(x[i]);				//нахождение максимального значения вектора невязки из массива x[n](оценка нормы)
		for (i = 0; i < n; i++)
		{
			if (norm < fabs(x[i]))
			{
				norm = fabs(x[i]);
			}
		}

		cout << "Норма: " << norm << endl;


		checkpoint++;
	}

}

int main()
{
	cout << "Gauss method....." << endl;
	setlocale(LC_ALL, "rus");

	int n;
	int checkpoint = 0;

	cout << "Введите размерность матрицы: " << endl;  //определение размеров матрицы
	cin >> n;

	double** a = new double* [n];
	for (int i = 0; i < n; i++)
	{
		a[i] = new double[n];
	}

	double** saveA = new double* [n];
	for (int i = 0; i < n; i++)
	{
		saveA[i] = new double[n];
	}
	
	double** b = new double* [n];
	for (int i = 0; i < n; i++)
	{
		b[i] = new double[n + 1];
	}

	double** b2 = new double* [n];					//b~
	for (int i = 0; i < n; i++)
	{
		b[i] = new double[n + 1];
	}

	double* x = new double[n];					  //выделение памяти под столбец-вектор

	inputM(a, b, n);
	inputV(a, b, n);

	//save A
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			saveA[i][j] = a[i][j];
		}
	}

	//place for x~~
	double* x2 = new double[n];

	double* solutions = new double[n];

	cout << "Ваша матрица имеет размерность: " << endl;
	outputM(a, n);

	Gauss(a, x, x2, b, n, checkpoint);

	//Ax~~ = Ax~
	//save A n+1 (changed)
	for (int i = 0; i < n; i++)
	{
		saveA[i][n] = x2[i];
	}
	Gauss(saveA, x, x2, b, n, checkpoint);


	cout << "Окончательно преобразованная матрица: " << endl;

	outputM(a, n);

	//погрешность
	int number = 0, i = 0;
	double max = fabs(saveA[i][n] - a[i][n]);
	for (int i = 0; i < n; i++)
	{
		if (fabs(saveA[i][n] - a[i][n]) > max)
		{
			max = fabs(saveA[i][n] - a[i][n]);
		}
	}

	i = 0;
	double max_down = fabs(a[i][n]);
	for (int i = 0; i < n; i++)
	{
		if (fabs(a[i][n]) > max)
		{
			max_down = fabs(a[i][n]);
		}
	}

	cout << "\nERROR = " << max / max_down << endl << endl;
	

	for (int number = 0; number < n; number++)
	{
		cout << "X" << (number + 1) << " = " << a[number][n] << endl;
	}
	cout << endl;

	return 0;
}