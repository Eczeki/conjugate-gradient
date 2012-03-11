/*
Author		:	Basiak Marcin, Pawe³ Nowak 
University	:	AGH Wydzial In¿ynierii Metali i Inoformatyki Przemyslowej
Topic		:	Conjugate Gradient Method
*/

#include<iostream>
#include<fstream>
#include <iomanip>

void error(const char * p1 , const char * p2="");												//Function, which prints error in console. 
double ** createMatrix(std::size_t const rows , std::size_t const columns);		//Function, which create Matrix from file.
void deleteMatrix(double ** matrix	, std::size_t size );										//Function, which delete Matrix.
bool isSymetric  (double ** matrix	, std::size_t size );										//Function, which check the matrix.
bool isPositiveDefinite(double ** matrix , std::size_t size,double x);							//Function, which check whether matrix is positive definte
double ** transpose(unsigned const rows , unsigned const columns ,  double const ** matrix);	//Function, which transpose matrix
double multiply_Ax_rows_vector(double ** A, unsigned rA, double ** x, unsigned  rx, unsigned cx, unsigned size);
void ConjugateGradientMethod(const char * from , const char * to);
double multiplyVectors(double ** v1, unsigned i1, double ** v2, unsigned i2, unsigned size);
double multiply_pT_A_p(double ** p, double ** A, unsigned k, unsigned size);

void StructuralTest(double **A, double * b, double ** x);

int main(int argc , char * argv[])
{
	if(argc != 3)	{	error("Wrong number of arguments...");		}

	//open input  file stream
	std::ifstream from(argv[1]);
	//open output file stream
	std::ofstream to(argv[2]);

	//**if(!from)	{	error("Can't open input file" , argv[1]);	}
	//**if(!to)		{	error("Can't open output file" , argv[2]);	}

	ConjugateGradientMethod(argv[1], argv[2]);

	from.close();
	to.close();
}

bool isPositiveDefinite(double ** matrix , std::size_t size , double x)
{
	/*	x^T * A * x>0 */
	return true;
}

void ConjugateGradientMethod(const char * from , const char * to)
{
	std::size_t size = 0;																//matrix[size][size]
	//**size << atoi(from);																//the first element in file is size of matrix
	size = 2;

	//**if(size<=0){	error("Wrong size of matrix in file...:",from);	}

	double ** A = 0;
	const std::size_t max_iter = 1000;		//max iteration
	const int min_R = 0.01;
	double ** r = 0;	//mo¿na tego te¿ nie robiæ na macierzy i nie pamiêtaæ wszystkich zmiennnych, tylko wykorzystaæ tablicê dwuelementow¹
	double ** p = 0;	//j.w.
	double ** x = 0;	//j.w.
	A = createMatrix(size, size);													//Create matrix A
	r = createMatrix(size, max_iter);		//rows are a vector, columns are a number of element	
	p = createMatrix(size, max_iter);		//rows are a vector, columns are a number of element
	x = createMatrix(size, max_iter);		//rows are a vector, columns are a number of element
	double * alpha = new double[max_iter];
	double * beta = new double[max_iter];

	double * b =  new double [size];

	//testing
	StructuralTest(A, b, x);

	//r0 = b - Ax0
	for(unsigned i = 0; i < size; i++){
		r[i][0] = b[i] - multiply_Ax_rows_vector(A, i, x, i, 0, size);
	}

	//p0 = r0
	for(unsigned i = 0; i < size; i++)
		p[i][0] = r[i][0];

	//algorithm
	unsigned k;
	bool flag = false;
	for(k = 0; k < max_iter; k++){
		if(abs(multiply_pT_A_p(p, A, k, size)) < 1e-12){ /*break;*/
			error("Too small denominator r!");
		}
		alpha[k] = multiplyVectors(r, k, r, k, size) / multiply_pT_A_p(p, A, k, size);
		for(unsigned i = 0; i < size; i++){
			x[i][k+1] = x[i][k] + alpha[k] * p[i][k];
			r[i][k+1] = r[i][k] - alpha[k] * multiply_Ax_rows_vector(A, i, p, i, k, size);
		}

		for(unsigned i = 0; i < size; i++){
			if( abs(r[i][k+1]) < min_R){ 
				flag = true;
				break;
			}
		}
		if (flag == true) break;

		beta[k] = multiplyVectors(r, k+1, r, k+1, size) / multiplyVectors(r, k, r, k, size);
		for(unsigned i = 0; i < size; i++){
			p[i][k+1] = r[i][k+1] + beta[k] * p[i][k];
		}

		std::cout << k << ".\t";
		for(unsigned i = 0; i < size; i++){
			std::cout << std::setprecision(5) << x[i][k+1] << "\t"; 
		}
		std::cout << "\n"; //**
		//if(k % 10 == 0)	system("pause"); 
	}
	//the result is vector x

	deleteMatrix(A, size);														//delete matrix A
}

void StructuralTest(double **A, double * b, double ** x){
	/*A[0][0] = 2;
	A[0][1] = 2;
	A[0][2] = 1;
	
	A[1][0] = 2;
	A[1][1] = 1;
	A[1][2] = 2;

	A[2][0] = 1;
	A[2][1] = 2;
	A[2][2] = 1;

	b[0] = 9;
	b[1] = 11;
	b[2] = 7;

	x[0][0] = 1000;
	x[1][0] = 100;
	x[2][0] = 1;*/

	/*A[0][0] = 3;
	A[0][1] = 3;
	A[0][2] = 4;
	
	A[1][0] = 1;
	A[1][1] = 2;
	A[1][2] = 3;

	A[2][0] = 3;
	A[2][1] = 1;
	A[2][2] = 2;

	b[0] = 38;
	b[1] = 25;
	b[2] = 20;

	x[0][0] = 1;
	x[1][0] = 1;
	x[2][0] = 1;*/

	A[0][0] = 4;
	A[0][1] = 1;
	
	A[1][0] = 1;
	A[1][1] = 3;

	b[0] = 1;
	b[1] = 2;

	x[0][0] = 2;
	x[1][0] = 1;
}

double multiply_pT_A_p(double ** p, double ** A, unsigned k, unsigned size){
	double  tmp = 0; 
	for(unsigned i=0; i < size; i++){
		tmp += p[i][k] * multiply_Ax_rows_vector(A, i, p, i, k, size);
	}
	return tmp;
}

double multiplyVectors(double ** v1, unsigned i1, double ** v2, unsigned i2, unsigned size){	//multiply two vectors
	double tmp = 0;
	for(unsigned i = 0; i<size; i++){
		tmp += v1[i][i1] * v2[i][i2];
	}
	return tmp;
}

double multiply_Ax_rows_vector(double ** A, unsigned rA, double ** x, unsigned  rx, unsigned cx, unsigned size){
	double tmp = 0;
	for(unsigned i = 0; i < size; i++){
		tmp += A[rA][i] * x[i][cx];
	}
	return tmp;
}


void error(const char * p1 , const char * p2)
{
	std::cerr << p1 << " " << p2 << "\n";
	std::exit(1);
}

double ** createMatrix(std::size_t const rows , std::size_t const columns  )
{
	double ** matrix = new double * [rows];
	for(std::size_t i = 0 ; i < rows ; ++i)
	{
		matrix[i]=new double[columns];
	}
	if(matrix==0){	error("Bad memory allocation...");	}
	else return matrix;
}

void deleteMatrix(double ** matrix , std::size_t size )
{
	for(std::size_t i=0 ; i<size ; ++i)
	{
		delete [] matrix[i];
	}
	delete [] matrix;
}

bool isSymetric(double ** matrix , std::size_t size)
{
	for(std::size_t i = 0 ; i < size ; ++i)
		for(std::size_t j = 0 ; j < size ; ++j)
			if(matrix[i][j] != matrix[j][i])
				return false;
	return true;
}

double ** transpose(unsigned const rows , unsigned const columns ,  double const ** matrix)
{
	double ** _transpose  = 0;
	_transpose = createMatrix(columns,rows);											//This function takes first:matrix,rows,columns but we must 
																					//swap rows and columns becouse we must allocate memory for 
																					//new matrix which is transpose.
	for(unsigned i = 0 ; i < rows ; ++i)
		for(unsigned j = 0 ; j < columns ; ++j)
			_transpose[j][i] = matrix[i][j];

	return _transpose;
}