/*
Author		:	Basiak Marcin, Pawe³ Nowak 
University	:	AGH Wydzial In¿ynierii Metali i Inoformatyki Przemyslowej
Topic		:	Conjugate Gradient Method

A		matrix[n][n] - symetryczna,rzeczywista,dodatnio okreœlona
Ax=b	gdzie x jest rozwi¹zaniem uk³adu 
		Mówimy, ¿e dwa niezerowe wektory u i v s¹ sprzê¿one (wzglêdem A) jeœli:
		u^T	* A * v == 0

Macierz symetryczna			-	macierz kwadratowa której wyrazy po³o¿one symetrycznie
								wzglêdem przek¹tnej g³ównej s¹ równe. 
								A[i][j]==A[j][i] czyli A^T=A.
								|2	1	3|
								|1	6	7|
								|3	7	9|

Macierz rzeczywista			-	Kiedy wszystkie elementy s¹ liczbami rzeczywistymi.

Macierz dodatnio okreœlona	-	nazywamy macierz A[n][n] która charaketyzuje siê:
								|W przypadku gdy A jest macierz¹ rzeczywist¹	|
								|A jest symetryczna i dla ka¿dego niezerowego	|
								|wektora x nale¿y do R^n zachodzi x^TAx>0		|
								|Co oznacza ¿e wszystkie wartoœci w³asne		|
								|macierzy A s¹ dodatnie							|

Metoda sprzê¿onych gradientów:
Rozwi¹zanie uk³adu równañ przedstawiamy jako rozwi¹zanie zagadnienia minimalizacji funkcji E(x), 
której gradient jest ró¿nic¹ pomiêdzy lew¹ a praw¹ stron¹ uk³adu równañ.

*/

#include<iostream>
#include<fstream>


void error(const char * p1 , const char * p2="");												//Function, which prints error in console. 
void createMatrix(double ** matrix , std::size_t const rows , std::size_t const columns);		//Function, which create Matrix from file.
void deleteMatrix(double ** matrix	, std::size_t size );										//Function, which delete Matrix.
bool isSymetric  (double ** matrix	, std::size_t size );										//Function, which check the matrix.
bool isPositiveDefinite(double ** matrix , std::size_t size,double x);							//Function, which check whether matrix is positive definte
double ** transpose(unsigned const rows , unsigned const columns ,  double const ** matrix);	//Function, which transpose matrix
double multiply_Ax_rows_vector(double ** A, unsigned rA, double ** x, unsigned  rx, unsigned cx, unsigned size);
void ConjugateGradientMethod(const char * from , const char * to);
double multiplyVectors(double ** v1, unsigned i1, double ** v2, unsigned i2, unsigned size);
double multiply_pT_A_p(double ** p, double ** A, unsigned k, unsigned size);

int main(int argc , char * argv[])
{
	if(argc != 3)	{	error("Wrong number of arguments...");		}

	//open input  file stream
	std::ifstream from(argv[1]);
	//open output file stream
	std::ofstream to(argv[2]);

	if(!from)	{	error("Can't open input file" , argv[1]);	}
	if(!to)		{	error("Can't open output file" , argv[2]);	}

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
	size << atoi(from);																//the first element in file is size of matrix
	
	if(size<=0){	error("Wrong size of matrix in file...:",from);	}

	double ** A = 0;
	const std::size_t max_iter = 100000;		//max iteration
	const int min_R = 0.0001;
	double ** r = 0;	//mo¿na tego te¿ nie robiæ na macierzy i nie pamiêtaæ wszystkich zmiennnych, tylko wykorzystaæ tablicê dwuelementow¹
	double ** p = 0;	//j.w.
	double ** x = 0;	//j.w.
	createMatrix(A, size, size);													//Create matrix A
	createMatrix(r, size, max_iter);		//rows are a vector, columns are a number of element	
	createMatrix(p, size, max_iter);		//rows are a vector, columns are a number of element
	createMatrix(x, size, max_iter);		//rows are a vector, columns are a number of element
	double * alpha = new double[size];
	double * beta = new double[size];

	double * b =  new double [size];

	//r0 = b - Ax0
	for(unsigned i = 0; i < size; i++){
		r[i][0] = b[i] - multiply_Ax_rows_vector(A, i, x, i, 0, size);
	}

	//p0 = r0
	for(unsigned i = 0; i < size; i++)
		r[i][0] = p[i][0];

	//algorithm
	for(unsigned k = 0; k < max_iter; k++){
		alpha[k] = multiplyVectors(r, k, r, k, size) / multiply_pT_A_p(p, A, k, size);
		for(unsigned i = 0; i < size; i++){
			x[i][k+1] = x[i][k] + alpha[k] * p[i][k];
			r[i][k+1] = r[i][k] - alpha[k] * multiply_Ax_rows_vector(A, i, p, i, k, size);
		}

		for(unsigned i = 0; i < size; i++){
			if( r[i][k+1] < min_R) break;
		}

		beta[k] = multiplyVectors(r, k, r, k, size);
		for(unsigned i = 0; i < size; i++){
			p[i][k+1] = r[i][k+1] + beta[k] * p[i][k];
		}
	}
	//the result is vector x
 
	deleteMatrix(A, size);														//delete matrix A
}

double multiply_pT_A_p(double ** p, double ** A, unsigned k, unsigned size){
	double  tmp; 
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


void error(const char * p1 , const char * p2/* = ""*/)
{
	std::cerr << p1 << " " << p2 << "\n";
	std::exit(1);
}

void createMatrix(double ** matrix , std::size_t const rows , std::size_t const columns  )
{
	matrix = new double * [rows];
	for(std::size_t i = 0 ; i < rows ; ++i)
	{
		matrix[i]=new double[columns];
	}
	if(matrix==0){	error("Bad memory allocation...");	}
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
	createMatrix(_transpose,columns,rows);											//This function takes first:matrix,rows,columns but we must 
																					//swap rows and columns becouse we must allocate memory for 
																					//new matrix which is transpose.
	for(unsigned i = 0 ; i < rows ; ++i)
		for(unsigned j = 0 ; j < columns ; ++j)
			_transpose[j][i] = matrix[i][j];

	return _transpose;
}