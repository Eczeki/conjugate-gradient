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

	double ** matrix=0;
	createMatrix(matrix, size, size);													//Create matrix A


	deleteMatrix(matrix,size);														//delete matrix A
}

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
	double ** _transpose;
	createMatrix(_transpose,columns,rows);											//This function takes first:matrix,rows,columns but we must 
																					//swap rows and columns becouse we must allocate memory for 
																					//new matrix which is transpose.
	for(unsigned i = 0 ; i < rows ; ++i)
		for(unsigned j = 0 ; j < columns ; ++j)
			_transpose[j][i] = matrix[i][j];

	return _transpose;
}