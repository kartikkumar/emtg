//EMTG Matrix class
//header-only library for matrices
//Jacob Englander 1/13/2013

#pragma once

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <exception>

using namespace std;

namespace EMTG { namespace math {

	enum MatrixType {standard, identity, Rxhat, Ryhat, Rzhat};

	template <class T> class Matrix
	{
	public:
		//default constructor
		Matrix();

		//constructor for standard matrices when size is known but there is no initialization data
		Matrix(const int& n_in, const int& m_in);

		//constructor for standard matrices when size is known and an array of data is supplied
		Matrix(const int& n_in, const int& m_in, const T* input_data);

		//constructor for standard matrices when size is known and a scalar is to be assigned to all values
		Matrix(const int& n_in, const int& m_in, const T& value);

		//constructor for special matrix types, all of which happen to be square
		Matrix(const int& n_in, const T& input_data, const MatrixType& type_in);

		//copy constructor
		Matrix(const Matrix& OtherMatrix);

		//destructor
		virtual ~Matrix();

		//********************************************************************
		//public methods. These methods are virtual so that other classes can overload them
		
		//value assignment functions
		virtual void assign_all(const T* input_data);
		virtual void assign_entry(const int& i, const int& j, const T& value);
		virtual void assign_zeros();
		virtual void assign_constant(const T& value);
		virtual void assign_row(const int i, const T* rowvalues);
		virtual void assign_row(const int i, const Matrix& TheRow);
		virtual void assign_column(const int j, const T* columnvalues);
		virtual void assign_column(const int j, const Matrix& TheColumn);
		virtual void clear();
		virtual void resize(const int& n_in, const int& m_in);
		virtual void construct_rotation_matrix(const T& angle);
		virtual void construct_identity_matrix();
		virtual Matrix& operator= (const T& scalar);
		virtual Matrix& operator= (const Matrix& OtherMatrix);
		

		//get methods
		virtual T getentry(const int i, const int j) const;
		virtual Matrix getrow(const int i) const;
		virtual Matrix getcolumn(const int j) const;
		virtual T operator() (const int i, const int j) const;
		virtual T operator() (const int i) const;
		virtual T& getreference (const int i, const int j);
		virtual T& operator() (const int i, const int j);
		virtual T& operator() (const int i);

		//print methods
		virtual void print_to_screen();
		virtual void print_to_file();
		virtual void print_to_file(std::string filename);

		//basic math operations

		//sign change
		virtual Matrix& operator- ();

		//matrix addition
		virtual Matrix operator+ (const Matrix& OtherMatrix);
		virtual Matrix& operator+= (const Matrix& OtherMatrix);
		virtual Matrix operator+ (const T& scalar);
		virtual Matrix& operator+= (const T& scalar);

		//matrix subtraction
		virtual Matrix operator- (const Matrix& OtherMatrix);
		virtual Matrix& operator-= (const Matrix& OtherMatrix);
		virtual Matrix operator- (const T& scalar);
		virtual Matrix& operator-= (const T& scalar);

		//right multiplication (left multiplication is done by using the other object's multiplier)
		virtual Matrix operator* (const Matrix& OtherMatrix);
		virtual Matrix operator* (const T& scalar);
		virtual Matrix& operator*= (const Matrix& OtherMatrix);
		virtual Matrix& operator*= (const T& scalar);
		
		//comparisons
		virtual bool operator== (const Matrix& OtherMatrix);
		virtual bool operator!= (const Matrix& OtherMatrix);
#ifdef _EMTG_thruth_table
		virtual Matrix<bool> operator<= (const Matrix& OtherMatrix);
		virtual Matrix<bool> operator>= (const Matrix& OtherMatrix);
		virtual Matrix<bool> operator< (const Matrix& OtherMatrix);
		virtual Matrix<bool> operator> (const Matrix& OtherMatrix);
#endif
		
		//multiplication(element-wise)
		virtual Matrix element_multiply (const Matrix& OtherMatrix);
		virtual Matrix& element_multiply_equals (const Matrix& OtherMatrix);

		//division (element-wise)
		virtual Matrix operator/ (const T& scalar);
		virtual Matrix element_divide (const Matrix& OtherMatrix);
		virtual Matrix& operator/= (const T& scalar);
		virtual Matrix& element_divide_equals (const Matrix& OtherMatrix);
		
		//special matrix math
		virtual Matrix transpose();
		virtual T determinant();
		virtual Matrix inverse();
		virtual Matrix horz_cat(const Matrix& OtherMatrix);
		virtual Matrix vert_cat(const Matrix& OtherMatrix);
		virtual Matrix reshape(const int new_n, const int new_m);
		/*
		virtual Matrix row_reduced_echelon_form();
		virtual int column_rank();
		virtual int row_rank();
		virtual Matrix principle_minors();
		virtual T trace();*/

		//vector math
		virtual T norm();
		virtual Matrix unitize();
		virtual T dot(const Matrix& OtherMatrix);
		virtual Matrix cross(const Matrix& OtherMatrix);
		virtual Matrix unitcross(const Matrix& OtherMatrix);
		virtual void cross_in_place(const Matrix& OtherMatrix, Matrix& TargetMatrix);
		

	protected:
		// fields
		int n; //number of rows
		int m; //number of columns
		MatrixType type;
		bool issquare;
		bool isvector;
		T* values;
	}; //end class definition: Matrix

	
	//***************************************************
	//method definitions
	template <class T> Matrix<T>::Matrix() : n(1), m(1), type(standard)
	{
		//default constructor constructs a 1x1 matrix initialized to zero

		issquare = true;
		isvector = true;

		values = new T[n*m];

		assign_zeros();
	}

	//constructor when size is known but there is no initialization data
	template <class T> Matrix<T>::Matrix(const int& n_in, const int& m_in) : n(n_in), m(m_in), type(standard)
	{

		issquare = (n == m ? true : false);
		isvector = issquare ? false : n == 1 ? true : m == 1 ? true : false;

		values = new T[n*m];

		assign_zeros();
	}

	//constructor when size is known and an array of data is supplied
	template <class T> Matrix<T>::Matrix(const int& n_in, const int& m_in, const T* input_data) : n(n_in), m(m_in), type(standard)
	{
		issquare = (n == m ? true : false);
		isvector = issquare ? false : n == 1 ? true : m == 1 ? true : false;

		values = new T[n*m];

		assign_all(input_data);
	}

	template <class T> Matrix<T>::Matrix(const int& n_in, const int& m_in, const T& value) : n(n_in), m(m_in), type(standard)
	{
		issquare = (n == m ? true : false);
		isvector = issquare ? false : n == 1 ? true : m == 1 ? true : false;

		values = new T[n*m];

		assign_constant(value);
	}

	
	template <class T> Matrix<T>::Matrix(const int& n_in, const T& input_data, const MatrixType& type_in) : n(n_in), m(n_in), type(type_in)
	{
		//all special matrix types are square
		issquare = true;
		isvector = issquare ? false : n == 1 ? true : m == 1 ? true : false;

		//allocate memory
		values = new T[n*n];

		//assign the matrix type
		type = type_in;

		switch (type)
		{
		case standard:
			{
				cout << "EXCEPTION: Do not use the special matrix type constructor to build a standard matrix!" << endl;
				throw 1711;
				break;
			}
		case identity:
			{
				construct_identity_matrix();
				break;
			}
		default: //currently this is all of the rotation matrices
			{
				construct_rotation_matrix(input_data);
				break;
			}
		}
	}

	//copy constructor
	template <class T> Matrix<T>::Matrix(const Matrix& OtherMatrix)
	{
		n = OtherMatrix.n;
		m = OtherMatrix.m;
		type = OtherMatrix.type;

		issquare = (n == m ? true : false);
		isvector = issquare ? false : n == 1 ? true : m == 1 ? true : false;

		values = new T[n*m];
		assign_all(OtherMatrix.values);
	}


	//***************************************************
	//destructor
	template <class T> Matrix<T>::~Matrix()
	{
		delete[] values;
	}

	//***************************************************
	//assignment functions
	template<class T> void Matrix<T>::assign_all(const T* input_data)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m + j] = input_data[i*m + j];
		}
	}

	template<class T> void Matrix<T>::assign_entry(const int& i, const int& j, const T& value)
	{
		if (i >= n || j >= m)
		{
			cout << "EXCEPTION: Invalid index in [" << i << ", " << j << "]" << endl;
			throw 1711;
		}
		values[i*m + j] = value;
	}

	template<class T> void Matrix<T>::assign_zeros()
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m + j] = 0;
		}
	}

	template <class T> void Matrix<T>::assign_constant(const T& value)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m + j] = value;
		}
	}

	template <class T> void Matrix<T>::assign_row(const int i, const T* rowvalues)
	{
		for (int j = 0; j < m; ++j)
			values[i*m+j] = rowvalues[j];
	}

	template <class T> void Matrix<T>::assign_row(const int i, const Matrix& TheRow)
	{
		if (!(TheRow.m == m && TheRow.n == 1))
		{
			cout << "EXCEPTION: Cannot assign row because dimension is incorrect." << endl;
			throw 1711;
		}

		for (int j = 0; j < m; ++j)
			values[i*m+j] = TheRow(j);
	}

	template <class T> void Matrix<T>::assign_column(const int j, const T* columnvalues)
	{
		for (int i = 0; i < n; ++i)
			values[i*m+j] = columnvalues[i];
	}

	template <class T> void Matrix<T>::assign_column(const int j, const Matrix& TheColumn)
	{
		if (!(TheColumn.n == n && TheColumn.m == 1))
		{
			cout << "EXCEPTION: Cannot assign column because dimension is incorrect." << endl;
			throw 1711;
		}

		for (int i = 0; i < n; ++i)
			values[i*m+j] = TheColumn(i);
	}

	template <class T> void Matrix<T>::construct_identity_matrix()
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				values[i*m+j] = (i == j) ? 1 : 0;
			}
		}
	}

	template <class T> void Matrix<T>::construct_rotation_matrix(const T& angle)
	{
		//pre-compute trig values
		T cangle = cos(angle);
		T sangle = sin(angle);

		switch (type)
		{
			case Rxhat:
			{
				assign_entry(0, 0,  1.0);
				assign_entry(0, 1,  0.0);
				assign_entry(0, 2,  0.0);
				assign_entry(1, 0,  0.0);
				assign_entry(1, 1,  cangle);
				assign_entry(1, 2, -sangle);
				assign_entry(2, 0,  0);
				assign_entry(2, 1,  sangle);
				assign_entry(2, 2,  cangle);
				break;
			}
			case Ryhat:
			{
				assign_entry(0, 0,  cangle);
				assign_entry(0, 1,  0.0);
				assign_entry(0, 2,  sangle);
				assign_entry(1, 0,  0.0);
				assign_entry(1, 1,  1.0);
				assign_entry(1, 2,  0.0);
				assign_entry(2, 0, -sangle);
				assign_entry(2, 1,  0.0);
				assign_entry(2, 2,  cangle);
				break;
			}
			case Rzhat:
			{
				assign_entry(0, 0,  cangle);
				assign_entry(0, 1, -sangle);
				assign_entry(0, 2,  0.0);
				assign_entry(1, 0,  sangle);
				assign_entry(1, 1,  cangle);
				assign_entry(1, 2,  0.0);
				assign_entry(2, 0,  0.0);
				assign_entry(2, 1,  0.0);
				assign_entry(2, 2,  1.0);
				break;
			}
		}
	}

	template<class T> void Matrix<T>::clear()
	{
		delete[] values;
		n = 1;
		m = 1;

		issquare = true;

		values = new T[n*m];
		assign_zeros();
	}
	
	template<class T> void Matrix<T>::resize(const int& n_in, const int& m_in)
	{
		delete[] values;
		n = n_in;
		m = m_in;

		issquare = (n == m ? true : false);

		values = new T[n*m];
		assign_zeros();
	}

	template<class T> Matrix<T>& Matrix<T>::operator= (const T& scalar)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m+j] = scalar;
		}

		return *this;
	}
	
	template<class T> Matrix<T>& Matrix<T>::operator= (const Matrix& OtherMatrix)
	{
		n = OtherMatrix.n;
		m = OtherMatrix.m;

		delete[] values;
		values = new T[n*m];
		type = OtherMatrix.type;

		assign_all(OtherMatrix.values);

		return *this;
	}

	//***************************************************
	//get methods
	template<class T> T Matrix<T>::getentry(const int i, const int j) const
	{
		if (i > n-1 || j > m-1 || i < 0 || j < 0)
		{
			cout << "EXCEPTION: Invalid index in [" << i << ", " << j << "]" << endl;
			throw 1711;
		}
		return values[i*m+j];
	}

	template<class T> Matrix<T> Matrix<T>::getrow(const int i) const
	{
		Matrix<T> TheRow(1, m);

		for (int j = 0; j < m; ++j)
			TheRow(0, j) = getentry(i, j);

		return TheRow;
	}

	template<class T> Matrix<T> Matrix<T>::getcolumn(const int j) const
	{
		Matrix<T> TheColumn(n, 1);

		for (int i = 0; i < n; ++i)
			TheColumn(i, 0) = getentry(i, j);

		return TheColumn;
	}

	template<class T> T Matrix<T>::operator() (const int i, const int j) const
	{
		return getentry(i, j);
	}

	template<class T> T Matrix<T>::operator() (const int i) const
	{
		if (!isvector)
		{
			cout << "EXCEPTION: Cannot access a non-vector matrix using only one index" << endl;
			throw 1711;
		}

		if (m == 1)
			return getentry(i, 0);
		else
			return getentry(0, i);
	}

	template<class T> T& Matrix<T>::getreference(const int i, const int j)
	{
		if (i > n-1 || j > m-1 || i < 0 || j < 0)
		{
			cout << "EXCEPTION: Invalid index in [" << i << ", " << j << "]" << endl;
			throw 1711;
		}
		return values[i*m+j];
	}
		
	template<class T> T& Matrix<T>::operator() (const int i, const int j)
	{
		return getreference(i, j);
	}

	template<class T> T& Matrix<T>::operator() (const int i)
	{
		if (!isvector)
		{
			cout << "EXCEPTION: Cannot access a non-vector matrix using only one index" << endl;
			throw 1711;
		}

		if (m == 1)
			return getreference(i, 0);
		else
			return getreference(0, i);
	}

	//***************************************************
	//print methods
	template<class T> void Matrix<T>::print_to_screen()
	{
		cout.width(15);
		cout.precision(20);

		for (int i = 0; i < n; ++i)
		{
			cout.precision(20);
			cout << values[i*m];
			for (int j = 1; j < m; ++j)
			{
				cout.width(25);
				cout.precision(20);
				cout << values[i*m+j];
			}
			cout << endl;
		}
	}

	template<class T> void Matrix<T>::print_to_file()
	{
		this->print_to_file("Matrix.txt");
	}

	template<class T> void Matrix<T>::print_to_file(std::string filename)
	{
		std::ofstream outputfile(filename.c_str(), ios::trunc);

		for (int i = 0; i < n; ++i)
		{
			outputfile.width(30);
			outputfile.precision(20);
			outputfile << values[i*m];
			for (int j = 1; j < m; ++j)
			{
				outputfile.width(30);
				outputfile.precision(20);
				outputfile << values[i*m+j];
			}
			outputfile << endl;
		}

		outputfile.close();
	}

	//***************************************************
	//basic math operations

	//sign change
	template<class T> Matrix<T>& Matrix<T>::operator- ()
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m+j] *= -1;
		}

		return *this;
	}

	//***************************************************
	//matrix addition
	template<class T> Matrix<T> Matrix<T>::operator+ (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are the same size
		if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
		{
			cout << "EXCEPTION: Matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<T> NewMatrix(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				NewMatrix.assign_entry(i, j, getentry(i,j) + OtherMatrix.getentry(i,j));
			}
		}

		return NewMatrix;
	}

	template<class T> Matrix<T> Matrix<T>::operator+ (const T& scalar)
	{
		Matrix<T> NewMatrix(n, m, values);

		NewMatrix += scalar;

		return NewMatrix;
	}

	template<class T> Matrix<T>& Matrix<T>::operator+= (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are the same size
		if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
		{
			cout << "EXCEPTION: Matrix sizes do not match" << endl;
			throw 1711;
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				values[i*m+j] += OtherMatrix.values[i*m+j];
			}
		}

		return *this;
	}

	template<class T> Matrix<T>& Matrix<T>::operator+= (const T& scalar)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				values[i*m+j] += scalar;
			}
		}

		return *this;
	}

	//***************************************************
	//matrix subtraction
	template<class T> Matrix<T> Matrix<T>::operator- (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are the same size
		if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
		{
			cout << "EXCEPTION: Matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<T> NewMatrix(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				NewMatrix.assign_entry(i, j, getentry(i,j) - OtherMatrix.getentry(i,j));
			}
		}

		return NewMatrix;
	}

	template<class T> Matrix<T> Matrix<T>::operator- (const T& scalar)
	{
		Matrix<T> NewMatrix(n, m, values);

		NewMatrix -= scalar;

		return NewMatrix;
	}

	template<class T> Matrix<T>& Matrix<T>::operator-= (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are the same size
		if (!(n == OtherMatrix.n) || !(m == OtherMatrix.m))
		{
			cout << "EXCEPTION: Matrix sizes do not match" << endl;
			throw 1711;
		}

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				values[i*m+j] -= OtherMatrix.values[i*m+j];
			}
		}

		return *this;
	}

	template<class T> Matrix<T>& Matrix<T>::operator-= (const T& scalar)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				values[i*m+j] -= scalar;
			}
		}

		return *this;
	}

	//right multiplication (left multiplication is done by using the other object's left multiplier
	template<class T> Matrix<T> Matrix<T>::operator* (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are compatible
		if (!(m == OtherMatrix.n))
		{
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<T> NewMatrix(n, OtherMatrix.m);

		for (int i = 0; i < NewMatrix.n; ++i)
		{
			for (int j = 0; j < NewMatrix.m; ++j)
			{
				for (int k = 0; k < m; ++k)
					NewMatrix(i,j) += values[i*m + k] * OtherMatrix.values[k*OtherMatrix.m + j];
			}
		}


		return NewMatrix;
	}

	template<class T> Matrix<T> Matrix<T>::operator* (const T& scalar)
	{
		Matrix<T> NewMatrix(n, m, values);

		NewMatrix *= scalar;

		return NewMatrix;
	}

	template<class T> Matrix<T>& Matrix<T>::operator*= (const Matrix& OtherMatrix)
	{
		Matrix<T> NewMatrix = *this * OtherMatrix;

		assign_all(NewMatrix.values);

		return *this;
	}
	
	template<class T> Matrix<T>& Matrix<T>::operator*= (const T& scalar)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m + j] *= scalar;
		}

		return *this;
	}

	//comparators
	template<class T> bool Matrix<T>::operator== (const Matrix& OtherMatrix)
	{
		if (!(m == OtherMatrix.m))
			return false;

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				if (!(values[i*m+j] == OtherMatrix.values[i*m+j]))
					return false;
			}
		}
		
		return true;
	}

	template<class T> bool Matrix<T>::operator!= (const Matrix& OtherMatrix)
	{
		return !(*this == OtherMatrix);
	}

#ifdef _EMTG_thruth_table
	template<class T> Matrix<bool> Matrix<T>::operator<= (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are compatible
		if (!(n == OtherMatrix.n && m == OtherMatrix.m))
		{
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<bool> TruthTable(n, m);
		
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				TruthTable(i,j) = values[i*m+j] <= OtherMatrix.values[i*m+j] ? true : false;
		}

		return TruthTable;
	}

	template<class T> Matrix<bool> Matrix<T>::operator>= (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are compatible
		if (!(n == OtherMatrix.n && m == OtherMatrix.m))
		{
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<bool> TruthTable(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				TruthTable(i,j) = values[i*m+j] >= OtherMatrix.values[i*m+j] ? true : false;
		}

		return TruthTable;
	}

	template<class T> Matrix<bool> Matrix<T>::operator< (const Matrix& OtherMatrix)
		{
		//first check to see if the matrices are compatible
		if (!(n == OtherMatrix.n && m == OtherMatrix.m))
		{
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<bool> TruthTable(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				TruthTable(i,j) = values[i*m+j] < OtherMatrix.values[i*m+j] ? true : false;
		}

		return TruthTable;
	}

	template<class T> Matrix<bool> Matrix<T>::operator> (const Matrix& OtherMatrix)
	{
		//first check to see if the matrices are compatible
		if (!(n == OtherMatrix.n && m == OtherMatrix.m))
		{
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}

		Matrix<bool> TruthTable(n, m);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				TruthTable(i,j) = values[i*m+j] > OtherMatrix.values[i*m+j] ? true : false;
		}

		return TruthTable;
	}
#endif //_EMTG_thruth_table

	//multiplication(element-wise)
	template<class T> Matrix<T> Matrix<T>::element_multiply (const Matrix& OtherMatrix)
	{
		//first check to make sure that the operation can be done
		//is the number of rows the same?
		if (n == OtherMatrix.n)
		{
			//is the number of columns the same?
			if (m == OtherMatrix.m)
			{
				//the number of rows and columns is the same, so operate on the entire matrix
				Matrix<T> NewMatrix(n,m);

				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] * OtherMatrix(i,j);
				}

				return NewMatrix;
			}
			else if (OtherMatrix.m == 1)
			{
				//we actually want to operate through by a column vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] * OtherMatrix.values[i];
				}

				return NewMatrix;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else if (m == OtherMatrix.m)
		{
			if (OtherMatrix.n == 1)
			{
				//we actually want to operate through by a row vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] * OtherMatrix.values[j];
				}

				return NewMatrix;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else
		{
			//can't do this
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}
	}
	
	template<class T> Matrix<T>& Matrix<T>::element_multiply_equals (const Matrix& OtherMatrix)
	{
		//first check to make sure that the operation can be done
		//is the number of rows the same?
		if (n == OtherMatrix.n)
		{
			//is the number of columns the same?
			if (m == OtherMatrix.m)
			{
				//the number of rows and columns is the same, so operate on the entire matrix

				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						values[i*m+j] *= OtherMatrix(i,j);
				}

				return *this;
			}
			else if (OtherMatrix.m == 1)
			{
				//we actually want to operate through by a column vector

				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						values[i*m+j] *= OtherMatrix.values[i];
				}
				
				return *this;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else if (m == OtherMatrix.m)
		{
			if (OtherMatrix.n == 1)
			{
				//we actually want to operate through by a row vector
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						values[i*m+j] *= OtherMatrix.values[j];
				}

				return *this;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else
		{
			//can't do this
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}
	}

	//division (element-wise)
	template<class T> Matrix<T> Matrix<T>::operator/ (const T& scalar)
	{
		Matrix<T> NewMatrix(n,m, values);

		NewMatrix /= scalar;

		return NewMatrix;
	}
	
	template<class T> Matrix<T> Matrix<T>::element_divide (const Matrix& OtherMatrix)
	{
		//first check to make sure that the operation can be done
		//is the number of rows the same?
		if (n == OtherMatrix.n)
		{
			//is the number of columns the same?
			if (m == OtherMatrix.m)
			{
				//the number of rows and columns is the same, so operate on the entire matrix
				Matrix<T> NewMatrix(n,m);

				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] / OtherMatrix(i,j);
				}

				return NewMatrix;
			}
			else if (OtherMatrix.m == 1)
			{
				//we actually want to operate through by a column vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] / OtherMatrix.values[i];
				}

				return NewMatrix;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else if (m == OtherMatrix.m)
		{
			if (OtherMatrix.n == 1)
			{
				//we actually want to operate through by a row vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] / OtherMatrix.values[j];
				}

				return NewMatrix;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else
		{
			//can't do this
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}
	}

	template<class T> Matrix<T>& Matrix<T>::operator/= (const T& scalar)
	{
		T inverse_scalar = 1 / scalar;

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				values[i*m+j] *= inverse_scalar;
		}

		return *this;
	}

	template<class T> Matrix<T>& Matrix<T>::element_divide_equals (const Matrix& OtherMatrix)
	{
		//first check to make sure that the operation can be done
		//is the number of rows the same?
		if (n == OtherMatrix.n)
		{
			//is the number of columns the same?
			if (m == OtherMatrix.m)
			{
				//the number of rows and columns is the same, so operate on the entire matrix
				Matrix<T> NewMatrix(n,m);

				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] /= OtherMatrix(i,j);
				}

				return *this;
			}
			else if (OtherMatrix.m == 1)
			{
				//we actually want to operate through by a column vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] /= OtherMatrix.values[i];
				}

				return *this;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else if (m == OtherMatrix.m)
		{
			if (OtherMatrix.n == 1)
			{
				//we actually want to operate through by a row vector
				Matrix<T> NewMatrix(n,m);
				
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < m; ++j)
						NewMatrix(i,j) = values[i*m+j] /= OtherMatrix.values[j];
				}

				return *this;
			}
			else
			{
				//can't do this
				cout << "EXCEPTION: matrix sizes do not match" << endl;
				throw 1711;
			}
		}
		else
		{
			//can't do this
			cout << "EXCEPTION: matrix sizes do not match" << endl;
			throw 1711;
		}
	}

	//*****************************************************
	//Special Matrix Math Functions
	template<class T> Matrix<T> Matrix<T>::transpose()
	{
		Matrix<T> NewMatrix(m, n);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
				NewMatrix(j,i) = getentry(i,j);
		}
		
		return NewMatrix;
	}

	template<class T> T Matrix<T>::determinant()
	{
		if (!issquare)
		{
			cout << "EXCEPTION: Cannot take the determinant of a non-square matrix!" << endl;
			throw 1711;
		}

		if (n == 1)
		{
			//the determinant of a 1x1 is trivial

			return values[0];
		}
		else if (n == 2)
		{
			//the determinant of a 2x2 is ad - bc

			return values[0] * values[3] - values[1] * values[2];
		}
		else if (n == 3)
		{
			//the determinant of a 3x3 is aei + bfg + cdh - ceg - bdi - afh

			return values[0]*values[4]*values[8] + values[1]*values[5]*values[6] + values[2]*values[3]*values[7] - values[2]*values[4]*values[6] - values[1]*values[3]*values[8] - values[0]*values[5]*values[7];
		}
		else
		{
			//determinants of larger matrices are much more complex. Perhaps someday I will link to LAPACK

			cout << "EXCEPTION: Determinants of matrices with n > 3 are not yet supported" << endl;
			throw 1711;
		}
	}

	template<class T> Matrix<T> Matrix<T>::inverse()
	{
		//can only invert square matrices
		if (!issquare)
		{
			cout << "EXCEPTTION: Non-square matrices cannot be inverted." << endl;
			throw 1711;
		}

		switch (type)
		{
			case standard:
			{
				//first check if the matrix is invertible, i.e is the determinant nonzero?
				T TheDeterminant = determinant();
				if (TheDeterminant == 0)
				{
					cout << "EXCEPTION: Cannot invert matrix because determinant is zero." << endl;
					throw 1711;
				}
				//create a matrix to hold the inverse
				Matrix<T> TheInverse(n, m);

				switch (n)
				{
					case 1:
					{
						TheInverse(0,0) = *values;
						break;
					}
					case 2:
					{
						TheInverse(0,0) = values[3];
						TheInverse(0,1) = -values[1];
						TheInverse(1,0) = -values[2];
						TheInverse(1,1) = values[0];
						TheInverse /= TheDeterminant;
						break;
					}
					case 3:
					{
						TheInverse(0,0) = values[4]*values[8] - values[5]*values[7];
						TheInverse(0,1) = values[5]*values[6] - values[3]*values[8];
						TheInverse(0,2) = values[3]*values[7] - values[4]*values[6];
						TheInverse(1,0) = values[2]*values[7] - values[1]*values[8];
						TheInverse(1,1) = values[0]*values[8] - values[2]*values[6];
						TheInverse(1,2) = values[6]*values[1] - values[0]*values[7];
						TheInverse(2,0) = values[1]*values[5] - values[2]*values[4];
						TheInverse(2,1) = values[2]*values[3] - values[0]*values[5];
						TheInverse(2,2) = values[0]*values[4] - values[1]*values[3];
						TheInverse /= TheDeterminant;
						break;
					}
					default:
					{
						//library currently only inverts up to 3x3
						cout << "EXCEPTION: Currently only up to a 3x3 matrix may be inverted." << endl;
						throw 1711;
						break;
					}
				}
				
				return TheInverse.transpose();
				break;
			}
			case identity:
			{
				return *this;
				break;
			}
			case Rxhat:
			{
				return transpose();
				break;
			}
			case Ryhat:
			{
				return transpose();
				break;
			}
			case Rzhat:
			{
				return transpose();
				break;
			}
		}

		return *this;
	}

	template<class T> Matrix<T> Matrix<T>::horz_cat(const Matrix& OtherMatrix)
	{
		if (!(n == OtherMatrix.n))
		{
			cout << "EXCEPTION: Cannot perform horizontal concatenation because vertical dimension does not match" << endl;
			throw 1711;
		}

		Matrix NewMatrix(n, m + OtherMatrix.m);
		
		for (int i = 0; i < n; ++i)
		{
			for (int j1 = 0; j1 < m; ++j1)
				NewMatrix(i, j1) = values[i*m+j1];

			for (int j2 = 0; j2 < OtherMatrix.m; ++j2)
				NewMatrix(i,j2 + m) = OtherMatrix(i,j2);
		}

		return NewMatrix;
	}

	template<class T> Matrix<T> Matrix<T>::vert_cat(const Matrix& OtherMatrix)
	{
		if (!(m == OtherMatrix.m))
		{
			cout << "EXCEPTION: Cannot perform vertical concatenation because horizontal dimension does not match" << endl;
			throw 1711;
		}

		Matrix NewMatrix(n + OtherMatrix.n, m);
		

		for (int j = 0; j < m; ++j)
		{
			for (int i1 = 0; i1 < n; ++i1)
				NewMatrix(i1, j) = values[i1*m+j];

			for (int i2 = 0; i2 < OtherMatrix.n; ++i2)
				NewMatrix(i2 + n,j) = OtherMatrix(i2,j);
		}

		return NewMatrix;
	}

	//reshape function; intended to work just like reshape in MATLAB
	template<class T> Matrix<T> Matrix<T>::reshape(const int new_n, const int new_m)
	{
		if (!(new_n*new_m == n*m))
		{
			cout << "EXCEPTION: Matrix sizes much match in a reshape!" << endl;
			throw 1711;
		}

		Matrix<T> NewMatrix(new_n, new_m, values);

		return NewMatrix;
	}
	
	//*************************************************
	//vector math
	template<class T> T Matrix<T>::norm()
	{
		return sqrt(dot(*this));		
	}

	template<class T> Matrix<T> Matrix<T>::unitize()
	{

		return (*this / this->norm());
	}

	template<class T> T Matrix<T>::dot(const Matrix& OtherMatrix)
	{
		T TheDot = 0;

		if (!isvector || !OtherMatrix.isvector)
		{
			cout << "EXCEPTION: Cannot take the norm of a non-vector" << endl;
			throw 1711;
		}

		int dim = n == 1 ? m : n;
		int Otherdim = OtherMatrix.n == 1 ? OtherMatrix.m : OtherMatrix.n;

		if (!(dim == Otherdim))
		{
			cout << "EXCEPTION: Can only take the dot product of equal length vectors" << endl;
			throw 1711;
		}

		for (int k = 0; k < dim; ++k)
			TheDot += values[k]*OtherMatrix.values[k];

		return TheDot;	
	}

	template<class T> Matrix<T> Matrix<T>::cross(const Matrix& OtherMatrix)
	{
		//cross product only works on 3-vectors
		if (!isvector || !OtherMatrix.isvector)
		{
			cout << "EXCEPTION: Cannot take the norm of a non-vector" << endl;
			throw 1711;
		}

		int dim = n == 1 ? m : n;
		int Otherdim = OtherMatrix.n == 1 ? OtherMatrix.m : OtherMatrix.n;

		if (!(dim == 3 && Otherdim == 3))
		{
			cout << "EXCEPTION: Can only take cross product of two 3-vectors" << endl;
			throw 1711;
		}

		Matrix TheCross(3,1);

		TheCross(0,0) = values[1]*OtherMatrix.values[2] - values[2]*OtherMatrix.values[1];
		TheCross(1,0) = values[2]*OtherMatrix.values[0] - values[0]*OtherMatrix.values[2];
		TheCross(2,0) = values[0]*OtherMatrix.values[1] - values[1]*OtherMatrix.values[0];

		return TheCross;
	}

	template<class T> Matrix<T> Matrix<T>::unitcross(const Matrix& OtherMatrix)
	{
		Matrix TheUnitCross = cross(OtherMatrix);

		return TheUnitCross / TheUnitCross.norm();
	}

	template<class T> void Matrix<T>::cross_in_place(const Matrix& OtherMatrix, Matrix& TargetMatrix)
	{
		int TargetDim = TargetMatrix.n == 1 ? TargetMatrix.m : TargetMatrix.n;
		if (!(TargetDim == 3))
		{
			cout << "EXCEPTION: Target matrix of a cross product must be a 3-vector" << endl;
			throw 1711;
		}

		TargetMatrix.values[0] = values[1]*OtherMatrix.values[2] - values[2]*OtherMatrix.values[1];
		TargetMatrix.values[1] = values[2]*OtherMatrix.values[0] - values[0]*OtherMatrix.values[2];
		TargetMatrix.values[2] = values[0]*OtherMatrix.values[1] - values[1]*OtherMatrix.values[0];
	}
}}// close namespace