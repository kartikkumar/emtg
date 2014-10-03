//basic 6x6 state transition matrix
//for use with Kepler propagators
//a header-only library

#include <iostream>
#include <fstream>

#include "STM.h"



namespace Kepler
{
	//constructors
	STM::STM(void)
	{
		this->clear();
	}

	STM::STM(double* input_data)
	{
		this->assign_all(input_data);
	}

	//destructor does not need to do anything
	STM::~STM(void) {}

	//methods
	//value assignment functions
	void STM::assign_all(const double* input_data)
	{
		for (int k = 0; k < 36; ++k)
			this->data[k] = input_data[k];
	}

	void STM::assign_entry(const int& i, const int& j, const double& value)
	{
		this->data[i*6+j] = value;
	}

	void STM::clear()
	{
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
				this->assign_entry(i,j, (i==j ? 1.0 : 0.0));
		}
	}

	STM& STM::operator= (const STM& OtherSTM)
	{
		for (int k = 0; k < 36; ++k)
			this->data[k] = OtherSTM.data[k];

		return *this;
	}

	//get methods
	double STM::getentry(const int& i, const int& j) const
	{
		return this->data[i*6+j];
	}

	double STM::operator() (const int& i, const int& j) const
	{
		return this->getentry(i, j);
	}

	double& STM::getreference (const int& i, const int& j)
	{
		return this->data[i*6+j];
	}

	double& STM::operator() (const int& i, const int& j)
	{
		return this->getreference(i, j);
	}

	//print methods
	void STM::print_to_screen()
	{
		std::cout.width(15);
		std::cout.precision(20);

		for (int i = 0; i < 6; ++i)
		{
			std::cout.precision(20);
			std::cout << this->data[i*6];
			for (int j = 1; j < 6; ++j)
			{
				std::cout.width(25);
				std::cout.precision(20);
				std::cout << this->data[i*6+j];
			}
			std::cout << std::endl;
		}
	}

	void STM::print_to_file()
	{
		this->print_to_file("STM_output.txt");
	}

	void STM::print_to_file(std::string filename)
	{
		std::ofstream outputfile(filename.c_str(), std::ios::trunc);

		for (int i = 0; i < 6; ++i)
		{
			outputfile.width(30);
			outputfile.precision(20);
			outputfile << this->data[i*6];
			for (int j = 1; j < 6; ++j)
			{
				outputfile.width(30);
				outputfile.precision(20);
				outputfile << this->data[i*6+j];
			}
			outputfile << std::endl;
		}

		outputfile.close();
	}

	//sign change
	STM& STM::operator- ()
	{
		for (int k = 0; k < 36; ++k)
			this->data[k] *= -1;

		return *this;
	}

	//right multiplication (left multiplication is done by using the other object's multiplier)
	STM STM::operator* (const STM& OtherSTM)
	{
		STM NewSTM;

		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				for (int k = 0; k < 6; ++k)
					NewSTM(i,j) += this->data[i*6 + k] * OtherSTM.data[k*6 + j];
			}
		}

		return NewSTM;
	}

	STM STM::operator* (const double& scalar)
	{
		STM NewSTM;

		NewSTM.assign_all(this->data);

		NewSTM *= scalar;

		return NewSTM;
	}

	void STM::state_vector_right_multiply (const double* initial_perturbation_vector, double* final_perturbation_vector)
	{
		for (int i = 0; i < 6; ++i)
		{
			final_perturbation_vector[i] = 0.0;
			for (int j = 0; j < 6; ++j)
				final_perturbation_vector[i] += this->getentry(i, j) * initial_perturbation_vector[j];
		}
	}

	STM& STM::operator*= (const STM& OtherSTM)
	{
		STM NewSTM = *this * OtherSTM;

		this->assign_all(NewSTM.data);

		return *this;
	}

	STM& STM::operator*= (const double& scalar)
	{
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
				this->data[i*6 + j] *= scalar;
		}

		return *this;
	}
		
	//comparisons
	bool STM::operator== (const STM& OtherSTM)
	{
		
		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
			{
				if (!(this->data[i*6+j] == OtherSTM.data[i*6+j]))
					return false;
			}
		}
		
		return true;
	}

	bool STM::operator!= (const STM& OtherSTM)
	{
		return !(*this == OtherSTM);
	}

	//special matrix math
	STM STM::transpose()
	{
		STM NewSTM;

		for (int i = 0; i < 6; ++i)
		{
			for (int j = 0; j < 6; ++j)
				NewSTM(j,i) = this->getentry(i,j);
		}
		
		return NewSTM;
	}

	STM STM::inverse()
	{
		//from Battin p458, top of page
		//this is annoying

		STM NewSTM;

//[  P33,  P43,  P53, -P03, -P13, -P23]
		NewSTM.assign_entry(0,0,this->getentry(3,3));
		NewSTM.assign_entry(0,1,this->getentry(4,3));
		NewSTM.assign_entry(0,2,this->getentry(5,3));
		NewSTM.assign_entry(0,3,-this->getentry(0,3));
		NewSTM.assign_entry(0,4,-this->getentry(1,3));
		NewSTM.assign_entry(0,5,-this->getentry(2,3));
//[  P34,  P44,  P54, -P04, -P14, -P24]
		NewSTM.assign_entry(1,0,this->getentry(3,4));
		NewSTM.assign_entry(1,1,this->getentry(4,4));
		NewSTM.assign_entry(1,2,this->getentry(5,4));
		NewSTM.assign_entry(1,3,-this->getentry(0,4));
		NewSTM.assign_entry(1,4,-this->getentry(1,4));
		NewSTM.assign_entry(1,5,-this->getentry(2,4));
//[  P35,  P45,  P55, -P05, -P15, -P25]
		NewSTM.assign_entry(2,0,this->getentry(3,5));
		NewSTM.assign_entry(2,1,this->getentry(4,5));
		NewSTM.assign_entry(2,2,this->getentry(5,5));
		NewSTM.assign_entry(2,3,-this->getentry(0,5));
		NewSTM.assign_entry(2,4,-this->getentry(1,5));
		NewSTM.assign_entry(2,5,-this->getentry(2,5));
//[ -P30, -P40, -P50,  P00,  P10,  P20]
		NewSTM.assign_entry(3,0,-this->getentry(3,0));
		NewSTM.assign_entry(3,1,-this->getentry(4,0));
		NewSTM.assign_entry(3,2,-this->getentry(5,0));
		NewSTM.assign_entry(3,3,this->getentry(0,0));
		NewSTM.assign_entry(3,4,this->getentry(1,0));
		NewSTM.assign_entry(3,5,this->getentry(2,0));
//[ -P31, -P41, -P51,  P01,  P11,  P21]
		NewSTM.assign_entry(4,0,-this->getentry(3,1));
		NewSTM.assign_entry(4,1,-this->getentry(4,1));
		NewSTM.assign_entry(4,2,-this->getentry(5,1));
		NewSTM.assign_entry(4,3,this->getentry(0,1));
		NewSTM.assign_entry(4,4,this->getentry(1,1));
		NewSTM.assign_entry(4,5,this->getentry(2,1));
//[ -P32, -P42, -P52,  P02,  P12,  P22]
		NewSTM.assign_entry(5,0,-this->getentry(3,2));
		NewSTM.assign_entry(5,1,-this->getentry(4,2));
		NewSTM.assign_entry(5,2,-this->getentry(5,2));
		NewSTM.assign_entry(5,3,this->getentry(0,2));
		NewSTM.assign_entry(5,4,this->getentry(1,2));
		NewSTM.assign_entry(5,5,this->getentry(2,2));
		
		return NewSTM;
	}
}//end namespace Kepler