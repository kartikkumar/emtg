//header file for basic 6x6 state transition matrix
//for use with Kepler propagators
//a header-only library

#ifndef KEPLER_STM_CLASS
#define KEPLER_STM_CLASS

#include <string>

namespace Kepler
{
	class STM
	{
	public:
		//only has the basic constructor/destructor
		STM(void);
		STM(double* input_data);

		//destructor
		virtual ~STM(void);

		//********************************************************************
		//public methods. These methods are virtual so that other classes can overload them
		
		//value assignment functions
		virtual void assign_all(const double* input_data);
		virtual void assign_entry(const int& i, const int& j, const double& value);
		virtual void clear();
		virtual STM& operator= (const STM& OtherSTM);

		//get methods
		virtual double getentry(const int& i, const int& j) const;
		virtual double operator() (const int& i, const int& j) const;
		virtual double& getreference (const int& i, const int& j);
		virtual double& operator() (const int& i, const int& j);

		//print methods
		virtual void print_to_screen();
		virtual void print_to_file();
		virtual void print_to_file(std::string filename);

		//basic math operations
		
		//sign change
		virtual STM& operator- ();

		//right multiplication (left multiplication is done by using the other object's multiplier)
		virtual STM operator* (const STM& OtherSTM);
		virtual STM operator* (const double& scalar);
		virtual void state_vector_right_multiply(const double* initial_perturbation_vector, double* final_perturbation_vector);
		virtual STM& operator*= (const STM& OtherSTM);
		virtual STM& operator*= (const double& scalar);
		
		//comparisons
		virtual bool operator== (const STM& OtherSTM);
		virtual bool operator!= (const STM& OtherSTM);

		//special matrix math
		virtual STM transpose();
		virtual STM inverse();

		//fields
	protected:
		double data[36];

	};//end class STM

}//end namespace Kepler

#endif //KEPLER_STM_CLASS