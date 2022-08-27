#include <iostream>
#include <math.h>
 
using namespace std;
 
class Complex
{
public:
    Complex(double real = 0.0,double image = 0.0){
        m_real = real; 
		m_image = image;
    }
 
    Complex operator+(const Complex &c){   	
        return Complex(this->m_real + c.m_real, this->m_image + c.m_image);
    }
    
    Complex operator-(const Complex &c){
    	return Complex(this->m_real - c.m_real, this->m_image - c.m_image);
	}
	
	Complex operator*(const Complex &c){
    	return Complex(this->m_real * c.m_real - this->m_image * c.m_image, 
		this->m_real * c.m_image + this->m_image * c.m_real);
	}
	
	Complex operator/(const Complex &c){
    	return Complex((this->m_real * c.m_real + this->m_image * c.m_image) / 
		(c.m_real * c.m_real + c.m_image * c.m_image), 
		(this->m_image * c.m_real - this->m_real * c.m_image) 
		/ (c.m_real * c.m_real + c.m_image * c.m_image));
	}
	
	Complex Sin(){
		return Complex(((exp(this->m_image) + exp(-this->m_image)) / 2) * sin(this->m_real), 
		((exp(this->m_image) - exp(-this->m_image)) / 2) * cos(this->m_real));
	}
	
	Complex Cos(){
		return Complex(((exp(this->m_image) + exp(-this->m_image)) / 2) * cos(this->m_real), 
		- ((exp(this->m_image) - exp(-this->m_image)) / 2) * sin(this->m_real));
	}
	
	double Modu(){
		return sqrt(this->m_real * this->m_real + this->m_image * this->m_image);
	}
	
	double Real(){
		return this->m_real;
	}
	
	double Image(){
		return this->m_image;
	} 
 
    friend ostream &operator<<(ostream &os, const Complex &c){
        os << c.m_real << " + " << c.m_image << "i";
        return os;
    }
 
    void print(Complex &c);

private:
    double m_real;
    double m_image;
};
 
