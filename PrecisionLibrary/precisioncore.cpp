/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2007-2021
 *                       Henrik Vestermark
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Henrik Vestermark Software License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :precisioncore.cpp
 * Module ID Nbr   :   
 * Description     :Arbitrary precision core functions for integer and floating
 *					point precision class
 * -----------------------------------------------------------------------------
 * Change Record   :   
 *
 * Version	Author/Date		Description of changes
 * -------  ---------------	----------------------------------------------------
 * 03.01	HVE/14-Aug-2021	Switch to new internal binary format iptype & fptype
 * 03.02	HVE/19-Oct-2021	Passing all internal testing
 * 03.03	HVE/19-Nov-2021	Fixed compiler bugs reported by GNU version on a mac
 * 03.04	HVE/20-Nov-2021 A few bugs fixed and change to avoid compiler warnings
 * 03.05	HVE/21-Nov-2021 More cleaning up and improvements
 * 03.06	HVE/22-Nov-2021	Fix a bug in _float_precision_fptoa() where fraction was incorret truncated
 *							Also fix a bug in the modf().
 * 03.07	HVE/26-Nov-2021 Added 1/sqrt(2) and sqrt(2) to the "constant" _float_table
 * 03.08	HVE/5-Dec-2021	Improved the algorithm for float_precion_ftoa() that improved performanen for a 100K digits varianble
 *							from more than 21000sec -> 300msec.
 * 03.09	HVE/6-Dec-2021	Fix an float_precision issue where numbers with small precision but high negative exponent was converted to 
 *							0E0 instead of a very small number. Also added gigatrunk splitting in float_precision_fptoa 
 * 03.10	HVE/8-Dec-2021	Change eptype to an intmax_t insead of int to raise the exponent limit to more than 300M digits
 * 03.11	HVE/10-Dec-2021	Added two more trunking levels for faster handling of digits exceedig 10-100M decimal digits
 * 03.12	HVE/11-Dec-2021	Allowed separaors ' symbols to be part of a string based number for atoip() and atofp(). atofp() can be called with a std::string or a char *
 * 03.13	HVE/25-Dec-2021	remove of double *a,*b with new and replaced it with vector<double> va, vb. 
 * 03.14	HVE/25-Dec-2021 Added _float_precision_schonhagen_strassen_umul() to add better multiplication for medium size multiplication of digits<6,000 the function was modified 
 *							from the int_precision counterpart
 * 03.15	HVE/26-Dec-2021	Change the FFT _int_precision_umul_fourier() and _float_precision_umul_fourier() to be able to handle digits in excess of 25E9 digits 
 * 03.16	HVE/29-Dec-2021 Fixed a few smaller bugs in _flot_precision_fptoa()
 * 03.17	HVE/3-Jan-2022	Special test version plus a fix in _umul_fourier() where the 4bit version was kicked in aftter 18M digits instead of >150M digits
 * 03.18	HVE/3-Jan-2022	A bug was found in the float_precision_schonhage-straasen_umul() and was change to just a call the _umul_fourier() as a temporary fix
 * 03.19	HVE/6-Jan-2022	Fixed the overflow bug in _int_precision_umul_fourier() & _float_precision_umul_fourier() that happens only for very large multiplications > 10M digits
 * 03.20	HVE/6-Jan-2022	Added _INVSQRT3 and SQRT3 as a build in constant
 * 03.21	HVE/8-Jan-2022	Added Multi threading in _int_multiplication_umul_fourier(), _float_precision_umul_fourier() and float_table(_PI)
 * 03.22	HVE/10-Jan-2022	Minor optimazation when multiply by power of 2. It is faster to do x.exponent(x.exponent()+'power of two');
 * 03.23	HVE/19-Jan-2022	Change name of schoonhage-strassen to just _umul_linear to followed the name convention. change _umul to _umul_school nd added a new _umul as the entry point 
 *							for all multiplication algorithm
 * 03.24	HVE/21-Jan-2022	Use string2number in the _float_precision_atofp for higher performance, particular when digits exceed 100,000+ digits
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VIP_[] = "@(#)precisioncore.cpp 03.24 -- Copyright (C) Henrik Vestermark";

#include <cstdint>
#include <ctime>
#include <cmath> 
#include <iostream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <algorithm>
#include <thread>

#include "iprecision.h"
#include "fprecision.h"

// #define HVE_DEBUG
#define HVE_THREAD				// Add Multi threading for _int_precison_umul_fourier(), _float_precison_umul_fourier() and PI calculations

///////////////////////////////////////////////
//
//
//    Integer Precision Core
//
//
///////////////////////////////////////////////


///////////////////////////////////////////////
//
//
//    Integer Precision Input and output operator
//
//
///////////////////////////////////////////////

std::ostream& operator<<( std::ostream& strm, const int_precision& d ) 
   { return strm << _int_precision_itoa(const_cast<int_precision *>(&d) ).c_str(); }

std::istream& operator>>( std::istream& strm, int_precision& d )
         { 
         char ch; std::string s;
         strm.get(ch);// strm >> ch; 
         while( ch == ' ' ) strm.get(ch);  // Ignore leading white space.
         if( ch == '+' || ch == '-' ) { s += ch; strm.get(ch); } else s += '+';  // Parse sign

         if( ch == '0' ) // Octal, Binary or Hexadecimal number
            {
            strm.get( ch );
            if( ch == 'x' || ch == 'X' ) // Parse Hexadecimal
               for( s += "0x"; (ch >= '0' && ch <= '9') || (ch >='a' && ch <= 'f') || (ch >= 'A' && ch <= 'F'); strm.get( ch ) ) s += ch;
            else
               if( ch == 'b' || ch == 'B' )  // Parse Binary
                  for( s += "0b"; ch >= '0' && ch <= '1'; strm.get( ch ) ) s += ch;
               else // Parse Octal
                  for( s += "0"; ch >= '0' && ch <= '7'; strm.get( ch ) ) s += ch;
            }
         else // Parse Decimal number
            for( ; ch >= '0' && ch <= '9'; strm.get(ch) /*strm >> ch*/ ) s += ch;

         strm.putback( ch );  // ch contains the first character not part of the number, so put it back
         if(!strm.fail() && s.length() >= 2 )  // Valid number has at least a length of 2 or higher
            d = int_precision( const_cast<char *>( s.c_str() ) );

         return strm;
         }


///////////////////////////////////////////////
//
//
//    Miscellaneous
//
//
///////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 		std::string _int_precision_itoa Convert number to ascii string
//	@return		std::string -	the converted number in ascii string format
//	@param		"a"	-	Number to convert to ascii
//	@param		"base"	-	base for conversion to ascii
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string
//   using base. Sign is only added if negative
//	  Base is default BASE_10 but can be anything from base 2..36
//
std::string _int_precision_itoa( int_precision *a, const int base )
   {
   return ( a->sign() < 0 ? "-": "" ) + _int_precision_itoa( a->pointer(), base );
   }



///////////////////////////////////////////////
//
//
//    Core Support Functions. 
//
//
///////////////////////////////////////////////

//
// Core functions
// The core functions all perform unsigned arithmetic un elements of the string class!
//    _int_precision_strip_leading_zeros	-- Strips non significant leading zeros
//    _int_precision_uadd_short			-- add a short digit [0..RADIX] to the string
//    _vector_reverse_binary			-- Reverse bit in the data buffer
//    _vector_fourier					-- Fourier transformn the data
//    _vector_real _fourier				 -- Convert n discrete double data into a fourier transform data set
//

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_reverse_binary
//	@return 	void	-	
//	@param   "data[]"	-	array of double complex number to permute
//	@param   "n"	-	number of element in data[]
//
//	@todo  
//
// Description:
//   Reverse binary permute
//   n must be a power of 2
//
/*
static void _int_reverse_binary( std::complex<double> data[], const size_t n )
   {
   size_t i, j, m;

   if( n <=2 ) return;

   for( j=1, i=1; i < n; i++ )
      {
      if( j > i ) 
         std::swap( data[ j - 1 ], data[ i - 1 ] );
 
      for( m = n >> 1; m >= 2 && j > m; m >>= 1 )
         j -= m;

      j += m;
      }
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_fourier do the fourier transformation
//	@return 	static void	-	
//	@param   "data[]"	-	complex<double> fourie data
//	@param   "n"	-	number of element in data (must be a power of 2)
//	@param   "isign"	-	transform in(1) or out(-1)
//
//	@todo
//
// Description:
//   Wk=exp(2* PI *i *j )  j=0..n/2-1
//   exp( k * i *j ) => exp( k * i * (j-1) + k * i ) => exp( t + o ) for each new step
//   exp( t + o ) => exp(t)-exp(t)*( 1 - cos(o) -isin(o) ) => exp(t)-exp(t)*(a-ib)
//   => exp(t)+exp(t)*(-a+ib) => exp(t)( 1 + (-a+b) )
//   sin(t+o)=sin(t)+[-a*sin(t)+b*cos(t)]
//   a=2sin^2(o/2), b=sin(o)
//   n must be a power of 2
//
/*
static void _int_fourier( std::complex<double> data[], const size_t n, const int isign )
   {
   double theta;
   std::complex<double> w, wp;
   size_t mh, m, r, j, i;

   _int_reverse_binary( data, n );

   for( m = 2; n >= m; m <<= 1 )
      {
      theta = isign * 2 * 3.14159265358979323846264 / m;
      wp = std::complex<double>( -2.0 * sin( 0.5 * theta ) * sin( 0.5 * theta ), sin( theta ) );
      w = std::complex<double> ( 1, 0 ); // exp(0) == exp( isign*2*PI*i/mmax * m-1 )
      mh = m >> 1;

      for( j = 0; j < mh; j++ )      // m/2 iteration
         {
         for( r = 0; r <= n - m; r += m )
            {
            std::complex<double> tempc;
            i = r + j;
            tempc = w * data[ i + mh ];              // u=data[i]; v=data[j]*w; data[i]=u+v;data[j]=u-v;
            data[ i + mh ] = data[ i ] - tempc;
            data[ i ] += tempc;
            }
      
         w =  w * wp + w;  // w = w(1+wp) ==> w *=1+wp;
         }
      }
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 			_int_real_fourier
//	@return 			static void	-	
//	@param   "data[]"	-	
//	@param   "n"	-	number of data element in data. n must be a power of 2)
//	@param   "isign"	-	Converting in(1) or out(-1)
//
//	@todo  
//
// Description:
//   Convert n discrete double data into a fourier transform data set
//   n must be a power of 2
//
/*
static void _int_real_fourier( double data[], const size_t n, const int isign )
   {
   size_t i;
   double theta, c1 = 0.5, c2;
   std::complex<double> w, wp, h1, h2;

   theta = 3.14159265358979323846264 / (double)( n >> 1 );
   if( isign == 1 )
      {
      c2 = -c1;
      _int_fourier( (std::complex<double> *)data, n >> 1, 1 );
      }
   else
      {
      c2 = c1;
      theta = -theta;
      }
   wp = std::complex<double>( -2.0 * sin( 0.5 * theta ) * sin( 0.5 * theta ), sin( theta ) );
   w = std::complex<double> ( 1 + wp.real(), wp.imag() );
   for( i = 1; i < (n>>2); i++ )
      {
      size_t i1, i2, i3, i4;
      std::complex<double> tc;

      i1 = i + i;
      i2 = i1 + 1;
      i3 = n + 1 - i2;
      i4 = i3 + 1;
      h1 = std::complex<double> ( c1 * ( data[i1] + data[i3] ), c1 * ( data[i2]-data[i4]));
      h2 = std::complex<double> ( -c2 * ( data[i2]+data[i4] ), c2 * ( data[i1]-data[i3]));
      tc = w * h2;
      data[i1]=h1.real()+tc.real();
      data[i2]=h1.imag()+tc.imag();
      data[i3]=h1.real() - tc.real();
      data[i4]=-h1.imag() + tc.imag();
      w *= ( std::complex<double>(1) + wp );
      }
   if( isign == 1 )
      {
      double t;
      data[0] = (t=data[0]) + data[1];
      data[1] = t - data[1];
      }
   else
      {
      double t;
      data[0]=c1*((t=data[0])+data[1]);
      data[1]=c1*(t-data[1]);
      _int_fourier( (std::complex<double> *)data, n>>1, -1 );
      }
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_reverse_binary
//	@return 	void	-	
//	@param   "data[]"	-	array of double complex number to permute
//	@param   "n"	-	number of element in data[]
//
//	@todo  
//
// Description:
//   Reverse binary permute
//   n must be a power of 2
//
static void _vector_reverse_binary(std::complex<double> data[], const size_t n)
	{
	size_t i, j, m;

	if (n <= 2) return;

	for (j = 1, i = 1; i < n; i++)
	{
		if (j > i)
		{
			//if (j - 1 >= n || i - 1 >= n)  // Debug
			//	j = j;// DEBUG Error
			std::swap(data[j - 1], data[i - 1]);
		}

		for (m = n >> 1; m >= 2 && j > m; m >>= 1)
			j -= m;

		j += m;
	}
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 	_int_fourier do the fourier transformation
//	@return 	static void	-	
//	@param   "data[]"	-	complex<double> fourie data
//	@param   "n"	-	number of element in data (must be a power of 2)
//	@param   "isign"	-	transform in(1) or out(-1)
//
//	@todo
//
// Description:
//   Wk=exp(2* PI *i *j )  j=0..n/2-1
//   exp( k * i *j ) => exp( k * i * (j-1) + k * i ) => exp( t + o ) for each new step
//   exp( t + o ) => exp(t)-exp(t)*( 1 - cos(o) -isin(o) ) => exp(t)-exp(t)*(a-ib)
//   => exp(t)+exp(t)*(-a+ib) => exp(t)( 1 + (-a+b) )
//   sin(t+o)=sin(t)+[-a*sin(t)+b*cos(t)]
//   a=2sin^2(o/2), b=sin(o)
//   n must be a power of 2
//
static void _vector_fourier(std::complex<double> data[], const size_t n, const int isign)
	{
	double theta;
	std::complex<double> w, wp;
	size_t mh, m, r, j, i;

	_vector_reverse_binary(data, n);

	for (m = 2; n >= m; m <<= 1)
		{
		theta = isign * 2 * 3.14159265358979323846264 / m;
		wp = std::complex<double>(-2.0 * sin(0.5 * theta) * sin(0.5 * theta), sin(theta));
		w = std::complex<double>(1, 0); // exp(0) == exp( isign*2*PI*i/mmax * m-1 )
		mh = m >> 1;

		for (j = 0; j < mh; j++)      // m/2 iteration
			{
			for (r = 0; r <= n - m; r += m)
				{
				std::complex<double> tempc;
				i = r + j;
				//if (i + mh >= n) //DEBUG
				//	i = i;  // DEBUG
				tempc = w * data[i + mh];              // u=data[i]; v=data[j]*w; data[i]=u+v;data[j]=u-v;
				data[i + mh] = data[i] - tempc;
				data[i] += tempc;
				}

			w = w * wp + w;  // w = w(1+wp) ==> w *=1+wp;
			}
		}
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/19/2005
//	@brief 			_int_real_fourier
//	@return 			static void	-	
//	@param   "data[]"	-	
//	@param   "n"	-	number of data element in data. n must be a power of 2)
//	@param   "isign"	-	Converting in(1) or out(-1)
//
//	@todo  
//
// Description:
//   Convert n discrete double data into a fourier transform data set
//   n must be a power of 2
//
static void _vector_real_fourier(std::vector<double>& data, const size_t n, const int isign)
	{
	size_t i;
	double theta, c1 = 0.5, c2;
	std::complex<double> w, wp, h1, h2;
	double *ptr = data.data() ;

	theta = 3.14159265358979323846264 / (double)(n >> 1);
	if (isign == 1)
		{
		c2 = -c1;
		_vector_fourier((std::complex<double> *)ptr, n>>1, 1);  // n>> 1 == complex numbers
		}
	else
		{
		c2 = c1;
		theta = -theta;
		}
	wp = std::complex<double>(-2.0 * sin(0.5 * theta) * sin(0.5 * theta), sin(theta));
	w = std::complex<double>(1 + wp.real(), wp.imag());
	for (i = 1; i < (n >> 2); i++)
		{
		size_t i1, i2, i3, i4;
		std::complex<double> tc;

		i1 = i + i;
		i2 = i1 + 1;
		i3 = n + 1 - i2;
		i4 = i3 + 1;
		h1 = std::complex<double>(c1 * (data[i1] + data[i3]), c1 * (data[i2] - data[i4]));
		h2 = std::complex<double>(-c2 * (data[i2] + data[i4]), c2 * (data[i1] - data[i3]));
		tc = w * h2;
		data[i1] = h1.real() + tc.real();
		data[i2] = h1.imag() + tc.imag();
		data[i3] = h1.real() - tc.real();
		data[i4] = -h1.imag() + tc.imag();
		w *= (std::complex<double>(1) + wp);
		}
	if (isign == 1)
		{
		double t;
		data[0] = (t = data[0]) + data[1];
		data[1] = t - data[1];
		}
	else
		{
		double t;
		data[0] = c1*((t = data[0]) + data[1]);
		data[1] = c1*(t - data[1]);
		_vector_fourier((std::complex<double> *)ptr, n>>1, -1);
		}
	}


///////////////////////////////////////////////
//
//
//    Core Functions. BINARY
//
// The core functions all perform unsigned arithmetic of elements of the vector<iptype> class!
//
//	  build_i_number					--	Build integer number as vector<iptype>
//    _int_precision_strip_leading_zeros	-- Strips non significant leading zeros
//    _int_precision_strip_trailing_zeros	-- Strips non significant trailing zeros
//	  _int_precision_clz				-- Count leading zeros in an iptype
//	  _int_precision_clz				-- Count leading zeros in an vector<iptype>
//	  _int_precision_ctz				-- Count trailing zeros in an iptype
//	  _int_precision_ctz				-- Count trailing zeros in an vector<iptype>
//	  _int_precision_csb				-- Bit position of the most significant bit in vector<iptype>
//    _int_precision_compare			-- Compare two strings for numeric order
//    _int_precision_uneg				-- Negate ones-complement of unsigned integer
//    _int_precision_uadd_short			-- Add a short iptype digit e.g. [0..2^64] to the number
//    _int_precision_uadd				-- Add two unsigned binary numbers
//    _int_precision_usub_short			-- Subtract a short iptype digit [0..2^64] from the number
//    _int_precision_usub				-- Subtract two unsigned binary numbers
//    _int_precision_umul_short			-- Multiply a iptype digit [0..2^64] to the number
//	  _int_precision_umul				-- Multiply two unsigned binary numbers using the most optimal multiplication algorithm
//    _int_precision_umul_school		-- Multiply two unsigned binary numbers using old fashion school algorithm
//    _int_precision_udiv_short			-- Divide a iptype digit [0..^64] up in the number
//    _int_precision_udiv				-- Divide two unsinged strings
//    _int_precision_urem				-- Remainder of dividing two unsigned numbers
//	  _int_precision_urem_short			-- Get the remainder an iptype digit[0.. ^ 64] up in the number
//	  _int_precision_shiftright			-- Shift right the vector<iptype> number
//	  _int_precision_shiftleft			-- Shift left the vector<iptype> number
//    _int_precision_itoa				-- Convert internal precision to BASE_10 string
//	  _int_precision_uand				-- And the binary numbers together
//	  _int_precision_uor				-- Or the binary numbers together
//	  _int_precision_xor				-- Xor the beinary numbers together
//	  _int_precision_unegate			-- Negate the binary number
//    _build_i_number					-- Build the binary number from string representation
//	
//	  _vector_reverse_binary			-- Reverse bit in the data buffer
//    _vector_fourier					-- Fourier transformn the data
//    _vector_real_fourier				-- Convert n discrete double data into a fourier transform data set
//	  _precision_umul64					-- Generic multiplication of two vector<iptype> or vector<fptype> numbers
//	  _int_precision_umul_fourier		-- Multiply two binary numbers using FFT
//	  _int_precision_umul_karatsuba		-- Multiply two unsigned binary numbers using karatsuba algorithm
//	  _int_precision_umul_linear		-- Multiply two unsigned binary numbers using schonhagen-strassen algorithm with linearconvolution
//
//
///////////////////////////////////////////////

// Use for various conversions to and from strings
static size_t _powerof10Table[20] = { 1,10,100,1000,10 * 1000,100 * 1000,1000 * 1000,10 * 1000000,100 * 1000000,1000 * 1000000,
									10000ull * 1000000,100000ull * 1000000,1000000ull * 1000000,10ull * 1000000 * 1000000,100ull * 1000000 * 1000000,
									1000ull * 1000000 * 1000000, 10000ull * 1000000 * 1000000, 100000ull * 1000000 * 1000000,1000000ull * 1000000 * 1000000,
									1000000ull * 1000000 * 1000000 * 10 };


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief		Build a binary repesentation of signle digit number
//	@return		std::vector<iptype>	- The integer precision string	
//	@param		"digit"		- The next digit to be added to the integer point number to convert. Is always positive
// @param		"base"		- The base of the digit being added
//
//	@todo 	
//
// Description:
//   Add a digit to the number being build for the integer precision number
//   The function dosnt create any leading significant zeros
//   To run it efficiently is is better to take advantages of iptype (64bit) instead of just one decimal,
//		binary,octal or hexdecimal digit at a time
//    
static inline std::vector<iptype> build_i_number(std::vector<iptype> &number, iptype digit, iptype base)
	{
	number = _int_precision_umul_short(&number, base);
	number = _int_precision_uadd_short(&number, digit);
	return number;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Dec/2021
//	@brief		Build a binary repesentation of a string of single digit number
//	@return		std::vector<iptype>	- The integer precision string
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect th binary number that correspond to the decimal representation of the number
//
static std::vector<iptype> decimal2number(std::vector<iptype>& number, std::string s, size_t start, size_t end)
	{
	size_t i=start, length=end-start;
	const size_t max_digits = 19;
	iptype pwr;

	for(i=start;length>=max_digits;length-=max_digits, i+=max_digits)
		{
		std::string s2 = s.substr(i, max_digits);
		uint64_t n = strtoull(s2.c_str(), NULL, BASE_10);
		pwr = _powerof10Table[max_digits];
		build_i_number(number, n, pwr);
		}
	
	if (length!= 0)
		{
		std::string s2 = s.substr(i, length);
		uint64_t n = strtoull(s2.c_str(), NULL, BASE_10);
		pwr = _powerof10Table[length];
		build_i_number(number, n, pwr);
		}

	return number;
	}

// Do  a single trunk and return it
// It is guarantee that there is trunk size of data
static std::vector<iptype> trunk2number(std::vector<iptype>& number, std::string s, size_t start)
	{
	size_t i;
	std::vector<iptype> sum(1,0);
	size_t thr = MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS;
	static std::vector<iptype> _trunkPowerof10;

	if ( _trunkPowerof10.size()==0) // is _trunkPowerof10 build or created
		{
		std::vector<iptype> p(1,_powerof10Table[MAX_DECIMAL_DIGITS]);
		for (i = MAX_TRUNK_SIZE, _trunkPowerof10.assign(1,1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 =_int_precision_umul(&_trunkPowerof10, &p);	// Odd
			if (i > 1) p = _int_precision_umul(&p, &p);					// square it
			}
		}

	decimal2number(sum, s, start, start + thr);
	number =_int_precision_umul(&number, &_trunkPowerof10);
	number =_int_precision_uadd(&number, &sum );

	return number;
	}

#if false
// Do a single kilo trunk and return it
// It is guarantee that there is trunk size of data
static std::vector<iptype> kilo2number(std::vector<iptype>& number, std::string s, size_t start)
	{
	size_t i;
	std::vector<iptype> sum(1, 0);
	size_t thr = MAX_KILOTRUNK_SIZE;
	static std::vector<iptype> _trunkPowerof10;

	if (_trunkPowerof10.size() == 0) // is _trunkPowerof10 build or created
		{
		std::vector<iptype> p(1, _powerof10Table[MAX_DECIMAL_DIGITS]);
		for (i = MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _trunkPowerof10.assign(1, 1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 = _int_precision_umul_fourier(&_trunkPowerof10, &p);	// Odd
			if (i > 1) p = _int_precision_umul_fourier(&p, &p);					// square it
			}
		}

	for (i=start;thr>0; --thr, i+=MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
		{
		trunk2number(sum, s, i);
		number = _int_precision_umul_fourier(&number, &_trunkPowerof10);
		number = _int_precision_uadd(&number, &sum);
		sum.assign(1, 0);
		}

	return number;
	}
#endif


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit decimal numbers
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//
/*static*/ std::vector<iptype> string2number(const std::string s, const size_t start, size_t end)
	{
	size_t i = end - start, j, radix_inx=0;
	std::string s2;
	std::vector<std::vector<iptype> > vn(0);
	std::vector<iptype> num(1);
	static std::vector<std::vector<iptype> > radix;

	vn.reserve(i/MAX_DECIMAL_DIGITS+16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > MAX_DECIMAL_DIGITS; i -= MAX_DECIMAL_DIGITS, end-= MAX_DECIMAL_DIGITS)
		{
		s2 = s.substr(end - MAX_DECIMAL_DIGITS, MAX_DECIMAL_DIGITS);
		num[0] = strtoull(s2.c_str(), NULL, BASE_10);
		vn.push_back(num);
		}
	s2 = s.substr(start, i);
	num[0] = strtoull(s2.c_str(), NULL, BASE_10);
	vn.push_back(num);

	// Step2 collected into higher binary values by reducing the vector with 2,3,...,n MAX_DECIMAL_DIGITS
	if (radix.size() == 0)
		{
		num[0] = _powerof10Table[MAX_DECIMAL_DIGITS];
		radix.push_back(num);
		}
	for (;vn.size() > 1; ++radix_inx )
		{
		if (radix_inx >= radix.size())
			{
			num = _int_precision_umul(&radix[radix_inx - 1], &radix[radix_inx - 1]); // replace by _int_precision_square_fourier() when ready
			radix.push_back(num);  
			}
		for (i = 0, j = 0; j < vn.size(); ++i, j += 2)
			{
			if (j + 1 < vn.size())
				{
				num = _int_precision_umul(&vn[j + 1], &radix[radix_inx]);
				vn[i] = _int_precision_uadd(&vn[j], &num);
				}
			else
				vn[i] = vn[j];
			}
		vn.resize(i);
		}

	return vn[0];
 	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		21/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit numbers in base  or base 16
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//	@param		"base"		- The base of the string number (either base 2 or base 16)
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//	Base is either base 2 or base 16
//
static std::vector<iptype> stringbase2number(const std::string s, const size_t start, size_t end, int base)
	{
	const uintmax_t baselength = base == BASE_2 ? 64 : 16;
	size_t i = end - start;
	std::string s2;
	std::vector<iptype> vn;
	iptype n;

	vn.reserve(i / baselength + 16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > baselength; i -= baselength, end -= baselength)
		{
		s2 = s.substr(end - baselength, baselength);
		n = strtoull(s2.c_str(), NULL, base);
		vn.push_back(n);
		}
	s2 = s.substr(start, i);
	n= strtoull(s2.c_str(), NULL, base);
	vn.push_back(n);

	return vn;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Jan/2022
//	@brief		Build a binary repesentation of a string of single digit numbers in base 8 (octal)
//	@return		std::vector<iptype>	- The integer precision number
//  @param		"s"			- The decimal string 
//	@param		"start"		- The start index of the string 
//  @param		"end"		- The End index of the string
//
//	@todo 	
//
// Description:
//	Build and collect the binary number that correspond to the decimal representation of the number
//	Currently Only base 8 is supproted
//
static std::vector<iptype> stringbase8_2number(const std::string s, const size_t start, size_t end)
	{
	size_t i = end - start, j, radix_inx = 0;
	std::string s2;
	std::vector<std::vector<iptype> > vn(0);
	std::vector<iptype> num(1);
	static std::vector<std::vector<iptype> > radix;

	vn.reserve(i / MAX_OCTAL_DIGITS + 16);
	// Step 1 partition the string into a binary vector with 1 binary digit in order of least to most significant 
	for (; i > MAX_OCTAL_DIGITS; i -= MAX_OCTAL_DIGITS, end -= MAX_OCTAL_DIGITS)
		{
		s2 = s.substr(end - MAX_OCTAL_DIGITS, MAX_OCTAL_DIGITS);
		num[0] = strtoull(s2.c_str(), NULL, BASE_8);
		vn.push_back(num);
		}
	s2 = s.substr(start, i);
	num[0] = strtoull(s2.c_str(), NULL, BASE_8);
	vn.push_back(num);

	// Step2 collected into higher binary values by reducing the vector with 2,3,...,n MAX_OCTAL_DIGITS
	if (radix.size() == 0)
		{
		num[0] = 01'000'000'000'000'000'000'000ull;
		radix.push_back(num);
		}
	for (; vn.size() > 1; ++radix_inx)
		{
		if (radix_inx >= radix.size())
			{
			num = _int_precision_umul(&radix[radix_inx - 1], &radix[radix_inx - 1]); // replace by _int_precision_square_fourier() when ready
			radix.push_back(num);
			}
		for (i = 0, j = 0; j < vn.size(); ++i, j += 2)
			{
			if (j + 1 < vn.size())
				{
				num = _int_precision_umul(&vn[j + 1], &radix[radix_inx]);
				vn[i] = _int_precision_uadd(&vn[j], &num);
				}
			else
				vn[i] = vn[j];
			}
		vn.resize(i);
		}

	return vn[0];
	}



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		_int_precision_strip_leading_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove leading nosignificant zeros of the binary number
//	This is from the start of the vector<iptype>
//
void _int_precision_strip_leading_zeros(std::vector<iptype> *s)
	{
	std::vector<iptype>::iterator pos;

	// Strip leading zeros
	for (pos = s->begin(); pos != s->end() && *pos == (iptype)0; ++pos);	// Find first not zero digit

	if (s->begin() != pos)
		s->erase(s->begin(), pos);
	if (s->empty())
		s->assign(1,(iptype)0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_strip_trailing_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove trailing nosignificant zeros of the binary number
//		this is from the top of the vector<iptype> 
//
void _int_precision_strip_trailing_zeros(std::vector<iptype> *s)
	{
	size_t i;
	std::vector<iptype>::reverse_iterator pos;

	// Strip leading zeros
	for (i = s->size() - 1, pos=s->rbegin(); i > 0 && *pos == (iptype)0; ++pos, --i);

	s->resize(i + 1);  // Keep at least one digit by default

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"a"	-	iptype operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary iptype 
//
size_t _int_precision_clz(const iptype a)
	{
	iptype x = a;
	size_t offset = 0;
	static const unsigned char lookup[16] = { 4,3,2,2,1,1,1,1,0,0,0,0,0,0 };

	if (sizeof(iptype) > 4 && x & 0xffffffff00000000u)
		x >>= 32; else offset += 32;

	if (sizeof(iptype) > 2 && x & 0xffff0000u)
		x >>= 16; else offset += 16;

	if (sizeof(iptype) > 1 && x & 0xff00u)
		x >>= 8; else offset += 8;

	if ( x & 0xf0u)
		x >>= 4; else offset += 4;
	offset += lookup[(unsigned char)x];
	return offset;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_int_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"mb"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary iptype 
//
size_t _int_precision_clz(const std::vector<iptype> &mb ) 
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = mb.size(); i > 0; --i)
		{
		cnt = _int_precision_clz(mb[i - 1]);
		tot_cnt += cnt;
		if (cnt != 64) break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Aug/2021
//	@brief 		_int_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"a"	-	iptype operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary iptype 
//	  iptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _int_precision_ctz(const iptype a)
	{
	iptype x = a;  // sizeof(iptype) can be 8, 4, 2, or 1 
	size_t cnt;
	if (x & 0x1) 
		cnt = 0;
	else
		{
		cnt = 1;
		if (sizeof(iptype) >= 8 && (x & 0xffffffffu) == 0)
			{
			x >>= 32; cnt += 32;
			}
		if (sizeof(iptype) >= 4 && (x & 0xffffu) == 0)
			{
			x >>= 16; cnt += 16;
			}
		if (sizeof(iptype) >= 2 && (x & 0xff) == 0)
			{
			x >>= 8; cnt += 8;
			}
		if ((x & 0xf) == 0)
			{
			x >>= 4; cnt += 4;
			}
		if ((x & 0x3) == 0)
			{
			x >>= 2; cnt += 2;
			}
		cnt -= x & 0x1;
		}
	
	return cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_int_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"mb"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary iptype 
//	  iptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _int_precision_ctz( const std::vector<iptype> &mb )
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = 0; i < mb.size(); ++i)
		{
		cnt = _int_precision_ctz(mb[i]);
		tot_cnt += cnt;
		if (cnt != 64) break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		_int_precision_csb
//	@return		size_t	-	the number of significant bits in vector<iptype>
//	@param		"a"	-	vector<iptype> operand
//
//	@todo
//
// Description:
//   Count number of significant bits in a vector<iptype> number
//
size_t _int_precision_csb(const std::vector<iptype> &a)
	{
	size_t bit_pos = 0, cnt;
	for (size_t i = a.size(); i > 0; --i)
		{
		if (a[i - 1] == 0) continue;
		cnt = _int_precision_clz(a[i - 1]);
		bit_pos = Bitsiptype - cnt;
		if (i - 1 > 0) 
			bit_pos += Bitsiptype * (i - 1);
		break;
		}
	return bit_pos;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		_int_precision_compare
//	@return 	int	-	The compare result. 0==equal, 1==s1>s2, -1==s1<s2
//	@param		"s1"	-	First operand to compare
//	@param		"s2"	-	Second operand to compare
//
//	@todo  
//
// Description:
//   Compare two unsigned vector<iptype> binary numbers 
//   and return 0 is equal, 1 if s1 > s2 otherwise -1
//   Optimized check length first and determine 1 or -1 if equal
//   compare the digits until a determination can be made.
//
int _int_precision_compare(const std::vector<iptype> *s1, const std::vector<iptype> *s2)
	{
	std::vector<iptype>::reverse_iterator s1_pos, s2_pos;
	int cmp; size_t i;

	if (s1->size() > s2->size())
		cmp = 1;
	else
		if (s1->size() < s2->size())
			cmp = -1;
		else
			{// Same size.
			s1_pos = const_cast<std::vector<iptype> *>(s1)->rbegin(); s2_pos = const_cast<std::vector<iptype> *>(s2)->rbegin();
			for (cmp=0, i = s1->size(); i > 0; --i, ++s1_pos, ++s2_pos)
				{
				if (*s1_pos == *s2_pos) continue;
				if (*s1_pos > *s2_pos) 
					{ cmp = 1; break; }
				else 
					{ cmp = -1; break; }
				//if (i == 0) break;
				}
			}

	return cmp;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		20/Sep/2021
//	@brief 		std::vector<_TY> _precision_uadd64
//	@return 	std::vector<_TY> - 	the result of the addition. 2 dimensional vector
//	@param      "a"		-	_TY operand a
//	@param      "b"	   -	_TY operand b   
//
//	@todo
//
// Description:
//	 Generic Addition function for two vector<iptype> or vector<fptype> numbers
//  Add two unsigned iptype numbers togeher and return the result as a vector<iptype> [2] or vctor<fptype> [2]:
//	 where [0] is the lower operand and [1] is the upper operand
//
template<class _TY> std::vector<_TY> _precision_uadd64(const _TY a, const _TY b)
	{
	std::vector<_TY> res(2);
	res[0] = a + b;
	res[1] = res[0] < a ? 1 : 0;  // Carry
	return res;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uadd_short
//	@return 	std::vector<iptype> - 	the result of the add
//	@param      "src1"	-	Source binary to add short number
//	@param      "d"	   -	Number to add.   
//
//	@todo
//
// Description:
//   Short Add: The digit d of iptype [0..2^64] for iptype=uint64_t is added to the unsigned binary vector 
//   Optimized 0 add or early out add is implemented
//
std::vector<iptype> _int_precision_uadd_short(const std::vector<iptype> *src1, const iptype d)
	{
	iptype carry;
	std::vector<iptype>::const_iterator s_pos;
	std::vector<iptype>::iterator d_pos;
	std::vector<iptype> des;

	if (d == 0)   // Zero add
		return *src1;

	carry = d;
	des = *const_cast<std::vector<iptype> *> (src1);		// Copy source to des1
	d_pos = des.begin();
	s_pos = src1->begin(); 

	for (; carry != 0 && d_pos != des.end(); ++s_pos, ++d_pos)
		{
		*d_pos += carry;
		if (*d_pos < *s_pos ) 
			carry = 1;  // Set Carry
		else carry = 0;
		}

	// Exhaust the smalles of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != des.end(); )
		{
		iptype tmp = *d_pos;
		*d_pos = tmp + carry;
		if (*d_pos < tmp) 
			carry = 1;  // Set Carry
		else carry = 0;
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des.push_back(1); 

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uadd
//	@return		std::vector<iptype>	-	the result of adding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Add two unsigned decimal strings
//   Optimized: Used early out add
//
std::vector<iptype> _int_precision_uadd(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	iptype carry = 0, tmp;
	std::vector<iptype> des1;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	if (src1->size() >= src2->size())
		{
		des1 = *src1;
		pos = src2->begin();  
		end = src2->end(); 
		}
	else
		{
		des1 = *src2;
		pos = src1->begin(); 
		end = src1->end(); 
		}
	d_pos = des1.begin();

	for (; pos != end; )
		{ // Adding element by element for the two numbers
		*d_pos = *pos + *d_pos + carry;
		carry = *d_pos < *pos ? 1 : 0;
		++pos;
		++d_pos;
		}

	// Exhaust the smallest of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != des1.end();  )
		{
		tmp = *d_pos;
		*d_pos = tmp + carry;
		carry = *d_pos < tmp ? 1 : 0; 
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des1.push_back(1);

	_int_precision_strip_trailing_zeros(&des1);

	return des1;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_usub_short
//	@return 	std::vector<iptype> - 	the result of the add
//	@param      "src1"	-	Source string to add short number
//	@param      "d"	   -	iptype Number to add.   
// @param		"result" - Indicated wrap around (1) or not (0)
//
//	@todo
//
// Description:
//   Short subtract: The iptype digit d [0..2^64] is subtracted from the unsigned vector 
//   if src1 < d result is set to -1 (wrap around) otherwise result is set to  0 (no wrap around)
//   Optimized for 0 subtract
std::vector<iptype> _int_precision_usub_short(int *result, const std::vector<iptype> *src1, const iptype d)
	{
	iptype r, borrow=0;
	std::vector<iptype>::const_iterator pos;
	std::vector<iptype> des1;

	if (d == 0) // Nothing to subtract
		{
		*result = 0;
		return *src1;
		}

	des1.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	pos = src1->begin();
	r = *pos - (d + borrow);
	borrow = *pos < (d + borrow) ? 1 : 0;
	des1.push_back(r);
	for (++pos; borrow>0 && pos != src1->end(); ++pos)
		{
		r = *pos - borrow;
		borrow = *pos < borrow ? 1 : 0;
		des1.push_back(r);
		}
	_int_precision_strip_trailing_zeros(&des1);

	*result = borrow > 0 ? -1 : 0;
	return des1;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_usub
//	@return 	std::vector<iptype>	-	the result of subtracting src2 from src1
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
// @param		"result" - Return indicate wrap around (-1) otherwise 0
//
//	@todo
//
// Description:
//   Subtract two unsigned decimal strings
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<iptype> _int_precision_usub(int *result, const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	iptype r, borrow=0;
	std::vector<iptype>::const_iterator pos1, pos2;
	std::vector<iptype> des1;

	if (src1->size() > src2->size())
		des1.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	else
		des1.reserve(src2->capacity());  // Reserver space to avoid time consuming reallocation
	pos1 = src1->begin();
	pos2 = src2->begin();

	for (; pos1 != src1->end() || pos2 != src2->end();)
		{
		if (pos1 != src1->end() && pos2 != src2->end())
			{
			r = *pos1 - (*pos2 + borrow);
			borrow = *pos1 < (*pos2 + borrow) ? 1 : 
					 *pos1==0? borrow : 0;      // if borrow was not paid then propagate it to next iptype subtraction
			++pos1; ++pos2;
			}
		else
			if ( pos1 != src1->end())
				{
				r = *pos1 - (borrow);
				borrow = *pos1 < borrow ? 1 : 0; 
				++pos1;
				}
			else
				{
				r = 0-(*pos2 + borrow);
				//borrow = 0 < (*pos2 + borrow) ? 1 : 0;
				borrow = 1;
				++pos2;
				}
		des1.push_back(r);
		}
	_int_precision_strip_trailing_zeros(&des1);

	*result = borrow>0 ? -1 : 0;
	return des1;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::vector<_TY> _precision_umul64
//	@return 	std::vector<_TY> - 	the result of the multiplication. 2 dimensional vector
//	@param      "a"		-	_TY operand a
//	@param      "b"	   -	_TY operand b   
//
//	@todo
//
// Description:
//	 Generic multiplication function for two vector<iptype> or vector<fptype> numbers
//  Multiply two unsigned iptype numbers togeher and return the result as a vector<iptype> [2]: where [0] is the lower operand and [1] is the upper operand
//
template<class _TY> inline std::vector<_TY> _precision_umul64(const _TY a, const _TY b)
	{
	const _TY mask = 0xffffffff;
	const unsigned int shift = 32;
	const _TY a0 = a & mask, a1 = a >> shift;
	const _TY b0 = b & mask, b1 = b >> shift;
	const _TY a0b0 = a0*b0, a0b1 = a0*b1, a1b0 = a1*b0, a1b1 = a1*b1;
	_TY carry0, carry1, mid;
	std::vector<_TY> res(2);
	mid = a0b1 + a1b0; carry1 = mid < a0b1 ? 1 : 0;
	res[0] = a0b0 + ( ( mid&mask ) << shift ); carry0 = res[0] < a0b0 ? 1 : 0;
	res[1] = (mid >> shift) + a1b1 + (carry1 << shift) + carry0; // no overflow
	return res;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Aug/2021
//	@brief 		std::string _int_precision_umul_short
//	@return 	std::string - 	the result of the short multiplication
//	@param      "src1"	-	Source string to multiply short number
//	@param      "d"	   -	Number to multiply   
//
//	@todo
//
// Description:
//   Short Add: The digit d [0..RADIX] is multiplied to the unsigned decimal string
//   Optimized Multiply with zero yields zero, Multiply with one return the original 
//   
//
std::vector<iptype> _int_precision_umul_short(const std::vector<iptype> *src1, const iptype d)
	{
	iptype carry=0;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype> des;
	std::vector<iptype> tmp(2);

	if (d == 0)  // Multiply by zero is zero.
		{
		des.push_back(0);
		return des;
		}

	if (d == 1)  // Multiply by one dont change the src1.
		{
		des = *const_cast<std::vector<iptype> *> (src1);
		_int_precision_strip_trailing_zeros(&des);
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation   
	pos = src1->begin();
	end = src1->end(); 

	for (; pos != end; ++pos)
		{
		tmp=_precision_umul64(*pos, d);
		tmp[0] += carry;
		carry = (tmp[0] < carry) ? tmp[1]+1 : tmp[1];
		des.push_back(tmp[0]);
		}

	if (carry != 0)
		des.push_back(carry);
	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		19/Jan/2022
//	@brief 		std::vector<iptype>  _int_precision_umul
//	@return		std::vector<iptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<iptype> using the most optimal multiplication algorithm based on operand sizes
//
std::vector<iptype> _int_precision_umul(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	std::vector<iptype> des;
	// Check for multiplication with 1 digit and use umul_short().
	if (src1->size() == 1)
		des = _int_precision_umul_short(src2, *src1->begin());
	else
		if (src2->size() == 1)
			des = _int_precision_umul_short(src1, *src2->begin());
		else
			if (src1->size() + src2->size() < 4000) // Use Schonhage-Strassen for multiplication
				des = _int_precision_umul_linear(src1, src2);
			else // Use FFT for multiplication
				des = _int_precision_umul_fourier(src1, src2);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021
//	@brief 		std::vector<iptype>  _int_precision_umul_school
//	@return		std::vector<iptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<iptype>.
// Not used anymore since the complexity is o(n^2)
//
std::vector<iptype> _int_precision_umul_school(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int disp;
	std::vector<iptype> des, tmp, offset;
	std::vector<iptype>::const_iterator pos2, pos2_end;

	pos2 = src2->begin();
	pos2_end = src2->end();

	des = _int_precision_umul_short(src1, *pos2);
	for (pos2++, disp = 1; pos2 != pos2_end; disp++, pos2++)
		{
		if (*pos2 != 0)
			{
			offset.push_back(0);
			tmp = _int_precision_umul_short(src1, *pos2);

			tmp.insert(tmp.begin(), disp, 0); // = offset + tmp;
			des = _int_precision_uadd(&des, &tmp);
			}
		}

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2double
//	@return		void	-	
//	@param		"dp"	-	pointer to array of doubles
//	@param		"d"		-	the binary iptype or fptype
//	@param	`   "bits"  -	Splitting bits (8 or 4)
//
//	@todo
//
// Description:
// convert an iptype or fptype into an array of doubles[]
//
template<class _TY> static size_t convertbinary2double(double *dp,  _TY d, const bool first, const int bits=8 )
	{
	const unsigned int mask = bits==8 ?0xff : 0xf;
	int k = 0;
	for (int i = sizeof(d) * 8 - bits; i >= 0; i -= bits)
		{
		_TY val = (d >> i) & mask;
		if (first == true && val == 0 && i != 0 && k==0) continue;
		k++;
		*dp++ = (double)(val);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2double
//	@return		void	-	
//	@param		"dp"	-	pointer to array of doubles
//	@param		"d"		-	the binary iptype or fptype
//	@param	`   "bits"  -	Splitting bits (8 or 4)
//
//	@todo
//
// Description:
// convert an iptype or fptype into an array of doubles[]
//
template<class _TY> static size_t convertbinary2double(std::vector<double>::iterator dp, _TY d, const bool first, const int bits = 8)
	{
	const unsigned int mask = bits == 8 ? 0xff : 0xf;
	int k = 0;
	for (int i = sizeof(d) * 8 - bits; i >= 0; i -= bits)
		{
		if (first == true && ((d >> i) & mask) == 0 && i != 0 && k == 0) continue;
		k++;
		*dp++ = (double)((d >> i) & mask);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertbinary2uint8
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned chars
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned bytes. Used in Schonhagen-Strassen
//
static inline size_t convertbinary2uint8(unsigned char *dp, iptype d, const bool first)
	{
	int k = 0;
	for (int i = Bitsiptype - 8; i >= 0; i -= 8)
		{
		if (first == true && ((d >> i) & 0xff) == 0 && i != 0 && k==0) continue;
		k++;
		*dp++ = (unsigned char)((d >> i) & 0xff);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Oct/2021
//	@brief 		void convertbinary2uint16
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned shorts (16bits)
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned shorts (16bit). Used in Schonhagen-Strassen
//
static inline size_t convertbinary2uint16(unsigned short *dp, iptype d, const bool first)
	{
	int k = 0;
	for (int i = Bitsiptype - 16; i >= 0; i -= 16)
		{
		if (first == true && ((d >> i) & 0xffff) == 0 && i != 0 && k == 0) continue;
		k++;
		*dp++ = (unsigned short)((d >> i) & 0xffff);
		}
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Oct/2021
//	@brief 		void convertbinary2Halfiptype
//	@return		void	-	
//	@param		"dp"	-	pointer to array of unsigned shorts (16bits)
//	@param		"d"		-	the binary iptype
//
//	@todo
//
// Description:
// convert an iptype into an array of unsigned shorts (16bit). Used in Schonhagen-Strassen
//
static inline size_t convertbinary2Halfiptype(iptype *dp, iptype d, const bool first)
	{
	const unsigned int HalfBitsiptype = Bitsiptype / 2;
	const iptype mask = (~(iptype)0) >> HalfBitsiptype;
	int k = 0;
	if (first == false || ((d >> HalfBitsiptype) & mask) != 0)
		{*dp++ = d >> HalfBitsiptype; ++k; }
	*dp++ = d & mask; k++;
	return k;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertdoublebinary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of doubles
//	@param		"bits"	-	Number of bits in *dp (either 8 or 4 bits)
//
//	@todo
//
// Description:
// convert an an array of doubles[] into a binary iptype
//
static iptype convertdouble2binary(double *dp, size_t maxinx, const double cy=0, const int bits=8)
	{
	const unsigned int mask = bits==8? 0xff: 0xf;
	iptype d=(unsigned char)cy;
	for (int i = 0; i < maxinx; ++i)
		{
		d <<= bits;
		d |= (unsigned char)(*dp++)&mask;
		if (bits == 4 && ++i < maxinx)
			{
			d <<= bits;
			d |= (unsigned char)(*dp++) & mask;
			}
		}
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertuint8_2binary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of unsigned int
//
//	@todo
//
// Description:
// convert an an array of unsigned int [] into a binary iptype. Used in Schonhagen-Strassen
//
static iptype convertuint8_2binary(unsigned int *dp, size_t maxinx, const unsigned int cy = 0)
	{
	iptype d = (unsigned char)cy;
	for (size_t i = maxinx; i > 0; --i)
		{
		d <<= 8;
		d |= (unsigned char)(dp[i-1]);
		}
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Aug/2021
//	@brief 		void convertuint8_2binary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of unsigned int
//
//	@todo
//
// Description:
// convert an an array of unsigned int [] into a binary iptype. Used in Schonhagen-Strassen
//
static iptype convertuint16_2binary(uint64_t *dp, size_t maxinx, const unsigned int cy = 0)
	{
	iptype d = (unsigned short)cy;
	for (size_t i = maxinx; i > 0; --i)
		{
		d <<= 16;
		d |= (unsigned short)(dp[i - 1]);
		}
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17/Aug/2021
//	@brief 		void convertHalfiptype2binary
//	@return		iptype	-	return the binary constructed number
//	@param		"dp"	-	pointer to array of half iptypes 
//
//	@todo
//
// Description:
// convert an an array of half iptypes [] into a binary iptype. Used in Schonhagen-Strassen
//
static inline iptype convertHalfiptype2binary( uintmax_t *dp, const size_t maxinx, const unsigned int cy = 0)
	{
	iptype d = (iptype)cy;
	if( maxinx > 1 )
		d |= dp[ 1 ];
	d <<= ( Bitsiptype / 2);
	d |= dp[0];
	return d;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_umul_fourier
//	@return 	std::vector<iptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//		Multiply two unsigned binary numbers
//		Optimized: Used FFT algorithm to performed the multiplication
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65  10^5 digits is not enough for arbitrary precsion so we are only using 1byte.	
//
std::vector<iptype> _int_precision_umul_fourier(const std::vector<iptype> *src1, const std::vector<iptype> *src2, int nbits)
	{
	std::vector<iptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits == 0 ? 8 : nbits;
	int radix = bits == 8 ? 256 : 16;
	size_t sz = sizeof(iptype);
	std::vector<double> va, vb;

	l1 = src1->size();
	l2 = src2->size();
	des.reserve(l1 + l2 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1 < l2 ? l2 : l1;  
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(iptype);  // Convert to byte 
	if (l > 8'000'000*sizeof(iptype) || bits == 4)
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2;  // use 2^4 instead of 2^8
		}
	for (n = 1; n < l; n <<= 1) ;
	n <<= 1;
	
#ifdef HVE_THREAD
	// Using parallel sections below speeds up the performance of the two calls to _int_real_Fourier() with a factor of 1.8 
	if (nbits == 0 || l1 + l2>10'000)
		{// Starting thread using lambda expressions
		// L1, l2, va, vb by reference since it is used after the thread has terminated
		std::thread first([&, n, bits]()
			{std::vector<iptype>::const_reverse_iterator pos, end;
			size_t i;
			va.resize(n);
			for (i = 0, pos = src1->rbegin(), end = src1->rend(); pos != end; ++pos)
				i += convertbinary2double(&va[i], *pos, i == 0, bits );
			l1 = i; // L1 now Number of bytes or nibbles
			_vector_real_fourier(va, n, 1); // FFT va
			});

		std::thread second([&, n, bits]()
			{std::vector<iptype>::const_reverse_iterator pos, end;
			size_t i;
			vb.resize(n);
			for (i = 0, pos = src2->rbegin(), end = src2->rend(); pos != end; ++pos)
				i += convertbinary2double(&vb[i], *pos, i == 0, bits);
			l2 = i; // L2 now Number of bytes or nibbles
			_vector_real_fourier(vb, n, 1); // FFT vb
			});

		first.join();
		second.join();
	}
	else
#endif
		{
		std::vector<iptype>::const_reverse_iterator pos, end;
		va.resize(n);
		vb.resize(n);
		// Start with most significant fptype e.g. src1[0]
		for (l1 = 0, pos = src1->rbegin(), end = src1->rend(); pos != end; ++pos)
			l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
		// L1 now Number of bytes or nibbles
		// Start with most significant fptype e.g. src2[0]
		for (l2 = 0, pos = src2->rbegin(), end = src2->rend(); pos != end; ++pos)
			l2 += convertbinary2double(&vb[l2], *pos, l2 == 0, bits);
		// L2 now number of bytes or nibbles
		_vector_real_fourier(va, n, 1); // FFT va
		_vector_real_fourier(vb, n, 1); // FFT vb
		}

	vb[0] *= va[0];
	vb[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t;
		vb[j] = (t = vb[j])*va[j] - vb[j + 1] * va[j + 1];
		vb[j + 1] = t*va[j + 1] + vb[j + 1] * va[j];
		}
	_vector_real_fourier(vb, n, -1);
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = vb[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		vb[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l2 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;   // Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		iptype num;
		size_t inx = l1 - sz*(l2 - l + 1);
		num = convertdouble2binary(&vb[inx], sizeof(iptype), 0, bits);
		des.push_back(num);
		}
	l2 = l1 % sz;   // Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)		// do the the last 64bit integers from b
		{
		fptype num;
		num = convertdouble2binary(&vb[0], l2, cy, bits);
		des.push_back(num);
		}
	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		17-Aug-2021
//	@brief 		std::vector<iptype> _int_precision_umul_karatsuba_
//	@return 	std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned binary numbers of vector<iptype>, using the karatsuba method
//	 Karatsuba is faster than umul_fourier with operands up to	 xx decimal digits whereafter _umul_fourier is faster
//   Notice when operands can fit into a 64bit integer we switch to native multiplications.
//
std::vector<iptype> _int_precision_umul_karatsuba(const std::vector<iptype> *lhs, const std::vector<iptype> *rhs)
	{
	std::vector<iptype> result, z0, z1, z2, z3;
	std::vector<iptype> lhshigh, lhslow, rhshigh, rhslow;
	int wrap;
	size_t half_length, length, l_length = lhs->size(), r_length = rhs->size(), tot_len;
	length = l_length < r_length ? r_length : l_length;
	tot_len = l_length + r_length;

	// Short cuts
	if ((l_length == 1 && lhs->front() == 0) || (r_length==1 && rhs->front() == 0) )
		{//lhs * 0 or rhs*0 is zero
		result=std::vector<iptype> (1,0);
		return result;
		}
	if (l_length == 1 && lhs->front() == 1 )
		{
		//rhs * 1 is rhs
		result = *rhs;
		return result;
		}
	if( r_length == 1 && rhs->front() == 1)
		{
		//lhs * 1 is lhs
		result = *lhs;
		return result;
		}
	if (l_length==1 && r_length==1)  // If max digits in lhs & rhs less than to fit into a 64 bit integer then do it the binary way
		{
		result = _precision_umul64(lhs->front(), rhs->front());
		if (result[0] == 0 || result[1] == 0)
			wrap = 0;
		_int_precision_strip_trailing_zeros(&result);
		return result;
		}

	// Splitting
	half_length = length >> 1;
	if (l_length <= half_length)
		{
		lhshigh.insert(lhshigh.begin(), 1, 0); lhslow.insert(lhslow.begin(),lhs->begin(), lhs->end() );
		}
	else if (l_length < length)
		{
		if (half_length >= l_length )
			lhshigh.insert(lhshigh.begin(), 1, 0);
		else
			lhshigh.insert(lhshigh.begin(),lhs->begin()+half_length,lhs->end());
		lhslow.insert(lhslow.begin(), lhs->begin(), lhs->begin() + half_length );
		}
	else
		{
		lhslow.insert(lhslow.begin(), lhs->begin(), lhs->begin() + half_length);
		lhshigh.insert(lhshigh.begin(), lhs->begin() + half_length, lhs->end());
		}

	if (r_length <= half_length)
		{
		rhshigh.insert(rhshigh.begin(), 1, 0); 
		rhslow.insert(rhslow.begin(),rhs->begin(), rhs->end() );
		}
	else if (r_length < length)
		{
		if (half_length >= r_length ) 
			rhshigh.insert(rhshigh.begin(), 1, 0);
		else 
			rhshigh.insert(rhshigh.begin(), rhs->begin()+half_length, rhs->end());
		rhslow.insert(rhslow.begin(), rhs->begin(), rhs->begin()+half_length);
		}
	else
		{
		rhslow.insert(rhslow.begin(), rhs->begin(), rhs->begin() + half_length);
		rhshigh.insert(rhshigh.begin(), rhs->begin() + half_length, rhs->end());
		}

	// Evaluation
	z0 = _int_precision_uadd(&lhshigh, &lhslow);
	z1 = _int_precision_uadd(&rhshigh, &rhslow);
	z2 = _int_precision_umul_karatsuba(&z0, &z1);
	z1 = _int_precision_umul_karatsuba(&lhslow, &rhslow);
	z0 = _int_precision_umul_karatsuba(&lhshigh, &rhshigh);
	z3 = _int_precision_uadd(&z0, &z1);
	z3 = _int_precision_usub(&wrap, &z2, &z3);

	// Recomposition
	z0 = _int_precision_ushiftleft( &z0, Bitsiptype*2*half_length);
	z3 = _int_precision_ushiftleft( &z3, Bitsiptype*half_length);
	z0 = _int_precision_uadd(&z0, &z1);
	result = _int_precision_uadd(&z0, &z3);
	return result;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22-Aug-2019
//	@brief 		std::vector<iptype> _int_precision_schonhage_strassen_linear_umul
//	@return		std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
//   Notice when operands can fit into a 64bit integer we switch to native multiplications.
//
/*
std::vector<iptype> _int_precision_schonhage_strassen_linear_umul1(std::vector<iptype> *lhs, std::vector<iptype> *rhs)
	{
	size_t i, j;
	size_t l_length = lhs->size(), r_length = rhs->size(), length=l_length+r_length;
	std::vector<iptype> res;
	std::vector<iptype>::reverse_iterator pos;
	std::vector<unsigned int> linearconvolution(sizeof(iptype)*(l_length + r_length), 0);  // initialize it with zero
	std::vector<unsigned char> ua(l_length * sizeof(iptype)), ub(r_length * sizeof(iptype));

	// Convert to unsigned byte from vector<iptype> and notice we stored in reverse order 
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into unsigned char from iptype by mapping each mBinary number  
	// ua[0]=an-1>>56, ua[1]=an-1>>48&0xff, ua[2]=an-1>>40&0xff,...ua[7]=an-1&0xff skipping leading zero byte in the most significat number
	// the function convertbinary2uint8() does this job per mBinary iptype number
	for (i = 0, pos = lhs->rbegin(); pos != lhs->rend(); ++pos)
		i += convertbinary2uint8( &ua[i], *pos, i == 0 );
	l_length = i;  // l_length now in bytes instead of iptype's
	for (j = 0, pos = rhs->rbegin(); pos != rhs->rend(); ++pos)
		j += convertbinary2uint8( &ub[j], *pos, j == 0 );
	r_length = j;  // l_length now in bytes instead of iptype's

	// do the linear convolution
	for (i = 0; i < r_length; ++i)
		for (j = 0; j < l_length; ++j)
			linearconvolution[i + j] += (unsigned int)ub[r_length - 1 - i] * (unsigned int)ua[l_length - 1 - j];
	
	res.reserve(length + 16);
	unsigned int nextCarry = 0;
	for (i = 0; i < l_length + r_length - 1; ++i )
		{
		linearconvolution[i] += nextCarry;
		nextCarry = linearconvolution[i] / 256;
		linearconvolution[i] %= 256;
		}
	if (nextCarry != 0)
		linearconvolution[i++] = nextCarry % 256;
	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte

	// Now convert then back into a vector<iptype> format. i is the number of bytes in the result
	j = ( i / sizeof(iptype) == 0) ? i : sizeof(iptype);   // Number of remaining 8bits digits in the most siginifcant iptype 
	if ( j > 0 )  // do the j bytes uptil 64bit integers first starting forward from linearconvolution
		{
		iptype num;
		num = convertuint8_2binary(&linearconvolution[0], j, 0 );
		res.push_back(num);
		}
	j = (i-1) / sizeof(iptype);
	for (size_t l = 1; l <= j; ++l)  // do the full 64bit integers first starting backwards from b
		{
		iptype num;
		size_t inx = sizeof(iptype)*l;
		num = convertuint8_2binary(&linearconvolution[inx], sizeof(iptype));
		res.push_back(num);
		}

	_int_precision_strip_trailing_zeros(&res);

	return res;
	}
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Oct/2021
//	@brief 		std::vector<iptype> _int_precision_schonhage_strassen_linear_umul
//	@return		std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
//   Notice when operands can fit into a 64bit integer we switch to native multiplications.
//
/*
std::vector<iptype> _int_precision_schonhage_strassen_linear_umul2(std::vector<iptype> *lhs, std::vector<iptype> *rhs, const unsigned bytes )
	{
	size_t i, j;
	size_t l_length = lhs->size(), r_length = rhs->size(), length = l_length + r_length;
	const unsigned chunksize = sizeof(iptype) / bytes;
	const uint64_t radix = (uint64_t)1 << bytes * 8;
	std::vector<iptype> res;
	std::vector<iptype>::reverse_iterator pos;
	std::vector<uint64_t> linearconvolution(sizeof(iptype)*(l_length + r_length), 0);  // initialize it with zero
	std::vector<unsigned short> ua(l_length * sizeof(iptype)), ub(r_length * sizeof(iptype));
	uint64_t nextCarry = 0;

	// Convert to unsigned byte from vector<iptype> and notice we stored in reverse order 
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into unsigned char from iptype by mapping each mBinary number  
	// ua[0]=an-1>>56, ua[1]=an-1>>48&0xff, ua[2]=an-1>>40&0xff,...ua[7]=an-1&0xff skipping leading zero byte in the most significat number
	// the function convertbinary2uint8() does this job per mBinary iptype number
	for (i = 0, pos = lhs->rbegin(); pos != lhs->rend(); ++pos)
		i += convertbinary2uint16(&ua[i], *pos, i == 0);
	l_length = i;  // l_length now in bytes instead of iptype's
	for (j = 0, pos = rhs->rbegin(); pos != rhs->rend(); ++pos)
		j += convertbinary2uint16(&ub[j], *pos, j == 0);
	r_length = j;  // l_length now in bytes instead of iptype's

	// do the linear convolution
	for (i = 0; i < r_length; ++i)
		for (j = 0; j < l_length; ++j)
			{
			uint64_t m = (uint64_t)ub[r_length - 1 - i] * (uint64_t)ua[l_length - 1 - j];
			linearconvolution[i + j] += m;
			if (linearconvolution[i + j] < m) // carry
				{
				// Propagate carry
				for (size_t k = i + j + 1; ; ++k)
					{
					linearconvolution[k] += 1;		// Add carry
					if (linearconvolution[k] >= 1) // Any carry?
						break;
					}
				}
			}

	res.reserve(length + 16);
	for (i = 0; i < l_length + r_length - 1; ++i)
		{
		linearconvolution[i] += nextCarry;
		nextCarry = linearconvolution[i] / radix;
		if (nextCarry != 0)
			nextCarry = nextCarry;
		linearconvolution[i] %= radix;
		}
	if (nextCarry != 0)
		linearconvolution[i++] = nextCarry % radix;

	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte

	// Now convert then back into a vector<iptype> format. i is the number of bytes in the result
	j = ( i / chunksize == 0) ? i : chunksize;   // Number of remaining 16bits digits in the most siginifcant iptype 
	if (j > 0)  // do the j bytes uptil 64bit integers first starting forward from linearconvolution
		{
		iptype num;
		num = convertuint16_2binary(&linearconvolution[0], j, 0);
		res.push_back(num);
		}
	j = (i - 1) / chunksize;
	for (size_t l = 1; l <= j; ++l)  // do the full 64bit integers first starting backwards from b
		{
		iptype num;
		size_t inx = chunksize*l;
		num = convertuint16_2binary(&linearconvolution[inx], chunksize);
		res.push_back(num);
		}

	_int_precision_strip_trailing_zeros(&res);

	return res;
	}
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		16/Oct/2021
//	@brief 		std::vector<iptype> _int_precision_schonhage_strassen_linear_umul
//	@return		std::vector<iptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	First unsigned source argument
//	@param		"rhs"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
//   Using the full advantages of Half bit size of iptype in the linear convolution
//
std::vector<iptype> _int_precision_umul_linear(const std::vector<iptype> *lhs, const std::vector<iptype> *rhs)
	{
	size_t i, j;
	size_t l_length = lhs->size(), r_length = rhs->size();
	const unsigned int HalfBitsiptype = Bitsiptype >> 1;  // same as / 2
	const iptype mask = (~(iptype)0) >> HalfBitsiptype;
	const uintmax_t radix = (uintmax_t)1 << HalfBitsiptype;
	std::vector<iptype> res;
	std::vector<iptype>::const_reverse_iterator pos, end;
	std::vector<uintmax_t> linearconvolution( 2 * (l_length + r_length), 0);  // initialize it with zero
	std::vector<iptype> ua( l_length * 2), ub( r_length * 2 );
	uintmax_t nextCarry = 0;

	// Convert to half iptype from vector<iptype> and notice we stored in reverse order 
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half iptype  from iptype by mapping each mBinary number into 2 half iptype numbers 
	// the function convertbinary2Halfiptype() does this job per mBinary iptype number
	for (i = 0, pos = lhs->rbegin(), end=lhs->rend(); pos != end; ++pos)
		i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
	l_length = i;  // l_length now in half iptype's  instead of iptype's
	for (j = 0, pos = rhs->rbegin(), end=rhs->rend(); pos != end; ++pos)
		j += convertbinary2Halfiptype(&ub[j], *pos, j == 0);
	r_length = j;  // l_length now in half iptype's instead of iptype's

	// do the linear convolution
	for (i = 0; i < r_length; ++i)
		for (j = 0; j < l_length; ++j)
			{
			uintmax_t m = (uintmax_t)ub[r_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
			linearconvolution[i + j] += m;
			if (linearconvolution[i + j] < m) // carry
				{
				// Propagate carry
				for (size_t k = i + j + 1; ; ++k)
					{
					linearconvolution[k] += radix;	// Add carry
					if (linearconvolution[k] >= 1)	// Continue adding carry?
						break;
					}
				}
			}

	res.reserve( r_length + l_length + 2 );
	for (i = 0; i < l_length + r_length - 1; ++i)
		{
		linearconvolution[i] += nextCarry;
		nextCarry = linearconvolution[i] >> HalfBitsiptype;  //  same as / radix;
		if (nextCarry != 0)
			nextCarry = nextCarry;
		linearconvolution[i] &= mask;  // same as %= radix;
		}
	if (nextCarry != 0)
		linearconvolution[i++] = nextCarry & mask; // same as % radix;

	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
	// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
	// do the full 64bit integers first starting from least significant HalfBitsiptype
	for ( j = 0; j < i; j+=2 ) 
		{
		iptype num;
		num = convertHalfiptype2binary(&linearconvolution[j], 2);
		res.push_back(num);
		}

	_int_precision_strip_trailing_zeros(&res);
	return res;
	}

// Short Division: The ptype digit d  is divide up into the unsigned ptype vector
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021
//	@brief 		std::vector<iptype>	_int_precision_udiv_short
//	@return 	std::vector<iptype>	- The result of the short division
//	@param      "src1"				- Source string to divide with the short number
//	@param      "d"					- Number to divide
// @param		"remaind"			- The remaind of the short division
//
//	@todo
//
// Description:
//   Short divide: The ptype digit d [0..2^32] is divided up in the unsigned vector<iptype> 
//	  Notice only up to int 32bit can be handle as short div.
//   Divide with zero throw an exception
//
std::vector<iptype> _int_precision_udiv_short(iptype *remaind, const std::vector<iptype> *src1, const iptype d)
	{
	const unsigned int shifts = 4 * sizeof(iptype);
	const iptype mask = (~((iptype)0) ) >> shifts;
	iptype ir;
	std::vector<iptype>::const_reverse_iterator s_pos, s_end;
	std::vector<iptype> des;

	if (d == 0)
		{
		throw int_precision::divide_by_zero();
		}

	if (d == 1)  // Divide by one dont change the src1.
		{
		des = *const_cast<std::vector<iptype> *>(src1);
		_int_precision_strip_trailing_zeros(&des);
		*remaind = 0;
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	s_pos = src1->rbegin();
	s_end = src1->rend();
	for ( ir=0; s_pos != s_end; ++s_pos)
		{
		iptype n, qh, ql;
		/*if (ir == 0)
			{// Shortcut when carry is zero
			des.push_back(*s_pos / d );
			ir = *s_pos % d;
			}
		else*/
			{
			n = *s_pos >> shifts;
			n |= ir << shifts;
			qh = n / d;	ir = n % d;
			n = *s_pos & mask;
			n |= ir << shifts; 
			ql = n / d;	ir = n % d;
			n = (qh << shifts) | ql;
			des.push_back(n);
			}
		}

	reverse(des.begin(), des.end());
	_int_precision_strip_trailing_zeros(&des);
	*remaind = ir;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_udiv
//	@return		std::vector<iptype>-	the result of disivison
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Divide two unsigned binary numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<iptype> _int_precision_udiv(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int plusdigit, plusbits, wrap, i;
	iptype underflow;
	std::vector<iptype> des, tmp, quotient, divisor;

	des.push_back(0);
	divisor = *const_cast<std::vector<iptype> *> (src1);
	if (src2->size() == 1 && ( src2->front() >> 32 ) == 0 ) // Make short div if denominator <= 32 bit integer.
		return _int_precision_udiv_short( &underflow, &divisor, src2->front());

	plusdigit = (int)divisor.size() - (int)src2->size();
	if (plusdigit < 0)  //src1 / src2 == 0
		return des;

	plusbits = (int)_int_precision_clz(src2->back())- (int)_int_precision_clz(divisor.back()) ;
	plusbits=plusdigit * Bitsiptype + plusbits;
	for(i=0; plusbits >= 1 ; ++i) 
		{
		tmp = _int_precision_ushiftleft(src2, plusbits);
		if (_int_precision_compare(&divisor, &tmp) < 0)
			{ // Too much reduce with one power of radix
			--plusbits; continue;
			}
		divisor = _int_precision_usub(&wrap, &divisor, &tmp);
		quotient.clear();
		quotient.insert(quotient.begin(), (plusbits / Bitsiptype ) + 1, 0);
		quotient[quotient.size() - 1] = (iptype)(1) << (plusbits % Bitsiptype );
		des = _int_precision_uadd(&des, &quotient);
		}

	for (wrap = 0; wrap == 0; )
		{
		divisor = _int_precision_usub(&wrap, &divisor, src2);
		if (wrap == 0) // src1 was indeed > src2
			des = _int_precision_uadd_short(&des, 1);
		}

	_int_precision_strip_trailing_zeros(&des);
	return des;
	}

// Short Remainder: The iptype digit d [1..2^64] is divide up into the unsigned vector<iptype> and the remaing is returned
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6-Aug-2021
//	@brief 		std::vector<iptype> _int_precision_urem_short
//	@return 	std::vector<iptype> - 	the result of the short remainder
//	@param      "src1"	-	Source string to divide with the short number
//	@param      "d"	   -	Number to divide
//
//	@todo
//
// Description:
//   Short remainder: The iptype digit d [0..2^64] is divided up in the unsigned vector<iptype> and the remaing is retuened
//   Divide with zero throw an exception
//   if d==1 then result == 0, for d==2,4,5,8,10 we only test the last few digits to get the result. This speed up rem for large integers with small rem value
//   since we dont have to run through every digits in src1
//
std::vector<iptype> _int_precision_urem_short(const std::vector<iptype> *src1, const iptype d)
	{
	const unsigned int shifts = 4 * sizeof(iptype);
	const iptype mask = (~((iptype)0)) >> shifts;
	iptype ir;
	std::vector<iptype>::const_reverse_iterator s_pos, s_end;
	std::vector<iptype> des;

	if (d == 0)
		{
		throw int_precision::divide_by_zero();
		}

	if (d == 1)  // Remainer is always 0 for d==1
		{
		des.push_back(0);
		return des;
		}

	// Short cut
	ir = *const_cast<std::vector<iptype> *>(src1)->begin();
	switch (d)
		{
		case 2: des.push_back(ir % 2); return des; break;
		case 4: des.push_back(ir % 4); return des; break;
		case 5: des.push_back(ir % 5); return des; break;
		case 8: des.push_back(ir % 8); return des; break;
		case 10: des.push_back(ir % 10); return des; break;
		default:;   // No Short cut
		}

	s_pos = src1->rbegin();
	s_end = src1->rend();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		iptype n;
		if (ir == 0)
			{// Shortcut when carry is zero
			ir = *s_pos % d;
			}
		else
			{
			n = *s_pos >> shifts;
			n += ir << shifts;
			ir = n % d;
			n = *s_pos & mask;
			n += ir << shifts;
			ir = n % d;
			}
		}

	des.push_back(ir);
	return des;
	 }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_urem
//	@return		std::vector<iptype>	-	the remaing result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Find the remainder when divide two unsigned vector<iptype> numbers
//   Optimized: Used early out add and multiplication w. zero
//

std::vector<iptype> _int_precision_urem(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	int wrap;
	std::vector<iptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 &&  (src2->front() >> 32) == 0 ) // Make short rem 
		{
		iptype rem;
		//des = _int_precision_urem_short( src1, src2->front());
		//return _int_precision_udiv_short(&underflow, &divisor, src2->front());
		_int_precision_udiv_short(&rem, src1, src2->front());
		des[0] = rem;
		return des;
		}

	tmp = _int_precision_udiv( src1, src2);
	tmp = _int_precision_umul( &tmp, src2);
	des = _int_precision_usub(&wrap, src1, &tmp);

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		23/Oct/2021
//	@brief 		std::vector<iptype> _int_precision_udivrem
//	@return		std::vector<iptype>	-	the quotient.  Result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
// @param		"r"		-	Pointer to remainder of the division
//
//	@todo
//
// Description:
//   Find both the remainder and the qoutient when divide two unsigned vector<iptype> numbers
//   
//
std::vector<iptype> _int_precision_udivrem(std::vector<iptype> *src1, std::vector<iptype> *src2, std::vector<iptype> *r )
	{
	int wrap;
	std::vector<iptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short rem 
		{
		iptype rem;
		//des = _int_precision_urem_short( src1, src2->front());	
		//return _int_precision_udiv_short(&underflow, &divisor, src2->front());
		des = _int_precision_udiv_short(&rem, src1, src2->front());
		tmp.assign(1, rem);
		*r = tmp;
		return des;
		}

	des = _int_precision_udiv(src1, src2);
	tmp = _int_precision_umul(&des, src2);
	tmp = _int_precision_usub(&wrap, src1, &tmp);

	_int_precision_strip_trailing_zeros(&des);
	*r = tmp;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_ushiftright
//	@return 	std::vector<iptype> - The negated number	
//	@param      "src1"	-	The number to be shifted
// @param		"shift"	-	The shift count
//
//	@todo
//
// Description:
// Implement the shift left operation src1 >> src2
//
std::vector<iptype> _int_precision_ushiftright(const std::vector<iptype> *src1, const size_t shift )
	{
	size_t shiftwidth = Bitsiptype;
	std::vector<iptype>::const_reverse_iterator pos, end;
	std::vector<iptype> des, c0(1,0);
	size_t discard, within;
	size_t i;
	iptype carry, mask;
	
	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsiptype
	if (_int_precision_compare(src1,&c0)==0)  // Short cut: a zero number zero shifting right is still zero.
		return *src1;
	if (shift == 0)  // Short cut: shift zero right does not change the number.
		return *src1;

	discard = shift / Bitsiptype;
	within = shift % Bitsiptype;
	shiftwidth -= within;
	mask = (~((iptype)0)) >> shiftwidth;  // check is resukltis still ok for shiftwidth == 64
	i = src1->size();
	for (carry=0, pos = src1->rbegin(), end=src1->rend(); pos != end && i>discard; --i, ++pos)
		{
		iptype n, nextcarry;
		n = *pos;
		nextcarry = n & mask;
		n >>= within;
		carry <<= shiftwidth;  // check shiftwidth==64 abd the resukt is stilkl ok
		n |= carry;
		carry = nextcarry;
		des.push_back(n);
		}
	reverse(des.begin(), des.end());// check for des.size()==0
	if (des.size() == 0)
		des.push_back(0);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_ushiftleft
//	@return 	std::vector<iptype> - The negated number	
//	@param      "src1"	-	The number to be shifted
// @param		"shift"	-	The shift count
//
//	@todo
//
// Description:
// Implement the shift left operation src1 << src2
//
std::vector<iptype> _int_precision_ushiftleft(const std::vector<iptype> *src1, const size_t shift)
	{
	const unsigned int shiftwidth = 8 * sizeof(iptype);
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype> des, c0(1,0);
	size_t adding, within;
	iptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % sizeof(iptype)*8
	if (_int_precision_compare(src1, &c0) == 0)  // Short cut: a zero number zero shifting left is still zero.
		return *const_cast<std::vector<iptype> *>(src1);
	if (shift == 0)  // Short cut: shift zero left does not change the number.
		return *const_cast<std::vector<iptype> *>(src1);

	adding = shift / Bitsiptype;
	within = shift % Bitsiptype;
	//for (; adding > 0; --adding)
	//	des.push_back(0);
	if( adding > 0 )
		des.insert(des.begin(), adding, 0);
	mask = (~((iptype)0)) >> within;  mask = ~mask;
	for (carry = 0, pos = src1->begin(), end=src1->end(); pos != end; ++pos)
		{
		iptype n, nextcarry;
		n = *pos;
		nextcarry = n & mask;
		n <<= within;
		n |= carry >> (shiftwidth - within);    // check shiftwidth 64 and witin == 0
		carry = nextcarry;
		des.push_back(n);
		}
	if (carry != 0)
		des.push_back(carry >> (shiftwidth - within)); // check shiftwidth 64 and witin == 0

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uneg
//	@return 	std::vector<iptype> - The negated number	
//	@param      "src"	-	The number to negate
//
//	@todo
//
// Description:
//   Negate one-complement the unsigned integer src
//
std::vector<iptype> _int_precision_uneg(const std::vector<iptype> *src)
	{
	std::vector<iptype>::iterator pos;
	std::vector<iptype> des;

	des = *const_cast<std::vector<iptype> *>(src);
	for (pos = des.begin(); pos != des.end(); ++pos)
		{
		*pos = ~*pos;
		}

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uand
//	@return		std::vector<iptype>	-	the result of anding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   And two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uand(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos;
	std::vector<iptype>::iterator d_end, d_pos;

	if (src1->size() >= src2->size()) // Making the sortest operand the result operand since that will be the maxium number of digits
		{
		des = *const_cast<std::vector<iptype> *>(src2);
		pos = src1->begin();
		d_end = des.end();
		}
	else
		{
		des = *const_cast<std::vector<iptype> *>(src1);
		pos = src2->begin();
		d_end = des.end();
		}
	d_pos = des.begin();

	for (; d_pos != d_end; ++pos, ++d_pos )
		{ // Anding element by element for the two numbers
		*d_pos &= *pos;
		}

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uor
//	@return		std::vector<iptype>	-	the result of oring src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Or two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uor(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	if (src1->size() >= src2->size()) // Making the sortest operand the result operand since that will be the maxium number of digits
		{
		des = *const_cast<std::vector<iptype> *>(src1);
		pos = src2->begin();
		end = src2->end();
		}
	else
		{
		des = *const_cast<std::vector<iptype> *>(src2);
		pos = src1->begin();
		end = src1->end();
		}
	d_pos = des.begin();

	for (; pos != end; ++pos, ++d_pos)
		{ // oring element by element for the two numbers
		*d_pos |= *pos;
		}

	//_int_precision_strip_leading_zeros(&des1);

	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_uxor
//	@return		std::vector<iptype>	-	the result of xoring src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Xor two unsigned vector<iptype> numbers
//   Optimized: Used early out and
//
std::vector<iptype> _int_precision_uxor(const std::vector<iptype> *src1, const std::vector<iptype> *src2)
	{
	std::vector<iptype> des;
	std::vector<iptype>::const_iterator pos, end;
	std::vector<iptype>::iterator d_pos;

	if (src1->size() >= src2->size()) // Making the sortest operand the result operand since that will be the maxium number of digits
		{
		des = *const_cast<std::vector<iptype> *>(src1);
		pos = src2->begin();
		end = src2->end();
		}
	else
		{
		des = *const_cast<std::vector<iptype> *>(src2);
		pos = src1->begin();
		end = src1->end();
		}
	d_pos = des.begin();

	for (; pos != end; ++pos, ++d_pos)
		{ // xoring element by element for the two numbers
			*d_pos ^= *pos;
		}

	//_int_precision_strip_leading_zeros(&des1);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		std::vector<iptype> _int_precision_unegate
//	@return		std::vector<iptype>	-	the result of negating src1 
//	@param		"src1"	-	First unsigned source argument
//
//	@todo
//
// Description:
//   Negate unsigned vector<iptype> number
//
std::vector<iptype> _int_precision_unegate(std::vector<iptype> *src1)
	{
	std::vector<iptype> des=*src1;
	std::vector<iptype>::iterator d_pos;

	for (d_pos = des.begin(); d_pos != des.end(); ++d_pos)
		{ // negating element of the number
		*d_pos = ~*d_pos;
		}

	return des;
	}

///////////////////////////////////////////////
//
//
//    To and from string conversion including power of 10 tables for up to 64bit unsigned integers
//		or higher power of tables for splitting numbers into trunks, kilotrunks, megatrunks or
//		gigatrunks.
//
//
///////////////////////////////////////////////

// Same table as above but converted to float_precisions. Need to be sure that defaul precision is >=20
static float_precision _fpPowerof10Table[20] = {float_precision(_powerof10Table[0]),float_precision(_powerof10Table[1]),float_precision(_powerof10Table[2]),
												float_precision(_powerof10Table[3]),float_precision(_powerof10Table[4]),float_precision(_powerof10Table[5]),
												float_precision(_powerof10Table[6]),float_precision(_powerof10Table[7]),float_precision(_powerof10Table[8]),
												float_precision(_powerof10Table[9]),float_precision(_powerof10Table[10]),float_precision(_powerof10Table[11]),
												float_precision(_powerof10Table[12]),float_precision(_powerof10Table[13]),float_precision(_powerof10Table[14]), 
												float_precision(_powerof10Table[15]),float_precision(_powerof10Table[16]),float_precision(_powerof10Table[17]),
												float_precision(_powerof10Table[18]),float_precision(_powerof10Table[19]) };


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		25/Oct/2021
//	@brief 		std::string check_digits
//	@return 	bool	-	Return true if digits is valid other throw an exception and return false
//	@param		"s"	-		String of digits to check
//	@param		"start"	-	Value to convert to ascii string based on RADIX
//	@param		"end"	-	RADIX value of conversion 
//	@param		"base"	-	Check for base (default BASE_10)
//
//	@todo  Add to do things	
//
// Description:
//		Check the ascii string for valid digits according to base.
//
static bool check_digits(const std::string s, const size_t start, const size_t end, const int base = BASE_10)
	{
	if (base <= BASE_10)
		{
		for (size_t i = start; i < end; i++)
			{
			if (s[i] < '0' || s[i] >= base + '0')
				{
				throw int_precision::bad_int_syntax();
				return false;
				}
			}
		}
	else
		for (size_t i = start; i < end; i++)
			{ // base 11..36
			if ((s[i] < '0' || s[i] > '9') && (tolower(s[i]) < 'a' || tolower(s[i]) > base -BASE_10 +'a' ) )
				{
				throw int_precision::bad_int_syntax();
				return false;
				}
			}
	return true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Dec/2021
//	@brief 		std::string remove_separators from string
//	@return 	std::string	-	Return the string without seperators '\''
//	@param		"s"	-		String of digits to remove separators from
//
//	@todo  Add to do things	
//
// Description:
//	Remove thousand seperators from string
//	A thousand separator is either a ' or a space  ' '
//
static std::string remove_separators(const std::string& s)
	{
	std::string r;
	std::string::const_iterator pos = s.begin();
	r.reserve(s.length() + 16);
	for (; pos != s.end(); ++pos)
		if (*pos != '\'' && *pos != ' ')
			r.push_back(*pos);

	return r;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		31/Oct/2021
//	@brief 		std::string uitostring10
//	@return 	static std::string	-	Return the ascii representation of number
//	@param		"value"	-	Value to convert to ascii string based on RADIX
//	@param		"minlength"	-	minlength of string. default to 0
//
//	@todo  Add to do things	
//
// Description:
//   This convert an uintmax_t (64bit) integer into a decimal string in base 10
//
static inline std::string uitostring10(const uintmax_t value, const unsigned minlength=0 )
	{
	std::string s;
	unsigned digit;
	uintmax_t uvalue=value;

	do
		{
		digit = (unsigned)(uvalue % BASE_10 );
		uvalue /= BASE_10;
		s.push_back((char)ICHARACTER10((unsigned char)digit));
	} while (uvalue > 0);

	while (s.length() < minlength)
		s.push_back('0');
	reverse(s.begin(), s.end());
	return s;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/3/2006
//	@brief 		std::string itostring
//	@return 	static std::string	-	Return the ascii representation of number
//	@param		"value"	-	Value to convert to ascii string based on RADIX
//	@param		"radix"	-	RADIX value of conversion 
//
//	@todo  Add to do things	
//
// Description:
//   This function replace Microsoft _itoa() to a generic function that return the
//   string representation of the number in Base Radix. 
//   Radix can be in the range from 2..256 (only 2..36 deliveres a readable string)
//   only if value is < 0 will it return with a leading sign
std::string itostring( const int value, const unsigned radix )
   {
   std::string s;
   unsigned digit;
   unsigned uvalue;
   
   if (radix < BASE_2 || radix > 36 )
      return s;  // Conversion not supported

   if( radix == BASE_10 && value < 0 )
      uvalue = -value;
   else
      uvalue = (unsigned)value;

   do 
      {
      digit = (unsigned) (uvalue % radix);
      uvalue /= radix;

      if( radix <= 36 )
         {
         // Convert to ascii and store
         if( digit < 10 )
            s.push_back( (char)ICHARACTER10( (unsigned char)digit ) );      
         else
            s.push_back( (char)( digit - 10 + 'a' ) );      
         }
      else
         { // Keep it 'binary' not readable string
         s.push_back( (unsigned char)digit );      
         }
   } while (uvalue > 0);

   if( radix == BASE_10 && value < 0 )
      s.push_back( '-' );
   reverse(s.begin(), s.end());
   return s;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		std::string _int_precision_atoip
//	@return 	vector<iptype>	-	The integer precision number
//	@param		"str"	-	The arbitrary precision string as a std::string
// @param		"sign"	-	Returned the sign as either +1 or -1
//
//	@todo  
//
// Description:
// Convert ascii string to vector<iptype> number
// A leading 0 is intepreted as a octal number
// a leading 0x is interpreted as a hexadecimal number 
// a leading 0b is interpreted as a binary number
// otherwise it's a decimal number.
//
std::vector<iptype> _int_precision_atoip(const std::string& str, int *sign)
	{
	std::vector<iptype> number;

	number = _int_precision_atoip(str.c_str(), sign);
	return number;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_atoi
//	@return 	string	-	The integer precision string
//	@param		"s"	-		The arbitrary precision string as char *
//	@param		"sign"	-	Returned the sign as either +1 or -1
//
//	@todo  
//
// Description:
// Convert ascii string to a binary number
// A leading 0 is intepreted as a octal number (BASE_8)
// a leading 0x is interpreted as a hexadecimal number  (BASE_16)
// a leading 0b is interpreted as a binary number (BASE_2)
// otherwise it's a decimal number.(BASE_10)
// The resulting number is a vector<iptype> type
//
std::vector<iptype> _int_precision_atoipold(const char *str, int *sign)
	{
	//extern std::string remove_separators(const std::string&); // Forward declaration
	std::string s(str);
	std::string::const_iterator pos;
	std::vector<iptype> number(1,0);
	size_t i; std::string::const_iterator pos_start;
	int base = BASE_10;
	size_t digits;
	size_t pwr = 10'000'000'000'000'000;  //1E16
	 
	s=remove_separators(s);
	number.reserve(s.size() + 16);
	*sign = +1;
	pos = s.begin();
	if (*pos == '+' || *pos == '-')
		{
		*sign = CHAR_SIGN(*pos);
		++pos;
		if (pos == s.end())
			{
			throw int_precision::bad_int_syntax();
			}
		}

	if (*pos == '0') // Octal, binary or hex representation
		{
		if (pos + 1 != s.end() && tolower(pos[1]) == 'x')
			{
			base = BASE_16; 
			digits = sizeof(iptype); 
			pwr = (uintmax_t)1u << 4*digits; pos += 2;
			}
		else
			if (pos + 1 != s.end() && tolower(pos[1]) == 'b')
				{// Collec Binary representation
				base = BASE_2; 
				digits = sizeof(iptype) * 4; 
				pwr = (uintmax_t)1u << digits; pos += 2;
				}
			else
				{ // Collect octal representation
				base = BASE_8; 
				digits = sizeof(iptype) * 2; 
				pwr = (uintmax_t)1u << 3*digits; 
				}
		}
	else
		{ // Collect decimal representation
		base = BASE_10; 
		digits = sizeof(iptype) * 2; // 16,8,4,2
		//pwr = digits == 16 ? 10000000000000000 : 100000000;
		pwr = 100; 
		if (digits > 2) pwr *= pwr;
		if (digits > 4) pwr *= pwr;
		if (digits > 8) pwr *= pwr;
		}

	check_digits(s, pos - s.begin(), s.size(), base);
	for (pos_start=pos, i=0; pos != s.end(); ++pos)
		{
		i++;
		if (i%digits == 0 )
			{
			std::string s2 = s.substr(pos_start-s.begin(), digits);
			uintmax_t n =  strtoull(s2.c_str(),NULL, base);
			build_i_number(number, n, pwr);
			pos_start = pos+1;
			}
		}

	if( i%digits != 0 )
		{
		std::string s2 = s.substr(pos_start - s.begin(), i%digits);
		int64_t n = strtoull(s2.c_str(),NULL,base);
		if(base==BASE_10) pwr=_powerof10Table[i%digits];
		else pwr= base << (i%digits-1);
		build_i_number(number, n, pwr);
		}

	if (number.size() == 1 && number[0] == 0 && *sign == -1)
		*sign = +1;

	return number;
	}
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Aug/2021
//	@brief 		std::vector<iptype> _int_precision_atoi
//	@return 	string	-	The integer precision string
//	@param		"s"	-		The arbitrary precision string as char *
//	@param		"sign"	-	Returned the sign as either +1 or -1
//
//	@todo  
//
// Description:
// Convert ascii string to a binary number
// A leading 0 is intepreted as a octal number (BASE_8)
// a leading 0x is interpreted as a hexadecimal number  (BASE_16)
// a leading 0b is interpreted as a binary number (BASE_2)
// otherwise it's a decimal number.(BASE_10)
// The resulting number is a vector<iptype> type
//
std::vector<iptype> _int_precision_atoip(const char *str, int *sign)
	{
	std::string s(str);
	std::string::const_iterator pos;
	std::vector<iptype> number(1, 0);
	std::string::const_iterator pos_start;
	int base = BASE_10;
	size_t digits;
	size_t pwr = 10'000'000'000'000'000;  //1E16

	s = remove_separators(s);
	number.reserve(s.size() + 16);
	*sign = +1;
	pos = s.begin();
	if (*pos == '+' || *pos == '-')
		{
		*sign = CHAR_SIGN(*pos);
		++pos;
		if (pos == s.end())
			throw int_precision::bad_int_syntax();
		}

	if (*pos == '0') // Octal, binary or hex representation
		{
		if (pos + 1 != s.end() && tolower(pos[1]) == 'x')
			{
			base = BASE_16;
			digits = sizeof(iptype);
			pwr = (uintmax_t)1u << 4 * digits; pos += 2;
			}
		else
			if (pos + 1 != s.end() && tolower(pos[1]) == 'b')
				{// Collec Binary representation
				base = BASE_2;
				digits = sizeof(iptype) * 4;
				pwr = (uintmax_t)1u << digits; pos += 2;
				}
			else
				{ // Collect octal representation
				base = BASE_8;
				digits = sizeof(iptype) * 2;
				pwr = (uintmax_t)1u << 3 * digits;
				}
		}

	check_digits(s, pos - s.begin(), s.size(), base);
	if(base == BASE_10)
		number = string2number(s, pos - s.begin(), s.size());
	else
		if(base==BASE_2 || base==BASE_16)
			number = stringbase2number(s, pos - s.begin(), s.size(),base);
		else
			{// BASE_8
			number = stringbase8_2number(s, pos - s.begin(), s.size() );
		/*	intmax_t i;
			for (pos_start = pos, i = 0; pos != s.end(); ++pos)
				{
				i++;
				if (i%digits == 0)
					{
					std::string s2 = s.substr(pos_start - s.begin(), digits);
					uintmax_t n = strtoull(s2.c_str(), NULL, base);
					build_i_number(number, n, pwr);
					pos_start = pos + 1;
					}
				}

			if (i%digits != 0)
				{
				std::string s2 = s.substr(pos_start - s.begin(), i%digits);
				int64_t n = strtoull(s2.c_str(), NULL, base);
				if (base == BASE_10) pwr = _powerof10Table[i%digits];
				else pwr = base << (i%digits - 1);
				build_i_number(number, n, pwr);
				}
				*/
			}

	if (number.size() == 1 && number[0] == 0 && *sign == -1)
		*sign = +1;

	return number;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Aug/2021
//	@brief 		std::vector<iptype>  _int_precision_itoa
//	@return 	std::string -	Return the inter precision as a string (no leading sign)
//	@param      "a"	-	the internal integer precision string
// @param		"base" -	the base to convert to. default is BASE_10
//
//	@todo 
//
// Description:
//   Convert int_precsion to ascii string using base
//   The string returned has no leading sign
//
std::string _int_precision_itoa(const std::vector<iptype>  *a, const int base )
	{
	std::vector<iptype> src, c0(1,0), cbase(1,base), tmp_rem;
	std::string s;

	src = *a;
	s.reserve(src.capacity());
	if (src.size() == 1 && base==BASE_10)  // Only one iptype binary number in vector. Convert it directly to string and return
		{
		return s = uitostring10((uintmax_t)src[0]);
		}
	if (base != BASE_10)
		{ // All other bases
		for (; src.size() > 1 || src[0] != 0;)
			{// Take one digit at  time. 
			tmp_rem = _int_precision_urem(&src, &cbase);
			src = _int_precision_udiv(&src, &cbase);
			if (base == BASE_16 && tmp_rem[0] >= BASE_10)
				s.push_back((char)((unsigned char)(tmp_rem[0] - BASE_10 + 'a')));
			else
				s.push_back((char)ICHARACTER10((unsigned char)tmp_rem[0]));
			}
		}
	else
		{// Base 10 only. Loop until we have a max of one iptype number left
		std::vector<iptype> cbasepower10_div(1, _powerof10Table[8]);// 1E8
		for (; src.size() > 1;)
			{ // Take 8 decimal digits at a time
			src = _int_precision_udivrem(&src, &cbasepower10_div, &tmp_rem);
			s = uitostring10((uintmax_t)tmp_rem[0], 8) + s;
			}
		if (src[0] != 0)
			{// Do the last iptype directly via "native" functions. yielding up to 18 digits, instead of a single digit
			s = uitostring10((uintmax_t)src[0]) + s;
			}
		}
	if (s.size() == 0)
		s.push_back((char)ICHARACTER10((unsigned char)0));
	if (base == BASE_2) s += "b0";
	if (base == BASE_8)	s += "0";
	if (base == BASE_16) s += "x0";
	if(base!=BASE_10)
		reverse(s.begin(), s.end());

	return s;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2012
//	@brief 		Calculate abs(x)
//	@return 	int_precision -	Return absolute value of x
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   int precision abs()
//    
//
int_precision abs(const int_precision& x)
	{
	int_precision i;
	
	if (x.sign() < 0)
		{
		i = x; i.change_sign();
		return i;
		}

	return x;
	}

///////////////////////////////////////////////
//
//
//    Miscellaneous function
//		ipow()
//		ipow_modulo()
//		iprime()
//		gcd()
//		lcm()
//
//
///////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		26/Aug/2007
//	@brief 		return the integer power of x^y
//	@return 	int_precision	-	The integer precision power of x^y
//	@param		"x"	-	The int precision x
//	@param		"y"	-	The int precision y. Max y is 2^32-1
//
//	@todo  
//
// Description:
// Return the integer power of x^y. For any pratical purpose the power y is restricted to 2^32-1
//
int_precision ipow( const int_precision& x, const int_precision& y )
   {
   int_precision p(x), r(1);

   for(uintmax_t n = y; n > 0; n >>= 1) 
        {
        if( ( n & 0x1 ) != 0 ) r *= p;  // Odd
        if( n > 1 ) p *= p;				// Square it				 
        }
   return r;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		24/Aug/2012
//	@brief 		return the integer power of x^y%z
//	@return 	int_precision	-	The integer precision power of x^y%z
//	@param		"x"	-	The int precision x
//	@param		"y"	-	The int precision y. Max y is 2^32-1
// @param		"z"	-	The int precision z.
//	@todo  
//
// Description:
// Return the integer power of x^y%z. For any pratical purose the power y is restricted to 2^32-1
//
int_precision ipow_modulo( const int_precision& x, const int_precision& y, const int_precision& z )
   {
   int_precision p(x), r(1);

   p%=z;
   for(int n = y; n > 0; n >>= 1) 
        {
        if( ( n & 0x1 ) != 0 ) { r *= p; r %= z; } // Odd
        p *= p;	p %= z;					 
        }
   return r;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Sep/2012
//	@brief 		Check a number for a prime
//	@return 	bool-	true is the integer is a prime number false otherwise
//	@param		"prime"	-	The int precision prime
//	@todo  
//
// Description:
// Return true if the integer prime is a prime number.
// All integers are of the form 30k + i for i = 0, 1, 2,...,29 and k an integer from 0..  However, 2 divides 0, 2, 4,...,28 and 3 divides 0, 3, 6,...,27 and 5 divides 0, 5, 10,...,25. 
// So all prime numbers are of the form 30k + i for i = 1, 7, 11, 13, 17, 19, 23, 29 (i.e. for i < 30 such that gcd(i,30) = 1). 
// Note that if i and 30 are not coprime, then 30k + i is divisible by a prime divisor of 30, namely 2, 3 or 5, and is therefore not a prime.
// 
//
bool iprime(const int_precision& prime)
	{
	int precheck[11] = { 10, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
	int primes[9] = { 8, 1, 7, 11, 13, 17, 19, 23, 29 };
	int_precision count, kp(30), mod;
	int i;

	for (i = 1; i <= precheck[0]; i++)
	if ((int)(prime % (int_precision)precheck[i]) == 0) return prime==int_precision(precheck[i]);

	for (; kp * kp < prime; kp += 30)   //Loop to divide the number by every number 6*count-1 and 6*count+1 and count < sqrt(i)
		{
		for (i = 1; i <= primes[0]; i++)
			{
			count = kp + primes[i];
			mod = prime % count;
			if (mod == 0)				// Statement to change the variable 'check' to 1 if the number gets divided
				return false;			//meaning its not prime
			}
		}

	return true;						// It is a prime
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		30/Aug/2021
//	@brief 		Greatest Common Divisor (gcd)
//	@return 	int_precision - gcd(a,b)
//	@param		"a"	-	first number - positive
// @param		"b"	-	second number - positive
//
//	@todo  
//
// Description:
// Return the greatest common divisor of the two numbers a & b.
// It used the Binary gcd method only using shifting and subtraction
// Changed to also handle negative arguments a,b;
//	Algorithm: while (b>0) { tmp = b; b = a%b; a = tmp; } return a;
//
int_precision gcd( const int_precision& a, const int_precision& b )
	{
	size_t shift;
	int_precision u, v, tmp, i0(0), i1(1);

	// GCD(0,v)==v; GCD(u,0)==0; GCD(0,0)==0
	if (a == i0) return b;
	if (b == i0) return a;
	u = a; v = b; 
	if(u < i0) u = -u; 
	if(v < i0) v = -v;
#if false
	while (b > i0)
		{
		tmp = v; v = u%v; u = tmp;
		}
	return u;
#else
		tmp = u | v;
		shift = tmp.ctz();
		u >>=u.ctz();
		do {
			v >>= v.ctz();
			if (u > v) {
				tmp = v; v = u; u = tmp;
				}
			v -= u;
		} while (v != i0);
		return u << shift;
	}
#endif


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Feb/2017
//	@brief 		Least Common Multiplier
//	@return 	int_precision - lcm(a,b)
//	@param		"a"	-	first number
// @param		"b"	-	second number
//
//	@todo  
//
// Description:
// Return the least common multiplier of the two numbers a & b.
// It used the Binary lcm method only using shifting, subtraction and one multiplication and one division
// 
//
int_precision lcm(const int_precision& a, const int_precision& b)
	{
	int_precision r(a), gcd_ab;

	gcd_ab = gcd(a, b);
	r /= gcd_ab;
	r *= b;
	return r;
	}



///////////////////////////////////////////////
//
//
//    End of Integer Precision Core
//
//
///////////////////////////////////////////////

///////////////////////////////////////////////
//
//
//    Floating point Precision Core
//
//
///////////////////////////////////////////////

class float_precision_ctrl float_precision_ctrl(PRECISION,ROUND_NEAR);

///////////////////////////////////////////////
//
//
//    Floating point Precision Input, Output operator
//
//
///////////////////////////////////////////////


std::ostream& operator<<( std::ostream& strm, const float_precision& d )
    { return strm << _float_precision_fptoa( const_cast<float_precision *>(&d) ).c_str();}

std::istream& operator>>( std::istream& strm, float_precision& d )
         { char ch; std::string s; int cnt, exp_cnt=0;
         strm >> ch;  while( ch == ' ' ) strm.get(ch);  // Ignore leading white space.
         if( ch == '+' || ch == '-' ) { s += ch; strm >> ch; } else s += '+';  // Parse sign
         for( cnt = 0; ch >= '0' && ch <= '9'; cnt++, strm >> ch ) s += ch;  // Parse integer part
         if( ch == '.' )  // Parse fraction
            for( s += '.', strm >> ch; ch >= '0' && ch <= '9'; cnt++, strm >> ch ) s += ch;   // Parse fraction part
         if( ch == 'e' || ch == 'E' )
            {
            s += 'e'; strm >> ch; if( ch == '+' || ch == '-' ) { s += ch; strm >> ch; } else s += '+';  // Parse Expo sign 
            for( exp_cnt =0; ch >= '0' && ch <= '9'; exp_cnt++, strm >> ch ) s += ch;  // Parse expo number
            }

         std::cin.putback( ch );  // ch contains the first character not part of the number, so put it back
         if( !strm.fail() && ( cnt > 0 || exp_cnt > 0 ) )  // Valid number 
            d = float_precision( const_cast<char *>( s.c_str() ), float_precision_ctrl.precision(), float_precision_ctrl.mode() );
         return strm;
         }


///////////////////////////////////////
//
// CONVERT BINARY FLOAT PRECISION to and from ascii representation
//
//   _float_precision_fptoa()		-- Convert from float_precision to string format
//	  _float_precision_fptoainteger()-- Convert from float_precision to integer string format. Obsolete?
//	  _float_precision_fptod()		-- Convert from float_precision to double
//   _float_precision_dtofp()		-- Convert from double to float_precision
//   _float_precision_atofp()		-- Convert from string format to float_precision format
//	  _float_precision_fptoip()		-- Convert from float_precision to int_precision
//	  static check_float_digits()	-- Check string for illegal number characers
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		25/Oct/2021
//	@brief 		std::string check_'float_digits
//	@return 	bool	-	Return true if digits is valid other throw an exception and return false
//	@param		"s"	-		String of digits to check
//	@param		"start"	-	Value to convert to ascii string based on RADIX
//	@param		"end"	-	RADIX value of conversion 
//
//	@todo  Add to do things	
//
// Description:
//		Check the ascii string for valid digits according to base.
//
static bool check_float_digits(const std::string& s, const size_t start, const size_t end )
	{
	for (size_t i = start; i < end; i++)
		{
		if (s[i] < '0' || s[i] > '9')
			{
			throw float_precision::bad_float_syntax();
			return false;
			}
		}
	return true;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert part of float_precision numbers into string (decimal representation)
//	@return		std::string -	The partly converted decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//	@param		"digits"	-	The remainf digis in the float_precision number 
//
//	@todo 	
//
// Description:
//   Convert partly a float_precision numbers into string (decimal representation)
//		A call to the function convert the next up to max_digits decimal number into a string
//		and return the string
//		Digits to be converted is inthe range 1..max_digits
//
static std::string number2Decimal(float_precision& fp, const int digits )
	{
	uintmax_t di;
	int min_width = digits;						// So many wanted digits. always >= 1

	if (min_width <= 0) min_width = 1;					// at least 1
	if (min_width > MAX_DECIMAL_DIGITS) min_width = MAX_DECIMAL_DIGITS;	// Take max 18 digits at a time
	fp *= _fpPowerof10Table[min_width];			// float_precision(_powerof10table[min_width]);
	di = (uintmax_t)fp.toFraction();  			// Get the next max digit decimal number										
	fp.precision(fp.precision() - min_width);	// Reduce precision with min_width decimal digits
	return uitostring10(di, min_width);			// Return min_width length string
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//

static std::string trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp); 
	static float_precision _trunkPowerof10(0, MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_trunkPowerof10.iszero())// is _trunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_TRUNK_SIZE, _trunkPowerof10 = float_precision(1); i > 0; i >>= 1 )	
			{// Build multiply factor for trunk size
			if (i & 1) _trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;					// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_TRUNK_SIZE; ++i)
		{
		str += number2Decimal(trunk, MAX_DECIMAL_DIGITS);
		}
	fp *= _trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a Megatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//

static std::string kilotrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _kilotrunkPowerof10(0, MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_kilotrunkPowerof10.iszero())// is _megatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _kilotrunkPowerof10 = float_precision(1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _kilotrunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_KILOTRUNK_SIZE; ++i)
		{
		str += trunk2Decimal(trunk);
		}
	fp *= _kilotrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Dec/2021
//	@brief 		Convert a Megatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string megatrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _megatrunkPowerof10(0, MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_megatrunkPowerof10.iszero())// is _megatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		for (i = MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE, _megatrunkPowerof10 = float_precision(1); i > 0; i >>= 1)
			{// Build multiply factor for trunk size
				if (i & 1) _megatrunkPowerof10 *= p;	// Odd
				if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGATRUNK_SIZE; ++i)
		{
		str += kilotrunk2Decimal(trunk);
		}
	fp *= _megatrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		14/Dec/2021
//	@brief 		Calculate the expoenent reduction for very small numbers.
//	@return		int -	The decimal 10 exponent reductions
//	@param		"f"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   The number f much be less than 1. e.g. exponent() < 0 
//	 This claulate the decimal exponent reducion from he base 2 internal exponent
//
static int exponent2reduction(float_precision& f)
	{
	int expo10reduction = 0;
	eptype expo = f.exponent();  // start with the power of 2 exponent
	static float_precision e1kfactor(0, 1'000 + 2);
	static float_precision e10kfactor(0, 10'000 + 2);
	static float_precision e100kfactor(0, 100'000 + 2);
	float_precision p(10, e1kfactor.precision());

	if (expo < 0)
		{
		expo10reduction = -(int)(expo / log2(10));
		if(expo10reduction!=0) 
			expo10reduction -= 1;  // Calculate the exponent decimal 10 reduction that ensure that 10^expo10reduction< 2^exponent
		expo = expo10reduction; // Now it is converted to a power of 10 exponent
		f.precision(f.precision() /*+ 5*/);  // Add extra guard bits
		//Need to biuld any factors?
		if (expo >= 1'000 && e1kfactor.iszero())
			{// Build the e10factor first time it is needed.
			e1kfactor = float_precision(1); p = float_precision(10);
			for (eptype i = 1'000; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e1kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}
		if (expo >= 10'000 && e10kfactor.iszero())
			{// Build the e10factor first time it is needed.
			p.precision(e10kfactor.precision());
			e10kfactor = float_precision(1); p = e1kfactor;
			for (eptype i = 10; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e10kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}

 		if (expo >= 100'000 && e100kfactor.iszero())
			{// Build the e10factor first time it is needed.
			p.precision(e100kfactor.precision());
			e100kfactor = float_precision(1); p = e10kfactor;
			for (eptype i = 10; i > 0; i >>= 1)
				{// Build multiply factor for trunk size
				if (i & 1) e100kfactor *= p; // Odd
				if (i > 1) p *= p;			// square it
				}
			}

		for (; expo >= 100'000; expo -= 100'000)
			f *= e100kfactor;

		for (; expo >= 10'000; expo -= 10'000)
			f *= e10kfactor;

		for (; expo >= 1'000; expo -= 1'000)
			f *= e1kfactor;
	
		for (; expo > 0; expo -= std::min((int)MAX_DECIMAL_DIGITS, (int)expo))
			f *= _fpPowerof10Table[std::min((int)MAX_DECIMAL_DIGITS, (int)expo)];

		f.precision(f.precision() /*- 5*/);  // reverse the extra guard bits
		}

	return expo10reduction;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Mega10trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string mega10trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _mega10trunkPowerof10(0, MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_mega10trunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_mega10trunkPowerof10 = float_precision(1);
		for (i = MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _mega10trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGA10TRUNK_SIZE; ++i)
		{
		str += megatrunk2Decimal(trunk);
		}
	fp *= _mega10trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Mega100trunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string mega100trunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _mega100trunkPowerof10(0, MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_mega100trunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_mega100trunkPowerof10 = float_precision(1);
		for (i = MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _mega100trunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_MEGA100TRUNK_SIZE; ++i)
		{
		str += mega10trunk2Decimal(trunk);
		}
	fp *= _mega100trunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Dec/2021
//	@brief 		Convert a Gigatrunk of float_precision numbers into string (decimal representation)
//	@return		std::string -	The trunk converted to a decimal floating point string	
//	@param		"fp"		-	float_precision number to partly convert
//
//	@todo 	
//
// Description:
//   Convert a trunk of a float_precision numbers into string (decimal representation)
//		A call to the function convert the next trunk size decimal number into a string
//		and return the trunk string
//		A trunk size is defined as threshold*max_digits number of decimal digits
//
static std::string gigatrunk2Decimal(float_precision& fp)
	{
	int i;
	std::string str;
	float_precision trunk(fp);
	static float_precision _gigatrunkPowerof10(0, MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE * MAX_TRUNK_SIZE * 20, ROUND_DOWN);

	if (_gigatrunkPowerof10.iszero())// is _gigatrunkPowerof10 build or created
		{
		float_precision p(_powerof10Table[MAX_DECIMAL_DIGITS], fp.precision(), fp.mode());
		_gigatrunkPowerof10 = float_precision(1);
		for (i = MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE; i > 0; i >>= 1)
			{// Build multiply factor for trunk size
			if (i & 1) _gigatrunkPowerof10 *= p;	// Odd
			if (i > 1) p *= p;						// square it
			}
		}

	// Set trunk precision to reduce precision to at least hold the needed number of digits we can ge out of the trunk size
	trunk.precision(MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);	// We overallocate the precision 20 versus 18 digits in the trunk
	str.reserve(MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE * MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE * 20);
	// Do it trunk size of thr (eshold) uintmax_t digits
	for (i = 0; i < MAX_GIGATRUNK_SIZE; ++i)
		{
		str += mega100trunk2Decimal(trunk);
		}
	fp *= _gigatrunkPowerof10;
	fp.toFraction();
	fp.precision(fp.precision() - MAX_GIGATRUNK_SIZE*MAX_MEGA100TRUNK_SIZE*MAX_MEGA10TRUNK_SIZE*MAX_MEGATRUNK_SIZE*MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS);
	return str;
	}

#if false
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Sep/2021
//	@brief 		Convert float_precision numbers into string (decimal representation)
//	@return		std::string -	The decimal floating point string	
//	@param		"a"			-	float_precision number to convert
//
//	@todo 	
//
// Description:
//   Convert float_precision numbers into string (decimal representation)
//   
//
std::string _float_precision_fptoaold(const float_precision *a, int maxdig )
	{
	float_precision fp(*a), fracp, intp, c0(0);
	int_precision ip;
	std::string str;
	int expo10, sign = 1;
	size_t found=0, i;
	const int max_digits = sizeof(intmax_t) >= 8 ? std::min(18,maxdig) : std::min(9, maxdig);
	float_precision ftable[19];

	fracp.precision(a->precision());
	fracp = modf(fp, &intp);
	if (intp.iszero() && fracp.iszero() )
		{str = "0E0"; return str;}
	if (intp.sign() < 0 || fracp.sign() < 0)
		sign = -1;
	intp = fabs(intp);
	ip = intp; 
	str += ip.toString();
	expo10 = (int)( str.size() - 1 );
	if (!fracp.iszero() || str.size() > 1 )
		str.insert(str.begin() + 1, '.'); 
	fracp = fabs(fracp);
	fracp.precision(fracp.precision() + 1);  // Add an extra guard bit
	str.reserve(fp.precision() + str.length() + 32);
	for (i = 0; i <= max_digits; ++i)  // build the ftable for faster access
		ftable[i] = float_precision(_powerof10Table[i], 20);
	i = 0;
	for (size_t len=fp.precision()+str.length(); 
		(str.length() < len||abs(fracp.exponent())>3) && !fracp.iszero(); ++i)
		{
		uintmax_t di;
		int min_width = found==0 ? max_digits : (int)(len - str.length());		// So many remaining digits. always >= 1

		if(min_width <=0 ) min_width=1;		// at least 1
		if(min_width > max_digits ) min_width=max_digits;	// Take max 18 digits at a time
		fracp *= ftable[min_width];			// float_precision(_powerof10table[min_width]);
		di=(uintmax_t)fracp.toFraction();  	// previous di = (uintmax_t)fracp; fracp -= float_precision(di);
		// Reduce precision with min_width
		fracp.precision(fracp.precision() - min_width);
		if (di != 0) found = 1;				// Indicate we have reach a no zero in the fraction
		str += uitostring10(di, min_width);
		if (i%1000==0&&i!=0)
			std::cout << "\t\tString length =" << str.length() << " Capacity=" << str.capacity() << std::endl;
		} 

	// str is on the form x.yyyyyyyyyyyy
	// Remove any digits exceeding the precision
	if (str.size() > fp.precision() + 2)
		str.erase(fp.precision() + 2);
	
	// Remove trailing zeros if any
	found = str.find_last_not_of('0');
	if (found != std::string::npos && str.size() > std::max(2, (int)found + 1))
		{
		str.erase(std::max(2, (int)found + 1));
		if(str[str.size() - 1] == '.') 
			str.erase(1);
		}

	// check and find negative exponent.
	if(str[0]=='0')
		{
		found = str.find_first_not_of('0', 2);
		if (found == std::string::npos)
			{
			str = "0E0"; return str;
			}
		else
			{
			expo10 = -(int)(found - 1);
			str.erase( 0, found );
			if(str.size()>1)
				str.insert(str.begin() + 1, '.' );
			}
		}

	// Do the exponent part. Add E[+-]exponent to back of str
	str += "E";
	str += itostring(expo10, BASE_10);
	str = (sign < 0 ? "-" : "") + str;
	return str;
	}

#endif

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Dec/2021
//	@brief 		Convert float_precision numbers into string (decimal representation)
//	@return		std::string -	The decimal floating point string	
//	@param		"a"			-	float_precision number to convert
//
//	@todo 	
//
// Description:
//   Convert float_precision numbers into string (decimal representation)
//   We do it in multiple steps.
//	 1)	Repeat doing a megatrunk size by extrating a megatrunk size from the fptoa number
//		in the same manner as step 2
//   2) Repeat doing a trunk size by extracting a trunk size from fptoa number 
//		Then repeat multiple group within the trunk size of max_digis number at the time
//		The most efficient trunk size has be found to be 50*max_digits decimal in a trunk size
//		For every loop we reduce the original fptoa precision  with the trunk size precision 
//	 3)	Take the remaning fpto number in the range 0..max_digits
//
std::string _float_precision_fptoa(const float_precision *a)
	{
	const size_t thr = MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS;
	const size_t kthr = thr * MAX_KILOTRUNK_SIZE;
	const size_t mthr = kthr * MAX_MEGATRUNK_SIZE;
	const size_t m10thr = mthr * MAX_MEGA10TRUNK_SIZE;
	const size_t m100thr = m10thr * MAX_MEGA100TRUNK_SIZE;
	const size_t gigathr = m100thr * MAX_GIGATRUNK_SIZE;
	float_precision fracp, intp;
	std::string str;
	int expo10, sign = 1;
	size_t found, len;
	
	if (a->iszero() ) 
		return std::string("0E0");
	sign = a->sign();
	fracp.precision(a->precision());		// ensure fracp and intp has the same precision as a
	intp.precision(a->precision());
	fracp = modf(*a, &intp);				// Separate the integer and fraction 
	intp = fabs(intp);
	str += _float_precision_fptoainteger(&intp);  // Convert integer part to string 
	expo10 = (int)(str.size() - 1);			// Extract Expo as the size of string - 1
	if (!fracp.iszero() || str.size() > 1)	// If fraction part then add "." to string
		str.insert(str.begin() + 1, '.');
	fracp = fabs(fracp);					// Remove any fraction sign if any
	//fracp.precision(fracp.precision()- (int)(fracp.exponent()/log2(10)));  // Adjust for large negative exponent
	fracp.precision(fracp.precision() + 2 );  // Add extra guard bits
	fracp.mode(ROUND_DOWN);

	// Check for large negative exponent to avoid generating leading zero that will be cut of anyway at the end
	if(intp.iszero())
		expo10 -= exponent2reduction(fracp);	// remove large negative exponent from fracp and adjust fracp accordingly
	
	str.reserve(fracp.precision() + str.length() + 32);  // Ensure enough room,for the string to avoid reallocating
	len = fracp.precision() -2 + str.length();// The expected len of the fraction before we have enough

	// Do it in Giga trunks size
	for (; (str.length() + gigathr < len) && !fracp.iszero(); )
		str += gigatrunk2Decimal(fracp);

	// Do it in 100M trunks size
	for (; (str.length() + m100thr < len) && !fracp.iszero(); )
		str += mega100trunk2Decimal(fracp);

	// Do it in 10M trunks size
	for (; (str.length() + m10thr < len) && !fracp.iszero(); )
		str += mega10trunk2Decimal(fracp);
//-------
	// Do it in mega trunks size
	for (; (str.length() + mthr < len) && !fracp.iszero(); )
		str += megatrunk2Decimal(fracp);
			
	// Do it in kilo trunks size
	for (; (str.length() + kthr < len) && !fracp.iszero(); )
		str += kilotrunk2Decimal(fracp);

	// Do it in trunks of treshold times (64bit numbers). threshold*max_digits
	for (; (str.length() + thr < len) && !fracp.iszero(); )
		str += trunk2Decimal(fracp);

	// Do it for the remaining less than a trunk size of threshold
	for (; (str.length() < len || abs(fracp.exponent()) > 3) && !fracp.iszero(); )
		str += number2Decimal(fracp, (int)(len - str.length()));

	// Remove trailing zeros if any
	/*found = str.find_last_not_of('0');
	if (found != std::string::npos && str.size() > std::max(2, (int)found + 1))
		{
		str.erase(std::max(2, (int)found + 1));
		if (str[str.size() - 1] == '.')
			str.erase(1);
		}*/

	// check and find negative exponent.
	if (str[0] == '0')
		{
		found = str.find_first_not_of('0', 2);
		if (found == std::string::npos)
			{// Everything is zero
			return std::string("0E0");
			}
		else
			{// Reduce exponent with the number of leading zeros found
			expo10 += -(int)(found - 1);
			str.erase(0, found);
			if (str.size()>1)
				str.insert(str.begin() + 1, '.');
			}
		}
	// str is on the form x.yyyyyyyyyyyy
	// Remove any digits exceeding the precision
	if (str.size() > a->precision() + 2)
		str.erase(a->precision() + 2);

	// Remove trailing zeros if any
	found = str.find_last_not_of('0');
	if (found != std::string::npos && str.size() > std::max(2, (int)found + 1))
		{
		str.erase(std::max(2, (int)found + 1));
		if (str[str.size() - 1] == '.')
			str.erase(1);
		}

	// Add E[+-]exponent to back of str and sign if negaive to the front
	str = (sign < 0 ? "-" : "") + str+"E"+ itostring(expo10, BASE_10);
	return str;
	}

   
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Convert double (IEE754) into a float_precision numbers 
//	@return 	float_precision - The converted float_precision number
//	@param		"d"				- The double IEEE754 number
// @param		"p"				- The number of significant digits
// @param		"m"				- Round mode
//
//	@todo 	
//
// Description:
//   Convert double numbers into float_precision binary form
//    Constructor for double floating point 
//
float_precision _float_precision_dtofp(double d, size_t p, enum round_mode m)
	{
	eptype expo;
	uintmax_t fpb;
	std::vector<fptype> mb(2, 0);
	float_precision fp(0, p, m);

	if (d == 0)
		return fp;
	if (d < 0)
		{
		fp.sign(-1); d = -d;
		}
	expo = 0;
	fpb = *(uintmax_t *)&d;
	expo = (fpb >> 52) & 0x7ff; 
	expo -= 1023; // unbiased the doubl exponnt
	fp.exponent(expo);
	fpb <<= 12;
	mb[0] = 1; mb[1] = fpb;
	fp.number(mb);

	return fp;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Convert double (IEE754) into a float_precision numbers 
//	@return 	float_precision - The converted float_precision number
//	@param		"d"				- The double IEEE754 number
// @param		"p"				- The number of significant digits
// @param		"m"				- Round mode
//
//	@todo 	
//
// Description:
//   Convert double numbers into float_precision binary form
//    Conversion from float_pecision to double  
//
double _float_precision_fptod( const float_precision *fp )
	{
	uintmax_t t=0, expo;
	double d;

	if (fp->size() == 1 && fp->index(0) == 0)
		return 0.0;
	expo = fp->exponent();
	expo += 1023;  // Biased double exponent format 
	if(fp->size()>1)
		t = (uintmax_t)fp->index(1);
	t >>= 12;
	t |= ( expo << 52 ) & 0x7fffffffffffffff;
	d = *(double *)&t;
	if (fp->sign() < 0)
		d = -d;

	return d;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		Convert a string decimal number into a float_precision number
//	@return 	std::string - The decimal floating point string	
//	@param		"str"		-	ascii string of floating point number to convert
// @param		"p"			- The precision of the number
// @param		"m"			- The round mode of the number
//
//	@todo 	
//
// Description:
//   Convert ascii string into a float_precision numbers 
//    The ascii float format is based on standard C notation
//	
//
float_precision _float_precision_atofp(const std::string& str, size_t p, enum round_mode m)
	{
	return _float_precision_atofp(str.c_str(), p, m);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		Convert a string decimal number into a float_precision number
//	@return 	std::string - The decimal floating point string	
//	@param		"str"		-	ascii string of floating point number to convert
// @param		"p"			- The precision of the number
// @param		"m"			- The round mode of the number
//
//	@todo 	
//
// Description:
//   Convert ascii string into a float_precision numbers 
//    The ascii float format is based on standard C notation
//		Using the string2number function for better performance
//
float_precision _float_precision_atofp(const char *str, size_t p, enum round_mode m)
	{
	int sign, sign_expo;
	int expo_e10=0, f_digits=0, i_digits=0;
	std::string::size_type i, nidx, idx;
	std::string s(str);
	std::string::iterator pos;
	std::vector<iptype> number(1,0), sum(1,0);
	bool ipart=false, fpart=false, epart=false;

	number.reserve(s.size() + 16);
	s=remove_separators(s);
	idx = 0;
	sign = +1;
	// Parse leading sign if any
	pos = s.begin();
	if (*pos == '+' || *pos == '-')  // 
		{
		sign = CHAR_SIGN(*pos);
		pos++;
		idx = 1;
		if (pos == s.end())
			throw float_precision::bad_int_syntax();
		}

	// Determine any significant, fraction sign or exponent sign
	nidx = s.find_first_of(".eE", idx);
	if (nidx == std::string::npos) // Only digits (INTEGER) if any
		{
		int_precision ip(str);			// Construct integer to int_precision
		float_precision fp(ip, p, m );	// Construct float_precision from int_precision 
		return fp;
		} // End of Integer parsing

	// Floating point number starts here
	// Pick up significant beteen idx and nidx 
	if (nidx > idx) // Number of digits before the . sign or exponent Ee
		{
		ipart = true;
		// Strip leading zeros
		for (i = idx; i != nidx; i++) if (s[i] != '0') break;
		if (check_float_digits(s, i, nidx))
			{
			i_digits += (int)( nidx - i );
			if(s[nidx]!='.')  // No fraction part so collect the number
				number = string2number(s, i, nidx); // Faster than decimal2number for large integer portion of the number
			}
		}

	// Floating point representation
	if (s[nidx] == '.') // Any fraction ?
		{
		idx = nidx + 1;                    // Find start of fraction
		nidx = s.find_first_of("eE", idx); // Find end of fraction
		if (nidx == std::string::npos)
			nidx = s.length();

		if (idx < nidx)
			fpart = true;
		// Remove trailing zero digits
		for (i = nidx - 1; i >= idx; i--, nidx--) if (s[i] != '0') break;
		if (check_float_digits(s, idx, nidx))
			{
			f_digits += (int)(nidx - idx);
			// collect both the integer and fraction portion of the number
			if (i_digits != 0)
				{
				if (i_digits < f_digits)
					{
					int j;
					for (i = idx - 1, j=i_digits; j>0; --i, --j )
						s[i] = s[i - 1];
					}
				else
					{
					s.erase(idx - 1, 1); --nidx; --idx;
					}
				}
			number = string2number(s, idx-i_digits, nidx);
			}
		nidx = s.find_first_of("eE", nidx);
		}

	if (nidx != std::string::npos && (s[nidx] == 'e' || s[nidx] == 'E'))
		{// Parse the exponent . which is max sizeof of int. so use regular int operations
		idx = nidx + 1;
		nidx = s.length();
		sign_expo = CHAR_SIGN('+');;
		if (idx < nidx && (s[idx] == '+' || s[idx] == '-'))
			{
			sign_expo = CHAR_SIGN(s[idx]);
			idx++;
			if (idx == nidx)
				{
				throw float_precision::bad_float_syntax();
				}	// Sign but no number
			}
		else
			if (idx >= nidx)
				{
				throw float_precision::bad_float_syntax();
				}  // E but no number

		if (idx < nidx)
			epart = true;
		if (check_float_digits(s, idx, nidx))
			{
			// Collect exponent using base 10
			for (i = idx; i < nidx; i++)
				{
				expo_e10 *= BASE_10;
				expo_e10 += s[i] - '0';
				}
			}
		if (sign_expo < 0)
			expo_e10 = -expo_e10;

		// for illegal number
		if (ipart == false && fpart == false )
			{
			throw float_precision::bad_float_syntax();
			}  // no number before the E or no number at all
		}

	// Build the float_precision number
	int_precision ip(number);
	float_precision fp(ip, p, m), fppow(0,p+2,m);
	expo_e10 -= f_digits;  
	
	// Build the correct power adjustment if needed
	if (expo_e10 != 0)
		{	
		ip = ipow(int_precision(10), int_precision(abs(expo_e10)));
		fppow = float_precision(ip, p, m);
		// Correct for expo
		if (expo_e10 > 0)
			fp *= fppow;
		else
			if (expo_e10 < 0)
				fp *= _float_precision_inverse(fppow);
		} 
	fp.sign(sign);
	return fp;
	}


	//	@author Henrik Vestermark (hve@hvks.com)
	//	@date		7/Sep/2021
	//	@brief 		Convert a string decimal number into a float_precision number
	//	@return 	std::string - The decimal floating point string	
	//	@param		"str"		-	ascii string of floating point number to convert
	// @param		"p"			- The precision of the number
	// @param		"m"			- The round mode of the number
	//
	//	@todo 	
	//
	// Description:
	//   Convert ascii string into a float_precision numbers 
	//    The ascii float format is based on standard C notation
	//	
	//
	float_precision _float_precision_atofp2(const char *str, size_t p, enum round_mode m)
	{
		int sign, sign_expo;
		int expo_e10 = 0, f_digit = 0;
		std::string::size_type i, nidx, idx;
		std::string s(str);
		std::string::iterator pos;
		std::vector<iptype> number(1, 0), sum(1, 0);
		bool ipart = false, fpart = false, epart = false;

		number.reserve(s.size() + 16);
		s = remove_separators(s);
		idx = 0;
		sign = +1;
		// Parse leading sign if any
		pos = s.begin();
		if (*pos == '+' || *pos == '-')  // 
		{
			sign = CHAR_SIGN(*pos);
			pos++;
			idx = 1;
			if (pos == s.end())
				throw float_precision::bad_int_syntax();
		}

		// Determine any significant, fraction sign or exponent sign
		nidx = s.find_first_of(".eE", idx);
		if (nidx == std::string::npos) // Only digits (INTEGER) if any
		{
			int_precision ip(str);			// Construct integer to int_precision
			float_precision fp(ip, p, m);	// Construct float_precision from int_precision 
			return fp;
		} // End of Integer parsing

		  // Floating point number starts here
		  // Pick up significant beteen idx and nidx 
		if (nidx > idx) // Number of digits before the . sign or exponent Ee
		{
			ipart = true;
			// Strip leading zeros
			for (i = idx; i != nidx; i++) if (s[i] != '0') break;
			if (check_float_digits(s, i, nidx))
			{
				//decimal2number(number, s, i, nidx);
				number = string2number(s, i, nidx); // Faster than decimal2number for large integer portion of the number
			}
		}

		// Floating point representation
		if (s[nidx] == '.') // Any fraction ?
			{
			idx = nidx + 1;                      // Find start of fraction
			nidx = s.find_first_of("eE", idx); // Find end of fraction
			if (nidx == std::string::npos)
				nidx = s.length();

			if (idx < nidx)
				fpart = true;
			// Remove trailing zero digits
			for (i = nidx - 1; i >= idx; i--, nidx--) if (s[i] != '0') break;
			if (check_float_digits(s, idx, nidx))
				{
				f_digit += (int)(nidx - idx);

				//for (; idx + MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS <= nidx; idx += MAX_KILOTRUNK_SIZE*MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
				//	kilo2number(number, s, idx);  // NOT Working

				for (; idx + MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS <= nidx; idx += MAX_TRUNK_SIZE*MAX_DECIMAL_DIGITS)
					trunk2number(number, s, idx);

				if (idx != nidx)
					decimal2number(number, s, idx, nidx);
				}
			nidx = s.find_first_of("eE", idx);
			}

		if (nidx != std::string::npos && (s[nidx] == 'e' || s[nidx] == 'E'))
			{// Parse the exponent . which is max sizeof of int. so use regular int operations
			idx = nidx + 1;
			nidx = s.length();
			sign_expo = CHAR_SIGN('+');;
			if (idx < nidx && (s[idx] == '+' || s[idx] == '-'))
			{
				sign_expo = CHAR_SIGN(s[idx]);
				idx++;
				if (idx == nidx)
				{
					throw float_precision::bad_float_syntax();
				}	// Sign but no number
			}
			else
				if (idx >= nidx)
				{
					throw float_precision::bad_float_syntax();
				}  // E but no number

			if (idx < nidx)
				epart = true;
			if (check_float_digits(s, idx, nidx))
			{
				// Collect exponent using base 10
				for (i = idx; i < nidx; i++)
				{
					expo_e10 *= BASE_10;
					expo_e10 += s[i] - '0';
				}
			}
			if (sign_expo < 0)
				expo_e10 = -expo_e10;

			// for illegal number
			if (ipart == false && fpart == false)
			{
				throw float_precision::bad_float_syntax();
			}  // no number before the E or no number at all
		}

		// Build the float_precision number
		int_precision ip(number);
		//std::string check = _int_precision_itoa(&ip);  // DEBUG
		float_precision fp(ip, p, m), fppow(0, p + 2, m);
		expo_e10 -= f_digit;

		// Build the correct power adjustment if needed
		if (expo_e10 != 0)
		{
			ip = ipow(int_precision(10), int_precision(abs(expo_e10)));
			fppow = float_precision(ip, p, m);
			// Correct for expo
			if (expo_e10 > 0)
				fp *= fppow;
			else
				if (expo_e10 < 0)
					fp *= _float_precision_inverse(fppow);
		}
		fp.sign(sign);
		return fp;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		Convert float_precision numbers into int_precision
//	@return		int_precision - The int_precision of the floating point
//	@param		"a"			- float_precision number to convert
//
//	@todo 	
//		Can posible be optimzed by 
//		by vector<iptype> x.insert(x.begin(), fp.x.begin(), fp.x.begin()+chunk+1)
//			ip >>= 64-within; or similar
//
// Description:
//   Convert float_precision numbers into int_precision
//		
//
int_precision _float_precision_fptoip(const float_precision *fp)
	{
	size_t i;
	int_precision ip(1);
	eptype expo = fp->exponent();
	if (expo < 0 || fp->iszero() == true )
		return int_precision(0);
	if (expo == 0) return int_precision( 1 * fp->sign() );
	// fecth expo bits from the fp number after the '.'
	size_t chunk = expo / Bitsfptype;
	int within = expo % Bitsfptype;
	for (i = 0; i < chunk; ++i)
		{
		ip <<= Bitsfptype;
		if( i+1 < fp->size() )
			ip += fp->index( i+1 );
		}
	if (within != 0)
		{
		ip <<= within;
		if (i + 1 < fp->size())
			ip += (fp->index(i + 1) >> (Bitsfptype - within));
		}

	ip.sign(fp->sign());  // Set sign
	return ip;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		Convert float_precision numbers into string (integer representation)
//	@return		std::string - The decimal floating point string	
//	@param		"a"			- float_precision number to convert
//
//	@todo 	
//
// Description:
//   Convert float_precision numbers into integer string (integer representation)
//
std::string _float_precision_fptoainteger(const float_precision *a)
	{
	int_precision ip = *a;
	return ip.toString();
	}

///////////////////////////////////////
//
// END CONVERT FLOAT PRECISION to and from ascii representation
//
///////////////////////////////////////

///////////////////////////////////////
//
// FLOATING POINT CORE FUNCTIONS. STRING
//
//   _float_precision_strip_trailing_zeros
//   _float_precision_rounding
//   _float_precision_uadd_short
//
//   Works Directly on the string class of the float number.
//		Left of the conversion to Binary. assest if this can be done differently
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Remove trailingnosignificant zeros from the string
//	@return 	nothing	
//	@param		"s"	-	digital string
//
//	@todo  
//
// Description:
//   Remove trailing nosignificant zeros
//
void _float_precision_strip_trailing_zeros( std::string *s )
	{
	std::string::reverse_iterator pos;
	int count;

	// Strip trailing zeros
	for( count = 0, pos = s->rbegin(); pos != s->rend() && FDIGIT( *pos ) == 0; pos++ )
         count++;
      
	s->erase( s->length() - count, count );
	if( s->length() == 0 )
		*s = FCHARACTER(0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Round the mantisaa to significant digits and rounding control
//	@return 	int - Return the exponent adjustment (0 or 1) 
//	@param		"m"	-	digital string
// @param		"sign"   - The sign of the number
// @param		"precision" - The digital precision
// @param		 "mode"   - Rounding mode 
//
//	@todo  
//
// Description:
//   Rounding control
//   Round the fraction to the number of precision based on the round mode 
//   Note that the mantissa number has ALWAYS been normalize prior to rounding
//   The mantissa NEVER contain a leading sign
//   Rounding Mode Positive numnber   Result    
//   Rounding to nearest              +�   
//   Rounding toward zero (Truncate)  Maximum, positive finite value   
//   Rounding up (toward +�)          +�   
//   Rounding down) (toward -�)       Maximum, positive finite value   
//
//   Rounding Mode Negative number    Result    
//   Rounding to nearest              -�   
//   Rounding toward zero (Truncate)  Maximum, negative finite value   
//   Rounding up (toward +�)          Maximum, negative finite value   
//   Rounding down) (toward -�)       -�   
//
int _float_precision_rounding( std::string *m, int sign, size_t precision, enum round_mode mode )
   {
   enum round_mode rm = mode;

   if( m->length() > precision )  // More digits than we need 
      {
      if( rm == ROUND_NEAR )
         {
         if( 2 * FDIGIT( (*m)[ precision ] ) >= F_RADIX )
            rm = ROUND_UP; //sign < 0 ? ROUND_DOWN : ROUND_UP;
         else
            rm = ROUND_DOWN; // sign < 0 ? ROUND_UP : ROUND_DOWN;
         }
      else
         if( rm == ROUND_UP && sign < 0 )
            rm = ROUND_DOWN;
         else
            if( rm == ROUND_DOWN && sign < 0 )
               rm = ROUND_UP;

      // Chuck excessive digits
      m->erase( (std::string::size_type)precision, m->length() - precision );

      if( rm == ROUND_UP ) 
         {
         size_t before;

         before = m->length();
         *m = _float_precision_uadd_short( m, 1 );
         if( m->length() > before )
            {
            if( m->length() > precision )
               m->erase( (std::string::size_type)precision, m->length() - precision );

            _float_precision_strip_trailing_zeros( m );            
            return 1;
            }
         }
      }

   _float_precision_strip_trailing_zeros( m );            
   return 0;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		add a short integer to a floating point string (mantissa)
//	@return 	std::string - Return the added string
//	@param		"src1"	-	The source string
// @param		"d"  - The number to add
//
//	@todo  
//
// Description:
//   Short float Add: The digit d [0..F_RADIX] is added to the unsigned fraction string
//   Optimized 0 add or early out add is implemented
//
std::string _float_precision_uadd_short( std::string *src1, unsigned int d )
   {
   unsigned short ireg;
   std::string::reverse_iterator r1_pos, rd_pos;
   std::string des1;

   des1.reserve( src1->capacity() );
   des1 = *src1;
   if( d > F_RADIX )
      {
      throw float_precision::out_of_range();
      }

   if( d == 0 )   // Zero add
      return des1;

   ireg = (unsigned short)( F_RADIX * d );
   rd_pos = des1.rbegin();
   r1_pos = src1->rbegin();
   
   for(; r1_pos != src1->rend(); r1_pos++, rd_pos++ )
      {
      ireg = (unsigned short)( FDIGIT( *r1_pos ) + FCARRY( ireg ) ); 
      *rd_pos = FCHARACTER( (unsigned char)FSINGLE( ireg ) );
      if( FCARRY( ireg ) == 0 ) // Early out add
         break;
      }

   if( FCARRY( ireg ) != 0 )  // Insert the carry in the front of the number
      des1.insert( (std::string::size_type)0, 1, FCHARACTER( (unsigned char)FCARRY( ireg ) ) );

   return des1;
   }

//////////////////////////////////////
//
//	END of CORE Functions. STRING
//
//////////////////////////////////////


///////////////////////////////////////
//
// FLOATING POINT CORE BINARY
//
//   _float_precision_clz						-- Count trailing zeros starting in fptype
//   _float_precision_clz						-- Count trailing zeros starting in vector<fptype>
//	 _float_precision_ctz						-- Count trailing zeros starting in fptype
//	 _float_precision_ctz						-- Count trailing zeros starting in vector<fptype>
//	 _float_precision_strip_leading_zeros		-- Strip leading significant zeros
//   _float_precision_strip_trailing_zeros		-- Strip trailing zeros
//   _float_precision_normalize					-- Normalize the float number 
//   _float_precision_rounding					-- Round the number to given precision
//   _float_precision_right_shift				-- >> shift a vector<fptype> number
//   _float_precision_left_shift				-- <<  shift a vecotr<fptype> number
//   _float_precision_compare					-- Compare to vector<fptype> numbers and determine <,==,>
//   _float_precision_uadd_short				-- Add a vector<fptype> with a ftype number
//   _float_precision_uadd						-- Add two vector<fptype> numbers
//   _float_precision_usub_short				-- Subtract a ftype number from a vector<fptype>
//   _float_precision_usub						-- Subtract two vector<fptype> numbers
//   _float_precision_umul_short				-- Multiply a vector<fptype> with an ftype number
//   _float_precision_umul						-- Multiply two vector<fptype> numbers
//   _float_precision_umul_school				-- Multiply two vector<fptype> numbers using schoolbook algorithm
//   _float_precision_umul_linear				-- Multiply two vector<fptype> numbers using Schonhagen Strassen with linear convolution
//   _float_precision_umul_fourier				-- Multiply two vector<fptype> numbers using FFT
//   _float_precision_udiv_short				-- Divide vector<fptype> with a fptype number
//   _float_precision_udiv						-- Divide two vector<fptype> numbers
//   _float_precision_urem_short				-- Rem vector<fptype> with fptype number
//	 _float_precision_urem						-- Rem two vector<fptype> numbers 
//
//   Works Directly on the string class of the float number
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"a"	-	fptype operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary fptype 
//
size_t _float_precision_clz(const fptype a)
	{
	fptype x = a;
	size_t offset = 0;
	static const unsigned char lookup[16] = { 4,3,2,2,1,1,1,1,0,0,0,0,0,0 };

	if (sizeof(fptype) > 4 && x & 0xffffffff00000000u)
		x >>= 32; else offset += 32;

	if (sizeof(fptype) > 2 && x & 0xffff0000u)
		x >>= 16; else offset += 16;

	if (sizeof(fptype) > 1 && x & 0xff00u)
		x >>= 8; else offset += 8;

	if (x & 0xf0u)
		x >>= 4; else offset += 4;
	offset += lookup[(unsigned char)x];
	return offset;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_float_precision_clz
//	@return		unsigned int	-	the count of leading zero bits in "a"
//	@param		"mb"	-	vector<fptype> operand
//
//	@todo
//
// Description:
//   Count leading nosignificant zeros of the binary fptype 
//
size_t _float_precision_clz(const std::vector<fptype> &mb)
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = 0; i < mb.size(); ++i)
		{
		cnt = _float_precision_clz(mb[i]);
		tot_cnt += cnt;
		if (cnt != Bitsfptype ) 
			break;
		}
	return tot_cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"a"	-	fptype operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary fptype 
//	  fptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _float_precision_ctz(const fptype a)
	{
	fptype x = a;  // sizeof(fptype) can be 8, 4, 2, or 1 
	size_t cnt;
	if (x & 0x1)
		cnt = 0;
	else
		{
		cnt = 1;
		if (sizeof(fptype) >= 8 && (x & 0xffffffffu) == 0)
			{
			x >>= 32; cnt += 32;
			}		
		if (sizeof(fptype) >= 4 && (x & 0xffffu) == 0)
			{
			x >>= 16; cnt += 16;
			}
		if (sizeof(fptype) >= 2 && (x & 0xff) == 0)
			{
			x >>= 8; cnt += 8;
			}
		if ((x & 0xf) == 0)
			{
			x >>= 4; cnt += 4;
			}
		if ((x & 0x3) == 0)
			{
			x >>= 2; cnt += 2;
			}
		cnt -= x & 0x1;
		}

	return cnt;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		_flot_precision_ctz
//	@return		unsigned int	-	the count of trailing zero bits in "a"
//	@param		"mb"	-	vector<fptype> operand
//
//	@todo
//
// Description:
//   Count trailing nosignificant zeros of the binary fptype 
//	  fptype bit word input to count zero bits on right
//   cnt will be the number of zero bits on the right,
//   so if a is 1101000 (base 2), then c will be 3
// NOTE: if 0 == a, then c = 64.
//
size_t _float_precision_ctz(const std::vector<fptype> &mb)
	{
	size_t tot_cnt = 0, cnt;
	for (size_t i = 0; i < mb.size(); ++i)
		{
		cnt = _float_precision_ctz(mb[i]);
		tot_cnt += cnt;
		if (cnt != Bitsfptype ) break;
		}
	return tot_cnt;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		_float_precision_strip_leading_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove leading nosignificant zeros of the binary number
//	This is from the start of the vector<fptype>
//
void _float_precision_strip_leading_zeros(std::vector<fptype> *s)
	{
	std::vector<fptype>::iterator pos;

	// Strip leading zeros
	for (pos = s->begin(); pos != s->end() && *pos == (fptype)0; ++pos);	// Find first not zero digit

	if (s->begin() != pos)
		s->erase(s->begin(), pos);
	if (s->empty())
		s->assign(1, (fptype)0);

	return;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/Sep/2021
//	@brief 		_float_precision_strip_trailing_zeros
//	@return		void	-	
//	@param		"s"	-	pointer to source operand
//
//	@todo
//
// Description:
//   Remove trailing nosignificant zeros of the binary number
//		this is from the top of the vector<fptype> 
//
void _float_precision_strip_trailing_zeros(std::vector<fptype> *s)
	{
	size_t i;
	std::vector<fptype>::reverse_iterator pos;

	// Strip leading zeros
	for (i = s->size() - 1, pos = s->rbegin(); i > 0 && *pos == (fptype)0; ++pos, --i);

	s->resize(i + 1);  // Keep at least one digit by default

	return;
	}
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Right shift a string number 
//	@return 	the result of the shift	
//	@param		"s"	-	digital string
// @param		"shift" - Number of digital shifts
//
//	@todo  
//
// Description:
//   Right shift number x byniary digits by inserting 0 bits in front of the number.
//	  Doing it from index 0..s->size()
//
std::vector<fptype> _float_precision_right_shift(const std::vector<fptype> *src, const size_t shift)
	{
	size_t shiftwidth = Bitsfptype, adding, within;
	std::vector<fptype>::const_iterator pos, end;
	std::vector<fptype> des, c0(1, 0);
	fptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsfptype
	if (_int_precision_compare(src, &c0) == 0)  // Short cut: a zero number zero shifting right is still zero.
		return *src;
	if (shift == 0)  // Short cut: shift zero right does not change the number.
		return *src;

	adding = shift / Bitsfptype;
	des.insert(des.begin(), adding, 0);

	within = shift % Bitsfptype;
	shiftwidth -= within;
	mask = (~((fptype)0)) >> shiftwidth;  // Potential issue if shiftwidth is 64. needs to be tested
	for (carry = 0, pos = src->begin(), end=src->end(); pos != end; ++pos)
		{
		fptype n, nextcarry;
		n = *pos;
		if (within != 0)
			{
			nextcarry = n & mask;
			n >>= within;
			carry <<= shiftwidth;
			n |= carry;
			carry = nextcarry;
			}
		des.push_back(n);
		}
	if (carry != 0)
		des.push_back(carry << shiftwidth);
	if (des.size() == 0)
		des.push_back(0);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Left shift a string number 
//	@return 	The result of the shift	
//	@param		"s"	-	digital string
// @param		"shift" - Number of digital binary shifts
//
//	@todo  
//
// Description:
//   Left shift number x bits. However bis that are shifted out of the first index ([0]) is discarded
//
std::vector<fptype> _float_precision_left_shift(const std::vector<fptype> *src, const size_t shift)
	{
	const unsigned int shiftwidth = Bitsfptype;
	std::vector<fptype>::const_reverse_iterator pos, end;
	std::vector<fptype> des, c0(1, 0);
	size_t discard, within;
	fptype carry, mask;

	// Determine how many full digit shift and the last shift (last shift = shift count % Bitsfptype
	if (_int_precision_compare(src, &c0) == 0)  // Short cut: a zero number zero shifting left is still zero.
		return *src;
	if (shift == 0)  // Short cut: shift zero left does not change the number.
		return *src;

	within = shift % Bitsfptype;
	mask = (~((fptype)0)) >> within;  mask = ~mask;
	discard	= shift / Bitsfptype;
	if (discard < src->size()) // Discard less  than size of src?
		{
		for (carry = 0, pos = src->rbegin(), end = src->rend(); pos != end - discard; ++pos)
			{
			fptype n, nextcarry;
			n = *pos;
			nextcarry = n & mask;
			n <<= within;
			n |= carry >> (shiftwidth - within);  // check is shiftwidth is 64 and witiin is zero and the result is still ok
			carry = nextcarry;
			des.push_back(n);
			}
		if (carry != 0&& (des.size()<src->size()-discard))
			des.push_back(carry >> (shiftwidth - within)); // check is shiftwidth is 64 and witiin is zero and the result is still ok
		reverse(des.begin(), des.end());
		}
	else
		des.push_back(0);  // Discarded more than size of src. Return 0
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/Sep/2021
//	@brief 		Normalize a floating point mantissa
//	@return 	int - Return the exponent adjustment factor due to normalization
//	@param		"m"	-	digital string
//
//	@todo  
//
// Description:
//   Normalize the mantissa
//   1) If a number does not have a leading digit != 0 then left shift until 
//   it has and adjust the exponent accordingly.
//   or 2) if a number does have more than one leading digits then right shift until it has ony one
//	  and adjust the exponent accordingly
//
int _float_precision_normalize(std::vector<fptype> *m)
	{
	int expo = 0;
	size_t shift;
	const size_t offset = Bitsfptype - 1;
	if (m->size() == 1 && *m == 0)  // m=0 is a special case
		return 0;
	shift = _float_precision_clz(*m);
	if (m->size()*Bitsfptype == shift) // All zeros (also a special case)
		{
		m->erase(m->begin()+1, m->end());  return 0;
		}
	if (shift != offset)
		{
		if (shift < offset)
			{
			*m = _float_precision_right_shift(m, offset - shift); 
			expo += (int)(offset - shift);
			}
		else
			{
			*m = _float_precision_left_shift(m, shift - offset);
			expo -= (int)( shift - offset );
			}
		}
	return expo;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Round the mantisaa to significant digits and rounding control
//	@return 	int - Return the exponent adjustment (0 or 1) 
//	@param		"m"	-	digital string
// @param		"sign"   - The sign of the number
// @param		"precision" - The digital precision
// @param		"mode"   - Rounding mode 
//
//	@todo  
//
// Description:
//   Rounding control
//   Round the fraction to the number of precision based on the round mode 
//   Note that the fptype number has ALWAYS been normalize prior to rounding
//   Rounding Mode Positive numnber   Result    
//   Rounding to nearest              +�   
//   Rounding toward zero (Truncate)  Maximum, positive finite value   
//   Rounding up (toward +�)          +�   
//   Rounding down) (toward -�)       Maximum, positive finite value   
//
//   Rounding Mode Negative number    Result    
//   Rounding to nearest              -�   
//   Rounding toward zero (Truncate)  Maximum, negative finite value   
//   Rounding up (toward +�)          Maximum, negative finite value   
//   Rounding down) (toward -�) 
//		1) first check if we need to do any rounding at all
//		2) If mode == ROUND_NEAR determine if we are doing ROUND_DOWn or ROUND_UP
//		3) Discard excesive fptype digits that is beyond the precision
//		4) Discard excessive bits that is beyond precision
//		5) if Round up add one to the last bit of the precision. 
//		6) if 5) carry a digit to the most significant bit then >> shift 1 bit and return 1 for expoenent adjustment 
//		7) otherwise retun 0 for expoenent adjustment
//
int _float_precision_rounding(std::vector<fptype> *m, const int sign, const size_t precision, const enum round_mode mode)
	{
	enum round_mode rm = mode;
	const size_t pbits = (int)(ceil(precision*log2(BASE_10))); // Number of precision bits needed
	const size_t n =  pbits / Bitsfptype + 1;
	const size_t bn = pbits % Bitsfptype;
	int extra;

	if ((m->size() - 1)*Bitsfptype + 1 > pbits)  // More digits than we need 
		{
		//fptype check0; 
		size_t n1, bn1;
		switch (rm)
			{	
			case ROUND_NEAR:
				n1 = (pbits+1) / Bitsfptype + 1;
				bn1 = (pbits+1) % Bitsfptype;
				//check0 = (fptype)(1) << (Bitsfptype - bn1);  // can be removed after debug
				if (n1 < m->size() && (*m)[n1] & ((fptype)(1) << (Bitsfptype - bn1)))  // is last bit set -=> do ROUND_UP
					rm = ROUND_UP;
				else
					rm = ROUND_DOWN;
				break;
			case ROUND_UP:
				if (sign < 0)
					rm = ROUND_DOWN;
				break;
			case ROUND_DOWN:
				if (sign < 0)
					rm = ROUND_UP;
				break;
			case ROUND_ZERO:
				if (sign > 0)
					rm = ROUND_DOWN;
				else
					rm = ROUND_UP;
				break;
			}

		// Now rm is either ROUND_DOWN or ROUND_UP

		// Discard excessive fptype digits
		extra = bn == 0 ? 0 : 1;
		if(m->size()>n+extra)
			m->erase(m->begin()+n+extra,m->end());

		// Discard excessive bits within the last fptype.
		if (n < m->size())
			{
			//fptype check1 = (~(fptype)0) << (Bitsfptype - bn); // Debug can be removed after debug
			(*m)[n] &= ((~(fptype)0) << (Bitsfptype - bn));
			}

		if (rm == ROUND_UP)
			{
			const fptype bitvalue = (fptype)(1) << (Bitsfptype - bn );
			if (bitvalue == 0 )
				m->erase(m->begin() + n, m->end() );   // if we need to add to bit 0 in the previous fptype index. Erase the unneeded fptype digit
			*m = _float_precision_uadd_short( m, bitvalue );  // Add the carry
			if ((*m)[0] > (fptype)1 )
				{ // a carry was added to the most significant bit (now 2 instead of 1)
				// dont do comparison like: != 1 since it can also be zero when doing adding
				// Shift everything one to the right 
				if ((*m)[0] == (fptype)3)
					extra = extra;
				*m = _float_precision_right_shift( m, 1 );
 				if (m->size() > n+1 )
					m->erase(m->begin()+n+1, m->end() );
				_float_precision_strip_trailing_zeros(m);
				return 1;
				}
			}
		}

	_float_precision_strip_trailing_zeros(m);
	return 0;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		Compare to floating point string (mantissa only)
//	@return 	int - Return the compared result. 0==same, 1==s1>s2 or -1==s1<s2
//	@param		"s1"	- First digital binary number
// @param		"s2"	- Second digital binary number
//
//	@todo  
//
// Description:
//   Compare two unsigned binary vector<fptype> numbers
//   and return 0 is equal, 1 if s1 > s2 otherwise -1
//
int _float_precision_compare(const std::vector<fptype> *s1, const std::vector<fptype> *s2)
	{
	std::vector<fptype>::const_iterator p1, p1_end, p2, p2_end;

	for (p1 = s1->begin(), p2 = s2->begin(), p1_end=s1->end(), p2_end=s2->end(); p1 != p1_end && p2 != p2_end; ++p1, ++p2)
		{
		if (*p1 > *p2 ) return 1;
		if (*p1 < *p2 ) return -1;
		}
	// The are still the same and one or both is exhausted
	if (p1 == p1_end && p2 == p2_end)
		return 0; // is the same

	for (; p1 != p1_end; ++p1 )
		{
		if (*p1 > 0) return 1;
		}
	for (; p2 != p2_end; ++p2 )
		{
		if (*p2 > 0) return -1;
		}

	return 0;  // Same 
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		9/Sep/2021
//	@brief 		add a short integer to a floating point vector<fptype>
//	@return 	std::vector<fptype>- Return the added number
//	@param		"src"	-	The source string
// @param		"d"  - The number to add
//
//	@todo  
//
// Description:
//   Short float Add: The digit d [0..2^64] is added to the unsigned fraction string
//   Optimized 0 add or early out add is implemented
//
std::vector<fptype> _float_precision_uadd_short(const std::vector<fptype> *src, const fptype d)
	{
	fptype carry;
	std::vector<fptype>::const_reverse_iterator s_pos;
	std::vector<fptype>::reverse_iterator d_pos, d_end;
	std::vector<fptype> des;

	if (d == 0)   // Zero add
		return *src;

	carry = d;
	des = *src;		// Copy source to des1
	d_pos = des.rbegin();
	d_end = des.rend();
	s_pos = src->rbegin();

	for (; carry != 0 && d_pos != d_end; ++s_pos, ++d_pos)
		{
		*d_pos += carry;
		if (*d_pos < *s_pos)
			carry = 1;  // Set Carry
		else carry = 0;
		}

	// Exhaust the smalles of the number, so only the carry can changes the uppper radix digits
	for (; carry != 0 && d_pos != d_end; )
		{
		fptype tmp = *d_pos;
		*d_pos = tmp + carry;
		if (*d_pos < tmp)
			carry = 1;  // Set Carry
		else carry = 0;
		++d_pos;
		}

	// No more carry or end of upper radix number. 
	if (carry != 0) // If carry add the carry as a extra digit to the front of the number
		des.insert(des.begin(),1,carry);

	_int_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Sep/2021
//	@brief 		std::vector<fptype> _int_precision_uadd
//	@return		std::vector<fptype>	-	the result of adding src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Add two unsigned vector<fptype> numbers.
//   Notice src1.size() doesnt necessary needs to be of the same size as src2.size()
//	  src1[0] is the most significant elements (only containing 1bit of information), src1[src1.size()-1] the least significant
//	  same for src2.
//   Optimized: Used early out add
//
std::vector<fptype> _float_precision_uadd(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype carry = 0;
	std::vector<fptype> des;
	std::vector<fptype>::const_iterator pos;
	size_t i,j;

	if (src1->size() >= src2->size())
		{
		des = *const_cast<std::vector<fptype> *>(src1);
		pos = src2->begin();
		i = src2->size();
		}
	else
		{
		des = *const_cast<std::vector<fptype> *>(src2);
		pos = src1->begin();
		i = src1->size();
		}

	for (j=i; i > 0; --i,--j )
		{ // Adding element by element for the two numbers starting with least significant vector<fptype> of the smallest number
		des[j-1] += pos[i-1] + carry;
		carry = des[j-1] < pos[i-1] ? 1 : 0;
		}

	//  The last (most significant elements) cant overflow since there is only 1 bit of significant for a normalized number
	if (carry != 0)
		carry = carry;  // Error

	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/Sep/2021
//	@brief 		subtract a short integer from a floating point vector<fptype>
//	@return 	std::vector<fptype> - Return the subtracted vector<fptype> number
// @param		"result" -  If the number wraps around (d > src ) then result=1 otherwise 0
//	@param		"src"	-	The source string
// @param		"d"  - The number to subtract
//
//	@todo  
//
// Description:
//   Short Subtract: The fptype digit d [0..2^64] is subtracted from a unsigned vector<fptype> number
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<fptype>_float_precision_usub_short(int *result, const std::vector<fptype> *src1, const fptype d)
	{
	fptype r, borrow=0;
	std::vector<fptype> ::const_reverse_iterator pos, end;
	std::vector<fptype> des;

	if (d == 0) // Nothing to subtract
		{
		*result = 0;
		return *src1;
		}

	pos = src1->rbegin();
	end = src1->rend();
	r = *pos - (d + borrow);
	borrow = *pos < (d + borrow) ? 1 : 0;
	des.push_back(r);
	for (++pos; borrow>0 && pos != end; ++pos)
		{
		r = *pos - borrow;
		borrow = *pos < borrow ? 1 : 0;
		des.push_back(r);
		}
	_float_precision_strip_trailing_zeros(&des);
	reverse(des.begin(), des.end());
	*result = borrow > 0 ? -1 : 0;
	return des;
	}



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		12/Sep/2021
//	@brief 		subtract two floating point string
//	@return 	std::vector<fptype> - Return the subtracted vector<fptype> number
// @param		"result" -  If the number wraps around (d > src1 ) then result=1 otherwise 0
//	@param		"src1"	-	The first source string
// @param		"src2"  - The second source string
//
//	@todo  
//
// Description:
//   Subtract two unsigned vector<fptype> numbers
//   src1 or src2 need not be of the same size
//   if src1 < src2 return -1 (wrap around) otherwise return 0 (no wrap around)
//
std::vector<fptype> _float_precision_usub(int *result, const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype r, borrow = 0;
	std::vector<fptype>::const_reverse_iterator pos1, pos2;
	std::vector<fptype> des;
	size_t s1len=src1->size(), s2len=src2->size(), icur=std::max(s1len, s2len);

	if (s1len > s2len)
		des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	else
		des.reserve(src2->capacity());  // Reserver space to avoid time consuming reallocation
	pos1 = src1->rbegin();
	pos2 = src2->rbegin();
	
	for ( ; icur > 0; --icur )
		{
		if (icur<=s1len && icur <= s2len )
			{
			r = *pos1 - (*pos2 + borrow);
			borrow = *pos1 < (*pos2 + borrow) ? 1 : 
					 *pos1 == 0 ? borrow : 0;      // if borrow was not paid then propagate it to next fptype subtraction
			++pos1; ++pos2;
			}
		else
			if (icur <= s1len )
				{
				r = *pos1 - (borrow);
				borrow = *pos1 < borrow ? 1 : 0;
				++pos1;
				}
			else
				{// icur <=s2len
				r = 0 - (*pos2 + borrow);
				borrow = 1;
				++pos2;
				}
		des.push_back(r);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	*result = borrow>0 ? -1 : 0;
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_umul_short
//	@return 	std::vector<fptype> - 	the result of the short multiplication
//	@param      "src1"	-	Source string to multiply short number
//	@param      "d"	   -	Number to multiply   
//
//	@todo
//
// Description:
//   Short Mul: The digit d [0..2^64] is multiplied to the unsigned vector<fptype> number
//   Optimized Multiply with zero yields zero, Multiply with one return the original.
//
std::vector<fptype> _float_precision_umul_short(const std::vector<fptype> *src1, const fptype d)
	{
	fptype carry = 0;
	std::vector<fptype>::const_reverse_iterator pos, end;
	std::vector<fptype> des;
	std::vector<fptype> tmp(2);

	if (d == 0)  // Multiply by zero is zero.
		{
		des.push_back(0);
		return des;
		}

	if (d == 1)  // Multiply by one dont change the src1.
		{
		des = *src1;
		_float_precision_strip_trailing_zeros(&des);
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation   
	pos = src1->rbegin();
	end = src1->rend();

	for (; pos != end; ++pos)
		{
		tmp = _precision_umul64(*pos, d);
		tmp[0] += carry;
		carry = (tmp[0] < carry) ? tmp[1] + 1 : tmp[1];
		des.push_back(tmp[0]);
		}

	if (carry != 0)
		des.push_back(carry);
	reverse(des.begin(), des.end()); 
	_float_precision_strip_trailing_zeros(&des);

	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>  _float_precision_umul
//	@return		std::vector<fptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<fptype>. Brach out to the different multiplication algorithm based on size of the operands for optimized performance
//
std::vector<fptype> _float_precision_umul(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	const size_t s1len = src1->size(), s2len = src2->size();
	std::vector<fptype> des;

	if (*src1->begin() == 0 || *src2->begin() == 0)  // zero can only arise if the number is zero and therefore no need to check the size
		{
		des.assign(1, 0);
		}
	else
		{  // Both s1[0] and s2[0] == 1 check size to dermine if it is 1 or any power of 2
		if (s1len == 1)
			{
			des = *src2;  // Mutiply with 1 or any true power of 2
			}
		else
			if (s2len == 1)
				{
				des = *src1; // Mutiply with 1 or any true power of 2
				}
			else
				{// src1 and src2 size > 1. 
				if (s2len == 2 || s1len + s2len < 20)	// Measured for best performance .
					{
					des = _float_precision_umul_school(src1, src2);
					}
				else
					if (s1len+s2len < 6000)
						{
						des = _float_precision_umul_linear(src1, src2);
						}
					else
						{
						des = _float_precision_umul_fourier(src1, src2);
						}
				}
		}
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>  _float_precision_umul_school
//	@return		std::vector<fptype> -	the result of multiplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
// Multiply two unsigned vector<fptype>.
//	Not used anymore since the complexity is o(n^2)
//
std::vector<fptype> _float_precision_umul_school(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	fptype carry;
	std::vector<fptype> des, tmp;
	std::vector<fptype>::const_reverse_iterator pos1, pos2;

	des.insert(des.begin(), src1->size() + src2->size(), 0);
	size_t i = src2->size(), j = src1->size(), k;
	tmp.assign(2, 0);
	for (pos2 = src2->rbegin(); i > 0; ++pos2, --i)
		{
		carry = 0; j = src1->size();
		for (pos1 = src1->rbegin(), k = i + j; j>0 && k >= i; ++pos1, --k, --j)
			{
			if (*pos2 == 1)
				{
				tmp[0] = *pos1; tmp[1] = 0;
				}
			else
				if (*pos1 == 1)
					{
					tmp[0] = *pos2; tmp[1] = 0;
					}
				else
					tmp = _precision_umul64(*pos1, *pos2);
			tmp[0] += carry;
			tmp[1] += tmp[0] < carry ? 1 : 0;
			des[k - 1] += tmp[0];
			carry = tmp[1];
			carry += des[k - 1] < tmp[0] ? 1 : 0;  //can't overflow by just adding 1 since tmp[1] will be sufficient less than maximum number
			}
		if (carry != 0)  // Note: no overflow for last carry since number can't be bigger than the sum of the two digits
			des[k - 1] = carry;
		}
	_float_precision_strip_leading_zeros(&des);
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_umul_fourier
//	@return 	std::vector<fptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//		Multiply two unsigned binary numbers
//		Optimized: Used FFT algorithm to performed the multiplication
//		Plus added multi threading speeding up the calculation
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65 Only 10^5 decimal digits is usualy not enough for arbitrary precsion so we are using a max of 1 byte.	
//
std::vector<fptype> _float_precision_umul_fourier(const std::vector<fptype> *src1, const std::vector<fptype> *src2, int nbits)
	{
	std::vector<fptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits==0?8:nbits;
	int radix = bits==8 ? 256: 16;
	size_t sz = sizeof(fptype);
	std::vector<double> va, vb;

	l1 = src1->size();
	l2 = src2->size();
	des.reserve(l1 + l2 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1 < l2 ? l2 : l1;
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(fptype);  // Convert to byte size 
	if (l > 8'000'000*sizeof(fptype) || bits==4)  
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2; // use 2^4 instead of 2^8
		}
	for (n = 1; n < l; n <<= 1)  ;
	n <<= 1;

#ifdef HVE_THREAD
	// Using parallel sections below speeds up the performance of the two calls to _int_real_Fourier() with a factor of 1.8 
	if (nbits == 0|| l1+l2>10'000)
		{// Starting thread using lambda expressions
		// L1, l2, va, vb by reference since it is used after the thread has terminated
		std::thread first( [&, n, bits]() 
			{std::vector<fptype>::const_iterator pos, end;
			size_t i;
			va.resize(n);
			for (i = 0, pos = src1->begin(), end = src1->end(); pos != end; ++pos)
				i += convertbinary2double(&va[i], *pos, i == 0, bits);
			l1 = i; // L1 now Number of bytes or nibbles
			_vector_real_fourier(/*std::ref*/(va), n, 1); // FFT va
			} );
		
		std::thread second([&, n, bits]() 
			{std::vector<fptype>::const_iterator pos, end;
			size_t i;
			vb.resize(n);
			for (i= 0, pos = src2->begin(), end = src2->end(); pos != end; ++pos)
				i += convertbinary2double(&vb[i], *pos, i == 0, bits);
			l2 = i; // L2 now Number of bytes or nibbles
			_vector_real_fourier(std::ref(vb), n, 1); // FFT vb
			});

		first.join();
		second.join();
		}
	else
#endif
		{
		std::vector<fptype>::const_iterator pos, end;
		va.resize(n);
		vb.resize(n);
		// Start with most significant fptype e.g. src1[0]
		for (l1 = 0, pos = src1->begin(), end = src1->end(); pos != end; ++pos)
			l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
		// L1 now Number of bytes or nibbles
		// Start with most significant fptype e.g. src2[0]
		for (l2 = 0, pos = src2->begin(), end = src2->end(); pos != end; ++pos)
			l2 += convertbinary2double(&vb[l2], *pos, l2 == 0, bits);
		// L2 now number of bytes or nibbles
 		_vector_real_fourier(va, n, 1); // FFT va
		_vector_real_fourier(vb, n, 1); // FFT vb
		}

	// Do the multiplication in the frequence domain
	vb[0] *= va[0];
	vb[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t;
		vb[j] = (t = vb[j])*va[j] - vb[j + 1] * va[j + 1];
		vb[j + 1] = t*va[j + 1] + vb[j + 1] * va[j];
		}

	_vector_real_fourier(vb, n, -1); // Reverse FFT vb
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = vb[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		vb[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l2 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;				// Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		fptype num;
		size_t inx = l1 - sz *(l2 - l + 1);
		num = convertdouble2binary(&vb[inx], sz, 0, bits);
		des.push_back(num);
		//std::cout << std::hex << num << std::endl;// DEBUG
		}
	l2 = l1 % sz;			// Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)	// do the the last 64bit integers from b
		{
		fptype num;
		num = convertdouble2binary(&vb[0], l2, cy,bits);
		des.push_back(num);
		//std::cout << std::hex << num  << std::endl;  // DEBUG
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	
	return des;
	}

	//	@author Henrik Vestermark (hve@hvks.com)
	//	@date		25/Dec/2021
	//	@brief 		std::vector<fptype> _float_precision_schonhage_strassen_linear_umul
	//	@return		std::vector<fptype>	-	the result of multiplying src1 and src2
	//	@param		"lhs"	-	First unsigned source argument
	//	@param		"rhs"	-	Second unsigned source argument
	//
	//	@todo
	//
	// Description:
	//   Multiply two unsigned decimal strings, using the Schonhage-Strassen (linear convolution) method
	//   Using the full advantages of Half bit size of fptype in the linear convolution
	//
	std::vector<fptype> _float_precision_umul_linear(const std::vector<fptype> *lhs, const std::vector<fptype> *rhs)
		{
		size_t i, j;
		size_t l_length = lhs->size(), r_length = rhs->size();
		const unsigned int HalfBitsfptype = Bitsfptype >> 1;  // same as / 2
		const fptype mask = (~(fptype)0) >> HalfBitsfptype;
		const uintmax_t radix = (uintmax_t)1 << HalfBitsfptype;
		std::vector<fptype> des;
		std::vector<fptype>::const_iterator pos, end;
		std::vector<uintmax_t> linearconvolution(2 * (l_length + r_length), 0);  // initialize it with zero
		std::vector<fptype> ua(l_length * 2), ub(r_length * 2);
		uintmax_t nextCarry = 0; 

		// Convert to half fptype from vector<fptype> and notice we dont stored in reverse order as we did for integers
		// by first converting lhs onto ua and then rhs into ub
		// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half fptype  from fptype by mapping each mBinary number into 2 half fptype numbers 
		// the function convertbinary2Halfiptype() does this job per mBinary fptype number
		for (i = 0, pos = lhs->begin(), end = lhs->end(); pos != end; ++pos)
			i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
		l_length = i;  // l_length now in half iptype's  instead of iptype's
		for (j = 0, pos = rhs->begin(), end = rhs->end(); pos != end; ++pos)
			j += convertbinary2Halfiptype(&ub[j], *pos, j == 0);
		r_length = j;  // l_length now in half iptype's instead of iptype's

		// do the linear convolution
		for (i = 0; i < r_length; ++i)
			for (j = 0; j < l_length; ++j)
				{
				uintmax_t m = (uintmax_t)ub[r_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
				linearconvolution[i + j] += m;
				if (linearconvolution[i + j] < m) // carry
					{
					// Propagate carry
#ifdef HVE_DEBUG
					//std::cout << "Propagating overflow at indx=" << i + j << std::endl;
#endif
					for (size_t k = i + j + 1; ; ++k)
						{
						linearconvolution[k] += radix;	// Add carry

						if (linearconvolution[k] >= radix) // Continue adding carry ?
							break;
#ifdef HVE_DEBUG
					//	else
					//		std::cout << "Propagating overflow at indx=" << k << " Radix=" << radix << " Val=" << linearconvolution[k] << std::endl;
#endif
						}
					}
				}

 		des.reserve(r_length + l_length + 2);
		for (i = 0; i < l_length + r_length - 1; ++i)
			{	
			linearconvolution[i] += nextCarry;
#ifdef HVE_DEBUG
	//		if (linearconvolution[i] < nextCarry)
	//			std::cout << "Final overflow at index=" << i  << " NextCarry=" << nextCarry << " Index value="<< linearconvolution[i] << std::endl;
#endif
			
			nextCarry = linearconvolution[i] >> HalfBitsfptype;  //  same as / radix;
		//	if (nextCarry != 0)
		//		nextCarry = nextCarry;
			linearconvolution[i] &= mask;  // same as %= radix;
			}
		if (nextCarry != 0)
			linearconvolution[i++] = nextCarry & mask; // same as % radix;

		//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
		// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
		// do the full 64bit integers first starting from least significant HalfBitsiptype
		for (j = 0; j < i; j += 2)
			{
			fptype num;
			num = convertHalfiptype2binary(&linearconvolution[j], 2);
			des.push_back(num);
			}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

///////////////////////////////////////////////////////////////////////////////////////
//
//
//	Specialized function for squaring using fourier and linear
//
//
///////////////////////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_umul2_fourier
//	@return 	std::vector<fptyp> -	the result of multplying src1 and src2
//	@param		"src1"	-	unsigned source argument
//
//	@todo
//
// Description:
//		squaring the unsigned binary number
//		Optimized: Used FFT algorithm to performed the multiplication
//		Plus added multi threading speeding up the calculation
//		Since we convert the binary numbers into float we have to ensure proper accuracy in the calculation.
//		In numerical recipies in C (2nd edition, chaper 20.6) they state that using double the equations that need to be fulfilled for accuracy
//		1byte binary:
//			log2(256^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^9 digits
//			16+30+4.9=50.9  which should be just Ok for 1 byte binary digits.
//		2byte binary:
//			log2(256^2^2)+log2(N)+safeguard*log2(log2(N))< digits in double which is 53. using safeguard as 1 we get for N=10^5 digits
//			32+16.6+4.05=52.65 Only 10^5 decimal digits is usualy not enough for arbitrary precsion so we are using a max of 1 byte.	
//
std::vector<fptype> _float_precision_umul2_fourier(const std::vector<fptype> *src, int nbits)
	{
	std::vector<fptype> des;
	size_t n, l, l1, l2, j;
	double cy;
	int bits = nbits == 0 ? 8 : nbits;
	int radix = bits == 8 ? 256 : 16;
	size_t sz = sizeof(fptype);
	std::vector<double> va;

	l1 = src->size();
	des.reserve(l1 + l1 + 16);  // Ensure enough space to hold the Multiplication result to avoid reallocation of des
	l = l1;
	// Since we split the 64bit numbers into chunk of 8bit to ensure we have enough accuray when using double 
	l *= sizeof(fptype);  // Convert to byte size 
	if (l > 8'000'000 * sizeof(fptype) || bits == 4)
		{
		bits = 4; l <<= 1; radix = 16; sz *= 2; // use 2^4 instead of 2^8
		}
	for (n = 1; n < l; n <<= 1);
	n <<= 1;

	std::vector<fptype>::const_iterator pos, end;
	va.resize(n);
	// Start with most significant fptype e.g. src1[0]
	for (l1 = 0, pos = src->begin(), end = src->end(); pos != end; ++pos)
		l1 += convertbinary2double(&va[l1], *pos, l1 == 0, bits);
	// L1 now Number of bytes or nibbles
	_vector_real_fourier(va, n, 1); // FFT va

	// Do the multiplication in the frequence domain
	va[0] *= va[0];
	va[1] *= va[1];
	for (j = 2; j < n; j += 2)
		{
		double t = va[j], u = va[j + 1];
		va[j] = t*t - u*u;
		va[j + 1] = 2*t*u;
		}

	_vector_real_fourier(va, n, -1); // Reverse FFT vb
	for (cy = 0, j = 0; j <= n - 1; ++j)
		{
		double t;
		t = va[n - 1 - j] / (n >> 1) + cy + 0.5;
		cy = (unsigned long)(t / radix);  // Byte Radix 2^8 or 2^4
		va[n - 1 - j] = t - cy * radix;
		}

	// Now collect then back into a vector<fptype> format
	l1 += l1 - 1;				// max number of bytes or nibbles plus the carry
	l2 = l1 / sz;				// Number of full 64bit digits
	for (l = l2; l > 0; --l)	// do the full 64bit integers first starting backwards from b
		{
		fptype num;
		size_t inx = l1 - sz *(l2 - l + 1);
		num = convertdouble2binary(&va[inx], sz, 0, bits);
		des.push_back(num);
		}
	l2 = l1 % sz;			// Number of remaing 8bits or 4bits digits
	if (l2>0 || cy != 0)	// do the the last 64bit integers from b
		{
		fptype num;
		num = convertdouble2binary(&va[0], l2, cy, bits);
		des.push_back(num);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);

	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		24/Jan/2022
//	@brief 		std::vector<fptype> _float_precision_schonhage_strassen_linear_umul
//	@return		std::vector<fptype>	-	the result of multiplying src1 and src2
//	@param		"lhs"	-	unsigned source argument
//
//	@todo
//
// Description:
//   square the the unsigned decimal string, using the Schonhage-Strassen (linear convolution) method
//   Using the full advantages of Half bit size of fptype in the linear convolution
//
std::vector<fptype> _float_precision_umul2_linear(const std::vector<fptype> *src )
	{
	size_t i, j;
	size_t l_length = src->size();
	const unsigned int HalfBitsfptype = Bitsfptype >> 1;  // same as / 2
	const fptype mask = (~(fptype)0) >> HalfBitsfptype;
	const uintmax_t radix = (uintmax_t)1 << HalfBitsfptype;
	std::vector<fptype> des;
	std::vector<fptype>::const_iterator pos, end;
	std::vector<uintmax_t> linearconvolution(2 * (l_length + l_length), 0);  // initialize it with zero
	std::vector<fptype> ua(l_length * 2);
	uintmax_t nextCarry = 0;

	// Convert to half fptype from vector<fptype> and notice we dont stored in reverse order as we did for integers
	// by first converting lhs onto ua and then rhs into ub
	// e.g. lhs=a0+a1*R+a2*R^2,...an-1*R^n-1, a0 can be subdivied into half fptype  from fptype by mapping each mBinary number into 2 half fptype numbers 
	// the function convertbinary2Halfiptype() does this job per mBinary fptype number
	for (i = 0, pos = src->begin(), end = src->end(); pos != end; ++pos)
		i += convertbinary2Halfiptype(&ua[i], *pos, i == 0);
	l_length = i;  // l_length now in half iptype's  instead of iptype's

	// do the linear convolution
	for (i = 0; i < l_length; ++i)
		for (j = 0; j < l_length; ++j)
			{
			uintmax_t m = (uintmax_t)ua[l_length - 1 - i] * (uintmax_t)ua[l_length - 1 - j];
			linearconvolution[i + j] += m;
			if (linearconvolution[i + j] < m) // carry
				{
				// Propagate carry
				for (size_t k = i + j + 1; ; ++k)
					{
					linearconvolution[k] += radix;	// Add carry
					if (linearconvolution[k] >= radix) // Continue adding carry ?
						break;
					}
				}
			}

	des.reserve(l_length + l_length + 2);
	for (i = 0; i < l_length + l_length - 1; ++i)
		{
		linearconvolution[i] += nextCarry;
		nextCarry = linearconvolution[i] >> HalfBitsfptype;  //  same as / radix;
		//	if (nextCarry != 0)
		//		nextCarry = nextCarry;
		linearconvolution[i] &= mask;  // same as %= radix;
		}
	if (nextCarry != 0)
		linearconvolution[i++] = nextCarry & mask; // same as % radix;

	//linearconvolution now holds the result with [0] as the most significant byte number and [i-1] as the least sinificant byte as Halfbitsiptype numbers
	// Now convert then back into a vector<iptype> format. i is the number of HalfBitsiptype's in the result
	// do the full 64bit integers first starting from least significant HalfBitsiptype
	for (j = 0; j < i; j += 2)
		{
		fptype num;
		num = convertHalfiptype2binary(&linearconvolution[j], 2);
		des.push_back(num);
		}
	reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	return des;
	}

///////////////////////////////////////////////////////////////////////////////////////
//
//
//	END squaring functions
//
//
///////////////////////////////////////////////////////////////////////////////////////

	 
// Short Division: The fptype digit d  is divide up into the unsigned fptype vector
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype>	_float_precision_udiv_short
//	@return 	std::vector<fptype>	- The result of the short division
//	@param      "src1"				- Source string to divide with the short number
//	@param      "d"					- Number to divide
// @param		"remaind"			- The remaind of the short division
//
//	@todo
//
// Description:
//   Short divide: The fptype digit d [0..2^32] is divided up in the unsigned vector<fptype> 
//	  Notice only up to int 32bit can be handle as short div.
//   Divide with zero throw an exception
//
std::vector<fptype> _float_precision_udiv_short(fptype *remaind, const std::vector<fptype> *src1, const fptype d)
	{
	const unsigned int shifts = 4 * sizeof(fptype);
	const fptype mask = (~((fptype)0)) >> shifts;
	fptype ir;
	std::vector<fptype>::const_iterator s_pos, s_end;
	std::vector<fptype> des;

	if (d == 0)
		{
		throw float_precision::divide_by_zero();
		}

	if (d == 1)  // Divide by one dont change the src1.
		{
		des = *const_cast<std::vector<fptype> *> (src1);
		_float_precision_strip_trailing_zeros(&des);
		*remaind = 0;
		return des;
		}

	des.reserve(src1->capacity());  // Reserver space to avoid time consuming reallocation
	s_pos = src1->begin();
	s_end = src1->end();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		fptype n, qh, ql;
		/*if (ir == 0)
		{// Shortcut when carry is zero
		des.push_back(*s_pos / d );
		ir = *s_pos % d;
		}
		else*/
		{
			n = *s_pos >> shifts;
			n |= ir << shifts;
			qh = n / d;	ir = n % d;
			n = *s_pos & mask;
			n |= ir << shifts;
			ql = n / d;	ir = n % d;
			n = (qh << shifts) | ql;
			des.push_back(n);
		}
		}

	//reverse(des.begin(), des.end());
	_float_precision_strip_trailing_zeros(&des);
	*remaind = ir;
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_udiv
//	@return		std::vector<fptype>-	the result of disivison
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Divide two unsigned binary numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<fptype> _float_precision_udiv(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	int plusdigit, plusbits, wrap, i;
	fptype underflow;
	std::vector<fptype> des, tmp, quotient, divisor;

	des.push_back(0);
	divisor = *const_cast<std::vector<fptype> *> (src1);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short div if denominator <= 32 bit integer.
		return _float_precision_udiv_short(&underflow, &divisor, src2->front());

	plusdigit = (int)divisor.size() - (int)src2->size();
	if (plusdigit < 0)  //src1 / src2 == 0
		return des;

	plusbits = (int)_float_precision_clz(src2->back()) - (int)_float_precision_clz(divisor.back());
	plusbits = plusdigit * Bitsfptype + plusbits;
	for (i = 0; plusbits >= 1; ++i)
		{
		tmp = _int_precision_ushiftleft(src2, plusbits);
		if (_int_precision_compare(&divisor, &tmp) < 0)
			{ // Too much reduce with one power of radix
			--plusbits; continue;
			}
		divisor = _int_precision_usub(&wrap, &divisor, &tmp);
		quotient.clear();
		quotient.insert(quotient.begin(), (plusbits / (Bitsfptype)) + 1, 0);
		quotient[quotient.size() - 1] = (fptype)(1) << (plusbits % (Bitsfptype));
		des = _int_precision_uadd(&des, &quotient);
		}

	for (wrap = 0; wrap == 0; )
		{
		divisor = _int_precision_usub(&wrap, &divisor, src2);
		if (wrap == 0) // src1 was indeed > src2
			des = _int_precision_uadd_short(&des, 1);
		}

	_int_precision_strip_trailing_zeros(&des);
	return des;
	}

// Short Remainder: The fptype digit d [1..2^64] is divide up into the unsigned vector<fptype> and the remaing is returned
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_urem_short
//	@return 	std::vector<fptype> - 	the result of the short remainder
//	@param      "src1"	-	Source string to divide with the short number
//	@param      "d"	   -	Number to divide
//
//	@todo
//
// Description:
//   Short remainder: The fptype digit d [0..2^64] is divided up in the unsigned vector<fptype> and the remaing is retuened
//   Divide with zero throw an exception
//   if d==1 then result == 0, for d==2,4,5,8,10 we only test the last few digits to get the result. This speed up rem for large integers with small rem value
//   since we dont have to run through every digits in src1
//
std::vector<fptype> _float_precision_urem_short(const std::vector<fptype> *src1, const fptype d)
	{
	const unsigned int shifts = 4 * sizeof(fptype);
	const fptype mask = (~((fptype)0)) >> shifts;
	fptype ir;
	std::vector<fptype>::const_reverse_iterator s_pos, s_end;
	std::vector<fptype> des;

	if (d == 0)
		{
		throw float_precision::divide_by_zero();
		}

	if (d == 1)  // Remainer is always 0 for d==1
		{
		des.push_back(0);
		return des;
		}

	// Short cut
	ir = *src1->begin();
	switch (d)
		{
		case 2: des.push_back(ir % 2); return des; break;
		case 4: des.push_back(ir % 4); return des; break;
		case 5: des.push_back(ir % 5); return des; break;
		case 8: des.push_back(ir % 8); return des; break;
		case 10: des.push_back(ir % 10); return des; break;
		default:;   // No Short cut
		}

	s_pos = src1->rbegin();
	s_end = src1->rend();
	for (ir = 0; s_pos != s_end; ++s_pos)
		{
		fptype n;
		if (ir == 0)
			{// Shortcut when carry is zero
			ir = *s_pos % d;
			}
		else
			{
			n = *s_pos >> shifts;
			n += ir << shifts;
			ir = n % d;
			n = *s_pos & mask;
			n += ir << shifts;
			ir = n % d;
			}
		}

	des.push_back(ir);
	return des;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		13/Sep/2021
//	@brief 		std::vector<fptype> _float_precision_urem
//	@return		std::vector<fptype>	-	the remaing result of divide src1 with src2
//	@param		"src1"	-	First unsigned source argument
//	@param		"src2"	-	Second unsigned source argument
//
//	@todo
//
// Description:
//   Find the remainder when divide two unsigned vector<fptype> numbers
//   Optimized: Used early out add and multiplication w. zero
//
std::vector<fptype> _float_precision_urem(const std::vector<fptype> *src1, const std::vector<fptype> *src2)
	{
	int wrap;
	std::vector<fptype> des, tmp;

	des.push_back(0);
	if (src2->size() == 1 && (src2->front() >> 32) == 0) // Make short rem 
		{
		fptype rem;
		_float_precision_udiv_short(&rem, src1, src2->front());
		des[0] = rem;
		return des;
		}

	tmp = _float_precision_udiv(src1, src2);
	tmp = _float_precision_umul(&tmp, src2);
	des = _float_precision_usub(&wrap, src1, &tmp);

	_float_precision_strip_trailing_zeros(&des);

	return des;
	}

//////////////////////////////////////
//
//	END of CORE Functions. BINARY
//
//////////////////////////////////////




///////////////////////////////////////
//
// END FLOATING POINT CORE FUNCTIONS
//
//
///////////////////////////////////////


///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//
///////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Oct/2021
//	@brief 		Calculate the inverse of a 
//	@return 	float_precision -	Return 1/a
//	@param      "a"	-	The float_precision number to inverse
//
//	@todo  
//
// Description:
//   Inverse of V
//   Using a Newton iterations Un = U(2-UV)
//   Always return the result with 2 digits higher precision that argument
//   _float_precision_inverse() return a interim result for a basic operation like /
//
float_precision _float_precision_inverse( const float_precision& a )
   {
   const size_t extra=3;
   size_t precision;
   eptype expo;
   double fu;
   float_precision r, u, v;
   const float_precision c1(1), c2(2);

   if (a.iszero() == true)
	{ throw float_precision::divide_by_zero(); }
   precision = a.precision();  
   v.precision( precision + extra );
   v = a;
   expo = v.exponent();
   v.exponent( 0 );
   r.precision( precision + extra ); // Do iteration using 3 digits higher precision
   u.precision( precision + extra );

 	// New. Get a initial guess using ordinary floating point
   fu = (double)v;
   fu = 1 / fu;
   u = float_precision( fu );
   
   // Now iterate using Netwon Un=U(2-UV)
   for(;;)
      {
      r = u * v;                 // UV
      r = c2-r;                  // 2-UV
      u *= r;                    // Un=U(2-UV)
	  r.precision(precision);
	  if (r == c1)
		  break;
	  r.precision(precision + extra);
      }

   u.exponent( u.exponent() - expo );
   u.mode( a.mode() );
   //fu = (double)u;			// DEBUG
   return u;
   }

// Float Precision support functions

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Oct/2021
//	@brief 		Calculate sqrt(x)
//	@return 	float_precision -	Return sqrt(a)
//	@param      "x"	-	The sqrt argument
//
//	@todo  
//
// Description:
//   sqrt(V)
//   Equivalent with the same standard C function call
//   Separate exponent. e.g. sqrt(V*2^x)=2^x/2*sqrt(V)
//   Un=0.5U(3-VU^2)
//   Then Un == 1/Sqrt(V). and sqrt(V) = VUn
// The functionhas been improved using Newton with iterative deepening creating
//	a speed up with a factor of 3 over the classic Newton method.
//
float_precision sqrt(const float_precision& x)
	{
	const unsigned int extra = 2;
	size_t precision, digits;
	eptype expo, expo_sq;
	double fv;
	float_precision r, u, v;
	const float_precision c1(1), c3(3), c05(0.5);

	if (x.iszero() || x == c1)  // Simple squareroot
		return x;

	if (x.sign() < 0)
		{ throw float_precision::domain_error(); }
	precision = x.precision();
	if (x == c3)
		return _float_table(_SQRT3, precision);
	v.precision(precision + extra);
	v = x;
	expo = v.exponent();
	expo_sq = expo / 2;
	v.exponent(expo - 2 * expo_sq);
	r.precision(precision + extra); // Do iteration using 2 digits higher precision
	u.precision(precision + extra);

	// Get a initial guess using ordinary floating point
	// New. Get a initial guess using ordinary floating point
	fv = v;				// Convert to double	
	fv = 1 / sqrt(fv);  // set the initial guess with at approx 16 correct digits
	u = float_precision(fv);
	// Now iterate using Netwon Un=0.5U(3-VU^2)
	for (digits = std::min((size_t)32, precision); ; digits = std::min(precision + extra, digits * 2))
		{
		// Increase precision by a factor of two for the working variable s r & u. 
		r.precision(digits);
		u.precision(digits);
		// Notice V is the original number to squareroot which has the full precision 
		// so we start by assigning it to r, rounding it to the precision of r
		r = v;						// V
		r *= u * u;					// VU^2
 		r = c3 - r;					// 3-VU^2
		r *= c05;					// (3-VU^2)/2
		u *= r;						// U=U(3-VU^2)/2
		if (digits == precision + extra) // Reach final iteration step in regards to precision
			{
			r.precision(precision+1);	// round to final precision
			if (r == c1)	// break if no improvement
				break;
			r.precision(precision + extra);
			}
		}

	u *= v;
	u.exponent(u.exponent() + expo_sq);
	// Round to same precision as argument and mrounding mode
	u.mode(x.mode());
	u.precision(precision);

	return u;
	}


////////////////////////////////////////////////////////////////////////////////////////
//
// FLOAT PRECISION
//    Universal Constants LN2, LN10, e, _SQRT2, _SQRT3, _INVSQRT2, _INVSQRT3 and PI
//
///////////////////////////////////////////////////////////////////////////////////////

// Spigot function for internal calculation of transcendental constants

#if false
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		29/Jan/2017 
//	@brief 		Calculate transcendetal constant of pi
//	@return		std::string -	return the constant as a standard std::string
//	@param		"digits"	-	Number of digits
// @param		"no_dig"	-	Number of digits calculated per loop
//
//	@todo
//
// Description:
//
// 64bit version of the spigot algorithm.
// Notice acc, a, g needs to be unsigned 64bit. 
// Emperisk for pi to 2^n digits, acc need to hold approx 2^(n+17) numbers. while a[] and g needs approx 2^(n+3) numbers
// a[] & g could potential be unsigned long (32bit) going to a max of 2^29 digit or 536millions digit of PI. but with 
// unsigned 64bit you can do "unlimited"
//
// This function is not used anymore and has been replace by the Brent-Salamin algorihm for PI
//
static std::string spigot_pi_64(const int digits, int no_dig = 4)
	{
	static unsigned long f_table[] = { 0, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000 };
	static unsigned long f2_table[] = { 0,  2,  20,  200,  2000,  20000,  200000,  2000000,  20000000 };
	const int TERMS = (10 * no_dig / 3 + 1);
	bool first_time = true;									// First time in loop flag
	bool overflow_flag = false;								// Overflow flag
	char buffer[32];
	std::string ss;											// The String that hold the calculated PI											// Timer
	long b, c;												// Loop counters
	int carry, no_carry = 0;								// Outer loop carrier, plus no of carroer adjustment counts
	unsigned long f, f2;									// New base 1 decimal digits at a time
	unsigned long dig_n = 0;								// dig_n holds the next no_dig digit to add
	unsigned long e = 0;									// Save previous 4 digits
	unsigned long long acc = 0, g = 0, tmp64;
	ss.reserve(digits + 16);								// Pre reserve the string size to be able to accumulate all digits plus 8
	if (no_dig > 8) no_dig = 8;								// ensure no_dig<=8
	if (no_dig < 1) no_dig = 1;								// Ensure no_dig>0
	c = (digits / no_dig + 1) * no_dig;						// Since we do collect PI in trunks of no_dig digit at a time we need to ensure digits is divisble by no_dig.
	if (no_dig == 1) c++;									// Extra guard digit for 1 digit at a time.
	c = (c / no_dig + 1) * TERMS;							// c ensure that the digits we seek is divisble by no_dig 
	f = f_table[no_dig];									// Load the initial f
	f2 = f2_table[no_dig];									// Load the initial f2

	unsigned long long  *a = new unsigned long long[c];		// Array of 4 digits decimals
															// b is the nominator previous base; c is the index
	for (; (b = c -= TERMS) > 0 && overflow_flag == false; first_time = false)
		{
		for (; --b > 0 && overflow_flag == false;)
			{
			if (acc > ULLONG_MAX / b) overflow_flag = true;		// Check for overflow
			acc *= b;											// Accumulator *= nom previous base
			tmp64 = f;
			if (first_time == true)								// Test for first run in the main loop
				tmp64 *= f2;									// First outer loop. a[b] is not yet initialized
			else
				tmp64 *= a[b];									// Non first outer loop. a[b] is initialized in the first loop
			if (acc > ULLONG_MAX - tmp64) overflow_flag = true;	// Check for overflow
			acc += tmp64;										// add it to accumulator
			g = b + b - 1;										// denominated previous base
			a[b] = acc % g;										// Update the accumulator
			acc /= g;											// save carry
			}
		dig_n = (unsigned long)(e + acc / f);					// Get previous no_dig digits. Could occasinaly be no_dig+1 digits in which case we have to propagate back the extra digit.
		carry = (unsigned)(dig_n / f);							// Check for extra carry that we need to propagate back into the current sum of PI digits
		dig_n %= f;												// Eliminate the extra carrier so now l contains no_dig digits to add to the string
																// Add the carrier to the existing number for PI calculate so far.
		if (carry > 0)
			{
			++no_carry;											// Keep count of how many carrier detect
			for (size_t i = ss.length(); carry > 0 && i > 0; --i)	// Loop and propagate back the extra carrier to the existing PI digits found so far
				{												// Never seen more than one loop here but it can handle multiple carry back propagation 
				int new_digit;
				new_digit = (ss[i - 1] - '0') + carry;			// Calculate new digit
				carry = new_digit / 10;							// Calculate new carry if any
				ss[i - 1] = new_digit % 10 + '0';				// Put the adjusted digit back in our PI digit list
				}
			}
		//return s = uitostring10((uintmax_t)src[0]);
		(void)snprintf(buffer, sizeof(buffer), "%0*lu", no_dig, dig_n);		// Print previous no_dig digits to buffer
		ss += std::string(buffer);								// Add it to PI string
		if (first_time == true)
			ss.insert(1, ".");									// add the decimal pointafter the first digit to create 3.14...
		acc = acc % f;											// save current no_dig digits and repeat loop
		e = (unsigned long)acc;
		}

	ss.erase(digits + 1);											// Remove the extra digits that we didnt requested but used as guard digits
	if (overflow_flag == true)
		ss = std::string("Overflow:") + ss;						// Set overflow in the return string
	delete [] a;												// Delete the a[];	
	return ss;													// Return Pi with the number of digits
	}

// End Spigot PI
#endif

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		29/Jan/2017
//	@brief 		Calculate transcendetal constant of exp(1)
//	@return		std::string -	return the constant as a standard std::string
//	@param		"digits"	-	Number of digits
//
//	@todo
//
// Description:
//
// Spigot algorithm for e
// From The computer Journal 1968 (A H J Sale) written in Algo 60 and ported with some modification 
// to c++
static std::string spigot_e( const size_t digits)
	{
	unsigned int m;
	unsigned int tmp, carry;
	double test = (digits+1) * log(10);
	bool first_time = true;
	unsigned int *coef;
	std::string ss("2.");
	ss.reserve(digits + 16);
	double xnew, xold;

	// Stirling approximation of m!~Sqrt(2*pi*digits)(digits/e)^digits.
	// Taken ln on both side you get: m*(Math.log((m)-1)+0.5*Math.log(2*Math.pi*m);
	// Use Newton method to find in less that 4-5 iteration
	for (xold = 5, xnew = 0; ; xold = xnew)
		{
		double  f = xold*(log(xold) - 1) + 0.5*log(2 * 3.141592653589793 * xold);
		double f1 = 0.5 / xold + log(xold);
		xnew = xold - (f - test) / f1;
		if ((int)ceil(xnew) == (int)ceil(xold))
			break;
		}
	m = (unsigned int)ceil(xnew);
	if (m < 5)
		m = 5;
	coef = new unsigned int[m + 1];

	for (size_t i = 1; i < digits; ++i, first_time = false)
		{
		carry = 0;
		for (unsigned int j = m; j >= 2; j--)
			{
			if (first_time == true)
				tmp = 10;
			else
				tmp = coef[j] * 10;
			tmp += carry;
			carry = tmp / (j);
			coef[j] = tmp % (j);
			}
		ss.append(1, (char)(carry + '0'));
		}
	delete [] coef;
	return ss;
	}

// End Spigot e

// Spigot LN(X/Y) where x and y are integers and x>0 && x>y as conditions

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		28/Jan/2017
//	@brief 		Calculate transcendetal constant of ln(x/y)
//	@return		std::string -	return the constant as a standard std::string
//	@param		"x"	-	The nominator of the number x
// @param		"y"	-	The Denominator of the number x
//	@param		"digits"	-	Number of digits
// @param		"no_dig"	-	Number of digits calculated per loop
//
//	@todo
//
// Description:
//
// 64 bit version of spigot algorithm for LN(x/y) fraction 
// It has automatic 64bit integer overflow detection in which case the result start with the string "Overflow...."
// A Column: x-1,x-1,x-1,...,x-1
// B Column: x,x,x,x,x,...,x
// Initialization values: (x-1)/(x(n+1))...
// The function is declare static since it only serve as a sub function for the function _float_table()
//
static std::string spigot_lnxy_64(const unsigned int x, const unsigned int y, const size_t digits, int no_dig = 1)
	{
	static unsigned long f_table[] = { 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000 };
	bool first_time = true;				// First iteration of the algorithm
	bool overflow_flag = false;			// 64bit integer overflow flag
	//char buffer[32];
	std::string ss;						// The std::string that holds the ln(x)
	size_t dig;
	unsigned int car;
	size_t no_terms;				// No of terms to complete as a function of digits
	unsigned long f;					// New base 1 decimal digits at a time
	unsigned long dig_n;				// dig_n holds the next no_dig digit to add
	unsigned long long carry;
	unsigned long long tmp_n, tmp_dn;
	ss.reserve(digits + 16);
	int factor;

	if (x < y || x < 1)
		{ throw float_precision::domain_error(); }

	if (no_dig > 8) no_dig = 8;			// Ensure no_dig<=8
	if (no_dig < 1) no_dig = 1;			// Ensure no_dig>0
	// Since we do it in trunks of no_dig digits at a time we need to ensure digits is divisble with no_dig.
	dig = (digits / no_dig + (digits%no_dig>0 ? 1 : 0)) * no_dig;
	dig += no_dig;						// Extra guard digits
										// Calculate the number of terms needed
	factor = (int)ceil(10 * log(0.5) / log((double)(x - y) / (double)x));
	no_terms = (unsigned int)(factor * dig / 3 + 3);
	// Allocate the needed accumulators
	unsigned long long *acc_n = new unsigned long long[no_terms + 1];
	unsigned long long *acc_dn = new unsigned long long[no_terms + 1];
	f = f_table[no_dig];				// Load the initial f
	carry = 0;							// Set carry to 0
	//Loop for each no_dig
	for (size_t i = 0; i <= dig && overflow_flag == false; i += first_time == true ? 1 : no_dig, first_time = false)
		{
		// Calculate new number of terms needed
		no_terms = (unsigned int)(factor * (dig-i) / 3 + 3);
		// Loop for each no_terms
		for (size_t j = no_terms; j>0 && overflow_flag == false; --j)
			{
			if (first_time == true)
				{// Calculate the initialize value
				tmp_dn = (j + 1) * x;
				tmp_n = (x - y);
				}
			else
				{
				tmp_n = acc_n[j];
				tmp_dn = acc_dn[j];
				}
			if (tmp_n > (ULLONG_MAX) / f)
				overflow_flag = true;
			tmp_n *= f;		// Scale it
			// Check for 64bit overflow. Not very likely 
			if (carry > 0 && tmp_dn > (ULLONG_MAX - tmp_n) / carry)
				overflow_flag = true;
			tmp_n += carry * tmp_dn;
			carry = (tmp_n / (x * tmp_dn));
			carry *= (x - y);
			acc_n[j] = tmp_n % (tmp_dn * x);
			acc_dn[j] = tmp_dn;
			}

		if (first_time == true)
			{
			tmp_n = (x - y) * f;
			if (carry > 0 && tmp_n > (ULLONG_MAX - carry * x))
				overflow_flag = true;

			acc_n[0] = (tmp_n + carry*x);
			acc_dn[0] = x;
			dig_n = (unsigned)(acc_n[0] / (f*acc_dn[0]));
			}
		else
			{
			if (acc_n[0] > (ULLONG_MAX - carry * acc_dn[0]) / f)
				overflow_flag = true;
			dig_n = (unsigned)((acc_n[0] * f + carry * acc_dn[0]) / (f*acc_dn[0]));
			}
		car = (unsigned)(dig_n / f);
		dig_n %= f;
		// Add the carry to the existing number for digits calculate so far.
		if (car > 0)
			{
			for (size_t j = ss.length(); car > 0 && j > 0; --j)
				{
				unsigned int dd;
				dd = (ss[j - 1] - '0') + car;
				car = dd / 10;
				ss[j - 1] = dd % 10 + '0';
				}
			}
		ss += uitostring10((uintmax_t)dig_n, first_time==true? 1:no_dig );  //need to use uitostring10 instead of snprintf for performance
		//(void)snprintf(buffer, sizeof(buffer), "%0*lu", first_time == true ? 1 : no_dig, dig_n);
		//ss += std::string(buffer);
		if (first_time == true)
			acc_n[0] %= f*acc_dn[0];
		else
			{
			acc_n[0] = acc_n[0] * f + carry *acc_dn[0];
			acc_n[0] %= f  * acc_dn[0];
			}
		carry = 0;
		}

	ss.insert(1, ".");// add a . come after the first digit to create 2.30...
	if (overflow_flag == false)
		ss.erase(digits + 1); // Remove the extra digits that we didnt requested.
	else
		ss = std::string("Overflow:") + ss;

	delete[] acc_n;
	delete[] acc_dn;
	return ss;
	}

// End Spigot LN(X/Y)

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		11/26/2021
//	@brief 		Lookup or generate "fixed" constant ln2, PI log10 etc
//	@return 	float_precision	-	return the new table lookup value
//	@param		"tt"	-	Which table type lookup is needed
//	@param		"precision"	-	Number of significant digits
//
//	@todo
//
// Description:
//   Dynamic tables for "fixed" constant like ln(2), ln(10), e, PI and 1/Sqrt(2), sqrt(2)
//   If a higher precision is requested we create it and return otherwise 
//   we just the "constant" at a higher precision which eventually will be
//   rounded to the destination variables precision 
//	Added 1/sqrt(2) and sqrt(2) as constants
//
float_precision _float_table( enum table_type tt, size_t precision )
   {
   static float_precision ln2( 0, 0, ROUND_NEAR );
   static float_precision ln10( 0, 0, ROUND_NEAR );
   static float_precision pi( 0, 0, ROUND_NEAR );
   static float_precision e( 0, 0, ROUND_NEAR);
   static float_precision invsqrt2(0, 0, ROUND_NEAR);
   static float_precision invsqrt3(0, 0, ROUND_NEAR);
   float_precision res(0, precision, ROUND_NEAR);

	switch( tt )
		{
		case _EXP1:
			if (e.precision() >= precision)
				res = e;
			else
				{// Using Spigot algorithm for exp(1) Calculation
				std::string ss;
				size_t prec = std::max((size_t)20U, precision + 2);
				e.precision(prec);
				ss = spigot_e( prec );				// The result as a string in BASE_10
				e = float_precision(ss, prec);		// Convert to float_precision
				e.precision(std::max((size_t)20U, precision));
				res = e;
				}
			break;
		case _LN2:
			if (ln2.precision() >= precision)
				res = ln2;
			else
				{ // Using Spigot algorithm for LN2 Calculation
				std::string ss;
				size_t prec = std::max((size_t)20U, precision+2);
				ln2.precision( prec );
				ss = spigot_lnxy_64(2, 1, prec, 4);	// The result as a string in BASE_10
				ln2 = float_precision(ss, prec);		// Convert to float_precision
				ln2.precision(std::max((size_t)20U, precision));
				res = ln2;								// Save the result
				}
			break;
		case _LN10:
			if( ln10.precision() >= precision )
				res = ln10;
			else
				{ // Using Spigot Algorithm for LN10. LN(10)=3*ln(2)+ln(10/8)
				std::string ss;
				size_t prec=std::max((size_t)20U, precision + 2);
				ln10.precision(prec);
				ss = spigot_lnxy_64(2, 1, prec, 4);		// The result as a string in BASE_10
				ln10 = float_precision(ss, prec);		// Convert to float_precision
				ln10 *= float_precision(3);
				ss = spigot_lnxy_64(10, 8, prec, 4);	// The result as a string in BASE_10
				ln10 += float_precision(ss, prec);		// Convert and add to float_precision ln10
				ln10.precision(std::max((size_t)20U, precision));
				res = ln10;								// Save the result
				}
			break;
		case _PI: 
			if( pi.precision() > precision )
			res = pi;
			else
				{ // Using Brent-Salamin method
				unsigned int loopcnt = 0;
				const size_t min_precision = precision + 5 + (int)(log10(precision)+0.5);
				const eptype limit = -(int)ceil((precision+2)*log2(10));  
				const float_precision c05(0.5);
				float_precision a(1, min_precision), b(2, min_precision), sum(0.5, min_precision);
				float_precision ak(0, min_precision), bk(0, min_precision), ck(1, min_precision);
				float_precision ab(0, min_precision), asq(0, min_precision );
				float_precision pow2(1, precision);
#ifdef HVE_DEBUG
				std::cout << "PI Min precision=" << min_precision << " limit=2^" << limit << std::endl;  // Debug
				float_precision dt(0, 7);  // Debug
				int tadd = clock();
#endif
				pi.precision(min_precision);
				b = _float_table(_INVSQRT2, min_precision ); //when ready  // c1 / sqrt(b);  
				for (; !ck.iszero() && ck.exponent() > limit; ++loopcnt)
				{
#ifdef HVE_DEBUG
					dt = ck; // DEBUG
					tadd = clock() - tadd;
					std::cout << "\tPI Iteration=" << loopcnt << " Error=" << dt << " Time=" << tadd / CLOCKS_PER_SEC << std::endl;  // Debug
					tadd = clock();
#endif

#ifdef HVE_THREAD
					if(min_precision > 100'000)
						{
						std::thread first([=, &ak, &asq]()
						{ak = c05*(a + b);
						asq = ak * ak; });
						std::thread second([=, &ab, &bk]()
						{ab = a * b;
						bk = sqrt(ab); });
						first.join();
						second.join();
						}
					else
#endif
						{
						ak = c05*(a + b);
						ab = a * b;
						bk = sqrt(ab);
						asq = ak * ak;
						}

					ck = asq - ab;
					pow2.exponent(pow2.exponent() + 1); //pow2 *= float_precision(2);
					sum -= pow2*ck;
					a = ak; b = bk;
					}
#ifdef HVE_DEBUG
				dt = ck; // DEBUG
				tadd = clock() - tadd;
				std::cout << "\tPI Final Iteration=" << loopcnt << " Error=" << dt << " Time=" << tadd / CLOCKS_PER_SEC << std::endl;  // Debug
#endif
				pi = asq / sum; // float_precision(2) * asq / sum;
				pi.exponent(pi.exponent() + 1);    // Faster way to multiply by 2
				res = pi;
				// Round and store it
				pi.precision(precision);
#ifdef HVE_DEBUG
				std::cout << "PI finish with prec=" << precision << std::endl;  // Debug
#endif
				}	
			break;
		case _INVSQRT2:
			if (invsqrt2.precision() > precision)
			res = invsqrt2;
			else
				{
				unsigned int loopcnt = 0;
				const unsigned int extra = 2+5;
				const float_precision c1(1), c3(3), c05(0.5);
				size_t digits;
				float_precision r;
	 
				digits = invsqrt2.precision();  // Get current precision
				if (invsqrt2.iszero() == true)  // First time, do initialization 
					{
					digits = 16;
					invsqrt2.precision(std::max(precision, digits));  // Ensure minimum as 16 decimal digits
					// New. Get a initial guess using ordinary floating point
					// set the initial guess with at approx 16 correct digits
					invsqrt2 = float_precision(1.0 / sqrt(2.0),digits); // Ensure same precision as standard IEEE754
					}
				precision = std::max(precision, digits); // Keep maxium precision already calculated
				invsqrt2.precision(precision);	
#ifdef HVE_DEBUG
				std::cout << "INVSQRT2 Max precision=" << precision+extra << " start Precision=" << digits*2 << std::endl;  // Debug
				int tadd = clock();
#endif
				// Now iterate using Netwon Un=0.5U(3-2U^2), where U=invsqrt2
				for (digits *= 2; ; digits = std::min(precision + extra, digits * 2), ++loopcnt)
					{
					// Increase precision by a factor of two for the working variable s r & u. 
					r.precision(digits);
					invsqrt2.precision(digits);
					// Notice V is the original number to squareroot which has the full precision 
					// so we start by assigning it to r, rounding it to the precision of r
					r = invsqrt2 * invsqrt2;	// U^2
					r.exponent(r.exponent() + 1);	// 2U^2
					r = c3 - r;					// 3-2U^2
					r *= c05;					// (3-2U^2)/2
					invsqrt2 *= r;				// U=U(3-2U^2)/2
					//std::cout << loopcnt << " Prec=" << digits << " Diff=" << ((r - c1).exponent()) / log2(10) << std::endl;  // DeBUG
#ifdef HVE_DEBUG
					tadd = clock() - tadd;
					std::cout << "\tINVSQRT2 Iteration=" << loopcnt << " Precision="<< digits << " Error exponent=" << ((r - c1).exponent()) / log2(10) << " Time=" << tadd / CLOCKS_PER_SEC << std::endl;  // Debug
					tadd = clock();
#endif
					if (digits == precision + extra) // Reach final iteration step in regards to precision
						{
  						r.precision(precision + 1);	// round to final precision
						if (r == c1)	// break if no improvement
							break;
						r.precision(precision + extra);
						//std::cout << loopcnt << " Prec=" << precision << " Diff=" << ((r - c1).exponent())/log2(10) << std::endl; // DEBUG
						}
					}

				// Round to same precision as argument
				invsqrt2.precision(precision);
				res = invsqrt2;
#ifdef HVE_DEBUG
				std::cout << "INVSQRT2 finish with precision=" << precision << std::endl;  // Debug
#endif
				}
			break;
		case _SQRT2:  // use 2/sqrt(2)=sqrt(2)
			res = _float_table(_INVSQRT2, precision); //float_precision(2) * _float_table(_INVSQRT2, precision);
			res.exponent(res.exponent() + 1);	// Faster way to multiply by 2
			break;
		case _INVSQRT3:
			if (invsqrt3.precision() > precision)
				res = invsqrt3;
			else
				{
				unsigned int loopcnt = 0;
				const unsigned int extra = 2 + 5;
				const float_precision c1(1), c3(3), c05(0.5);
				size_t digits;
				float_precision r;

				digits = invsqrt3.precision();  // Get current precision
				if (invsqrt3.iszero() == true)  // First time, do initialization 
					{
					digits = 16;
					invsqrt3.precision(std::max(precision, digits));  // Ensure minimum as 16 decimal digits
																	  // New. Get a initial guess using ordinary floating point
																	  // set the initial guess with at approx 16 correct digits
					invsqrt3 = float_precision(1.0 / sqrt(3.0), digits); // Ensure same precision as standard IEEE754
					}
				precision = std::max(precision, digits); // Keep maxium precision already calculated
				invsqrt3.precision(precision);
				// Now iterate using Netwon Un=0.5U(3-2U^2), where U=invsqrt2
				for (digits *= 2; ; digits = std::min(precision + extra, digits * 2), ++loopcnt)
					{
					// Increase precision by a factor of two for the working variable s r & u. 
					r.precision(digits);
					invsqrt3.precision(digits);
					// Notice V is the original number to squareroot which has the full precision 
					// so we start by assigning it to r, rounding it to the precision of r
					r = c3;						// 3
					r *= invsqrt3 * invsqrt3;	// 3U^2
					r = c3 - r;					// 3-3U^2
					r *= c05;					// (3-3U^2)/2
					invsqrt3 *= r;				// U=U(3-3U^2)/2
					if (digits == precision + extra) // Reach final iteration step in regards to precision
						{
						r.precision(precision + 1);	// round to final precision
						if (r == c1)	// break if no improvement
							break;
						r.precision(precision + extra);
						}
					}

				// Round to same precision as argument
				invsqrt3.precision(precision);
				res = invsqrt3;
				}
			break;
		case _SQRT3:  // use 3/sqrt(3)=sqrt(3)
			res = float_precision(3) * _float_table(_INVSQRT3, precision);
			break;
		}

   return res;
   }

///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//    Exp(), Log(), Log10(), Nroot()
//
///////////////////////////////////////




// Experimental new exp() using sinh() and sqrt()
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Aug/2013
//	@brief 		Calculate exp(x)
//	@return 	float_precision -	Return exp(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a the identity that exp(x)=sinh(x)+sqrt(1+sinh(x)^2)
//	  This has proven to be faster than the standard taylor series for exp()
//   exp(x) == 1 + x + x^2/2!+x^3/3!+....
//	  A test i smade for x is an integer in which case we do pow(e,x) which is more than 400 times faster
//    sine e is calculated using the spigot algorithm using pure 64bit integer arithmetic
//
float_precision exp( const float_precision& x )
   {
   size_t precision;
   float_precision v;
   const float_precision c1(1);

   precision = x.precision()+2;  
   v.precision( precision );
   v = x;
   if( v.sign() < 0 )
      v.change_sign();

   if(floor(v)==v)  // v is an Integer
	  {// use the 100 times faster shortcut exp(v)=exp(1)^v
	  v = _float_table(_EXP1, precision );
	  v = pow( v, abs(x) );
	  }
   else
	  {
      v=sinh(v);
      v.precision( 2 * precision );  // Double the precision to avoid loss of significant when performaing 1+v*v
      v=v+sqrt(c1+v*v);
	  v.precision( precision );
      }

   if( x.sign() < 0 )
      v = _float_precision_inverse( v );
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );v.mode( x.mode() );
   v.precision( x.precision() );  

   return v;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/Oct/2021
//	@brief 		Calculate log(x)
//	@return 	float_precision -	Return log(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	 Simplified and improve with the new binary storage
//   Use a taylor series until their is no more change in the result
//   Equivalent with the same standard C function call
//   ln(x) == 2( z + z^3/3 + z^5/5 ...
//   z = (x-1)/(x+1)
//
float_precision log( const float_precision& x )
	{
	size_t precision;
	eptype expo;
	int j, k, no_reduction;
	double zd;
	float_precision res, r, z, z2;
	const float_precision c1(1);

	if( x <= float_precision(0) ) 
		{ throw float_precision::domain_error(); }

	precision = x.precision() + 2;  
	z.precision( precision ); // Do calc at 2 higher precision to allow correct rounding of result
	z = x;
	expo = z.exponent();		// Get original exponent
	z.exponent( 0 );			// Set exponent to zero getting a z number beween [1..2)

	// Check for augument reduction and increase precision if necessary
	zd = (double)z;
	no_reduction = (int)ceil(log( log(zd) / log(1.1)) / log(2));
	no_reduction = std::max(no_reduction, 0);
	precision += no_reduction;

	z.precision( precision ); // adjust precision to allow correct rounding of result
	r.precision( precision ); 
	z2.precision( precision );
	res.precision( precision );
   
	// In order to get a fast Taylor series result we need to get the fraction closer to one
	// The fraction part is [1...1.1) (base 10) at this point
	// Repeat a series of no_reduction square root 
	for( k = 0; k < no_reduction; ++k )
		z = sqrt(z);

	// number now at [1...1.1). Setup the iteration
	z = ( z - c1 ) / ( z + c1 );
	z2 = z * z;
	res = z;
	// Iterate using taylor series ln(x) == 2( z + z^3/3 + z^5/5 ... )
	for( j=3;;j+=2 )
		{
		z *= z2;
		r = z/float_precision(j);
		if( res + r == res )
			break;
		res += r;
		}

	// Adjust result from the reduction by multiply it with 2^(k+1)
	res *= float_precision( pow( 2.0, (double)( k + 1 ) ) );
	if( expo != 0 )  // Adjust for original exponent
		{// Ln(x^y) = Ln(x) + Ln(2^y) = Ln(x) + y * ln(2) 
		res += float_precision(expo) * _float_table(_LN2, precision + 1);
		}

	// Round to same precision as argument and rounding mode
	res.mode( x.mode() );
	res.precision( x.precision() );  

	return res;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate log10(x)
//	@return 	float_precision -	Return log10(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Log10. Use the equation log10(x)=log(x)/log(10)
//   Equivalent with the same standard C function call
//
float_precision log10( const float_precision& x )
	{
	size_t precision = x.precision();  
	float_precision res( 0, precision + 1 );

	if( x <= float_precision(0) ) 
		{ throw float_precision::domain_error(); }

	res = x;
	res = log( res ) / _float_table( _LN10, precision + 1 );
   
	// Round to same precision as argument and rounding mode
	res.mode( x.mode() );
	res.precision( x.precision() );  

	return res;
	}

///////////////////////////////////////
//
// FLOAT PRECISION FUNCTIONS
//    Special functions: pow(), fmod(), floor(), ceil(), modf(), fabs(), ldexp(), frexp()
//
///////////////////////////////////////



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		4/Oct/2021
//	@brief 		Calculate pow(x,y)
//	@return 	float_precision -	Return pow(x,y)
//	@param      "x"	- The argument
// @param      "y" - The power argument
//
//	@todo  
//
// Description:
//   x^y == exp( y * ln( x ) ) ); in general, however if y is an integer then we use the ipow() algorithm instead.
//   Update to use the new method .toInteger()
// 
float_precision pow( const float_precision& x, const float_precision& y )
   {
   float_precision i, res(1);
   eptype expo;
   bool yinteger=false;

   // add two extra guard digits to avoid loss of precision when performing  exp( y * ln(x) ) )
   res.precision( x.precision()+2 );  

   expo = y.exponent();
   if( expo >= 0 )
      {
      i.precision( y.precision() );
      i = y;
	  i.toInteger(); // now i is the integer part of y. 
	  // Check that y is a true integer, with a max range of a 32 bit integer
	  if( y == i && i <= float_precision( INT_MAX ) )
		  yinteger = true;
      }
   
   if( yinteger == false ) // y is not an integer so do x^y= exp^(y*log(x)) the regular way
      {
	  res = x;
      res = log( res ) * y;
      res= exp( res );
      }
   else
		{ // raise to the power of y when y is an integer. Use optimized method.
		int sign = i.sign();
		if( sign < 0 )
			i.change_sign();

		int ie = (int)i;  // raise to the power of ie which is at max a standard 32bit integer
		float_precision p( x );

		for( int n = ie; n > 0; n >>= 1 ) 
			{
			if( ( n & 0x1 ) != 0 ) 
				res *= p;  // Odd
			p *= p;						 
			}
		if( sign < 0 )
			res = _float_precision_inverse( res );
		 }
   
   res.precision(x.precision());
   return res;
   }
 
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		27/Sep/2021
//	@brief 		Calculate fmod(x,y)
//	@return 	float_precision -	Return fmod(x,y)
//	@param      "x"	- The argument
// @param      "y"   - The argument
//
//	@todo  
//
// Description:
//   float precision. fmod remainder of x/y
//   Equivalent with the standard C function fmod
//   x = i * y + f or f = x - i * y; and i = integer(x/y)
//	  Revised to use the method.toInteger() for faster calculation
//
float_precision fmod( const float_precision& x, const float_precision& y )
   {
   float_precision i, f;
   
   f.precision( x.precision() );
   i.precision( x.precision() );
   i = x / y;
   if( i.exponent() < 0 )
		f = x;
   else
		{
		i.toInteger();
		f = x - i * y;
		}	

   return f;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Oct/2021
//	@brief 		Calculate floor(x)
//	@return 	float_precision -	Return floor(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision floor
//   Equivalent with the same standard C function floor()
//	  Rounds x downward, returning the largest integral value that is not greater than x.
//
float_precision floor( const float_precision& x )
	{
	float_precision f(0, x.precision() );
	const float_precision c1(1);

	if( x.exponent() < 0 ) // is number less than |1|
		{
		if( x.sign() < 0 )
			f = -c1;
		}
	else
		{
		f = x;
		f.toInteger();
		if (f.sign() < 0 && (x - f).iszero() == false )
			f -= c1;
		}
 
	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate ceil(x)
//	@return 	float_precision -	Return ceil(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision ceil
//   Equivalent with the same standard C function ceil()
//   Rounds x upward, returning the smallest integral value that is not less than x.
//
float_precision ceil( const float_precision& x )
	{
	float_precision f(0, x.precision() );
	const float_precision c1(1);

	if( x.exponent() < 0 ) // is number less than |1|
		{
		if( x.sign() > 0 )
			f = c1;
		}
	else
		{
		f = x;
		f.toInteger();
		if (f.sign() > 0 && (x - f).iszero() == false)
			f += c1;
		}

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate trunc(x)
//	@return 	float_precision -	Return trunc(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision trunc
//   Equivalent with the same standard C function trunc()
//   Rounds x towards zero.
//
float_precision trunc(const float_precision& x)
	{
	float_precision f(0, x.precision());

	if (x.exponent() < 0) // is number less than |1|
		{
		if (x.sign() > 0)
			f = float_precision(0);
		}
	else
		{
		f = x;
		f.toInteger();
		}

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		Calculate round(x)
//	@return 	float_precision -	Return round(x)
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision round
//   Equivalent with the same standard C function round()
//
float_precision round(const float_precision& x)
	{
	float_precision f(x);
	const float_precision c05(0.5);

	f += (f.sign() < 0 ? -c05 : c05 );
	f = trunc(f);

	return f;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		30/Sep/2021
//	@brief 		Split a float number into integer part and fraction
//	@return 	float_precision -	Fraction part of x
//	@param      "x"	- The argument
// @param      "intptr" - Float_precision pointer to integer part of x
//
//	@todo  
//
// Description:
//   Float Precision fmod
//   Split a Floating point number into integer part and fraction
//   Equivalent with the same standard C function call
//	  Use a modified version of the function fmod(x,c1)
//	  Converted to binary version
//   Notice that if x < 0 then BOTH the return value and intptr becomes negative
//	  even when inptr==0
//
float_precision modf( const float_precision& x, float_precision *intptr )
	{
	float_precision i, f;
	eptype expo;

	f.precision(x.precision());
	i.precision(x.precision());
	intptr->precision(x.precision());
	f = x;
	expo = f.exponent();
	if (expo < 0)
		{
		i = float_precision(0);
		i.sign(x.sign());
		}
	else
		{
		i=f.toFraction();
		//i.toInteger();
		//f = x - i;
		}
	*intptr = i;
	return f;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2012
//	@brief 		Calculate abs(x)
//	@return 	float_precision -	Return absolute value of x
//	@param      "x"	- The argument
//
//	@todo  
//
// Description:
//   Float Precision abs()
//   Equivalent with the same standard C function call fabs()
//
float_precision abs( const float_precision& x )
   {
   float_precision f(0, x.precision() );

   f = x;
   if( f.sign() < 0 )
      f.change_sign();

   return f;
   }



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate ldexp((x)
//	@return 	float_precision -	Return ldexp(x)
//	@param      "x"	- The argument
// @param      "exp" - exponent argument
//
//	@todo  
//
// Description:
//   The ldexp function returns the value of x * 2^exp
//
float_precision ldexp( const float_precision& x, eptype exp )
   {
   if( exp == 0 )
      return x;
   if( exp > 0 && exp <= 63 )
      return x * float_precision( 1U << exp );

   return x * pow( float_precision( 2 ), float_precision( exp ) );
   }



//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		Calculate frexp(x,expptr)
//	@return 	float_precision -	Return mantissa part of number x
//	@param      "x"	- The argument
// @param      "expptr" - Pointer to the exponent part of number
//
//	@todo  
//
// Description:
//   The frexp()
//   The frexp function breaks down the floating-point value (x) into a mantissa (m) and an exponent (n), 
//   such that the absolute value of m is greater than or equal to 1/RADIX and less than RADIX, and x = m*Radix^n. 
//   The integer exponent n is stored at the location pointed to by expptr. 
//
float_precision frexp( float_precision& x, eptype *expptr )
   {
   if( x.iszero() )
      *expptr = 0;

   *expptr = x.exponent();
   x.exponent( 0 );

   return x;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		12/Nov/2015
//	@brief 		Calculate n root of x  (x^(1/n)
//	@return 	float_precision -	Return nroot(a)
//	@param      "x"	-	The nroot argument
//
//	@todo  
//
// Description:
//   nroot(V)
//   The nth root of x^(1/n) No Equivalent standard C function call
//   Seperate exponent. e.g. nroot(V*10^x)=10^x/2*nroot(V)
//   Un=U(1/n)(n+1-VU^n)
//   Then Un == 1/nroot(V). and nroot(V) = 1/Un
// The function has been improved using Newton with iterative deepening creating
//	a speed up with a factor of 3 over the classic Newton method.
// This is a much much faster option instead of the traditional pow() function x^y
// and that is why it has been added as a separate function
//
float_precision nroot(const float_precision& x, unsigned int n)
	{
	const size_t extra = 2;
	size_t precision;
	size_t digits;
	eptype expo, expo_sq;
	double fv;
	float_precision r, u, v, tmp, fn(n);
	const float_precision c1(1);

	if (x.iszero() || x == c1 || n == 1)
		return x;
	if (x.sign() < 0)
		{ throw float_precision::domain_error(); }
	precision = x.precision();
	v.precision(precision + extra);
	v = x;
	expo = v.exponent();
	expo_sq = expo / 2;
	v.exponent(expo - 2 * expo_sq);
	r.precision(precision + extra); // Do iteration using 2 digits higher precision
	u.precision(precision + extra);
	fn.precision(precision + extra);

	// Get a initial guess using ordinary floating point
	fv = v;		
	if (expo - 2 * expo_sq > 0)
		fv *= (double)F_RADIX;
	else
		if (expo - 2 * expo_sq < 0)
			fv /= (double)F_RADIX;
	fv = pow(fv, 1.0 / n);  // set the initial guess with at approx 16 correct digits
	fv = 1 / fv;

	tmp.precision(precision+1);
	u = float_precision(fv);
	fn = 1 / fn;
	// Now iterate using Netwon  Un=U*(-VU^n+(n+1))/n
	for (digits = std::min((size_t)32, precision); ; digits = std::min(precision + extra, digits * 2))
		{
		// Increase precision by a factor of two
		r.precision(digits);
		u.precision(digits);
		float_precision p(u);
		float_precision res(1, digits);
		// DO U^N
		for (int i = n; i > 0; i >>= 1)
			{
			if ((i & 0x1) != 0)
				res *= p;  // Odd
			if (i>1)
				p *= p;
			}
		// Notice V is the original number to nroot which has the full precision 
		// so we start by assigning it to r, rounding it to the precision of r
		r = v;							// V
		r *= res;						// VU^n
		r = float_precision(n + 1) - r; // (n+1)-VU^n
		r *= fn;						// (-VU^n+(n+1))/n
		u *= r;							// U=U*(-VU^n+(n+1))/n
		if (digits == precision + extra) // Reach final iteration step in regards to precision
			{
			tmp = r;			// round to final precision
			if (tmp == c1)		// break if no improvement
				break;
			}
		}

	u = 1 / u;			// n root of u is now 1/u;
	u.exponent(u.exponent() + expo_sq);
	// Round to same precision as argument and rounding mode
	u.mode(x.mode());
	u.precision(precision);

	return u;
	}


///////////////////////////////////////
//
// TRIGONOMETRIC FUNCTIONS
//   atan()
//   atan2()
//   asin()
//   acos()
//   sin()
//   cos()
//   tan()
//
///////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		atan
//	@return		float_precision	-	return atan(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the taylot series. ArcTan(x) = x - x^3/3 + x^5/5 ...
//   However first reduce x to abs(x)< 0.5 to improve taylor series
//   using the identity. ArcTan(x)=2*ArcTan(x/(1+sqrt(1+x^2)))
//   We actually dynamically adjust the argument reduction factor by applying
//   more with higher precision numbers.
//
float_precision atan( const float_precision& x )
   {
   size_t precision;
   int j, k;
   double zd, dlimit;
   float_precision r, u, v, v2;
   const float_precision c1(1);

   precision = x.precision()+2;  
   v.precision( precision );
   v = x;

   // Check for augument reduction and increase precision if necessary
   zd = log(precision ) / log(BASE_10);
   zd *= zd;
   j=(int)zd; if(j>32) j=32;
   
   // Lets just do one reduction because that quarantee us that it is less than 1
   // and we can then use standard IEEE754 to calculate the needed argument reduction.
   zd = (double)abs( v / ( c1 + sqrt( c1 + v * v ) ) );
   // Calculate the number of reduction needed 
   dlimit=0.5/pow(2.0,j); // Only estimated target reduction limit based on x/(1+sqrt(1+x*x)->x/2 for small x
   for( j=1; zd > dlimit; j++ )
       zd=zd/(1.0+sqrt(1.0+zd*zd));

   // Adjust the precision
   if (j > 0)
	   precision += j / 4; // PADJUST(j / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );

   // Transform the solution to ArcTan(x)=2*ArcTan(x/(1+sqrt(1+x^2)))
   for( k=1; j>0; k *= 2, j-- )
        v = v / ( c1 + sqrt( c1 + v * v ) );

   v2 = v * v;
   r = v;
   u = v;
   // Now iterate using taylor expansion
   for( j=3;; j+=2 )
      {
      v *= v2;
      v.change_sign();
      r = v / float_precision(j);;
      if( u + r == u )
         break;
      u += r;
      }

   u *= float_precision( k );

   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   return u;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		atan2
//	@return 	float_precision	-	return the angle (in radians) from the X axis to a point (y,x).
// @param		"y"   -  float_precision y-axis
//	@param      "x"	-	float_precision x-axis
//
//	@todo 
//
// Description:
//   use atan() to calculate atan2()
//
float_precision atan2( const float_precision& y, const float_precision& x )
   {
   size_t precision;
   float_precision u;
   const float_precision c0(0), c05(0.5);

   if( x.iszero() && y.iszero() )
      return c0;

   precision = x.precision()+2;  
   u.precision( precision );
   if( x.iszero() )
      {
      u = _float_table( _PI, precision );
      if( y < c0 )
         u *= -c05;
      else
         u *= c05;
      }
   else
      if( y.iszero() )
         {
         if( x < c0 )
            u = _float_table( _PI, precision );
         else
            u = c0;
         }
      else
         {
         u = atan( y / x );
         if( x < c0  && y < c0 )
            u -= _float_table( _PI, precision );

         if( x < c0 && y >= c0 )
            u += _float_table( _PI, precision );
         }

	// Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   return u;
   }



   
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		8/26/2013
//	@brief 		Calculate asin(x)
//	@return 	float_precision -	Return asin(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a taylor series until their is no more change in the result
//   asin(x) == x + x^3/(2*3)+(1*3)x^5/(2*4*5)+(1*3*5)x^7/(2*4*6*7)....
//   Use argument reduction via the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
//	  This function replace the other function using newton iteration. Taylor series is significant
//    faster e.g 50% for 10 digits, 3x for 100 digits and 5x for 1000 digits.
//
float_precision asin( const float_precision& x )
   {
   size_t precision;
   int k, j, sign;
   double zd, dlimit;
   float_precision r, u, v, v2, sqrt2, lc, uc;
   const float_precision c1(1), c2(2);

    if( x > c1 || x < -c1 )
      { throw float_precision::domain_error(); }

   precision = x.precision() + 2;  
   v.precision( precision );
   v = x;

   sign = v.sign();
   if( sign < 0 )
      v.change_sign();

   // Check for augument reduction and increase precision if necessary
   zd = log(precision) / log(BASE_10); 
   zd *= zd;
   j=(int)zd; if(j>16) j=16;
   zd=v;
   // Find the argument reduction factor
   for( dlimit=1.0; j > 0; j-- )
       {
       dlimit/=sqrt(2.0)* sqrt( 1.0 + sqrt( 1.0 - dlimit * dlimit ) );
       if( dlimit < zd ) break;
       }
   // j is the number of argument reduction
    // Adjust the precision
   if (j > 0)
	   precision += j / 4; // PADJUST(j / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );
   sqrt2.precision( precision );
   lc.precision( precision );
   uc.precision( precision );
  
   // Now use the identity arcsin(x)=2arcsin(x)/(sqrt(2)+sqrt(1-x*x))
  // until argument is less than dlimit
  // Reduce the argument to below 0.5 to make the newton run faster
   sqrt2=c2;				// Ensure correct number of digits
   sqrt2=sqrt( sqrt2 );	// Now calculate sqrt2 with precision digits
   for( k=0; j > 0; k++, j-- )
      v /= sqrt2 * sqrt( c1 + sqrt( c1 - v * v ) );
  
   v2 = v * v;
   r = v;
   u = v;
   // Now iterate using taylor expansion
   for( unsigned int m=3;; m+=2 )
      {
      if( j < 65536 ) 
          {
            uc = float_precision ( ( m - 2 ) * ( m - 2 ) );
            lc = float_precision( m * m - m ); 
        }
      else 
          {
          uc = float_precision( m - 2 ) * float_precision( m - 2 );
          lc = float_precision( m - 1 ) * float_precision( m );
        }
      v = uc * v2 / lc;
      r *= v;
      if( u + r == u )
         break;
      u += r;
      }

   if( k > 0 )
       u *= float_precision( 1 << k );

   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   if( sign < 0 )
      u.change_sign();

   return u;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		acos
//	@return 	float_precision	-	return acos(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//    Use Arccos(x)=PI/2 - Arcsin(x) or ArcCos(x)=0.5*PI-ArcSin(x)
//
float_precision acos( const float_precision& x )
   {
   size_t precision;
   float_precision y;
   const float_precision c1(1);
   
      if( x > c1 || x < -c1 )
      { throw float_precision::domain_error(); }
      
   precision = x.precision();  
   y = _float_table( _PI, precision );
   y *= float_precision( 0.5 );
   y -= asin( x );

   // Round to same precision as argument and rounding mode
   y.mode( x.mode() );
   y.precision( precision );  

   return y;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		sin
//	@return 	float_precision	-	return sin(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the taylor series. Sin(x) = x - x^3/3! + x^5/5! ...
//   1) However first reduce x to between 0..2*PI 
//   2) Then reduced further to between 0..PI using sin(x+PI)=-Sin(x)
//   3) Finally reduced it to below 0.5/3^reduction factor, using the trisection identity
//         sin(3x)=3*sin(x)-4*sin(x)^3
//   4) Then do the taylor. 
//   The argument reduction is used to reduced the number of taylor iteration 
//   and to minimize round off erros and calculation time
//
float_precision sin( const float_precision& x )
   {
   size_t precision;
   int k, sign, j;
   double zd;
   float_precision r, u, v, v2, de(0);
   const float_precision c1(1), c2(2), c3(3), c4(4);

   precision = x.precision() + 2;  
   // Check for augument reduction and increase precision if necessary
   zd = 2 * ( log(precision ) / log(BASE_10) ); 
   j=(int)zd; if(j>1 && j<5) j--; if(j>8) j=8;
   // Adjust the precision
   if (j > 0)
	   precision += j / 4; // PADJUST(j / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );

   v = x;
   sign = v.sign();
   if( sign < 0 )
      v.change_sign();
   
   // Check that argument is larger than 2*PI and reduce it if needed. 
   // No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
   if( v > float_precision( 2*3.14159265 ) )
      {
      // Reduce argument to between 0..2PI
      u = _float_table( _PI, precision );
      u *= c2;
      if( abs( v ) > u )
         {
         r = v / u; 
         (void)modf( r, &r ); 
         v -= r * u;
         }
      if( v < float_precision( 0 ) )
         v += u;
	  }   
   
   // Reduced it further to between 0..PI
   // However avoid calculating PI is not needed.
   // No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
   if( v > float_precision( 3.14159265 ) )
      {
	  u = _float_table( _PI, precision ); // We dont need to worry that we called it a second time since it will be cached
	  if( v > u )
		  { v -= u; sign *= -1; }
      }

   // Now use the trisection identity sin(3x)=3*sin(x)-4*sin(x)^3
   // until argument is less than 0.5/3^j  Where J is the number of reduction factor based on the needed precision of the argument.
   v2= v * float_precision( 2 * pow( 3.0, j ) );
   for( k = 0, r = c1; v2 > r; k++ )
      r *= c3;
   v /= r;

   v2 = v * v;
   r = v;
   u = v;

   // Now iterate using taylor expansion
   for( unsigned int m=3;; m+=2 )
      {
      de += float_precision( 4 * m - 6 ); // Avoid the multiplication in float_precision. 
      v = v2 / de;  // de is only 20 digits standard preecision but since v2 is precision it will be extended to v2 precision //( j<USHRT_MAX? float_precision( m * (m-1) ) : float_precision(m) * float_precision(m-1) );
      r *= v;
      r.change_sign();
      if( u + r == u )
         break;
      u += r;
      }

   for( ; k > 0 ; k-- )
      u *= ( c3 - c4 * u * u );

   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   if( sign < 0 )
      u.change_sign();

   return u;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		cos
//	@return 	float_precision	-	return cos(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the taylor series. Cos(x) = 1 - x^2/2! + x^4/4! - x^6/6! ...
//   1) However first reduce x to between 0..2*PI
//   2) Then reduced it further to between 0..PI using cos(x)=Cos(2PI-x) for x >= PI
//   3) Now use the trisection identity cos(3x)=-3*cos(x)+4*cos(x)^3
//      until argument is less than 0.5/3^argument reduction
//   4) Finally use Taylor 
//
float_precision cos( const float_precision& x )
   {
   size_t precision;
   int k, j;
   double zd;
   float_precision r, u, v, v2, de(0);
   const float_precision c05(0.5), c1(1), c2(2), c3(3), c4(4);

   precision = x.precision() + 2;  
   // Check for augument reduction and increase precision if necessary
   zd = 2 * ( log(precision ) / log(BASE_10) );
   j=(int)zd; if(j>1 && j<5) j--; if(j>8) j=8;
   // Adjust the precision
   if (j > 0)
	   precision += j / 4; // PADJUST(j / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );

   v = x;
 
   // Check that argument is larger than 2*PI and reduce it if needed. 
   // No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
   if( abs( v ) > float_precision( 2*3.14159265 ) )
      {  // Reduce argument to between 0..2P
	  u = _float_table( _PI, precision );
      u *= c2;
      if( abs( v ) > u )
         {
         r = v / u; 
         (void)modf( r, &r ); 
         v -= r * u;
         }
      if( v < float_precision( 0 ) )
         v += u;
      }

   // Reduced it further to between 0..PI. 
   // However avoid calculating PI is not needed.
   // No need for high perecision. we just need to figure out if we need to Calculate PI with a higher precision
   if( abs( v ) > float_precision( 3.14159265 ) )
      {
      r = _float_table( _PI, precision );
      if( v > r )
         v = r * c2 - v;
      }

   // Now use the trisection identity cos(3x)=-3*cos(x)+4*cos(x)^3
   // until argument is less than 0.5/3^j  Where J is the number of reduction factor based on the needed precision of the argument.
   v2= abs( v * float_precision( 2 * pow( 3.0, j ) ) );
   for( k = 0, r = c1; v2 > r; k++ )
      r *= c3;
   v /= r;

   v2 = v * v;
   r = c1;
   u = r;
   // Now iterate using taylor expansion
   for( unsigned int m=2;; m+=2 )
      {
      de += float_precision( 4 * m - 6 ); // Avoid the multiplication in float_precision. 
      v = v2 / de;  // de is only 20 digits standard preecision but since v2 is precision is will be extended to v2 precision //( m<USHRT_MAX? float_precision( m * (m-1) ) : float_precision(m) * float_precision(m-1) );
      r *= v;
      r.change_sign();
      if( u + r == u )
         break;
      u += r;
      }

   for( ; k > 0 ; k-- )
      u *= ( c4 * u * u - c3 );
 
   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   return u;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		tan
//	@return 	float_precision	-	return tan(x)
//	@param      "x"	-	float_precision argument
//
//	@todo 
//
// Description:
//   Use the identity tan(x)=Sin(x)/Sqrt(1-Sin(x)^2)
//   1) However first reduce x to between 0..2*PI
//   2) Use taylor
//
float_precision tan( const float_precision& x )
   {
   size_t precision;
   float_precision u, r, v, p;
   const float_precision c1(1), c2(2), c3(3), c05(0.5);

   precision = x.precision() + 2;  
   u.precision( precision );
   v.precision( precision );
   p.precision( precision );
   v = x;
  
   // Check that argument is larger than 2*PI and reduce it if needed. 
   p = _float_table( _PI, precision );
   u = c2 * p;
   if( abs( v ) > u )
      {
      r = v / u; 
      (void)modf( r, &r ); 
      v -= r * u;
      }
   if( v < float_precision( 0 ) )
      v += u;
    
   p *= c05;
   if( v == p || v ==  p * c3 )
      { throw float_precision::domain_error(); }

   u = sin( v ); 
   if( v < p || v > p * c3 ) 
      u /= sqrt( c1 - u * u );
   else
      u /= -sqrt( c1 - u * u );
   
   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   return u;
   }


///////////////////////////////////////
//
// END TRIGONOMETRIC FUNCTIONS
//
///////////////////////////////////////

///////////////////////////////////////
//
// Hyperbolic FUNCTIONS
//   sinh()
//   cosh()
//   tanh()
//	  arcsinh()
//	  arccosh()
//	  arctanh()
//
///////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013
//	@brief 		Calculate Sinh(x)
//	@return 	float_precision -	Return Sinh(x)
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a taylor series until their is no more change in the result
//   sinh(x) == x + x^3/3!+x^5/5!+....
//   Use argument reduction via sinh(3x)=sinh(x)(3+4sinh^2(x))	
//
float_precision sinh( const float_precision& x )
   {
   size_t precision;
   int k, j, sign;
   double zd, dlimit;
   float_precision r, u, v, v2, de(0);
   const float_precision c1(1), c3(3), c4(4);

   precision = x.precision() + 2;  
   v.precision( precision );
   v = x;
   r.precision( precision ); 
   u.precision( precision );
   v2.precision( precision );

   sign = v.sign();
   if( sign < 0 )
      v.change_sign();

   // Check for augument reduction and increase precision if necessary
   zd = log(precision) / log(BASE_10); // zd=PLOG10( precision );
   zd *= zd;
   j=(int)zd; j -= 1; if(j<0) j=0;  if(j>16) j=16;
   dlimit=pow( 3.0, j );
   // Now use the trisection identity sinh(3x)=sinh(x)(3+4Sinh^2(x))
   // until argument is less than 0.5 * (1/3)^j
   j = int( 2.0 * dlimit );
   v2= v * float_precision( j );
   for( k = 0, r = c1; v2 > r; k++ )
      r *= c3;
   // Adjust the precision
   if (k > 0)
	   precision += j / 4; // PADJUST(k / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );
  
   v /= r;
   v2 = v * v;
   r = v;
   u = v;
   // Now iterate using taylor expansion
   for( unsigned int m=3;; m+=2 )
      {
      de += float_precision( 4 * m - 6 ); // Avoid the multiplication in float_precision. 
      v = v2 / de;  // de is only 20 digits standard preecision but since v2 is precision is will be extended to v2 precision //( j<USHRT_MAX? float_precision( m * (m-1) ) : float_precision(m) * float_precision(m-1) );
      r *= v;
      if( u + r == u )
         break;
      u += r;
      }

   for( ; k > 0 ; k-- )
      u *= ( c3 + c4 * u * u );

   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   if( sign < 0 )
      u.change_sign();

   return u;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013
//	@brief 		Calculate Cosh(x)
//	@return 	float_precision -	Return Cosh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//   Use a taylor series until their is no more change in the result
//   cosh(x) == 1 + x^2/2!+x^4/4!+....
//   Use argument reduction via cosh(3x)=cosh(x)(4cosh^2(x)-3)	
//
float_precision cosh( const float_precision& x )
   {
   size_t precision;
   int k, j, sign;
   double zd, dlimit;
   float_precision r, u, v, v2, de(0);
   const float_precision c1(1), c3(3), c4(4);

   precision = x.precision() + 2;  
   v.precision( precision );
   v = x;
   r.precision( precision ); 
   u.precision( precision );
   v2.precision( precision );

   sign = v.sign();
   if( sign < 0 )
      v.change_sign();  // cosh(-x) = cosh(x)

    // Check for augument reduction and increase precision if necessary
   zd = log(precision) / log(BASE_10); // zd=PLOG10( precision );
   zd *= zd;
   j=(int)zd; j -= 1; if(j<0) j=0;  if(j>16) j=16;
   dlimit=pow( 3.0, j );
   // Now use the trisection identity cosh(3x)=cosh(x)(4cosh^2(xx)-3)
   // until argument is less than 0.5 * (1/3)^j
   j = (int)(  2.0 * dlimit );
   v2= v * float_precision( j );
   for( k = 0, r = c1; v2 > r; k++ )
      r *= c3;
   // Adjust the precision
   if (k > 0)
	   precision += j / 4; // PADJUST(k / 4);
   r.precision( precision );
   u.precision( precision );
   v.precision( precision );
   v2.precision( precision );
   
   v /= r;
   v2 = v * v;
   r = c1;
   u = r;
   // Now iterate using taylor expansion
   for( unsigned int m=2;; m+=2 )
      {
      de += float_precision( 4 * m - 6 ); // Avoid the multiplication of float_precision.  
      v = v2 / de;  // de is only 20 digits standard preecision but since v2 is precision is will be extended to v2 precision //( j<USHRT_MAX? float_precision(m * (m-1) ) : float_precision(m) * float_precision(m-1) );
      r *= v;
      if( u + r == u )
         break;
      u += r;
      }

   for( ; k > 0 ; k-- )
      u *= ( c4 * u * u -c3 );
      
   // Round to same precision as argument and rounding mode
   u.mode( x.mode() );
   u.precision( x.precision() );  

   return u;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/21/2013
//	@brief 		Calculate Tanh(x)
//	@return 	float_precision -	Return Tanh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	tanh = ( exp(x) - exp(-x) ) / ( exp( x) + exp(-x) )=(e^(2x)-1/(e^(2x)+1)
// 
//
float_precision tanh( const float_precision& x )
   {
   float_precision v, v2;
   const float_precision c1(1);

   v.precision( x.precision() + 1 );
   v2.precision( x.precision() + 1 );
   v = x;
   v = exp( v );
   v2= v * v;
   v = (v2-c1)/(v2+c1);

   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  

   return v;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcSinh(x)
//	@return 	float_precision -	Return ArcSinh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcSinh=Ln(x+Sqrt(x^2+1))
// 
float_precision asinh( const float_precision& x )
   {
   float_precision v;
   const float_precision c1(1);

   v.precision( x.precision() + 1 );
   v = x;
   v = log(v+sqrt(v*v+c1));
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  

   return v;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcCosh(x)
//	@return 	float_precision -	Return ArcCosh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcCosh=Ln(x+Sqrt(x^2-1))
// 
//
float_precision acosh( const float_precision& x )
   {
   float_precision v;
   const float_precision c1(1);

   if( x < c1 )
      { throw float_precision::domain_error(); }
   
   v.precision( x.precision() + 1 );
   v = x;
   v = log(v+sqrt(v*v-c1));
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  

   return v;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2013
//	@brief 		Calculate ArcTanh(x)
//	@return 	float_precision -	Return ArcTanh()
//	@param      "x"	-	   The argument
//
//	@todo  
//
// Description:
//	ArcTanh=0.5*Ln((1+x)/(1-x))
// 
//
float_precision atanh( const float_precision& x )
   {
   float_precision v;
   const float_precision c05(0.5), c1(1);

   if( x >= c1 || x <= -c1 )
      { throw float_precision::domain_error(); }

   v.precision( x.precision() + 1 );
   v = x;
   v = c05*log((c1+v)/(c1-v));
   
   // Round to same precision as argument and rounding mode
   v.mode( x.mode() );
   v.precision( x.precision() );  

   return v;
   }

///////////////////////////////////////
//
// END TRIGONOMETRIC FUNCTIONS
//
///////////////////////////////////////


///////////////////////////////////////
//
// FAST integer division and remaining using Floating point arithmetic
//
///////////////////////////////////////

#if false

int_precision _int_precision_fastdiv( const int_precision &s1, const int_precision &s2 )
	{
	size_t ss;
	int_precision r2;
	float_precision f1, f2, rf;

	ss = std::max(s1.size(), s2.size());
	ss = (int)(ceil(Bitsiptype*ss / log2(BASE_10)));
	f1.precision( ss+2);
	f2.precision( ss+2);
	rf.precision( ss+2 );
	f1=float_precision(s1, ss+ 2);
	f2=float_precision( s2, ss+ 2);
	rf=f1/f2;
	r2 = (int_precision)rf; 
	return r2;
	}

int_precision _int_precision_fastrem( const int_precision &s1, const int_precision &s2 )
	{
	size_t ss;
	int_precision r2;
	float_precision f1, f2, rf;

	ss = std::max(s1.size(), s2.size());
	ss = (int)(ceil(Bitsiptype*ss / log2(BASE_10)));
	f1.precision( ss+2);
	f2.precision( ss+2);
	rf.precision( ss+2 );
	f1=float_precision(s1, ss+ 2);
	f2=float_precision( s2, ss+ 2);
	rf=f1/f2;
	r2 = (int_precision)rf; 
	r2=s1-s2*r2;
	return r2;
	}
#endif

///////////////////////////////////////
//
// FLOATING POINT FUNCTIONS
//
///////////////////////////////////////