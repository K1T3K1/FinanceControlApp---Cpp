#ifndef INC_FPRECISION
#define INC_FPRECISION

/*
 *******************************************************************************
 *
 *
 *                       Copyright (c) 2002-2021
 *                       Henrik Vestermark
 *                       Denmark
 *
 *                       All Rights Reserved
 *
 *   This source file is subject to the terms and conditions of the
 *   Future Team Software License Agreement which restricts the manner
 *   in which it may be used.
 *   Mail: hve@hvks.com
 *
 *******************************************************************************
*/

/*
 *******************************************************************************
 *
 *
 * Module name     :   fprecision.h
 * Module ID Nbr   :
 * Description     :   Arbitrary floating poiint precision class
 * --------------------------------------------------------------------------
 * Change Record   :
 *
 * Version	Author/Date		Description of changes
 * -------  -----------		----------------------
 * 03.01	HVE/3-Oct-2021	Switch to from Decimal String to Binary fptype float numbers. The previous chnage record has been removed to simplify
 * 03.02	HVE/19-Oct-2021	Passed all float testing 
 * 03.03	HVE/19-Nov-2021 Fixed complier bugs reported n GNU version 14 running on a Mac
 * 03.04	HVE/20-Nov-2021 A few bugs fixed and change to avoid compiler warnings
 * 03.05	HVE/22-Nov-2021	Minor change and improvement
 * 03.06	HVE/2-Dec-2021	Fix a bug in _float_precision_left_shift() that did not shift correctly if shift count was >= 64
 * 03.07	HVE/8-Dec-2021	Change eptype to an intmax_t insead of int to raise the exponent limit to more than 300M digits
 * 03.08	HVE/10-Dec-2021	Added two more trunking levels for faster handling of digits exceedig 10-100M decimal digits
 * 03.09	HVE/11-Dec-2021	atofp() can be called with a std::string or a char *
 * 03.10	HVE/25-Dec-2021 Added _float_precision_schonhage_strassen_linear_umul() to add better multiplication for medium size multiplication of digits<6,000 the function was modified
 *							from the int_precision counterpart
 * 03.11	HVE/6-Jan-2022	Added _INVSQRT3 and SQRT3 as a build in constant
 *
 * End of Change Record
 * --------------------------------------------------------------------------
*/

/* define version string */
static char _VF_[] = "@(#)fprecision.h 03.11 -- Copyright (C) Henrik Vestermark";

#include <algorithm>
#include "iprecision.h"

typedef uintmax_t fptype;	// The default size of the internal binary vector type unsigned 64bit
							// option is unsigned int  32bit, unsigned short 16bi or unsigned char 8bit
							// However performance will suffer if fptype < 64bit or the maximum natural integer size 
typedef intmax_t eptype;	// The defalut size for the exponent. int (64bit) should give a range from -2^64+1 to 2^64-1

const unsigned int Bitsfptype = sizeof(fptype) * 8; // Const use throughtout the source which is the number of bits the fptype can hold.

// The four different rounding modes
// # ROUND_NEAR  Rounded result is the closest to the infinitely precise result.
// # ROUND_DOWN  Rounded result is close to but no greater than the infinitely precise result.
// # ROUND_UP    Rounded result is close to but no less than the infinitely precise result.
// # ROUND ZERO  Rounded result is towards zero.
enum round_mode { ROUND_NEAR, ROUND_UP, ROUND_DOWN, ROUND_ZERO };

// The build in constant!
enum table_type { _LN2, _LN10, _PI, _EXP1, _SQRT2, _INVSQRT2, _SQRT3, _INVSQRT3 };

// Default precision of 20 Radix digits if not specified
static const size_t PRECISION = 20;

// Default max number of digits to convert natively to decimal 
static const int MAX_DECIMAL_DIGITS = sizeof(uintmax_t) >= 8 ? 19 : 9;
static const int MAX_OCTAL_DIGITS = sizeof(uintmax_t) >= 8 ? 21 : 10;
static const int MAX_BINARY_DIGITS = sizeof(uintmax_t) >= 8 ? 64 : 32;
static const int MAX_HEX_DIGITS = sizeof(uintmax_t) >= 8 ? 16 : 8;
static const int MAX_TRUNK_SIZE = 50;
static const int MAX_KILOTRUNK_SIZE = 10;
static const int MAX_MEGATRUNK_SIZE = 10;
static const int MAX_MEGA10TRUNK_SIZE = 10;
static const int MAX_MEGA100TRUNK_SIZE = 10;
static const int MAX_GIGATRUNK_SIZE = 10;

// Float_precision radix. Can be either BASE 2, BASE_10, BASE 16 or BASE_256
static const int F_RADIX = BASE_10;		// Obsolete

inline unsigned char FDIGIT(char x)				{ return (unsigned char)( x - '0'); }
inline unsigned char FCHARACTER( char x )		{ return (unsigned char)( x + '0'); }

inline int FCARRY( unsigned int x )				{ return (int)( x / BASE_10 ); }
inline int FSINGLE( unsigned int x )			{ return (int)( x % BASE_10 ); }

//
// @class float_precision_ctrl
// @author Henrik Vestermark (hve@hvks.com)
// @date  2/15/2006
// @version 1.0
// @brief  Arbitrary float precision control class
//
// @todo
//
//// Float Precision control class
//  This keep track of the global settings of default precision and round mode.
//  Everytime a new float_precision constructor is invoked it takes the default
//  precision and round mode from this float_precision_ctrl class. Unless a precision and/or
//  rounding mode has explicit been specified.
//  Default precision is the manifest constant PRECISION
//  Default rounding mode is ROUND_NEAR (round to nearest)
//
class float_precision_ctrl {
   enum round_mode   mRmode;  // Global Rounding mode. Default Round Nearest
   size_t      mPrec;   // Global Number of decimals in mantissa. Default PRECISION.

   public:
      // Constructor
      float_precision_ctrl( unsigned int p=PRECISION, enum round_mode rm=ROUND_NEAR ): mRmode(rm), mPrec(p) {}

      // Coordinate functions
      enum round_mode mode() const              { return mRmode; }
      enum round_mode mode( enum round_mode m ) { return( mRmode = m ); }
      size_t precision() const					{ return mPrec>0 ? mPrec : PRECISION; }
      size_t precision( unsigned int p )        { mPrec = p > 0 ? p : PRECISION; return mPrec; }
   };

extern class float_precision_ctrl float_precision_ctrl;

class float_precision;

// Arithmetic + Binary and Unary
template <class _Ty> inline float_precision operator+( float_precision&, const _Ty& );
template <class _Ty> inline float_precision operator+( const _Ty&, const float_precision& );
inline float_precision operator+( int_precision&, float_precision& );					// Override int_precision - other type in iprecision.h
//inline float_precision operator+( float_precision&,int_precision&);					// Override int_precision - other type in iprecision.h
inline float_precision operator+( const float_precision& );								// Unary

// Arithmetic - Binary and Unary
template <class _Ty> inline float_precision operator-(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator-(const _Ty&, const float_precision&);
inline float_precision operator-(int_precision&, float_precision&);						// Override int_precision - other type in iprecision.h
inline float_precision operator-( const float_precision& );								// Unary

// Arithmetic * Binary
template <class _Ty> inline float_precision operator*(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator*(const _Ty&, const float_precision&);
inline float_precision operator*(int_precision&, float_precision&);						// Override int_precision * other type in iprecision.h

// Arithmetic / Binary
template <class _Ty> inline float_precision operator/(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator/(const _Ty&, const float_precision&);
inline float_precision operator/(int_precision&, float_precision&);						// Override int_precision / other type in iprecision.h

// Arithmetic % Binary
template <class _Ty> inline float_precision operator%(float_precision&, const _Ty&);
template <class _Ty> inline float_precision operator%(const _Ty&, const float_precision&);
inline float_precision operator%(int_precision&, float_precision&);						// Override int_precision % other type in iprecision.h
																						
// Boolean Comparision Operators
inline bool operator> ( const float_precision&, const float_precision& );
inline bool operator< ( const float_precision&, const float_precision& );
inline bool operator==( const float_precision&, const float_precision& );
inline bool operator!=( const float_precision&, const float_precision& );
inline bool operator>=( const float_precision&, const float_precision& );
inline bool operator<=( const float_precision&, const float_precision& );

// Precision Floating point functions equivalent with the std C functions
extern float_precision modf( const float_precision&, float_precision * );
extern float_precision fmod( const float_precision&, const float_precision& );
extern float_precision floor( const float_precision& );
extern float_precision ceil( const float_precision& );
extern float_precision trunc(const float_precision&);
extern float_precision round(const float_precision&);
extern float_precision fabs( const float_precision& );  // Obsolete. replaced by overloaded abs(). But here for backward compatitbility
extern float_precision sqrt( const float_precision& );
extern float_precision log10( const float_precision& );
extern float_precision log( const float_precision& );
extern float_precision exp( const float_precision& );
extern float_precision pow( const float_precision&, const float_precision& );
extern float_precision frexp( float_precision&, eptype * );
extern float_precision ldexp( const float_precision&, eptype );
extern float_precision abs( const float_precision& );

// Trigonometric functions
extern float_precision atan( const float_precision& );
extern float_precision atan2( const float_precision&, const float_precision& );
extern float_precision asin( const float_precision& );
extern float_precision acos( const float_precision& );
extern float_precision sin( const float_precision& );
extern float_precision cos( const float_precision& );
extern float_precision tan( const float_precision& );

// Hyperbolic functions
extern float_precision sinh( const float_precision& );
extern float_precision cosh( const float_precision& );
extern float_precision tanh( const float_precision& );
extern float_precision asinh( const float_precision& );
extern float_precision acosh( const float_precision& );
extern float_precision atanh( const float_precision& );

// Miscelanneous support function
extern float_precision nroot(const float_precision&, unsigned int);

// Support functions. Works on float_precision
float_precision _float_precision_inverse( const float_precision& );
float_precision _float_table( enum table_type, size_t );
// Binary version 
float_precision _float_precision_atofp(const char *, size_t, enum round_mode);
float_precision _float_precision_atofp(const std::string&, size_t, enum round_mode);
//float_precision _float_precision_atofp2(const char *, size_t, enum round_mode);  // DEBUG
float_precision _float_precision_dtofp(double, size_t, enum round_mode);
double _float_precision_fptod(const float_precision *);
std::string _float_precision_fptoa(const float_precision *);
//std::string _float_precision_fptoaold(const float_precision *, int = 18);  // DEBUG

int_precision  _float_precision_fptoip(const float_precision *);
std::string _float_precision_fptoainteger(const float_precision *);

// Core Supporting functions. Works directly on string class. NEEDS Maybe tO BE REOMOVED after conversion
int _float_precision_rounding( std::string *, int, size_t, enum round_mode );
void _float_precision_strip_trailing_zeros( std::string * );
std::string _float_precision_uadd_short( std::string *, unsigned int );

// Core Supporting functions. Works directly on vector<fptype> class
size_t _float_precision_clz(const fptype);
size_t _float_precision_clz(const std::vector<fptype> &);
size_t _float_precision_ctz(const fptype);
size_t _float_precision_ctz(const std::vector<fptype> &);
void _float_precision_strip_leading_zeros(std::vector<fptype> *);
void _float_precision_strip_trailing_zeros(std::vector<fptype> *);
int _float_precision_normalize(std::vector<fptype> *);
int _float_precision_rounding(std::vector<fptype> *, const int, const size_t, const enum round_mode);
std::vector<fptype> _float_precision_right_shift(const std::vector<fptype> *, const size_t);
std::vector<fptype> _float_precision_left_shift(const std::vector<fptype> *, const size_t);
int _float_precision_compare(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_uadd_short(const std::vector<fptype> *, const fptype );
std::vector<fptype> _float_precision_uadd(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_usub_short(int *, const std::vector<fptype> *, const fptype);
std::vector<fptype> _float_precision_usub(int *, const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_umul_short(const std::vector<fptype> *, const fptype );
std::vector<fptype> _float_precision_umul(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_umul_school(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_umul_fourier(const std::vector<fptype> *, const std::vector<fptype> *, int=8);
std::vector<fptype> _float_precision_umul_linear(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_umul2_fourier(const std::vector<fptype> *, int = 8);
std::vector<fptype> _float_precision_umul2_linear(const std::vector<fptype> * );
std::vector<fptype> _float_precision_udiv_short(fptype *, const std::vector<fptype> *, const fptype);
std::vector<fptype> _float_precision_udiv(const std::vector<fptype> *, const std::vector<fptype> *);
std::vector<fptype> _float_precision_urem_short(const std::vector<fptype> *, const fptype);
std::vector<fptype> _float_precision_urem(const std::vector<fptype> *, const std::vector<fptype> *);

#if false
int_precision _int_precision_fastdiv(const int_precision &, const int_precision &);
int_precision _int_precision_fastrem(const int_precision &, const int_precision &);
#endif


//
// @class float_precision
// @author Henrik Vestermark (hve@hvks.com)
// @date  1/24/2005
// @version 1.0
// @brief  Arbitrary float precision class
//
// @todo
//
//// Float Precision class
//  An Arbitrary float always has the format [sign][digit][.[digit]*][E[sign][digits]+] where sign is either '+' or '-'
//  And is always stored in normalized mode after an operation or conversion
//  The length or the representation is always >= 2
//  A null string is considered as an error and an exception is thrown
//  Floating Point Numbers is stored in BASE 2^64 (Radix R=2^64). 
//	
//  Also number is always strip for leading nosignificant zeros
//
class float_precision {
   enum round_mode   mRmode;	// Rounding mode. Default Round Nearest
   size_t			 mPrec;		// Number of decimals in mantissa. Default 20, We make a shot cut by assuming the number of digits can't exceed 2^32-1 on 32bit or 2^64-1 on 64bit system
   eptype            mExpo;		// Exponent as a power of 2 as in IEEE 754. We make a short cut here and use the eptype to hold the exponent
								// the exponent. This will allow us exponent in the range from -RADIX^2^32 to RADIX^2^32. Which should be enough
   int				 mSign;		// The sign +1 for "+"a and -1 for "-". Notice in version 2+ the sign has been separated frm the mNumber, same as for int_precision
   std::vector<fptype> mBinary; // The binary vector of fptype that holds the float number. Per definition the vector when the constructor is invoked will always be initialized to zero if no argument is provided.
								// The fraction point is always after the first digit and is implied. mBinary[0] holds the most significant part of the number.
								// e.g. R=2^64. Number=mBinary[0]*R^0+mBinary[1]*R^-1+mBInary[2]*R^-2,...mBinary[n-1]*R^-(n-1) etc.
   //size_t			mPrecBits;	// Precision in Bits. Experimental use only
   
   public:
      // Constructors
	  float_precision();											// When initialized with no parameters
      float_precision( char, size_t, enum round_mode );				// When initialized through a char
      float_precision( unsigned char, size_t, enum round_mode );	// When initialized through a unsigned char
      float_precision( short, size_t, enum round_mode );			// When initialized through a short
      float_precision( unsigned short, size_t, enum round_mode );	// When initialized through a unsigned short
      float_precision( int, size_t, enum round_mode );				// When initialized through a int
      float_precision( unsigned int, size_t, enum round_mode );		// When initialized through a unsigned int
      float_precision( long, size_t, enum round_mode );				// When initialized through a long
      float_precision( unsigned long, size_t, enum round_mode );	// When initialized through a unsigned long
	  float_precision( int64_t, size_t, enum round_mode);			// When initialized through a int64_t
	  float_precision( uint64_t, size_t, enum round_mode);			// When initialized through a uint64_t
	  float_precision( float, size_t, enum round_mode);				// When initialized through a double
      float_precision( double, size_t, enum round_mode );			// When initialized through a double
      float_precision( const char *, size_t, enum round_mode );		// When initialized through a char string
	  float_precision( const std::string&, size_t, enum round_mode);// When initialized through a std::string
	  float_precision( const std::vector<fptype>& v) : mRmode(float_precision_ctrl.mode()), mPrec(float_precision_ctrl.precision()), mSign(1), mExpo(0), mBinary(v) {} // When initialized through a std::vector<iptype>. Notice sign will be 1 since vector<iptype> is unsigned
	  float_precision( const float_precision& s) : mRmode(s.mRmode), mPrec(s.mPrec), mExpo(s.mExpo), mSign(s.mSign), mBinary(s.mBinary) {}  // When initialized through another float_precision
      float_precision( const int_precision&, size_t, enum round_mode );

      // Coordinate functions
	  enum round_mode mode() const;
	  enum round_mode mode(enum round_mode m); 
	  eptype exponent() const;
	  eptype exponent(eptype e);
	  int sign() const;
	  int sign(int s);
	  int change_sign();
	  size_t precision() const;
	  size_t precision(size_t p);
	  std::vector<fptype> number();
	  std::vector<fptype> number(std::vector<fptype> &mb);
	  std::vector<fptype> *pointer();
	  fptype index(size_t inx)	const;
	  size_t size() const;
	  bool iszero() const;
      float_precision epsilon();				// Return Beta^(1-t)

	  float_precision assign(const float_precision& a);
	  	  
	  // Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
	  std::string toString() const; 
	  std::string toFixed(int);
	  std::string toPrecision(int);
	  std::string toExponential(int);
	  float_precision toInteger();
	  float_precision toFraction();

	  // Implicit/explicit conversion operators
      operator char() const;
      operator short() const;
      operator int() const;
      operator long() const;
      operator unsigned char() const;
      operator unsigned short() const;
      operator unsigned int() const;
      operator unsigned long() const;
	  operator long long() const;
	  operator unsigned long long() const;
	  operator float() const;
      operator double() const;
      operator int_precision() const;

      // Essential operators
      float_precision& operator= ( const float_precision& );
      float_precision& operator+=( const float_precision& );
      float_precision& operator-=( const float_precision& );
      float_precision& operator*=( const float_precision& );
      float_precision& operator/=( const float_precision& );
	  float_precision& operator%=( const float_precision& );

      // Specialization
      friend std::ostream& operator<<( std::ostream& strm, const float_precision& d );
      friend std::istream& operator>>( std::istream& strm, float_precision& d );

      // Exception class
      class bad_int_syntax		{};
      class bad_float_syntax	{};
      class out_of_range		{};
      class divide_by_zero		{};
      class domain_error		{};
      class base_error			{};
   };


//////////////////////////////////////////////////////////
//
// FLOAT PRECISION CONSTRUCTORS
//
//////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		27/Nov/2021
//	@brief 		constructor for float_precision
//	@return 	nothing
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with empty parameter
//
inline float_precision::float_precision() 
	{
	mRmode = float_precision_ctrl.mode();
	mPrec = float_precision_ctrl.precision();
	//mPrecBits = (size_t)(log(mPrec) + 0.5);
	mExpo = 0;
	mBinary.assign(1, 0);	// Default set to 0
	mSign = +1;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		5/Sep/2021
//	@brief 		Char constructor for float_precision
//	@return 		nothing
//	@param      "c"	-	char digit
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with a character
//
inline float_precision::float_precision( const char c, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	if (c < 0) { mSign = -1; }
	mBinary.assign( 1, abs(c) ); 
	if (c != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Unsigned Char constructor for float_precision
//	@return 		nothing
//	@param      "c"	-	Integer unsigned char digit
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initilize with a character
//
inline float_precision::float_precision( const unsigned char c, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	mBinary.assign( 1, c );
	if (c != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


// 	@author Henrik Vestermark (hve@hvks.com)
// 	@date		6/Sep/2021
// 	@brief 		Short constructor for float_precision
// 	@return 		nothing
// 	@param      "i"	-	Integer number
// 	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
// 	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
// 
// 	@todo
// 
//  Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( short i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	if (i < 0 ) mSign = -1;
	mBinary.assign( 1, abs(i) );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Unsigned Short constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( unsigned short i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	mBinary.assign( 1, i );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Integer constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( int i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	if (i < 0) mSign = -1;
	mBinary.assign( 1, abs(i) );
	if(i!=0)
		mExpo+= _float_precision_normalize(&mBinary);		
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Unsigned Integer constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( unsigned int i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	mBinary.assign( 1, i );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		long constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( long i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	if (i < 0) mSign = -1;
	mBinary.assign( 1, abs(i) );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Unsigned Long constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( unsigned long i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	mBinary.assign( 1, i );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		int64_t 64bit constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( int64_t i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	if (i < 0) mSign = -1;
	mBinary.assign( 1, abs(i) );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Sep/2021
//	@brief 		Usigned 64bit constructor for float_precision
//	@return 		nothing
//	@param      "i"	-	64bit Integer number
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate and initialize with integer
//
inline float_precision::float_precision( uint64_t i, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	mExpo = 0;
	mRmode = m;
	mPrec = p;
	mSign = +1;
	mBinary.resize(2, 0);
	mBinary.assign( 1, i );
	if (i != 0)
		mExpo += _float_precision_normalize(&mBinary);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		string constructor for float_precision
//	@return 		nothing
//	@param      "str"	-	Floating point number as a string
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate input and convert to internal representation
//  Always add sign if not specified
//  Only use core base functions to create multi precision numbers
//  The float can be any integer or decimal float representation
//
inline float_precision::float_precision( const char *str, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   if( str == NULL || *str == '\0' )
      { throw bad_int_syntax(); return; }

   mRmode = m;  
   mPrec = p;   
   *this = _float_precision_atofp( str, p, m );
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  Aug-19-2016
//	@brief 		std::string constructor for float_precision
//	@return 		nothing
//	@param      "str"	-	Floating point number as a std::string
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor
//  Validate input and convert to internal representation
//  Always add sign if not specified
//  Only use core base functions to create multi precision numbers
//  The float can be any integer or decimal float representation
//

inline float_precision::float_precision(const std::string& str, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	if (str.empty())
	{ throw bad_int_syntax(); return; }

	mRmode = m;
	mPrec = p;
	*this = _float_precision_atofp( str.c_str(), p, m);
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		float constructor for float_precision
//	@return 		nothing
//	@param      "d"	-	Floating point number in IEE754 format
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor for float
//  Validate input and convert to internal representation
//  The float can be any integer or decimal float representation
//	  Convert the float to double and call _float_precision_dtofp
//
inline float_precision::float_precision(float f, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode())
	{
	mRmode = m;
	mPrec = p;
	*this = _float_precision_dtofp((double)f, p, m);   
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		7/Sep/2021
//	@brief 		double constructor for float_precision
//	@return 		nothing
//	@param      "d"	-	Floating point number in IEE754 format
//	@param      "p"	-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"	-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor for double
//  Validate input and convert to internal representation
//  The float can be any integer or decimal float representation
//
inline float_precision::float_precision( double d, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   mRmode = m;
   mPrec = p;
   *this = _float_precision_dtofp( d, p, m );
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		30/Sep/2021
//	@brief 		int_precision constructor for float_precision
//	@return 	nothing
//	@param      "ip"	-	arbitrary integer precision
//	@param      "p"		-	Number of precision (default float_precision_ctrl.precision())
//	@param      "m"		-	rounding mode (default float_precision_ctrl.mode())
//
//	@todo
//
// Description:
//  Constructor for int_precision to float_precision
//
inline float_precision::float_precision( const int_precision& ip, size_t p = float_precision_ctrl.precision(), enum round_mode m = float_precision_ctrl.mode() )
   {
   std::string s;
   const std::vector<iptype> inumber(ip.number());

   mSign = ip.sign();
   for (size_t i = inumber.size(); i > 0; --i)
	   mBinary.push_back((fptype)inumber[i - 1]);
   mRmode = m;
   mPrec = p;
   mExpo =_float_precision_normalize(&mBinary);
   mExpo += (int)( ( inumber.size() - 1) * Bitsfptype );
   mExpo += _float_precision_rounding(&mBinary, ip.sign(), mPrec, mRmode);
   }

//////////////////////////////////////////////////////////
//
// END FLOAT PRECISION CONSTRUCTORS
//
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
//
//
//    Implict/Explicit conversions to base types int, short, long, char
//
//
//////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		float_precision::operator double
//	@return 	double -
//
//	@todo  Add to do things
//
// Description:
//   conversion from float_precision to double operator
//
inline float_precision::operator double() const
   {// Conversion to double
   return _float_precision_fptod( this );
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator float
//	@return 	float
//
//	@todo  Add to do things
//
// Description:
//   Conversion from float_precision to float
//
inline float_precision::operator float() const
   {// Conversion to float
   return (float)(double)*this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Oct/2021
//	@brief 		float_precision::operator long long
//	@return 	long long
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to long long (64bit)
//
inline float_precision::operator long long () const
	{// Conversion to long long
	return (long long)(unsigned long long)*this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator long
//	@return 	long
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to long
//
inline float_precision::operator long() const
   {// Conversion to long
   return (long)(unsigned long long)*this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator int
//	@return 	int
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to int
//
inline float_precision::operator int() const
   {// Conversion to int
   return (int)(unsigned long long)*this;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator short
//	@return 	short
//
//	@todo  Add to do things
//
// Description:
//  Conversion float_precision to short
//
inline float_precision::operator short() const
   {// Conversion to short
   return (short)(unsigned long long)*this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator char
//	@return 	char
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to char
//
inline float_precision::operator char() const
   {// Conversion to char
   return (char)(unsigned long long)*this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Oct/2021
//	@brief 		float_precision::operator unsigned long long
//	@return 	unsigned long long
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to unsigned long long (64bit)
//
inline float_precision::operator unsigned long long() const
	{// Conversion to unsigned long long
	unsigned long long ull;
	unsigned shift;
	if (mExpo < 0)	{ return ull=0; }
	if (mExpo == 0) { return ull=(unsigned long long)mBinary[0]; }  // Return either 1 or zero
	if (mExpo > 0 && mBinary.size() == 1)  // true power of 2
		{
		ull = mBinary[0];
		if (mExpo < 64) 
			ull <<= mExpo;
		else
			ull = 0;  // Overflow all lower 4bits is zero
		return ull;
		}
	//All other cases
	size_t n = (mExpo - 1) / Bitsfptype + 1;
	if(mBinary.size()<n-1) 
		{
		return ull = 0;	// number exceed ull max number
		}
	// Take the lowest 64bit of the number
	if (mBinary.size() > n)
		ull = (unsigned long long)mBinary[n];
	else
		ull = 0;
	shift = mExpo%Bitsfptype;
	if (shift == 0)
		return ull;
	ull >>= Bitsfptype - shift;
	ull |= (unsigned long long)(mBinary[n - 1]) << shift;
	return ull;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator unsigned long
//	@return 	unsigned long
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to unsigned long
//
inline float_precision::operator unsigned long() const
   {// Conversion to unsigned long
   return (unsigned long)(unsigned long long)*this;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator unsigned int
//	@return 	unsigned int
//
//	@todo  Add to do things
//
// Description:
//  Conversion float_precision to unsigned int
//
inline float_precision::operator unsigned int() const
   {// Conversion to unsigned int
   return (unsigned int)(unsigned long long)*this;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator unsigned short
//	@return 	unsigned short
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to unsigned short
//
inline float_precision::operator unsigned short() const
   {// Conversion to unsigned short
   return (unsigned short)(unsigned long long)*this;
   }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		3/12/2006
//	@brief 		float_precision::operator unsigned char
//	@return 	unsigned char
//
//	@todo  Add to do things
//
// Description:
//  Conversion from float_precision to unsigned char
//
inline float_precision::operator unsigned char() const
   {// Conversion to unsigned char
   return (unsigned char)(unsigned long long)*this;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Sep/2021
//	@brief 		float_precision::operator int_precision
//	@return 	int_precision
//
//	@todo  Add to do things
//
// Description:
//  Conversion float_precision to int_precision
//
inline float_precision::operator int_precision() const
    {// Conversion to int_precision
    return _float_precision_fptoip( this );
    }

////////////////////////////////////////////////////////////////////////////////////
//
// Class Methods
//
///////////////////////////////////////////////////////////////////////////////////

// Coordinate functions

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::mode
//	@return 	return a copy of the current rounding mode
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mRmode
//
inline enum round_mode float_precision::mode() const 
	{ return mRmode; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::mode
//	@param		"m"	-	New rounding mode
//	@return 	Set and return a copy of the new rounding mode
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mRmode
//
inline enum round_mode float_precision::mode(enum round_mode m) { return(mRmode = m); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::exponent
//	@return 	return a copy of the current exponent
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mExpo
//
inline eptype float_precision::exponent() const
	{ return mExpo; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::exponent
//	@param		"e"	-	New exponent
//	@return 	Set and return a copy of the new exponent
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mExpo
//	Notice that exponent is a power of base 2, 2^expo
//
inline eptype float_precision::exponent(eptype e) 
	{ return(mExpo = e); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::sign
//	@return 	return a copy of the sign
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mSign
//
inline int float_precision::sign() const 
	{ return mSign; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::sign
//	@param		"mSign"	-	New mSign
//	@return 	Set and return a copy of the new sign
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mSign
//
inline int float_precision::sign(int s) 
	{ /*if(s<0 && mBinary.size() == 1 && mBinary[0] == 0 ) return mSign; */ return (mSign = s); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::change_sign
//	@return 	Change and return a copy of the sign
//
//	@todo  Add to do things
//
// Description:
//  Change sign and Return a copy of mSign
//
inline int float_precision::change_sign()
	{ mSign *= -1;  return mSign; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::precision
//	@return 	return a copy of the current precision
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mPrec
//
inline size_t float_precision::precision() const
	{ return mPrec; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::precision
//	@return 	Set and return a new copy of the decimal precision
//
//	@todo  Add to do things
//
// Description:
//  Set and return a copy of the new decimal precision
//	If new precision is less than previous precision the number is convered to new precision 
//	using the current rounding mode for the float_precision object
//
inline size_t float_precision::precision(size_t p) 
	{
	std::string m;
	mPrec = p > 0 ? p : float_precision_ctrl.precision();
	//mPrecBits = (size_t)ceil((mPrec*log2(BASE_10) + 0.5));
	mExpo += _float_precision_rounding(&mBinary, mSign, mPrec, mRmode);
	return mPrec;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::number
//	@return 	return a copy of mBinary
//
//	@todo  Add to do things
//
// Description:
//  Return a copy of mBinary
//
inline std::vector<fptype> float_precision::number() 
	{ return mBinary; }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::number
//	@param		"mBinary"	-	New mBinary vector<fptype>
//	@return 	Set and return a copy of mBinary
//
//	@todo  Add to do things
//
// Description:
//  Set and Return a copy of mBinary
//
inline std::vector<fptype> float_precision::number(std::vector<fptype> &mb) 
	{ return mBinary = mb; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		2/Nov/2021
//	@brief 		float_precision::pointer
//	@return 	return a pointer to mBinary
//
//	@todo  Add to do things
//
// Description:
//  Return a pointer to mBinary
//
inline std::vector<fptype> *float_precision::pointer()
	{ return &mBinary; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/201
//	@brief 		float_precision::index
//	@return		The the index  of mBinary[inx]
//	@param		"inx"	the index in the range [0..mBinary.size()]
//
//	@todo  
//
// Description:
//   Return the index of the mBinary vector<fptype>
//
inline fptype float_precision::index(size_t inx)	const
	{ return mBinary[inx]; }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/201
//	@brief 		float_precision::size
//	@return		The size of mBinary
//
//	@todo  
//
// Description:
//   Return the size of the mBinary vector<fptype>
//
inline size_t float_precision::size() const 
	{ return mBinary.size(); }

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/201
//	@brief 		float_precision::iszero
//	@return		Boolean true is float_precision object is zero otherwise false
//
//	@todo  
//
// Description:
//   compare with 0 and return the boolean comparision value
//
inline bool float_precision::iszero() const 
	{ return mBinary.size() == 1 && mBinary[0] == 0; }// Notice both +0 and -0 is allowed and return true
																							
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/Oct/2021
//	@brief 		float_precision::epsilon  return the epsilon such that 1.0+epsilon!=1.0
//	@return		float_precision	Return the epsilon for the given Radix and precision
//
//	@todo  
//
// Description:
//   Return the epsilon: The number where 1.0+epsilon!=1.0
//	 This function was rewritten to use B^1-<precision> code with floating arguments instead of below standard pow()
//   that was way to time consuming to execute in favor of this algorithm
//   res = pow(float_precision(F_RADIX, mPrec), float_precision(1-(int)mPrec));  // beta^1-t
//
inline float_precision float_precision::epsilon()
	{
	float_precision res(1);
	res.exponent(1 - (int)(ceil(mPrec*log2(BASE_10)))); // Number of precision bits needed
	return res;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::assign  assign the right hand side to the left hand side
//	@param		"a"	-	the right hand side of the assignment
//	@return 	The assigned floating point object
//
//	@todo
//		TBD
//
// Description:
//  assign the right hand side to the left hand side. Notice the left hand side field mRmode and mPrec is preserved
//	same as the = operator
//
inline float_precision float_precision::assign(const float_precision& a)
	{
	mRmode = a.mRmode;
	mPrec = a.mPrec;
	mExpo = a.mExpo;
	mBinary = a.mBinary;
	mSign = a.mSign; 
	return *this; }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::toInteger()  reduce float precision to a integer and return the discared fraction part
//	@return 	the fraction part as a float_precision object
//
//	@todo
//		TBD
//
// Description:
//  reduce float precision to a integer and return the discared fraction part
//
inline float_precision float_precision::toInteger()
	{// Trunck towards zero the floating number to its integer part and return the remaing fraction part
	float_precision frac = *this;
	if (mExpo < 0)
		{ mBinary.assign(1, 0); mExpo = 0;  return frac; }
	if (mExpo == 0) 
		{ mBinary.assign(1, mBinary[0]);  frac.mBinary[0] = 0; frac.mExpo += _float_precision_normalize(&frac.mBinary); return frac; }
	// fecth expo bits from the fp number after the '.'
	size_t n = mExpo / Bitsfptype + 1;
	unsigned bn = mExpo % Bitsfptype;
	// Discard excessive fptype digits
	int extra = bn == 0 ? 0 : 1;
	if (mBinary.size()> n + extra)
		mBinary.erase(mBinary.begin() + n + extra, mBinary.end());
	// Discard excessive bits within the last fptype.
	if (n < mBinary.size())
		{
		mBinary[n] &= ((~(fptype)0) << (Bitsfptype - bn));
		}
	// Adjust the fraction part
	frac.mBinary = _float_precision_left_shift(&frac.mBinary, frac.mExpo);
	if (frac.mBinary.size() == 0) frac.mBinary.push_back(0); else frac.mBinary[0] = 0;
	frac.mExpo = _float_precision_normalize(&frac.mBinary);
	_float_precision_strip_trailing_zeros(&frac.mBinary);
	return frac;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date		22/Nov/2021
//	@brief 		float_precision::toFraction  reduce float precision to a fraction and return the discared integer part
//	@return 	the integer part as a float_precision object
//
//	@todo
//		TBD
//
// Description:
//  reduce float precision to a fraction and return the discared integer part
//
inline float_precision float_precision::toFraction() 
	{// Return the integer as a float_precision and discard the inteer portion on this
	const fptype mask = ~(fptype)0;
	float_precision integer = *this;
	if (mExpo < 0) 
		{ integer.mBinary.assign(1, 0); integer.mExpo = 0; integer.mSign = 1;  return integer; }
	if (mExpo == 0) 
		{ integer.mBinary.assign(1, mBinary[0]); mBinary[0] = 0; mExpo += _float_precision_normalize(&mBinary); return integer; }
	// Fecth expo bits from the fp number after the '.' to adjust the integer portion
	const size_t n = integer.mExpo / Bitsfptype + 1;
	const unsigned bn = integer.mExpo % Bitsfptype;  // mExpo always>0 
	// Discard excessive fptype digits for the integer
	const int extra = bn == 0 ? 0 : 1;
	if (integer.mBinary.size()> n + extra)
		integer.mBinary.erase(integer.mBinary.begin() + n + extra, integer.mBinary.end());
	// Discard excessive bits within the last fptype.
	if (n < integer.mBinary.size())
		{
		integer.mBinary[n] &= (mask << (Bitsfptype - bn));
		}
	// Adjust the fraction part
	mBinary = _float_precision_left_shift(&mBinary, mExpo);
	if (mBinary.size() == 0) mBinary.push_back(0); else
		mBinary[0] = 0;

	mExpo = _float_precision_normalize(&mBinary);
	_float_precision_strip_trailing_zeros(&mBinary);
	return integer;
	}

//////////////////////////////////////////////////////////
//
// BEGIN FLOAT PRECISION FORMATTING FUNCTIONS
//
//////////////////////////////////////////////////////////

// Conversion methods. Safer and less ambiguios than overloading implicit/explivit conversion operators
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Oct/2021
//	@brief 		float_precision::toString()  float precision toString method
//	@return 	the string representation of the float_precision object
//
//	@todo
//		TBD 
//
// Description:
//  return the string value value of the float_precision object
//
inline std::string float_precision::toString() const 
	{ return _float_precision_fptoa(this); }

// .toFixed()
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Oct/2021
//	@brief 		float_precision::toFixed(fix)  float precision toFixed method
//	@return 	the string representation of the fixed presentation
//	@param		"fix"	-	Number of fixed decimal
//
//	@todo
//		remove the need to do a string rounding and adding 
//
// Description:
//  return the string value value of the fixed precision of the float_precision number
//	  same functionality as the javascript .toFixed() method	
//
inline std::string float_precision::toFixed(int fix = 0)
{
	std::string ss;
	int sign, expo;
	size_t inx;
	float_precision fp = *this;

	if (fix < 0) fix = 0;
	ss = this->toString();					// Now we have it in exponetial form and in Base 10 
	sign = this->mSign;						// get sign
	if (sign<0)
		ss.erase(0, 1);						// Erase sign
	inx = ss.find("E");						// Find start of Exponent
	expo = atoi(ss.substr(inx + 1).c_str());  // Get exponent value
	ss.erase(inx, std::string::npos);		// Erase exponent value from string
	ss.erase(1, 1);							// Erase .
	inx = 1;								// Where dot should be inserted			
	if (expo > 0)
	{
		ss.insert(ss.length(), expo, '0'); inx += expo;
	}		// Trailing with zeros
	else if (expo < 0) ss.insert(0, -expo, '0');				// Padd with leading zeros
	_float_precision_rounding(&ss, sign, inx + fix, this->mRmode);  // Round to fix. NOTICE string rounding 
	if (ss.length()<(unsigned)(inx + fix))
		ss.insert(ss.length(), inx + fix - ss.length(), '0');		// Add trailing zeros
	if (ss.length()>inx) ss.insert(inx, 1, '.');					// Insert fraction unless it after the last digit
	if (sign < 0 && ss != "0") ss.insert(0, 1, '-');		// Add sign if negative and not 0
	return ss;								// Return formatted representation of number
}

// .toPrecision()
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Oct/2021
//	@brief 		float_precision::toPrecision(fix)  float precision toPrecision method
//	@return 	the string representation of the Precision presentation
//	@param		"fix"	-	Number of fixed decimal
//
//	@todo
//		remove the need to do a string rounding and adding 
//
// Description:
//  return the string value value of the precision of the float_precision number
//	  same functionality as the javascript .toPrecision() method	
//
inline std::string float_precision::toPrecision(int fix = 1)
{
	std::string ss;
	int sign, expo;
	size_t inx, shf;

	if (fix <= 1) fix = 1;
	sign = this->mSign;
	ss = this->toString();					// Now we have it in exponetial form and in Base 10 with leading sign
	if (sign<0)
		ss.erase(0, 1);						// Erase sign
	inx = ss.find("E");						// Find start of Exponent
	expo = atoi(ss.substr(inx + 1).c_str());// Get exponent value
	ss.erase(inx, std::string::npos);		// Erase exponent value from string
	ss.erase(1, 1);							// Erase .
	inx = 1;								// Where dot should be inserted			
	_float_precision_rounding(&ss, sign, fix, this->mRmode);  // Round to fix. old fashion string rounding
	if (expo >= 0)
	{
		if ((unsigned)fix> ss.length()) ss.insert(ss.length(), fix - ss.length(), '0');  // Trailing with zeros, so we have fix decimals
		if (expo > 0 && inx < ss.length())	 // Adjust the decimal sign as long as we can accomodate all digits
		{
			shf = ss.length() - inx;
			if (shf >(unsigned)expo) shf = expo;
			expo -= (int)shf;
			inx += shf;
		}
	}
	else
	{  // Expo < 0
		if (expo<0 && (unsigned)fix > ss.length())	// Room for adding leading zeros to accomodate the exponen which is negative
		{
			shf = (unsigned)fix - ss.length();
			if (shf > (unsigned)-expo) shf = -expo;
			expo += (int)shf; ss.insert(0, shf, '0');
		}
	}
	ss = ss.substr(0, inx) + ((unsigned)fix > inx ? "." : "") + ss.substr(inx, fix);
	if (sign < 0) ss.insert(0, 1, '-');		// Add sign if negative
	if (expo != 0) { ss += "E"; ss += (expo < 0 ? "-" : ""); ss += itostring(abs(expo), BASE_10); }
	return ss;
}

// .toExponential()
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		10/Oct/2021
//	@brief 		float_precision::toExponential(fix)  float precision toExponential method
//	@return 	the string representation of the toExponential presentation
//	@param		"fix"	-	Number of fixed decimal
//
//	@todo
//		remove the need to do a string rounding and adding 
//
// Description:
//  return the string value value of the toExponential precision of the float_precision number
//	  same functionality as the javascript .toExponential() method	
//
inline std::string float_precision::toExponential(int fix = 0)
{
	std::string ss;
	int sign, expo;
	size_t inx;

	if (fix < 0) fix = 0;
	sign = this->mSign;
	ss = this->toString();					// Now we have it in exponetial form and in Base 10 with leading sign
	if (sign<0)
		ss.erase(0, 1);						// Erase sign
	inx = ss.find("E");						// Find start of Exponent
	expo = atoi(ss.substr(inx + 1).c_str());// Get exponent value
	ss.erase(inx, std::string::npos);		// Erase exponent value from string
	ss.erase(1, 1);							// Erase .
	inx = 1;								// Where dot should be inserted			
	_float_precision_rounding(&ss, sign, fix, this->mRmode);  // Round to fix. Old fashion string rounding
	if ((unsigned)fix> ss.length()) ss.insert(ss.length(), fix - ss.length(), '0');  // Trailing with zeros, so we have fix decimals
	ss = ss.substr(0, inx) + ((unsigned)fix > inx ? "." : "") + ss.substr(inx, fix);
	if (sign < 0) ss.insert(0, 1, '-');		// Add sign if negative
	if (expo != 0) { ss += "E"; ss += (expo < 0 ? "-" : ""); ss += itostring(abs(expo), BASE_10); }

	return ss;
}

//////////////////////////////////////////////////////////
//
// END FLOAT PRECISION FORMATTING FUNCTIONS
//
//////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////
//
// FLOAT PRECISION OPERATORS
//
//////////////////////////////////////////////////////////


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	Assign float precision numbers
//	@return 	float_precision&	-
//	@param   "a"	-	float precsion number to assign
//
//	@todo
//
// Description:
//  Assign operator
//  Round it to precision and mode of the left hand side
//  Only the exponent and mantissa is assigned
//  Mode and precision is not affected by the assignment.
//
inline float_precision& float_precision::operator=( const float_precision& a )
   {
   mExpo = a.mExpo;
   mSign = a.mSign;
   mBinary = a.mBinary;
   if( _float_precision_rounding( &mBinary, mSign, mPrec, mRmode ) != 0 )  // Round back to left hand side precision
      mExpo++;

   return *this;
   }


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	+= float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number to assign
//
//	@todo    Still missing code for x += a where add make sense. fx. if a is so small it does
//           not affect the result within the given precision is should ignored. same is true
//           if x is insignififcant comapre to a the just assign a to x
//
// Description:
//  The essential += operator
//  1) Align to same exponent
//  2) Align to same precision
//  3) Add Mantissa
//  4) Add carry to exponent
//  4) Normalize
//  5) Rounding to precission
//  Early out algorithm. i.e.
//     - x+=0 return x
//     - x+=a wher x is 0 return a
//
inline float_precision& float_precision::operator+=( const float_precision& a )
	{
	int sign, sign1, sign2, wrap;
	eptype expo_max;
	int digits_max;
	size_t precision_max;
	std::vector<fptype> s, s1, s2;

	if( a.iszero() )  // Add zero
		return *this;
	if( iszero() )      // Add a (not zero) to *this (is zero) Same as *this = a;
		return *this = a;

	// extract sign and unsigned portion of number
	sign1 = a.mSign;
	s1 = a.mBinary;		// Extract Mantissa
	sign2 = mSign;
	s2 = mBinary;		// Extract Mantissa
	expo_max = std::max( mExpo, a.mExpo );
	precision_max = std::max( mPrec, a.mPrec );

	// Check if add makes sense. Still missing

	// Right shift (padd leading zeros) to the smallest number
	if( a.mExpo != expo_max )
		s1 =_float_precision_right_shift( &s1, expo_max - a.mExpo );
	if( mExpo != expo_max )
		s2 = _float_precision_right_shift( &s2, expo_max - mExpo );
	// Now s1 and s2 is aligned to the same exponent. The biggest of the two

	// Round to same precision
	if( _float_precision_rounding( &s1, sign1, precision_max, a.mode() ) != 0 ) // If carry when rounding up then one right shift
		_float_precision_right_shift( &s1, 1 );
	if( _float_precision_rounding( &s2, sign2, precision_max, mRmode ) != 0 ) // If carry when rounding up then one right shift
		_float_precision_right_shift( &s2, 1 );
	
	digits_max = (int)std::max(s1.size(), s2.size());

	if( sign1 == sign2 )
		{
		s = _float_precision_uadd( &s1, &s2 );
		if( s.size() > digits_max ) // One more digit
			expo_max++;
		sign = sign1;
		}
	else
		{
		int cmp = _float_precision_compare( &s1, &s2 );
		if( cmp > 0 ) // Since we subctract less the wrap indicater need not to be checked
			{
 			s = _float_precision_usub( &wrap, &s1, &s2 );
			sign = sign1;
			}
		else
			if( cmp < 0 )
				{
				s = _float_precision_usub( &wrap, &s2, &s1 );
				sign = sign2;
				}
			else
				{  // Result zero
				sign = 1;
				s.push_back(0);
				expo_max = 0;
				}
		}

	expo_max += _float_precision_normalize( &s );            // Normalize the number
	if( _float_precision_rounding( &s, sign, mPrec, mRmode ) != 0 )  // Round back left hand side precision
		expo_max++;

	mSign = sign;
	mBinary = s;
	mExpo = expo_max;
	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	-= float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number to assign
//
//	@todo
//
// Description:
//  The essential -= operator
//  n = n - a is the same as n = n + (-a). so change sign and use the += operator instead
//
inline float_precision& float_precision::operator-=( const float_precision& a )
	{
	float_precision b;

	b.precision( a.precision() );
	b.mode( a.mode() );
	b = a;
	b.change_sign();
	*this += b;

	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	*= float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number to assign
//
//	@todo
//
// Description:
//  The essential *= operator
//  1) Multiply mantissa
//  2) Add exponent
//  3) Normalize
//  4) Rounding to precision
//
inline float_precision& float_precision::operator*=( const float_precision& a )
	{
	eptype expo_res;
	int sign, sign1, sign2;
	std::vector<fptype> s, s1, s2;

	// extract sign and unsigned portion of number
	sign1 = a.mSign;
	s1 = a.mBinary;
	sign2 = mSign;
	s2 = mBinary;
	sign = sign1 * sign2;
	s = _float_precision_umul(&s1, &s2);
	expo_res = mExpo + a.mExpo;
	if( s.size() -1 > s1.size() + s2.size() -2 ) // A carry
		expo_res++;
	expo_res += _float_precision_normalize( &s );            // Normalize the number
	if( _float_precision_rounding( &s, sign, mPrec, mRmode ) != 0 )  // Round back left hand side precision
		expo_res++;

	if( sign == -1 && s.size() == 1 && s[0] == 0 )  // Avoid -0 as result +0 is right
		sign = 1; // Change sign
	if( s.size() == 1 && s[0] == 0 )  // Result 0 clear exponent
		expo_res = 0;

	mSign = sign;
	mExpo = expo_res;
	mBinary = s;
	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/27/2006
//	@brief 	/= float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number to assign
//
//	@todo
//
// Description:
//  The essential /= operator
//  We do a /= b as a *= (1/b)
// Bug
//  1/27/2006 Inverse was always done with the precision of a instead of the Max precision of both this & a
//
inline float_precision& float_precision::operator/=( const float_precision& a )
	{
	float_precision c;

	if (this->iszero() ) // If divisor is zero the result is zero
		return *this;

	c.precision( a.precision() );
	if( a.precision() < mPrec )
		c.precision( mPrec );
	c = a;
	*this *= _float_precision_inverse(c);
	return *this;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  24-Mar-2021
//	@brief 	%= float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number to assign
//
//	@todo
//
// Description:
//  The essential /= operator
//  We do a %= b as fmod(a,b)
//
inline float_precision& float_precision::operator%=(const float_precision& a)
	{
	float_precision c;
	
	if (this->iszero()) // If divisor is zero the result is zero
		return *this;

	c.precision(a.precision());
	if (a.precision() < mPrec)
		c.precision(mPrec);
	c = a;
	*this = fmod( *this, c );
	return *this;
	}


//////////////////////////////////////////////////////////////////////////
//
//
// Mixed Mode arithmetic Unary and Binary operators +, -
//
//
//////////////////////////////////////////////////////////////////////////

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  10/Aug/2014
//	@brief 			operator+
//	@return 	float_precision	-	return addition of lhs + rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for float_precision + float_precision
//		Specialization for float_precision to avoid ambigous overload
//
/*
inline float_precision operator+( float_precision& lhs, float_precision& rhs)
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c += lhs;
	}
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  3/19/2006
//	@brief 			operator+
//	@return 	float_precision	-	return addition of lhs + rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Add operator for float_precision + <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator+( float_precision& lhs, const _Ty& rhs )
	{
	float_precision c(rhs);

	if( lhs.precision() > c.precision() )
		c.precision( lhs.precision() );

	return c += lhs;
	}


//  @author Henrik Vestermark (hve@hvks.com)
//  @date  3/19/2006
//  @version 1.0
//	@brief 			operator+
//	@return 	float_precision	-	return addition of lhs + rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Add operator for <any other type> + float_precision
//
template <class _Ty> inline float_precision operator+( const _Ty& lhs, const float_precision& rhs )
	{
	float_precision c(lhs);

	if( rhs.precision() > c.precision() )
		c.precision( rhs.precision() );

	return c += rhs;
	}

//  @author Henrik Vestermark (hve@hvks.com)
//  @date  3/19/2006
//  @version 1.0
//	@brief 			operator+
//	@return 	float_precision	-	return addition of lhs + rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Add operator for int_precision + float_precision
//
inline float_precision operator+( int_precision& lhs, float_precision& rhs )
	{
	float_precision c(lhs);

	if( rhs.precision() > c.precision() )
		c.precision( rhs.precision() );

	return c += rhs;
	}


//  @author Henrik Vestermark (hve@hvks.com)
//  @date  Aug/9/2014
//  @version 1.0
//	@brief 			operator+
//	@return 	float_precision	-	return addition of lhs + rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Add operator for int_precision + float_precision
//
inline float_precision operator+( float_precision& lhs, int_precision& rhs )
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c += lhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	Unary + float precision number
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number
//
//	@todo
//
// Description:
//  Unary add. Do nothing and return a
//
inline float_precision operator+( const float_precision& a )
	{
	// Otherwise do nothing
	return a;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator-
//	@return 	float_precision	-	return addition of lhs - rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Sub operator for float_precision - <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator-( float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d -= c;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator-
//	@return 	float_precision	-	return addition of lhs - rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Sub operator for  <any other type> - float_precision
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator-( const _Ty& lhs, const float_precision& rhs )
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d -= c;
	}

//  @author Henrik Vestermark (hve@hvks.com)
//  @date  7/29/2014
//  @version 1.0
//	@brief 			operator-
//	@return 	float_precision	-	return addition of lhs - rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Add operator for int_precision - float_precision
//
inline float_precision operator-(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c -= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	- float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	first float precsion number
//	@param   "b"	-	second float precsion number to subtract
//
//	@todo
//
// Description:
//  Binary subtract two float_precision numbers
//  Implenting using the essential -= operator
//
/*
inline float_precision operator-( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision );
   c = a;
   c -= b;

   return c;
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	Unary - float precision number
//	@return 	the resulting float_precision number
//	@param   "a"	-	float precsion number
//
//	@todo
//
// Description:
//  Unary hypen Just change sign
//
inline float_precision operator-( const float_precision& a )
	{
	float_precision b;

	b.precision( a.precision() );
	b = a;
	b.change_sign();

	return b;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator-
//	@return 	float_precision	*	return addition of lhs * rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Mul operator for float_precision * <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator*(float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs);

	if (lhs.precision() > c.precision())
		c.precision(lhs.precision());

	return c *= lhs;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator-
//	@return 	float_precision	*	return addition of lhs * rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Mul operator for float_precision * <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator*( const _Ty&lhs, const float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c *= rhs;
	}

//  @author Henrik Vestermark (hve@hvks.com)
//  @date  7/29/2014
//  @version 1.0
//	@brief 			operator*
//	@return 	float_precision	-	return addition of lhs * rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Add operator for int_precision * float_precision
//
inline float_precision operator*(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c *= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	* float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	first float precsion number
//	@param   "b"	-	second float precsion number to multiply
//
//	@todo
//
// Description:
//  Binary multiplying two float_precision numbers
//  Implenting using the essential *= operator
//
/*
inline float_precision operator*( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision );

   c = a;
   c *= b;

   return c;
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator/
//	@return 	float_precision	/	return addition of lhs / rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Div operator for float_precision / <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator/(float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d /= c;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  7/29/2014
//	@brief 			operator/
//	@return 	float_precision	-	return addition of lhs / rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  Div operator for  <any other type> / float_precision
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator/(const _Ty& lhs, const float_precision& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d /= c;
	}

//  @author Henrik Vestermark (hve@hvks.com)
//  @date  7/29/2014
//  @version 1.0
//	@brief 			operator/
//	@return 	float_precision	-	return division of lhs / rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  Div operator for int_precision / float_precision
//
inline float_precision operator/(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c /= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date  1/21/2005
//	@brief 	/ float precision numbers
//	@return 	the resulting float_precision number
//	@param   "a"	-	first float precsion number
//	@param   "b"	-	second float precsion number to divide
//
//	@todo
//
// Description:
//  Binary divide two float_precision numbers
//  Implenting using the essential /= operator
//
/*
inline float_precision operator/( const float_precision& a, const float_precision& b )
   {
   unsigned int precision;
   float_precision c;

   precision = a.precision();
   if( precision < b.precision() )
      precision = b.precision();

   c.precision( precision + 1 );

   c = a;
   c /= b;

   return c;
   }
*/

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  23-Mar-2021
//	@brief 			operator%
//	@return 	float_precision	%	return addition of lhs % rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  % operator for float_precision % <any other type>
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator%(float_precision& lhs, const _Ty& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d %= c;
	}

//	@author Henrik Vestermark (hve@hvks.com)
//	@date  23-Mar-2021
//	@brief 			operator%
//	@return 	float_precision	-	return addition of lhs % rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
//	@todo  Add to do things
//
// Description:
//  % operator for  <any other type> % float_precision
//  no const on the lhs parameter to prevent ambigous overload
//
template <class _Ty> inline float_precision operator%(const _Ty& lhs, const float_precision& rhs)
	{
	float_precision c(rhs), d(lhs);

	if (d.precision() < c.precision())
		d.precision(c.precision());

	return d %= c;
	}

//  @author Henrik Vestermark (hve@hvks.com)
//  @date  23-mar-2021
//  @version 1.0
//	@brief 			operator%
//	@return 	float_precision	-	return modulo of lhs % rhs
//	@param   "lhs"	-	First operand
//	@param   "rhs"	-	Second operand
//
// Description:
//  % operator for int_precision % float_precision
//
inline float_precision operator%(int_precision& lhs, float_precision& rhs)
	{
	float_precision c(lhs);

	if (rhs.precision() > c.precision())
		c.precision(rhs.precision());

	return c %= rhs;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		28/Sep/2021
//	@brief 		== float precision numberss
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  If both operands has the same sign, mantissa length and same exponent
//  and if the mantissa is identical then it's the same.
//  Precsion and rounding mode does not affect the comparison
//
inline bool operator==( const float_precision& a, const float_precision& b )
	{
	if (a.iszero() && b.iszero()) // if both zero then ignore sign. e.g. -0.0==+0.0
		return true;
	if( a.sign() != b.sign() )  
		return false; // sign differes then return false
	if (a.exponent() != b.exponent())
		return false; // exponent differes then return false
	if( _float_precision_compare(const_cast<float_precision&>(a).pointer(), const_cast<float_precision&>(b).pointer()) != 0 )
		return false;
	return true;
	} 


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		28/Sep/2021
//	@brief 		compare < float precision numbers
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  1) Test for both operand is zero and return false if condition is meet
//  2) If signs differs then return the boolean result based on that
//  3) Now if same sign and one operand is zero then return the boolean result
//  4) If same sign and not zero check the exponent
//  5) If same sign and same exponent then check the mantissa for boolean result
//  Precsion and rounding mode does not affect the comparison
//
inline bool operator<( const float_precision& a, const float_precision& b )
	{
	int sign1, sign2, cmp;
	bool zero1, zero2;

	zero1 = a.iszero(); 
	zero2 = b.iszero(); 
	if( zero1 == true && zero2 == true )  // Both zero
		return false;

	sign1 = a.sign();
	sign2 = b.sign();
	// Different signs
	if( sign1 < sign2 )
		return true;
	if( sign1 > sign2 )
		return false;

	// Now a &  b has the same sign
	if( zero1 == true )   // If a is zero and a & b has the same sign and b is not zero then a < b
		return true;
	if( zero2 == true )   // If b is zero and a & b has the same sign and a is not zero then a > b
		return false;

	// Same sign and not zero . Check exponent
	if( a.exponent() < b.exponent() )
		return sign1 > 0 ? true : false;
	if( a.exponent() > b.exponent() )
		return sign1 > 0 ? false: true;

	// Same sign & same exponent. Check mantissa
	cmp = _float_precision_compare(const_cast<float_precision&>(a).pointer(), const_cast<float_precision&>(b).pointer());
	if( cmp < 0 && sign1 == 1 )
		return true;
	else
		if( cmp > 0 && sign1 == -1 )
			return true;

	return false;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		!= float precision numberss
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  implemented negating the == comparison
//
inline bool operator!=( const float_precision& a, const float_precision& b )
	{
	return b == a ? false : true;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		> float precision numberss
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  Implemented using the equality a>b => b<a
//
inline bool operator>( const float_precision& a, const float_precision& b )
	{
	return b < a ? true : false;
	}


//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		<= float precision numberss
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  Implemented using the equality a<=b => not b<a
//
inline bool operator<=( const float_precision& a, const float_precision& b )
	{
	return b < a ? false : true;
	}


// Boolean >=
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		1/21/2005
//	@brief 		>= float precision numberss
//	@return 	the boolean result
//	@param		"a"	-	first float precsion number
//	@param		"b"	-	second float precsion number
//
//	@todo
//
// Description:
//  Implemented using the equality a>=b => not a<b
//
inline bool operator>=( const float_precision& a, const float_precision& b )
	{
	return a < b ? false: true;
	}


//////////////////////////////////////////////////////////
//
// END FLOAT PRECISION OPERATORS
//
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//
// BEGIN FLOAT PRECISION FUNCTIONS
//
//////////////////////////////////////////////////////////

// absolute()
//
//	@author Henrik Vestermark (hve@hvks.com)
//	@date		6/25/2014
//	@brief 		abs(a) float precision numbers
//	@return 	the absoluet value
//	@param		 "a"	-	first float precsion number
//
//	@todo
//
// Description:
//  return the absolute value of the float_precision number
//
inline float_precision fabs( const float_precision& a )
	{
	return abs(a);
	}

//////////////////////////////////////////////////////////
//
// END FLOAT PRECISION FUNCTIONS
//
//////////////////////////////////////////////////////////

#endif
