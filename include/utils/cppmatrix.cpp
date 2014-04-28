#include "cppmatrix.h"

#ifndef WIN32
	template<>
	mxClassID matrix<double>::typeToID()
	{
		return mxDOUBLE_CLASS;
	}
	template<>
	mxClassID matrix<float>::typeToID()
	{
		return mxSINGLE_CLASS;
	}
	template<>
	mxClassID matrix<unsigned char>::typeToID()
	{
		return mxUINT8_CLASS;
	}
	template<>
	mxClassID matrix<signed char>::typeToID()
	{
		return mxINT8_CLASS;
	}
	template<>
	mxClassID matrix<unsigned int>::typeToID()
	{
		ASSERT(sizeof(unsigned int)==4);
		return mxUINT32_CLASS;
	}
	template<>
	mxClassID matrix<signed int>::typeToID()
	{
		ASSERT(sizeof(signed int)==4);
		return mxINT32_CLASS;
	}
	template<>
	mxClassID matrix<unsigned short>::typeToID()
	{
		ASSERT(sizeof(unsigned short)==4);
		return mxUINT16_CLASS;
	}
	template<>
	mxClassID matrix<signed short>::typeToID()
	{
		ASSERT(sizeof(signed short)==4);
		return mxINT16_CLASS;
	}
    template<>
	mxClassID matrix<int64_t>::typeToID()
	{
		return mxINT64_CLASS;
	}
#endif

