/*
    Common preprocessor macros not tied to any specific component.
*/

#ifndef __NUMFORT_H__
#define __NUMFORT_H__

#define __IDENTITY(x) x
#define __APPEND(s1,s2) __IDENTITY(__IDENTITY(s1)_)s2
#define __APPEND2(s1,s2,s3) __IDENTITY(__APPEND(s1,s2)_)s3

#endif
