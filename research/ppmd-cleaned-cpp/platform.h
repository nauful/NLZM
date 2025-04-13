#pragma once

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;
typedef unsigned long long u64;

typedef signed char i8;
typedef signed short i16;
typedef signed int i32;
typedef signed long long i64;

typedef signed char s8;
typedef signed short s16;
typedef signed int s32;
typedef signed long long s64;

typedef float f32;
typedef double f64;

#define _CRT_SECURE_NO_WARNINGS

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <Windows.h>
#include <shellapi.h>

#include <math.h>
#include <stdio.h>

#define PLATFORM_OFFSETOF(type, member) ((size_t) &((type *)0)->member)
#define PLATFORM_COUNTOF(ARRAY) (sizeof(ARRAY) / sizeof((ARRAY)[0]))
#define PLATFORM_MIN(x, y) ((x) < (y) ? (x) : (y))
#define PLATFORM_MAX(x, y) ((x) < (y) ? (y) : (x))
#define PLATFORM_CLAMP(x, a, b) ((x) < (a) ? (a) : ((x) > (b) ? (b) : (x)))

#define PLATFORM_ASSERT(x) { if (!(x)) { __debugbreak(); } }

#ifdef _DEBUG
# define DEBUG
#endif

#ifdef DEBUG
# define PLATFORM_ASSERT_DEBUG(x) PLATFORM_ASSERT(x)
#else
# define PLATFORM_ASSERT_DEBUG(x) (x)
#endif

#define FPATH_MAX_LEN 1024

#ifndef M_PI
# define M_PI 3.14159265358979323846
#endif

#ifndef M_2PI
# define M_2PI (2.0 * M_PI)
#endif

#ifndef M_HPI
# define M_HPI (M_PI / 2.0)
#endif
