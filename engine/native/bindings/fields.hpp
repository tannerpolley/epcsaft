// Field lists for bound names equal to their C++ names: EPCSAFT_READONLY(T, a, b) expands to
// .def_ro("a", &T::a).def_ro("b", &T::b), EPCSAFT_READWRITE to .def_rw, and EPCSAFT_VALUES(E, X)
// to .value("X", E::X). Renamed fields, docstrings and other policies stay explicit nanobind calls.
#pragma once

#define EPCSAFT_READONLY(T, ...) EPCSAFT_EACH(EPCSAFT_READONLY_ONE, T, __VA_ARGS__)
#define EPCSAFT_READWRITE(T, ...) EPCSAFT_EACH(EPCSAFT_READWRITE_ONE, T, __VA_ARGS__)
#define EPCSAFT_VALUES(E, ...) EPCSAFT_EACH(EPCSAFT_VALUE_ONE, E, __VA_ARGS__)
#define EPCSAFT_READONLY_ONE(T, name) .def_ro(#name, &T::name)
#define EPCSAFT_READWRITE_ONE(T, name) .def_rw(#name, &T::name)
#define EPCSAFT_VALUE_ONE(E, name) .value(#name, E::name)

// EPCSAFT_EACH(f, T, names...) is f(T, name) for each of 1 to 40 names; a longer list does not compile.
#define EPCSAFT_EACH(f, T, ...) EPCSAFT_CAT(EPCSAFT_EACH_, EPCSAFT_COUNT(__VA_ARGS__))(f, T, __VA_ARGS__)
#define EPCSAFT_CAT(a, b) EPCSAFT_CAT_(a, b)
#define EPCSAFT_CAT_(a, b) a##b
#define EPCSAFT_COUNT(...)                                                                                             \
    EPCSAFT_COUNT_(__VA_ARGS__, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20,    \
                   19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, )
#define EPCSAFT_COUNT_(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, \
                       _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39, _40,  \
                       n, ...)                                                                                         \
    n
#define EPCSAFT_EACH_1(f, T, a) f(T, a)
#define EPCSAFT_EACH_2(f, T, a, ...) f(T, a) EPCSAFT_EACH_1(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_3(f, T, a, ...) f(T, a) EPCSAFT_EACH_2(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_4(f, T, a, ...) f(T, a) EPCSAFT_EACH_3(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_5(f, T, a, ...) f(T, a) EPCSAFT_EACH_4(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_6(f, T, a, ...) f(T, a) EPCSAFT_EACH_5(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_7(f, T, a, ...) f(T, a) EPCSAFT_EACH_6(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_8(f, T, a, ...) f(T, a) EPCSAFT_EACH_7(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_9(f, T, a, ...) f(T, a) EPCSAFT_EACH_8(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_10(f, T, a, ...) f(T, a) EPCSAFT_EACH_9(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_11(f, T, a, ...) f(T, a) EPCSAFT_EACH_10(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_12(f, T, a, ...) f(T, a) EPCSAFT_EACH_11(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_13(f, T, a, ...) f(T, a) EPCSAFT_EACH_12(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_14(f, T, a, ...) f(T, a) EPCSAFT_EACH_13(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_15(f, T, a, ...) f(T, a) EPCSAFT_EACH_14(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_16(f, T, a, ...) f(T, a) EPCSAFT_EACH_15(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_17(f, T, a, ...) f(T, a) EPCSAFT_EACH_16(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_18(f, T, a, ...) f(T, a) EPCSAFT_EACH_17(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_19(f, T, a, ...) f(T, a) EPCSAFT_EACH_18(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_20(f, T, a, ...) f(T, a) EPCSAFT_EACH_19(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_21(f, T, a, ...) f(T, a) EPCSAFT_EACH_20(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_22(f, T, a, ...) f(T, a) EPCSAFT_EACH_21(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_23(f, T, a, ...) f(T, a) EPCSAFT_EACH_22(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_24(f, T, a, ...) f(T, a) EPCSAFT_EACH_23(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_25(f, T, a, ...) f(T, a) EPCSAFT_EACH_24(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_26(f, T, a, ...) f(T, a) EPCSAFT_EACH_25(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_27(f, T, a, ...) f(T, a) EPCSAFT_EACH_26(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_28(f, T, a, ...) f(T, a) EPCSAFT_EACH_27(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_29(f, T, a, ...) f(T, a) EPCSAFT_EACH_28(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_30(f, T, a, ...) f(T, a) EPCSAFT_EACH_29(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_31(f, T, a, ...) f(T, a) EPCSAFT_EACH_30(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_32(f, T, a, ...) f(T, a) EPCSAFT_EACH_31(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_33(f, T, a, ...) f(T, a) EPCSAFT_EACH_32(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_34(f, T, a, ...) f(T, a) EPCSAFT_EACH_33(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_35(f, T, a, ...) f(T, a) EPCSAFT_EACH_34(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_36(f, T, a, ...) f(T, a) EPCSAFT_EACH_35(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_37(f, T, a, ...) f(T, a) EPCSAFT_EACH_36(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_38(f, T, a, ...) f(T, a) EPCSAFT_EACH_37(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_39(f, T, a, ...) f(T, a) EPCSAFT_EACH_38(f, T, __VA_ARGS__)
#define EPCSAFT_EACH_40(f, T, a, ...) f(T, a) EPCSAFT_EACH_39(f, T, __VA_ARGS__)
