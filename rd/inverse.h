//
// Created by Wu Zhenyu on 04/01/2022.
//

#ifndef RD_INVERSE_H
#define RD_INVERSE_H

#include <cblas.h>
#include <lapacke.h>
lapack_int mat_inv(double*, unsigned, double, double, int, int);

lapack_int matFac(double*, unsigned);

#endif //RD_INVERSE_H
