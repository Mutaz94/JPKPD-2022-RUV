// Source MD5: 2ddbf0a4c987922000b80089bf7b4db4


// FIXED:
// No fixed parameters.

// INCLUDES:


// NAMESPACES:

// MODEL HEADER FILES:
#include "mrgsolv.h"
#include "modelheader.h"

//INCLUDE databox functions:
#include "databox_cpp.h"

//USING plugins

// GLOBAL CODE BLOCK:
// GLOBAL VARS FROM BLOCKS & TYPEDEFS:
typedef double capture;
namespace {
  int i;
}

// GLOBAL START USER CODE:
#define eta(n) ETA(n)
#define eps(n) EPS(n)
double Cl, Vc, Q, Vp, Ka;
double Cc; 

// DEFS:
#define __INITFUN___ _model_base__cpp_main__
#define __ODEFUN___ _model_base__cpp_ode__
#define __TABLECODE___ _model_base__cpp_table__
#define __CONFIGFUN___ _model_base__cpp_config__
#define __REGISTERFUN___ R_init_base_cpp
#define _nEQ 3
#define _nPAR 5
#define depot_0 _A_0_[0]
#define cent_0 _A_0_[1]
#define periph_0 _A_0_[2]
#define depot _A_[0]
#define cent _A_[1]
#define periph _A_[2]
#define dxdt_depot _DADT_[0]
#define dxdt_cent _DADT_[1]
#define dxdt_periph _DADT_[2]
#define th1 _THETA_[0]
#define th2 _THETA_[1]
#define th3 _THETA_[2]
#define th4 _THETA_[3]
#define th5 _THETA_[4]
