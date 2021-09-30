[ prob ]
Project :Investigating the contribution of residual unexplained variability components in NLME approach
Model: Base model
- Two compartment model with linear elimination
- Individual parameters are log-normally distributed
- 5% proportional error model

[ set ] delta=0.1, end=50

[ cmt ] depot cent periph
[ theta ] @name th 
5 35 50 50 0.7
[ omega ]  0.09 0.09 0.09 0.09 0.09
[ sigma ]  0.0025 0.0

[ global ]
#include "include/math.h"
#define eta(n) ETA(n)
#define eps(n) EPS(n)

double Cl, Vc, Q, Vp, Ka;
double Cc; 
[ main ]

Cl = th1 * exp(eta(1));
Vc = th2 * exp(eta(2));
Q  = th3 * exp(eta(3));
Vp = th4 * exp(eta(4));
Ka = th5 * exp(eta(5));

[ ode ]

dxdt_depot = - Ka * depot;
dxdt_cent  =   Ka * depot - (Cl/Vc) * cent - (Q/Vc) * cent + (Q/Vp) * periph;
dxdt_periph = (Q/Vc) * cent - (Q/Vp) * periph;

[ error ]
Cc = (cent/Vc) * (1 + eps(1)) + eps(2);

int i = 0; 
while (Cc < 0) {
        if (++i > 200) {
                mrg::report("Model gave up, There are negative values associated with Cc");
                break;
        }
        simeps();
        Cc = (cent/Vc) * (1 + eps(1)) + eps(2); 
}

[ capture ] Cc 
