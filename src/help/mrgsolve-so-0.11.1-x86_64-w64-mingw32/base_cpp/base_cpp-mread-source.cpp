// Source MD5: 2ddbf0a4c987922000b80089bf7b4db4

#include "base_cpp-mread-header.h"

// PREAMBLE CODE BLOCK:
__BEGIN_config__
__END_config__

// MAIN CODE BLOCK:
__BEGIN_main__
Cl = th1 * exp(eta(1));
Vc = th2 * exp(eta(2));
Q  = th3 * exp(eta(3));
Vp = th4 * exp(eta(4));
Ka = th5 * exp(eta(5));
__END_main__

// DIFFERENTIAL EQUATIONS:
__BEGIN_ode__
dxdt_depot = - Ka * depot;
dxdt_cent  = Ka * depot - (Cl/Vc) * cent - (Q/Vc) * cent + (Q/Vp) * periph;
dxdt_periph = (Q/Vc) * cent - (Q/Vp) * periph;
__END_ode__

// TABLE CODE BLOCK:
__BEGIN_table__
Cc = (cent/Vc) * (1 + eps(1)) + eps(2);
i = 0; 
while (Cc < 0) {
        if (++i > 100) {
                mrg::report("Model gave up, negative values associated with Cc");
                break;
        }
        simeps();
        Cc = (cent/Vc) * (1 + eps(1)) + eps(2); 
}
_capture_[0] = Cc;
__END_table__
