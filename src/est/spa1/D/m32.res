Thu Sep 30 02:54:51 CDT 2021
$PROB template control stream
;-----------------------------------------------------------------------
; Project: 	Investigating the contribution of residual unexplained
; 	   	variability in nonlinear mixed-effect approach
; Model: 	Two-compartment model with linear elimination
; Estim:	First-order conditional est. with interaction
; Author: 	Mutaz M. Jaber <jaber038@umn.edu>
; Date created: 9/7/2021
; Date modified: 9/7/2021
;-----------------------------------------------------------------------
$INPUT ID TIME DV AMT MDV EVID
$DATA ../../../../data/spa1/D/dat32.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER NSIG=2
$PK
ET1 = EXP(ETA(1)*THETA(6))
ET2 = EXP(ETA(2)*THETA(7))
ET3 = EXP(ETA(3)*THETA(8))
ET4 = EXP(ETA(4)*THETA(9))
ET5 = EXP(ETA(5)*THETA(10))

CL = 5.0 * THETA(1) * ET1
V2 = 35  * THETA(2) * ET2
Q  = 50  * THETA(3) * ET3
V3 = 50  * THETA(4) * ET4
KA = 0.7 * THETA(5) * ET5
SC = V2
$ERROR
CVERR = 0.05
W = THETA(11)*F*CVERR

Y 	= F + W*ERR(1)

$THETA
(0,1) ; CL
(0,1) ; V2
(0,1) ; Q
(0,1) ; V3
(0,1) ; KA
(0,1) ; IIVCL
(0,1) ; IIVV2
(0,1) ; IIVQ
(0,1) ; IIVV3
(0,1) ; IIVKA
(0,1) ; CVPropErr

$OMEGA  (0.09 FIX)x5
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       30 SEP 2021
Days until program expires : 199
1NONLINEAR MIXED EFFECTS MODEL PROGRAM (NONMEM) VERSION 7.5.0
 ORIGINALLY DEVELOPED BY STUART BEAL, LEWIS SHEINER, AND ALISON BOECKMANN
 CURRENT DEVELOPERS ARE ROBERT BAUER, ICON DEVELOPMENT SOLUTIONS,
 AND ALISON BOECKMANN. IMPLEMENTATION, EFFICIENCY, AND STANDARDIZATION
 PERFORMED BY NOUS INFOSYSTEMS.

 PROBLEM NO.:         1
 template control stream
0DATA CHECKOUT RUN:              NO
 DATA SET LOCATED ON UNIT NO.:    2
 THIS UNIT TO BE REWOUND:        NO
 NO. OF DATA RECS IN DATA SET:      600
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      500
 TOT. NO. OF INDIVIDUALS:      100
0LENGTH OF THETA:  11
0DEFAULT THETA BOUNDARY TEST OMITTED:    NO
0OMEGA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   5
0DEFAULT OMEGA BOUNDARY TEST OMITTED:    NO
0SIGMA HAS SIMPLE DIAGONAL FORM WITH DIMENSION:   1
0DEFAULT SIGMA BOUNDARY TEST OMITTED:    NO
0INITIAL ESTIMATE OF THETA:
 LOWER BOUND    INITIAL EST    UPPER BOUND
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
  0.0000E+00     0.1000E+01     0.1000E+07
0INITIAL ESTIMATE OF OMEGA:
 0.9000E-01
 0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
 0.0000E+00   0.0000E+00   0.0000E+00   0.0000E+00   0.9000E-01
0OMEGA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0INITIAL ESTIMATE OF SIGMA:
 0.1000E+01
0SIGMA CONSTRAINED TO BE THIS INITIAL ESTIMATE
0COVARIANCE STEP OMITTED:        NO
 EIGENVLS. PRINTED:              NO
 SPECIAL COMPUTATION:            NO
 COMPRESSED FORMAT:              NO
 GRADIENT METHOD USED:     NOSLOW
 SIGDIGITS ETAHAT (SIGLO):                  -1
 SIGDIGITS GRADIENTS (SIGL):                -1
 EXCLUDE COV FOR FOCE (NOFCOV):              NO
 Cholesky Transposition of R Matrix (CHOLROFF):0
 KNUTHSUMOFF:                                -1
 RESUME COV ANALYSIS (RESUME):               NO
 SIR SAMPLE SIZE (SIRSAMPLE):
 NON-LINEARLY TRANSFORM THETAS DURING COV (THBND): 1
 PRECONDTIONING CYCLES (PRECOND):        0
 PRECONDTIONING TYPES (PRECONDS):        TOS
 FORCED PRECONDTIONING CYCLES (PFCOND):0
 PRECONDTIONING TYPE (PRETYPE):        0
 FORCED POS. DEFINITE SETTING DURING PRECONDITIONING: (FPOSDEF):0
 SIMPLE POS. DEFINITE SETTING: (POSDEF):-1
1DOUBLE PRECISION PREDPP VERSION 7.5.0

 TWO COMPARTMENT MODEL WITH FIRST-ORDER ABSORPTION (ADVAN4)
0MAXIMUM NO. OF BASIC PK PARAMETERS:   5
0BASIC PK PARAMETERS (AFTER TRANSLATION):
   BASIC PK PARAMETER NO.  1: ELIMINATION RATE (K)
   BASIC PK PARAMETER NO.  2: CENTRAL-TO-PERIPH. RATE (K23)
   BASIC PK PARAMETER NO.  3: PERIPH.-TO-CENTRAL RATE (K32)
   BASIC PK PARAMETER NO.  5: ABSORPTION RATE (KA)
 TRANSLATOR WILL CONVERT PARAMETERS
 CL, V2, Q, V3 TO K, K23, K32 (TRANS4)
0COMPARTMENT ATTRIBUTES
 COMPT. NO.   FUNCTION   INITIAL    ON/OFF      DOSE      DEFAULT    DEFAULT
                         STATUS     ALLOWED    ALLOWED    FOR DOSE   FOR OBS.
    1         DEPOT        OFF        YES        YES        YES        NO
    2         CENTRAL      ON         NO         YES        NO         YES
    3         PERIPH.      ON         NO         YES        NO         NO
    4         OUTPUT       OFF        YES        NO         NO         NO
1
 ADDITIONAL PK PARAMETERS - ASSIGNMENT OF ROWS IN GG
 COMPT. NO.                             INDICES
              SCALE      BIOAVAIL.   ZERO-ORDER  ZERO-ORDER  ABSORB
                         FRACTION    RATE        DURATION    LAG
    1            *           *           *           *           *
    2            6           *           *           *           *
    3            *           *           *           *           *
    4            *           -           -           -           -
             - PARAMETER IS NOT ALLOWED FOR THIS MODEL
             * PARAMETER IS NOT SUPPLIED BY PK SUBROUTINE;
               WILL DEFAULT TO ONE IF APPLICABLE
0DATA ITEM INDICES USED BY PRED ARE:
   EVENT ID DATA ITEM IS DATA ITEM NO.:      6
   TIME DATA ITEM IS DATA ITEM NO.:          2
   DOSE AMOUNT DATA ITEM IS DATA ITEM NO.:   4

0PK SUBROUTINE CALLED WITH EVERY EVENT RECORD.
 PK SUBROUTINE NOT CALLED AT NONEVENT (ADDITIONAL OR LAGGED) DOSE TIMES.
0ERROR SUBROUTINE CALLED WITH EVERY EVENT RECORD.
1


 #TBLN:      1
 #METH: First Order Conditional Estimation with Interaction

 ESTIMATION STEP OMITTED:                 NO
 ANALYSIS TYPE:                           POPULATION
 NUMBER OF SADDLE POINT RESET ITERATIONS:      0
 GRADIENT METHOD USED:               NOSLOW
 CONDITIONAL ESTIMATES USED:              YES
 CENTERED ETA:                            NO
 EPS-ETA INTERACTION:                     YES
 LAPLACIAN OBJ. FUNC.:                    NO
 NO. OF FUNCT. EVALS. ALLOWED:            10000
 NO. OF SIG. FIGURES REQUIRED:            2
 INTERMEDIATE PRINTOUT:                   YES
 ESTIMATE OUTPUT TO MSF:                  NO
 ABORT WITH PRED EXIT CODE 1:             NO
 IND. OBJ. FUNC. VALUES SORTED:           NO
 NUMERICAL DERIVATIVE
       FILE REQUEST (NUMDER):               NONE
 MAP (ETAHAT) ESTIMATION METHOD (OPTMAP):   0
 ETA HESSIAN EVALUATION METHOD (ETADER):    0
 INITIAL ETA FOR MAP ESTIMATION (MCETA):    0
 SIGDIGITS FOR MAP ESTIMATION (SIGLO):      100
 GRADIENT SIGDIGITS OF
       FIXED EFFECTS PARAMETERS (SIGL):     100
 NOPRIOR SETTING (NOPRIOR):                 0
 NOCOV SETTING (NOCOV):                     OFF
 DERCONT SETTING (DERCONT):                 OFF
 FINAL ETA RE-EVALUATION (FNLETA):          1
 EXCLUDE NON-INFLUENTIAL (NON-INFL.) ETAS
       IN SHRINKAGE (ETASTYPE):             NO
 NON-INFL. ETA CORRECTION (NONINFETA):      0
 RAW OUTPUT FILE (FILE): m32.ext
 EXCLUDE TITLE (NOTITLE):                   NO
 EXCLUDE COLUMN LABELS (NOLABEL):           NO
 FORMAT FOR ADDITIONAL FILES (FORMAT):      S1PE12.5
 PARAMETER ORDER FOR OUTPUTS (ORDER):       TSOL
 KNUTHSUMOFF:                               0
 INCLUDE LNTWOPI:                           NO
 INCLUDE CONSTANT TERM TO PRIOR (PRIORC):   NO
 INCLUDE CONSTANT TERM TO OMEGA (ETA) (OLNTWOPI):NO
 ADDITIONAL CONVERGENCE TEST (CTYPE=4)?:    NO
 EM OR BAYESIAN METHOD USED:                 NONE


 THE FOLLOWING LABELS ARE EQUIVALENT
 PRED=PREDI
 RES=RESI
 WRES=WRESI
 IWRS=IWRESI
 IPRD=IPREDI
 IRS=IRESI

 MONITORING OF SEARCH:


0ITERATION NO.:    0    OBJECTIVE VALUE:   16719.8212656428        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2169E+02  1.9720E+02 -7.3488E+01  1.5059E+01  1.8696E+02 -1.2979E+03 -6.6584E+02 -4.6298E+01 -1.3022E+03 -2.6482E+02
            -3.3743E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -723.127634019630        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2537E+00  9.8060E-01  1.0738E+00  1.8452E+00  1.2505E+00  2.2167E+00  1.3416E+00  9.1601E-01  1.9679E+00  9.9991E-01
             1.3592E+01
 PARAMETER:  3.2607E-01  8.0408E-02  1.7122E-01  7.1257E-01  3.2351E-01  8.9601E-01  3.9389E-01  1.2268E-02  7.7699E-01  9.9908E-02
             2.7095E+00
 GRADIENT:  -5.9987E+01  1.3250E+01 -1.8122E+01  2.2459E+01 -2.2062E+00  6.1112E+01  1.7484E+00  5.9895E+00  1.0205E+01  2.4616E+00
             2.7025E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -750.823298426501        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3130E+00  8.8943E-01  1.8066E+00  2.0874E+00  4.1087E+00  1.9107E+00  4.2904E+00  4.1170E-01  2.0396E+00  3.9794E+00
             1.2613E+01
 PARAMETER:  3.7231E-01 -1.7170E-02  6.9144E-01  8.3591E-01  1.5131E+00  7.4749E-01  1.5564E+00 -7.8746E-01  8.1275E-01  1.4811E+00
             2.6347E+00
 GRADIENT:  -2.9496E+01  1.8969E+01  5.3638E+00  4.5871E+01 -5.4229E+00  1.3876E+01  1.1378E+01 -5.7329E-02  3.0641E+01  3.0878E-01
             2.2890E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -758.105999281857        NO. OF FUNC. EVALS.: 124
 CUMULATIVE NO. OF FUNC. EVALS.:      279             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3191E+00  8.6190E-01  1.7746E+00  2.0256E+00  4.6177E+00  1.8967E+00  3.6461E+00  3.4012E-01  2.0282E+00  6.5696E+00
             1.2338E+01
 PARAMETER:  3.7692E-01 -4.8616E-02  6.7360E-01  8.0587E-01  1.6299E+00  7.4010E-01  1.3937E+00 -9.7847E-01  8.0713E-01  1.9824E+00
             2.6127E+00
 GRADIENT:  -1.9733E+01  1.6507E+01  5.8051E+00  4.4296E+01 -4.3095E+00  1.4685E+01  9.5039E+00 -7.3948E-02  2.8418E+01  2.3368E+00
             2.1587E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -830.974814827278        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  8.3632E-01  4.9194E-02  2.1206E-01  1.1282E+00  4.8647E+00  1.6812E+00  3.0081E-01  1.3780E+00  1.2786E+00  6.2695E+00
             8.1864E+00
 PARAMETER: -7.8749E-02 -2.9120E+00 -1.4509E+00  2.2058E-01  1.6820E+00  6.1949E-01 -1.1013E+00  4.2063E-01  3.4580E-01  1.9357E+00
             2.2025E+00
 GRADIENT:  -1.5669E+01  9.3009E+00 -4.0374E+01  1.0392E+02 -7.3491E+01  3.6676E+01  1.4141E-03 -1.4254E+00  4.7652E+01  4.3265E+01
            -9.5232E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -832.505914324467        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  7.3229E-01  2.4432E-02  1.3996E-01  9.2947E-01  5.0267E+00  1.5888E+00  2.5419E-01  1.8339E+00  1.1213E+00  6.0885E+00
             8.0379E+00
 PARAMETER: -2.1158E-01 -3.6119E+00 -1.8664E+00  2.6863E-02  1.7148E+00  5.6298E-01 -1.2697E+00  7.0642E-01  2.1453E-01  1.9064E+00
             2.1842E+00
 GRADIENT:  -5.2402E+01  5.8524E+00 -6.9653E+01  1.1115E+02 -8.8947E+01  2.5973E+01  2.0483E-04  1.7852E+01  4.2387E+01  4.2979E+01
            -1.2871E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -910.474838606438        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  5.4727E-01  1.0000E-02  4.5308E-02  4.1269E-01  8.2087E+00  1.4920E+00  3.2056E-01  1.0377E+00  3.4465E-01  4.1293E+00
             8.8592E+00
 PARAMETER: -5.0281E-01 -8.3553E+00 -2.9943E+00 -7.8505E-01  2.2052E+00  5.0015E-01 -1.0377E+00  1.3701E-01 -9.6522E-01  1.5181E+00
             2.2815E+00
 GRADIENT:  -2.7306E+01  0.0000E+00  8.7669E+01 -1.1490E+02 -8.3358E+00  3.2296E+01  1.5280E-04 -8.3376E+00 -4.7085E-03  5.4517E+00
             1.4248E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -916.714859749734        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  5.6227E-01  1.0000E-02  4.2733E-02  4.2031E-01  8.4869E+00  1.3585E+00  3.0597E-01  1.0868E+00  3.7782E-01  3.9749E+00
             8.8807E+00
 PARAMETER: -4.7577E-01 -8.4297E+00 -3.0528E+00 -7.6677E-01  2.2385E+00  4.0637E-01 -1.0843E+00  1.8319E-01 -8.7333E-01  1.4800E+00
             2.2839E+00
 GRADIENT:  -1.5092E+00  0.0000E+00 -1.0345E+00 -3.6576E+00 -5.2179E+00  3.7228E+00  9.6577E-05  1.8226E+00  1.2932E+00  5.0900E+00
             3.6909E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -921.100279183683        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  5.7903E-01  1.0000E-02  4.6404E-02  4.4508E-01  9.1806E+00  1.3440E+00  3.0363E-01  1.0792E+00  3.7431E-01  3.6237E+00
             8.8855E+00
 PARAMETER: -4.4640E-01 -8.5587E+00 -2.9704E+00 -7.0950E-01  2.3171E+00  3.9565E-01 -1.0919E+00  1.7625E-01 -8.8266E-01  1.3875E+00
             2.2844E+00
 GRADIENT:   4.9296E-02  0.0000E+00  1.1978E+00 -2.0758E+00 -1.0317E+00 -8.6527E-01  2.4625E-05 -1.5336E+00  9.4917E-01  4.4908E-02
             2.5093E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -921.336674287560        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  5.8423E-01  1.0000E-02  4.7993E-02  4.5524E-01  1.2571E+01  1.3480E+00  2.9045E-01  1.1562E+00  2.8857E-01  2.6220E+00
             8.8454E+00
 PARAMETER: -4.3746E-01 -9.4350E+00 -2.9367E+00 -6.8692E-01  2.6314E+00  3.9864E-01 -1.1363E+00  2.4515E-01 -1.1428E+00  1.0639E+00
             2.2799E+00
 GRADIENT:   1.4800E-01  0.0000E+00 -1.7343E+00  1.3226E+00  6.0631E-02 -2.2691E-01  1.6624E-05 -1.2564E-01  5.1330E-01 -4.2607E-03
            -1.3924E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -921.337111453549        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1354
 NPARAMETR:  5.8998E-01  1.0000E-02  4.9521E-02  4.6467E-01  1.4416E+01  1.3494E+00  2.8660E-01  1.1836E+00  2.5887E-01  2.2740E+00
             8.8493E+00
 PARAMETER: -4.2766E-01 -9.7929E+00 -2.9054E+00 -6.6644E-01  2.7683E+00  3.9967E-01 -1.1497E+00  2.6859E-01 -1.2514E+00  9.2153E-01
             2.2803E+00
 GRADIENT:   1.1990E-01  0.0000E+00 -1.0360E+00  6.9742E-01  5.9742E-02 -1.5294E-01  1.5630E-05  4.1697E-02  4.1138E-01 -2.6314E-03
            -8.4938E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -921.368714829981        NO. OF FUNC. EVALS.: 108
 CUMULATIVE NO. OF FUNC. EVALS.:     1462
 NPARAMETR:  5.8694E-01  1.0000E-02  4.9093E-02  4.6285E-01  1.3775E+01  1.3470E+00  2.8295E-01  1.1852E+00  2.4228E-01  6.8643E+00
             8.8290E+00
 PARAMETER: -4.3283E-01 -9.7929E+00 -2.9140E+00 -6.7035E-01  2.7229E+00  3.9790E-01 -1.1625E+00  2.6991E-01 -1.3177E+00  2.0263E+00
             2.2780E+00
 GRADIENT:   3.7079E+01  0.0000E+00  5.3102E+01  2.6257E+01  6.2863E-02  5.1991E+00  3.4377E-05  3.9707E-01  5.3962E-01  3.7847E-02
             1.9652E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -921.436304316580        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1622
 NPARAMETR:  5.9348E-01  1.0000E-02  5.0002E-02  4.6632E-01  1.2012E+01  1.3541E+00  2.6729E-01  1.2213E+00  1.5177E-01  5.9048E+00
             8.8870E+00
 PARAMETER: -4.2175E-01 -9.7929E+00 -2.8957E+00 -6.6288E-01  2.5859E+00  4.0312E-01 -1.2194E+00  2.9988E-01 -1.7854E+00  1.8758E+00
             2.2846E+00
 GRADIENT:   2.6795E+00  0.0000E+00  2.0339E+00 -5.8925E+00 -2.8296E-02  8.7463E-01  1.0163E-05  5.3900E-01  1.0596E-01  6.2868E-02
             3.2378E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -921.487957619129        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:     1788
 NPARAMETR:  5.8784E-01  1.0000E-02  4.8985E-02  4.6129E-01  1.0892E+01  1.3502E+00  2.3707E-01  1.2159E+00  9.2405E-02  4.9522E+00
             8.8440E+00
 PARAMETER: -4.3129E-01 -9.7929E+00 -2.9162E+00 -6.7373E-01  2.4881E+00  4.0026E-01 -1.3394E+00  2.9552E-01 -2.2816E+00  1.6998E+00
             2.2797E+00
 GRADIENT:  -8.9453E-01  0.0000E+00 -2.6868E+00  2.7545E+00 -8.6741E-03 -3.9376E-01  7.6676E-06 -7.5936E-01  1.1533E-02 -9.3410E-03
            -2.3253E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -921.498279066703        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1970             RESET HESSIAN, TYPE I
 NPARAMETR:  5.8915E-01  1.0000E-02  4.9114E-02  4.6207E-01  1.1092E+01  1.3516E+00  2.0368E-01  1.2258E+00  3.5882E-02  5.1211E+00
             8.8600E+00
 PARAMETER: -4.2908E-01 -9.7929E+00 -2.9136E+00 -6.7203E-01  2.5063E+00  4.0132E-01 -1.4912E+00  3.0357E-01 -3.2275E+00  1.7334E+00
             2.2815E+00
 GRADIENT:   3.9699E+01  0.0000E+00  5.3956E+01  2.1952E+01  5.6246E-02  6.0597E+00  1.8806E-05  5.6726E-01  1.6216E-02  4.7491E-03
             2.2171E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -921.500490426837        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2147
 NPARAMETR:  5.8852E-01  1.0000E-02  4.9132E-02  4.6118E-01  1.1085E+01  1.3516E+00  1.8896E-01  1.2279E+00  1.0000E-02  5.2557E+00
             8.8653E+00
 PARAMETER: -4.2896E-01 -9.7929E+00 -2.9136E+00 -6.7269E-01  2.5112E+00  4.0152E-01 -1.5583E+00  3.0609E-01 -5.0074E+00  1.7419E+00
             2.2817E+00
 GRADIENT:   4.0816E-01  0.0000E+00 -1.1313E-01  8.1951E-01  4.2952E-03  2.4440E-02  6.3722E-06  3.6507E-02  0.0000E+00 -3.3201E-03
            -1.5266E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2147
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.2088E-03 -1.4717E-05  7.4213E-03 -4.0163E-04 -1.2038E-03
 SE:             2.8737E-02  2.0790E-05  2.3892E-02  3.0059E-04  8.4304E-04
 N:                     100         100         100         100         100

 P VAL.:         8.0193E-01  4.7902E-01  7.5610E-01  1.8151E-01  1.5332E-01

 ETASHRINKSD(%)  3.7259E+00  9.9930E+01  1.9957E+01  9.8993E+01  9.7176E+01
 ETASHRINKVR(%)  7.3130E+00  1.0000E+02  3.5932E+01  9.9990E+01  9.9920E+01
 EBVSHRINKSD(%)  4.0830E+00  9.9925E+01  1.9975E+01  9.9017E+01  9.7833E+01
 EBVSHRINKVR(%)  7.9993E+00  1.0000E+02  3.5960E+01  9.9990E+01  9.9953E+01
 RELATIVEINF(%)  3.4720E+00  2.0120E-06  5.5762E-01  5.6368E-05  5.0227E-03
 EPSSHRINKSD(%)  1.0626E+01
 EPSSHRINKVR(%)  2.0123E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -921.50049042683725     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2.5619572221645512     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -921.500       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.89E-01  1.00E-02  4.91E-02  4.62E-01  1.11E+01  1.35E+00  1.90E-01  1.23E+00  1.00E-02  5.16E+00  8.86E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        9.00E-02
 
 ETA2
+        0.00E+00  9.00E-02
 
 ETA3
+        0.00E+00  0.00E+00  9.00E-02
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  9.00E-02
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  9.00E-02
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+        1.00E+00
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+        3.00E-01
 
 ETA2
+        0.00E+00  3.00E-01
 
 ETA3
+        0.00E+00  0.00E+00  3.00E-01
 
 ETA4
+        0.00E+00  0.00E+00  0.00E+00  3.00E-01
 
 ETA5
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  3.00E-01
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+        1.00E+00
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     R MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.58E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.81E+03  0.00E+00  1.86E+05
 
 TH 4
+       -7.13E+02  0.00E+00 -2.73E+04  4.78E+03
 
 TH 5
+        3.10E-01  0.00E+00 -2.79E+00  2.60E-01  1.85E+00
 
 TH 6
+        4.90E+00  0.00E+00  1.14E+02 -3.25E+01 -8.39E-04  9.07E+01
 
 TH 7
+       -6.50E-03  0.00E+00 -1.53E-03  2.91E-04 -4.30E-05  7.23E-04  2.11E-03
 
 TH 8
+        1.57E+00  0.00E+00 -1.78E+02 -3.76E+01  3.54E-03  1.63E+00 -2.27E-04  4.75E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        2.10E-02  0.00E+00 -1.06E-01  3.01E-02 -5.87E+00  1.37E-02 -1.14E-04 -1.49E-02  0.00E+00  1.74E-02
 
 TH11
+       -1.97E+01  0.00E+00  1.87E+02 -2.12E+01  2.57E+00  1.15E+00  4.83E-05  3.49E+00  0.00E+00 -4.40E-03  6.53E+00
 
 OM11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM12
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         .........
 
 OM13
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... .........
 
 OM14
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... .........
 
 OM15
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... .........
 
 OM22
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... .........
 
 OM23
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... .........
 
 OM24
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... .........
 
1

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 OM25
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... .........
 
 OM33
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM34
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM35
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM44
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
 
 OM45
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        .........
 
 OM55
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... .........
 
 SG11
+       ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
         ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... ......... .........
        ......... ......... .........
 
 Elapsed finaloutput time in seconds:     0.01
 #CPUT: Total CPU Time in Seconds,       50.487
Stop Time:
Thu Sep 30 02:55:43 CDT 2021
