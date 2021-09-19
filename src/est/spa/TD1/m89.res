Sat Sep 18 14:20:18 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat89.csv ignore=@
$SUBR ADVAN4 TRANS4
$EST MET=1 NOABORT MAX=10000 PRINT=5 INTER
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

$OMEGA  0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX 0.09 FIX
$SIGMA  1 FIX ;        [P]
$COVARIANCE UNCOND

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       18 SEP 2021
Days until program expires : 211
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
 NO. OF DATA RECS IN DATA SET:      500
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      400
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
 NO. OF SIG. FIGURES REQUIRED:            3
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
 RAW OUTPUT FILE (FILE): m89.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1643.12974510537        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7877E+01 -7.3715E+01 -2.1104E+01 -6.2235E+01  6.8356E+01  3.7816E+01 -5.7097E+00 -1.7414E+00  9.8617E+00 -1.5771E+01
            -2.1592E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.92741656122        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9432E-01  1.0457E+00  9.7510E-01  1.0284E+00  9.6278E-01  8.9272E-01  1.0181E+00  1.0154E+00  9.3785E-01  1.0377E+00
             1.0505E+00
 PARAMETER:  9.4306E-02  1.4469E-01  7.4784E-02  1.2805E-01  6.2068E-02 -1.3480E-02  1.1792E-01  1.1533E-01  3.5838E-02  1.3703E-01
             1.4923E-01
 GRADIENT:   2.7377E+01  2.4875E+00 -5.9368E+00  1.4572E+01  1.2121E+01 -3.2068E+00 -2.9117E+00  3.9671E-01  1.1862E+00  1.6657E-02
            -2.1332E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1650.41708960328        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.9274E-01  9.3906E-01  8.9484E-01  1.0960E+00  8.6504E-01  9.0785E-01  1.2353E+00  8.7551E-01  8.1957E-01  9.2723E-01
             1.0659E+00
 PARAMETER:  9.2713E-02  3.7122E-02 -1.1107E-02  1.9168E-01 -4.4985E-02  3.3207E-03  3.1132E-01 -3.2946E-02 -9.8971E-02  2.4441E-02
             1.6384E-01
 GRADIENT:   2.0675E+01  1.6639E+01 -2.6327E+00  3.1431E+01  3.8803E+00  3.1287E+00 -5.3748E-01  1.3814E+00 -6.8999E+00  1.8240E-01
             4.9071E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1650.94150638734        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.8624E-01  9.2378E-01  6.4201E-01  1.0642E+00  7.1847E-01  9.0339E-01  1.2726E+00  5.1033E-01  8.2414E-01  7.7505E-01
             1.0488E+00
 PARAMETER:  8.6148E-02  2.0715E-02 -3.4315E-01  1.6219E-01 -2.3063E-01 -1.5990E-03  3.4103E-01 -5.7269E-01 -9.3412E-02 -1.5483E-01
             1.4768E-01
 GRADIENT:  -3.0099E+00  4.5476E+00 -1.0989E+01  1.4243E+01  1.0433E+01 -7.8913E-02  2.1936E+00  2.5011E+00 -2.4348E+00  3.7302E+00
            -1.3793E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1651.13026763453        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      327
 NPARAMETR:  9.8566E-01  8.1823E-01  5.1361E-01  1.0914E+00  5.9277E-01  9.0266E-01  1.3897E+00  3.0129E-01  7.8824E-01  6.0982E-01
             1.0534E+00
 PARAMETER:  8.5555E-02 -1.0062E-01 -5.6629E-01  1.8748E-01 -4.2295E-01 -2.4058E-03  4.2908E-01 -1.0997E+00 -1.3796E-01 -3.9459E-01
             1.5199E-01
 GRADIENT:  -4.0923E+01  3.3218E+00 -1.0717E+01  1.0492E+00  3.0360E-01 -3.3352E+00 -8.2640E-02  1.4830E+00 -1.5907E+00  2.8098E+00
            -2.6517E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.59373361344        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      505
 NPARAMETR:  1.0011E+00  6.9010E-01  6.0888E-01  1.1822E+00  6.1062E-01  9.1041E-01  1.6153E+00  2.7681E-01  7.6148E-01  6.9499E-01
             1.0642E+00
 PARAMETER:  1.0107E-01 -2.7092E-01 -3.9613E-01  2.6739E-01 -3.9328E-01  6.1354E-03  5.7950E-01 -1.1844E+00 -1.7249E-01 -2.6386E-01
             1.6227E-01
 GRADIENT:   4.1634E+00  1.7681E+00  2.2027E-01  4.7926E-01 -9.9083E-01  8.1310E-01  5.2758E-01  4.6136E-02 -6.7879E-01  5.1775E-01
             1.4675E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.61645041723        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      683
 NPARAMETR:  9.9950E-01  6.6594E-01  6.0530E-01  1.1935E+00  6.0125E-01  9.0872E-01  1.6548E+00  2.6298E-01  7.5738E-01  6.8861E-01
             1.0619E+00
 PARAMETER:  9.9503E-02 -3.0656E-01 -4.0203E-01  2.7691E-01 -4.0875E-01  4.2787E-03  6.0365E-01 -1.2357E+00 -1.7789E-01 -2.7308E-01
             1.6006E-01
 GRADIENT:   3.3964E-01  2.8864E-01 -1.0904E-01  9.6315E-01  2.5158E-01  6.7698E-02  5.8323E-02 -6.1775E-02 -1.8569E-01 -6.3869E-02
             1.4962E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1652.66779443875        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      866
 NPARAMETR:  9.9742E-01  5.9769E-01  6.7908E-01  1.2456E+00  6.2236E-01  9.0731E-01  1.7941E+00  4.5591E-01  7.4029E-01  7.2821E-01
             1.0583E+00
 PARAMETER:  9.7421E-02 -4.1468E-01 -2.8702E-01  3.1958E-01 -3.7424E-01  2.7326E-03  6.8449E-01 -6.8546E-01 -2.0072E-01 -2.1717E-01
             1.5671E-01
 GRADIENT:  -6.6625E-01  7.4029E-01  9.6758E-01 -1.3629E+00 -3.5833E+00 -1.3145E-01 -3.6269E-02  1.9675E-01  2.8445E-01  7.0053E-01
            -1.6828E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1652.85610457750        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1046
 NPARAMETR:  9.9250E-01  4.8149E-01  9.0475E-01  1.3456E+00  7.1183E-01  9.0395E-01  2.1055E+00  7.7218E-01  7.1475E-01  8.2515E-01
             1.0522E+00
 PARAMETER:  9.2468E-02 -6.3086E-01 -9.3224E-05  3.9681E-01 -2.3992E-01 -9.8185E-04  8.4454E-01 -1.5853E-01 -2.3583E-01 -9.2188E-02
             1.5091E-01
 GRADIENT:  -1.4981E+00  1.7527E+00  6.1783E-01  4.3679E-01  4.7603E-01 -3.6714E-01  8.1654E-01  7.2244E-02  1.2462E+00 -9.3228E-01
            -7.7103E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1653.07133014516        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  9.8954E-01  3.4576E-01  1.0600E+00  1.4387E+00  7.4751E-01  9.0255E-01  2.5666E+00  9.4528E-01  6.8638E-01  8.8201E-01
             1.0515E+00
 PARAMETER:  8.9486E-02 -9.6201E-01  1.5822E-01  4.6371E-01 -1.9101E-01 -2.5316E-03  1.0426E+00  4.3729E-02 -2.7633E-01 -2.5553E-02
             1.5025E-01
 GRADIENT:  -1.1808E-01  7.1238E-01  2.3003E+00  1.3013E+00 -2.5915E+00 -7.2567E-02 -1.6176E-01 -2.7367E-01  5.0731E-01 -4.4176E-01
            -4.5208E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1653.09092264966        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1404
 NPARAMETR:  9.8893E-01  3.1317E-01  1.0763E+00  1.4578E+00  7.4794E-01  9.0238E-01  2.7204E+00  9.6685E-01  6.8001E-01  8.8894E-01
             1.0524E+00
 PARAMETER:  8.8867E-02 -1.0610E+00  1.7354E-01  4.7691E-01 -1.9043E-01 -2.7203E-03  1.1008E+00  6.6284E-02 -2.8565E-01 -1.7721E-02
             1.5105E-01
 GRADIENT:   1.4877E-02 -1.0150E-02 -3.3046E-02 -5.9137E-02  2.6717E-02  2.3925E-03 -8.2299E-04  1.4382E-02  6.6746E-03  8.2307E-03
             2.6807E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1653.09092507321        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1584             RESET HESSIAN, TYPE I
 NPARAMETR:  9.8893E-01  3.1320E-01  1.0763E+00  1.4578E+00  7.4794E-01  9.0238E-01  2.7205E+00  9.6666E-01  6.7999E-01  8.8890E-01
             1.0524E+00
 PARAMETER:  8.8866E-02 -1.0609E+00  1.7355E-01  4.7690E-01 -1.9043E-01 -2.7231E-03  1.1008E+00  6.6091E-02 -2.8568E-01 -1.7766E-02
             1.5105E-01
 GRADIENT:   3.3801E+01  3.9885E+00  4.9566E-01  5.3704E+01  1.5790E+00  2.2868E+00  3.6217E+00  2.7314E-02  1.5273E+00  6.4344E-02
             1.1996E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1653.09092507321        NO. OF FUNC. EVALS.:  59
 CUMULATIVE NO. OF FUNC. EVALS.:     1643
 NPARAMETR:  9.8893E-01  3.1320E-01  1.0763E+00  1.4578E+00  7.4794E-01  9.0238E-01  2.7205E+00  9.6666E-01  6.7999E-01  8.8890E-01
             1.0524E+00
 PARAMETER:  8.8866E-02 -1.0609E+00  1.7355E-01  4.7690E-01 -1.9043E-01 -2.7231E-03  1.1008E+00  6.6091E-02 -2.8568E-01 -1.7766E-02
             1.5105E-01
 GRADIENT:  -3.8983E-04 -5.0741E-03  2.1886E-03  1.2574E-02  8.5649E-03  1.5642E-03  1.1014E-03 -1.4248E-03 -2.6099E-04 -2.1008E-03
            -1.7849E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1643
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0088E-03  3.1068E-02 -3.0781E-02 -2.6568E-02 -1.8397E-02
 SE:             2.9812E-02  1.7548E-02  1.6487E-02  2.4488E-02  2.0491E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7300E-01  7.6650E-02  6.1894E-02  2.7795E-01  3.6928E-01

 ETASHRINKSD(%)  1.2755E-01  4.1213E+01  4.4768E+01  1.7963E+01  3.1353E+01
 ETASHRINKVR(%)  2.5493E-01  6.5441E+01  6.9494E+01  3.2700E+01  5.2875E+01
 EBVSHRINKSD(%)  5.7258E-01  4.7201E+01  4.6678E+01  1.5153E+01  2.7623E+01
 EBVSHRINKVR(%)  1.1419E+00  7.2122E+01  7.1568E+01  2.8009E+01  4.7616E+01
 RELATIVEINF(%)  9.7984E+01  4.3440E+00  4.0766E+00  1.2841E+01  6.9047E+00
 EPSSHRINKSD(%)  4.4735E+01
 EPSSHRINKVR(%)  6.9458E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1653.0909250732063     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -917.94009850946816     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.28
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1653.091       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.89E-01  3.13E-01  1.08E+00  1.46E+00  7.48E-01  9.02E-01  2.72E+00  9.67E-01  6.80E-01  8.89E-01  1.05E+00
 


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
+        1.38E+03
 
 TH 2
+       -2.63E+01  5.00E+02
 
 TH 3
+        8.80E+00  8.44E+01  3.10E+02
 
 TH 4
+       -8.40E+00  4.23E+02 -9.32E+01  8.41E+02
 
 TH 5
+        1.41E+00 -2.94E+02 -5.55E+02  4.10E+01  1.34E+03
 
 TH 6
+       -6.58E+00 -3.78E+00  8.13E-01 -2.53E+00 -3.50E+00  2.40E+02
 
 TH 7
+        9.84E-01  3.45E+01 -9.84E-01 -6.19E+00  1.96E+00 -3.23E-01  7.81E+00
 
 TH 8
+        2.58E+00  6.34E-01 -5.23E+01 -8.57E+00  7.09E+00 -3.02E+00  4.13E-01  4.15E+01
 
 TH 9
+        5.30E+00 -2.95E+01  3.49E+00 -1.29E+01  9.84E+00 -6.09E-02  2.85E+00 -3.26E+00  2.63E+02
 
 TH10
+       -4.92E+00  1.15E+01 -1.17E+01 -2.52E+01 -7.85E+01  3.57E+00  2.23E+00  2.16E+01 -2.41E+00  7.43E+01
 
 TH11
+       -9.78E+00 -3.44E+00 -6.99E+00 -8.61E+00 -8.82E-01  1.56E+00  3.30E-01  5.11E+00  1.31E+01  8.00E+00  1.92E+02
 
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
 #CPUT: Total CPU Time in Seconds,       25.870
Stop Time:
Sat Sep 18 14:20:46 CDT 2021
