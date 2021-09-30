Wed Sep 29 13:36:41 CDT 2021
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
$DATA ../../../../data/spa/A3/dat52.csv ignore=@
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
Current Date:       29 SEP 2021
Days until program expires : 200
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
 RAW OUTPUT FILE (FILE): m52.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -25.3506324495773        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4029E+02  3.0259E+01  1.4437E+02 -1.5803E+02  4.6408E+01  3.3805E+01 -4.4396E+01 -2.5249E+01 -1.4773E+02 -1.1948E+02
            -2.8842E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1045.50830454761        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7448E-01  8.9428E-01  8.9485E-01  1.1533E+00  8.8409E-01  9.2329E-01  1.0212E+00  8.5115E-01  1.4518E+00  1.2309E+00
             1.8126E+00
 PARAMETER:  7.4148E-02 -1.1734E-02 -1.1099E-02  2.4266E-01 -2.3200E-02  2.0193E-02  1.2101E-01 -6.1161E-02  4.7277E-01  3.0775E-01
             6.9477E-01
 GRADIENT:   1.0698E+02  2.5004E+01  3.2943E+01  2.6209E+01 -1.4094E+01 -6.4231E+00  2.6439E+00  9.1561E+00  3.3612E+01 -1.0156E+01
            -7.2108E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1134.91519709698        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7745E-01  5.2522E-01  5.1288E-01  1.3969E+00  5.2560E-01  9.2822E-01  3.0580E+00  1.7834E-02  1.1490E+00  6.4474E-01
             2.0908E+00
 PARAMETER:  7.7187E-02 -5.4393E-01 -5.6771E-01  4.3428E-01 -5.4321E-01  2.5509E-02  1.2178E+00 -3.9266E+00  2.3886E-01 -3.3891E-01
             8.3754E-01
 GRADIENT:   6.8273E+01  4.0548E+01 -8.7190E+01  1.4768E+02  1.8058E+02 -6.0662E+00  7.3141E+01  6.7479E-03  1.6952E+01 -2.3742E+01
            -4.5926E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1266.76293278673        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.9400E-01  3.7437E-01  3.0490E-01  1.3024E+00  3.1574E-01  9.9088E-01  1.1415E+00  1.0000E-02  1.0020E+00  4.7173E-01
             3.3565E+00
 PARAMETER:  9.3986E-02 -8.8252E-01 -1.0878E+00  3.6425E-01 -1.0528E+00  9.0838E-02  2.3235E-01 -6.9965E+00  1.0201E-01 -6.5135E-01
             1.3109E+00
 GRADIENT:   6.6283E-01  2.2536E+01  4.1697E+01  2.9760E+01 -1.9805E+01  1.7818E+01  6.4718E-01  0.0000E+00 -1.9597E+01  1.5815E+00
            -4.9630E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1268.87583817035        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      348
 NPARAMETR:  1.0036E+00  3.5994E-01  2.9020E-01  1.2966E+00  3.0773E-01  9.5605E-01  1.0014E+00  1.0000E-02  1.0639E+00  4.6581E-01
             3.4736E+00
 PARAMETER:  1.0355E-01 -9.2181E-01 -1.1372E+00  3.5977E-01 -1.0785E+00  5.5058E-02  1.0139E-01 -7.4925E+00  1.6194E-01 -6.6398E-01
             1.3452E+00
 GRADIENT:  -8.0604E+00  7.8657E+00  3.2300E+00  3.5400E+00 -1.1810E+01  4.5762E+00  8.4798E-01  0.0000E+00 -7.1748E+00  3.3561E+00
             9.5712E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1270.94024693550        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  9.9818E-01  2.4001E-01  3.4492E-01  1.3888E+00  3.2630E-01  9.3916E-01  1.0652E+00  1.0000E-02  1.0357E+00  4.3163E-01
             3.4946E+00
 PARAMETER:  9.8180E-02 -1.3271E+00 -9.6443E-01  4.2844E-01 -1.0199E+00  3.7234E-02  1.6321E-01 -8.3706E+00  1.3510E-01 -7.4019E-01
             1.3512E+00
 GRADIENT:  -8.0470E+00  5.3095E+00  1.3432E+01  3.2500E+00 -2.3846E+01  2.2341E+00 -1.2212E-01  0.0000E+00 -2.2421E+00  4.0394E-01
             1.5522E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1272.41281909871        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  9.9110E-01  1.2758E-01  4.6389E-01  1.5145E+00  4.0247E-01  9.2287E-01  8.0652E-01  1.0000E-02  9.6243E-01  3.2955E-01
             3.6125E+00
 PARAMETER:  9.1059E-02 -1.9590E+00 -6.6810E-01  5.1506E-01 -8.1014E-01  1.9735E-02 -1.1503E-01 -1.0116E+01  6.1710E-02 -1.0100E+00
             1.3844E+00
 GRADIENT:  -4.9805E-01  5.5853E-01 -3.0668E+00  1.1038E+00  6.0477E+00  2.4206E+00  2.7963E-02  0.0000E+00  1.0373E+00  4.2194E-01
             4.4879E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1272.64877781469        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      875
 NPARAMETR:  9.8837E-01  7.6111E-02  4.4857E-01  1.5248E+00  3.8727E-01  9.1446E-01  6.2142E-01  1.0000E-02  9.4889E-01  2.3418E-01
             3.6299E+00
 PARAMETER:  8.8305E-02 -2.4756E+00 -7.0169E-01  5.2185E-01 -8.4863E-01  1.0582E-02 -3.7575E-01 -1.1707E+01  4.7536E-02 -1.3517E+00
             1.3892E+00
 GRADIENT:   5.0444E-01  3.6863E-01  9.9207E-02  6.1833E-01  2.9536E-01 -3.1263E-02  1.0705E-03  0.0000E+00  7.1321E-02 -5.4773E-01
             7.5782E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1273.34889393670        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1055
 NPARAMETR:  9.8379E-01  1.0195E-02  3.7799E-01  1.5012E+00  3.3061E-01  9.2471E-01  1.3051E-01  1.0000E-02  9.7255E-01  3.9269E-01
             3.5216E+00
 PARAMETER:  8.3658E-02 -4.4859E+00 -8.7289E-01  5.0624E-01 -1.0068E+00  2.1729E-02 -1.9363E+00 -1.6495E+01  7.2168E-02 -8.3473E-01
             1.3589E+00
 GRADIENT:  -3.8592E+00  5.4555E-02  2.3833E+00  3.1684E+00 -4.7199E+00  1.2213E+00  2.2397E-06  0.0000E+00  2.6295E-01 -1.6880E-01
             2.3496E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1273.36652057096        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  9.8507E-01  1.0000E-02  3.8160E-01  1.5013E+00  3.3354E-01  9.2108E-01  1.2389E-01  1.0000E-02  9.6967E-01  4.0183E-01
             3.5099E+00
 PARAMETER:  8.4957E-02 -4.5700E+00 -8.6339E-01  5.0634E-01 -9.9799E-01  1.7790E-02 -1.9884E+00 -1.6606E+01  6.9196E-02 -8.1173E-01
             1.3556E+00
 GRADIENT:   1.8657E-01  0.0000E+00 -1.0039E+00 -8.6154E-01  1.5694E+00  4.1913E-02  2.6509E-06  0.0000E+00  9.2742E-02  1.0475E-01
             4.1600E-01

0ITERATION NO.:   46    OBJECTIVE VALUE:  -1273.36652057096        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1254
 NPARAMETR:  9.8507E-01  1.0000E-02  3.8160E-01  1.5013E+00  3.3354E-01  9.2108E-01  1.2389E-01  1.0000E-02  9.6967E-01  4.0183E-01
             3.5099E+00
 PARAMETER:  8.4957E-02 -4.5700E+00 -8.6339E-01  5.0634E-01 -9.9799E-01  1.7790E-02 -1.9884E+00 -1.6606E+01  6.9196E-02 -8.1173E-01
             1.3556E+00
 GRADIENT:   1.8657E-01  0.0000E+00 -1.0039E+00 -8.6154E-01  1.5694E+00  4.1913E-02  2.6509E-06  0.0000E+00  9.2742E-02  1.0475E-01
             4.1600E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1254
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.2448E-04 -3.1893E-05  1.2792E-04 -1.4668E-02 -9.1118E-04
 SE:             2.8344E-02  1.5612E-05  2.2031E-04  2.5913E-02  1.2485E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9368E-01  4.1062E-02  5.6148E-01  5.7135E-01  9.4182E-01

 ETASHRINKSD(%)  5.0438E+00  9.9948E+01  9.9262E+01  1.3188E+01  5.8173E+01
 ETASHRINKVR(%)  9.8331E+00  1.0000E+02  9.9995E+01  2.4637E+01  8.2505E+01
 EBVSHRINKSD(%)  4.6360E+00  9.9950E+01  9.9207E+01  1.2370E+01  5.7976E+01
 EBVSHRINKVR(%)  9.0570E+00  1.0000E+02  9.9994E+01  2.3210E+01  8.2340E+01
 RELATIVEINF(%)  6.5992E+01  1.9669E-06  1.4004E-04  1.7713E+01  2.6578E-01
 EPSSHRINKSD(%)  2.4367E+01
 EPSSHRINKVR(%)  4.2797E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1273.3665205709615     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -538.21569400722331     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1273.367       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.00E-02  3.82E-01  1.50E+00  3.34E-01  9.21E-01  1.24E-01  1.00E-02  9.70E-01  4.02E-01  3.51E+00
 


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
+        1.23E+03
 
 TH 2
+        0.00E+00  7.46E+02
 
 TH 3
+       -1.12E+02  0.00E+00  4.56E+03
 
 TH 4
+       -5.15E+01  0.00E+00 -1.76E+02  4.04E+02
 
 TH 5
+        2.71E+02  0.00E+00 -6.58E+03 -2.60E+02  1.05E+04
 
 TH 6
+       -5.52E-01  0.00E+00  3.36E+01 -1.61E+01 -5.04E+00  1.93E+02
 
 TH 7
+        4.84E-02  0.00E+00  2.07E-04 -9.70E-03 -8.86E-03  3.15E-02  1.32E-02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.46E+01  0.00E+00  3.92E+01 -1.64E+01  6.37E+01  7.56E-01  5.56E-02  0.00E+00  1.25E+02
 
 TH10
+       -1.29E+01  0.00E+00 -1.20E+02 -6.97E+00  2.10E+02 -3.42E-01 -2.81E-03  0.00E+00  8.63E-01  4.29E+01
 
 TH11
+       -2.02E+01  0.00E+00 -8.11E+00 -5.14E+00  2.47E+00  4.14E+00  1.56E-05  0.00E+00  8.05E+00  1.83E+01  3.08E+01
 
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
 #CPUT: Total CPU Time in Seconds,       21.252
Stop Time:
Wed Sep 29 13:37:04 CDT 2021
