Thu Sep 30 10:06:50 CDT 2021
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
$DATA ../../../../data/spa2/D/dat94.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   38082.7964344841        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0269E+03  8.0171E+02  1.1238E+02  7.0857E+02 -1.6967E+02 -3.2813E+03 -1.7774E+03 -1.0935E+02 -2.2621E+03 -8.4439E+02
            -7.1735E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -334.136481830447        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.1459E+00  1.2331E+00  9.7826E-01  1.4595E+00  1.1381E+00  2.4921E+00  1.3415E+00  9.6252E-01  1.2217E+00  8.7586E-01
             1.4227E+01
 PARAMETER:  2.3621E-01  3.0952E-01  7.8016E-02  4.7812E-01  2.2940E-01  1.0131E+00  3.9382E-01  6.1796E-02  3.0024E-01 -3.2548E-02
             2.7551E+00
 GRADIENT:  -1.6690E+01  3.1521E+01 -4.3135E+00  6.2986E+01  2.6129E+00  7.3070E+01 -2.6663E+01  4.7287E+00 -3.4980E+01  8.9421E+00
            -8.8398E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -406.983265634456        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1327E+00  2.7874E+00  1.1319E+00  8.4420E-01  1.8538E+00  2.7938E+00  4.4135E+00  5.4660E-01  2.9305E+00  3.3585E-02
             1.3991E+01
 PARAMETER:  2.2462E-01  1.1251E+00  2.2390E-01 -6.9361E-02  7.1723E-01  1.1274E+00  1.5847E+00 -5.0403E-01  1.1752E+00 -3.2937E+00
             2.7384E+00
 GRADIENT:  -2.0317E+01  1.8421E+01 -1.9311E+01  1.0799E+01  1.1405E+01  5.4536E+01  2.4441E+01  1.9129E-01  2.7650E+01  9.4289E-03
             8.3494E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -451.192764610025        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.2196E+00  1.2354E+00  2.5561E+00  1.4083E+00  1.4263E+00  2.3535E+00  5.6442E+00  1.2933E+00  1.3944E+00  3.7378E-01
             1.3544E+01
 PARAMETER:  2.9849E-01  3.1136E-01  1.0385E+00  4.4236E-01  4.5511E-01  9.5590E-01  1.8306E+00  3.5717E-01  4.3244E-01 -8.8409E-01
             2.7059E+00
 GRADIENT:   1.0624E+01  7.4313E+00 -3.6399E+00 -2.7169E+00 -1.0052E+01  4.8805E+00  2.8518E+01  2.2933E+00  1.2167E+01  1.5144E+00
             4.0811E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -452.185185594000        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.1933E+00  1.1815E+00  3.3135E+00  1.3887E+00  1.6651E+00  2.3478E+00  5.3836E+00  1.1681E+00  1.2291E+00  3.7261E-01
             1.3432E+01
 PARAMETER:  2.7674E-01  2.6675E-01  1.2980E+00  4.2840E-01  6.0990E-01  9.5348E-01  1.7834E+00  2.5539E-01  3.0625E-01 -8.8722E-01
             2.6976E+00
 GRADIENT:   4.1503E+00  1.5000E+00 -2.3453E+00 -1.6882E+00 -1.9046E+00  5.1431E+00  1.0335E+01  1.0791E+00  5.7830E+00  1.2721E+00
             3.7740E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -466.277752515970        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  1.1832E+00  8.7339E-01  7.6268E+00  1.7010E+00  2.1681E+00  2.4988E+00  7.4537E+00  1.4333E+00  1.1784E+00  2.1704E-01
             1.3368E+01
 PARAMETER:  2.6822E-01 -3.5375E-02  2.1317E+00  6.3123E-01  8.7385E-01  1.0158E+00  2.1087E+00  4.5996E-01  2.6414E-01 -1.4277E+00
             2.6929E+00
 GRADIENT:  -1.7352E+00  8.0486E+00 -1.3340E+00  5.1167E+00  3.1597E+00  1.3623E+01  4.3236E+00  3.5654E-01 -4.8758E+00  2.7962E-01
            -1.5914E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -471.679793136391        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  1.1875E+00  3.8210E-01  3.3852E+01  1.9564E+00  2.2714E+00  2.3661E+00  8.4465E+00  5.4077E+00  1.4647E+00  9.0669E-02
             1.3367E+01
 PARAMETER:  2.7183E-01 -8.6206E-01  3.6220E+00  7.7111E-01  9.2039E-01  9.6124E-01  2.2338E+00  1.7878E+00  4.8164E-01 -2.3005E+00
             2.6928E+00
 GRADIENT:  -1.7872E-01  3.3525E-01  7.4153E-02 -2.9119E+00 -9.0969E-01  7.8838E-01  1.4329E+00  2.9876E-01 -6.9291E-01  5.2019E-02
            -1.0529E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -476.010579966029        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  1.1933E+00  3.2329E-01  3.9702E+01  1.9883E+00  2.3610E+00  2.3576E+00  1.1018E+01  9.4221E-01  1.4769E+00  1.0565E-02
             1.3612E+01
 PARAMETER:  2.7672E-01 -1.0292E+00  3.7814E+00  7.8729E-01  9.5908E-01  9.5766E-01  2.4995E+00  4.0477E-02  4.8996E-01 -4.4502E+00
             2.7109E+00
 GRADIENT:   2.5693E+01  7.8125E+00  2.3265E-01 -1.4512E+01  2.4244E+00  4.4089E+01  2.2170E+02  6.1702E-03  5.0520E+00  7.1260E-04
             3.3064E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -479.903937489673        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      922             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1464E+00  2.8794E-01  1.8538E+01  2.0625E+00  2.3242E+00  2.1870E+00  9.9349E+00  1.4390E-01  1.4683E+00  1.1166E-02
             1.3715E+01
 PARAMETER:  2.3660E-01 -1.1450E+00  3.0198E+00  8.2391E-01  9.4339E-01  8.8253E-01  2.3960E+00 -1.8387E+00  4.8413E-01 -4.3949E+00
             2.7185E+00
 GRADIENT:   7.7875E+00  5.3043E+00  7.7404E-02  1.3371E+01  1.9738E+00  2.2597E+01  1.8330E+02  8.5675E-04  1.1584E+00  7.4855E-04
             3.1428E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -479.909633844336        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  1.1442E+00  2.8882E-01  1.6282E+01  2.0570E+00  2.2929E+00  2.1803E+00  9.8529E+00  1.0000E-02  1.4715E+00  1.1109E-02
             1.3753E+01
 PARAMETER:  2.3474E-01 -1.1419E+00  2.8900E+00  8.2126E-01  9.2982E-01  8.7945E-01  2.3878E+00 -5.0737E+00  4.8630E-01 -4.4000E+00
             2.7212E+00
 GRADIENT:  -1.5331E+00  1.9143E+00 -1.2928E-02 -6.0397E+00  1.3920E+00  3.6275E+00  2.2920E+01  0.0000E+00  1.6569E-01  7.0150E-04
             2.3633E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -479.976383912828        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.1442E+00  2.8388E-01  1.6211E+01  2.0752E+00  2.2562E+00  2.1728E+00  9.9106E+00  1.0000E-02  1.4707E+00  1.0000E-02
             1.3829E+01
 PARAMETER:  2.3468E-01 -1.1592E+00  2.8857E+00  8.3007E-01  9.1366E-01  8.7602E-01  2.3936E+00 -5.3721E+00  4.8574E-01 -4.5172E+00
             2.7268E+00
 GRADIENT:  -2.7253E+00  2.0484E+00  9.7233E-02 -2.5271E+00 -1.9917E-01  3.0203E+00  2.2891E+01  0.0000E+00 -5.3499E-01  2.2295E-04
             5.6848E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -480.054909867323        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1393
 NPARAMETR:  1.1482E+00  2.7463E-01  1.5635E+01  2.0875E+00  2.2609E+00  2.1677E+00  9.9330E+00  1.0000E-02  1.4785E+00  1.0000E-02
             1.3863E+01
 PARAMETER:  2.3581E-01 -1.1934E+00  2.8532E+00  8.3588E-01  9.1594E-01  8.7446E-01  2.3936E+00 -5.3721E+00  4.8949E-01 -4.5380E+00
             2.7275E+00
 GRADIENT:  -2.4350E+00 -1.1551E+02  1.5388E-02 -3.4067E-01  7.2362E-02  8.0196E+01 -2.9996E+01  0.0000E+00 -4.6618E-01  4.3831E-05
            -4.3861E+01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1393
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0610E-02  6.5224E-02  5.6168E-06 -9.3115E-02  4.3105E-07
 SE:             2.7810E-02  2.0347E-02  2.4941E-06  1.6550E-02  3.1197E-05
 N:                     100         100         100         100         100

 P VAL.:         4.5865E-01  1.3478E-03  2.4317E-02  1.8447E-08  9.8898E-01

 ETASHRINKSD(%)  6.8315E+00  3.1835E+01  9.9992E+01  4.4557E+01  9.9895E+01
 ETASHRINKVR(%)  1.3196E+01  5.3536E+01  1.0000E+02  6.9261E+01  1.0000E+02
 EBVSHRINKSD(%)  8.3406E+00  3.6583E+01  9.9980E+01  3.5689E+01  9.9801E+01
 EBVSHRINKVR(%)  1.5986E+01  5.9783E+01  1.0000E+02  5.8641E+01  1.0000E+02
 RELATIVEINF(%)  8.2659E+01  2.5752E+01  6.0951E-07  2.2456E+01  5.6257E-05
 EPSSHRINKSD(%)  5.5360E+00
 EPSSHRINKVR(%)  1.0766E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -480.05490986732269     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       622.67132997828435     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    33.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -480.055       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.15E+00  2.74E-01  1.57E+01  2.09E+00  2.26E+00  2.17E+00  9.91E+00  1.00E-02  1.48E+00  1.00E-02  1.38E+01
 


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
+        4.78E+04
 
 TH 2
+       -3.94E+04  3.27E+04
 
 TH 3
+       -1.28E-02  1.14E-01  3.18E-03
 
 TH 4
+       -7.44E+03  6.16E+03  1.13E-02  1.23E+03
 
 TH 5
+        8.69E-02 -2.72E+01 -1.04E-01 -6.95E+00  8.15E+00
 
 TH 6
+       -1.54E+02  1.24E+02 -1.07E-02  1.07E+03 -9.00E+02  1.93E+03
 
 TH 7
+        4.28E+00  3.91E+00 -2.34E-04 -1.83E+00 -1.84E-01  2.38E+00  6.88E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.79E+04 -1.48E+04 -3.17E-02 -2.81E+03  2.29E+00 -5.93E+01  2.97E-01  0.00E+00  6.73E+03
 
 TH10
+        1.12E-01  3.68E-03  2.70E-03  2.54E-02  4.35E-03  9.69E-05 -1.97E-03  0.00E+00  9.24E-03  2.26E+01
 
 TH11
+        1.11E+01 -2.88E+00  4.44E-03 -5.04E+00 -2.69E-01  5.72E+00 -6.31E-01  0.00E+00  2.09E+00 -1.01E-03  4.83E+00
 
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
 #CPUT: Total CPU Time in Seconds,       46.694
Stop Time:
Thu Sep 30 10:07:38 CDT 2021
