Wed Sep 29 11:10:31 CDT 2021
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
$DATA ../../../../data/spa/B/dat36.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1668.85749862822        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1138E+02 -3.3482E+01 -1.3753E+01  2.8802E+01  8.8323E+01  4.5151E+01 -4.4821E+00 -9.1840E-01  3.3993E+01 -5.0251E+01
             3.4389E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1682.03215510380        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  1.0247E+00  1.0583E+00  9.0234E-01  1.0122E+00  9.6908E-01  1.1078E+00  1.0394E+00  9.9376E-01  8.1682E-01  1.3255E+00
             8.6719E-01
 PARAMETER:  1.2443E-01  1.5666E-01 -2.7659E-03  1.1209E-01  6.8591E-02  2.0239E-01  1.3861E-01  9.3743E-02 -1.0234E-01  3.8183E-01
            -4.2492E-02
 GRADIENT:   1.7243E+01 -5.9600E+00 -2.7724E+01  2.4819E+01  2.2864E+01  1.2934E+01 -5.1388E+00  1.0299E+01 -2.8049E+00  9.3299E+00
            -1.5763E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1683.71212938761        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0035E+00  8.2452E-01  8.7759E-01  1.1734E+00  8.6235E-01  1.0922E+00  1.4232E+00  4.9956E-01  6.9753E-01  1.3325E+00
             8.6186E-01
 PARAMETER:  1.0345E-01 -9.2952E-02 -3.0572E-02  2.5990E-01 -4.8095E-02  1.8818E-01  4.5292E-01 -5.9403E-01 -2.6021E-01  3.8709E-01
            -4.8664E-02
 GRADIENT:  -2.0939E+01  2.2582E+01 -4.0828E+01  8.8202E+01  3.4618E+01  7.7718E+00 -5.8098E-01  3.8332E+00 -6.8998E+00  1.3460E+01
            -1.4688E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1687.76142935830        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      550
 NPARAMETR:  1.0163E+00  7.2246E-01  8.2577E-01  1.1997E+00  7.6049E-01  1.0713E+00  1.5541E+00  2.4618E-01  7.0634E-01  1.1889E+00
             8.9393E-01
 PARAMETER:  1.1619E-01 -2.2509E-01 -9.1438E-02  2.8209E-01 -1.7379E-01  1.6888E-01  5.4088E-01 -1.3017E+00 -2.4766E-01  2.7299E-01
            -1.2131E-02
 GRADIENT:   3.1825E+00  8.2669E+00 -6.6508E-01  1.7071E+01 -7.0068E+00  4.3259E-01  5.5213E-01  6.8000E-01  4.1899E-01  3.0168E+00
             2.2835E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1688.09664435319        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      708
 NPARAMETR:  1.0174E+00  6.8126E-01  8.2738E-01  1.2123E+00  7.5390E-01  1.0734E+00  1.6284E+00  4.7309E-02  6.9542E-01  1.1800E+00
             8.8952E-01
 PARAMETER:  1.1723E-01 -2.8381E-01 -8.9487E-02  2.9248E-01 -1.8250E-01  1.7084E-01  5.8763E-01 -2.9511E+00 -2.6325E-01  2.6549E-01
            -1.7076E-02
 GRADIENT:   6.3375E+02  5.7563E+01  8.7744E+00  4.2367E+02  2.3350E+01  1.4379E+02  2.6593E+01  4.2805E-02  1.6244E+01  3.7792E+00
             5.4364E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1688.10756580237        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      886
 NPARAMETR:  1.0151E+00  6.8248E-01  8.2339E-01  1.2136E+00  7.5282E-01  1.0704E+00  1.6144E+00  4.3317E-02  6.9762E-01  1.1865E+00
             8.8885E-01
 PARAMETER:  1.1496E-01 -2.8202E-01 -9.4330E-02  2.9358E-01 -1.8393E-01  1.6800E-01  5.7899E-01 -3.0392E+00 -2.6008E-01  2.7103E-01
            -1.7826E-02
 GRADIENT:   1.5572E+00  1.4909E-01 -4.0847E-01  2.3524E+00  5.8416E-01  2.9396E-01 -4.8795E-02  1.6301E-02  1.1819E-01  4.4603E-01
            -2.4476E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1688.13197944633        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1064
 NPARAMETR:  1.0157E+00  6.8490E-01  8.1777E-01  1.2094E+00  7.5012E-01  1.0723E+00  1.6201E+00  1.3048E-02  6.9658E-01  1.1778E+00
             8.8902E-01
 PARAMETER:  1.1560E-01 -2.7848E-01 -1.0118E-01  2.9009E-01 -1.8752E-01  1.6985E-01  5.8249E-01 -4.2392E+00 -2.6157E-01  2.6362E-01
            -1.7633E-02
 GRADIENT:   2.6673E+00 -8.8820E-01  3.0237E-01 -2.5228E+00 -2.6337E-01  1.0082E+00  3.1820E-01  1.6050E-03  1.8378E-01  1.7126E-01
             9.7309E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1688.13900308827        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0158E+00  6.8916E-01  8.1581E-01  1.2077E+00  7.5089E-01  1.0724E+00  1.6147E+00  1.0000E-02  6.9681E-01  1.1768E+00
             8.8900E-01
 PARAMETER:  1.1563E-01 -2.7228E-01 -1.0358E-01  2.8870E-01 -1.8650E-01  1.6986E-01  5.7917E-01 -6.2061E+00 -2.6125E-01  2.6281E-01
            -1.7660E-02
 GRADIENT:   2.5885E+00 -3.7484E-01 -9.5510E-01  1.7760E-03  8.5707E-01  9.7842E-01  3.7461E-01  0.0000E+00  1.1956E-01  2.6882E-01
             1.2571E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1688.14476645749        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1435
 NPARAMETR:  1.0161E+00  6.9192E-01  8.1458E-01  1.2054E+00  7.5148E-01  1.0730E+00  1.6097E+00  1.0000E-02  6.9794E-01  1.1762E+00
             8.8912E-01
 PARAMETER:  1.1599E-01 -2.6829E-01 -1.0509E-01  2.8682E-01 -1.8572E-01  1.7048E-01  5.7607E-01 -6.2061E+00 -2.5962E-01  2.6227E-01
            -1.7527E-02
 GRADIENT:   3.2016E+00 -9.0481E-01 -1.3465E+00 -1.1992E+00  1.3107E+00  1.2067E+00  4.0318E-01  0.0000E+00  3.0188E-01  3.6593E-01
             2.2848E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1688.14577738563        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1492
 NPARAMETR:  1.0160E+00  6.9258E-01  8.1543E-01  1.2058E+00  7.5122E-01  1.0724E+00  1.6081E+00  1.0000E-02  6.9760E-01  1.1755E+00
             8.8896E-01
 PARAMETER:  1.1587E-01 -2.6733E-01 -1.0403E-01  2.8717E-01 -1.8605E-01  1.6994E-01  5.7504E-01 -6.2061E+00 -2.6011E-01  2.6173E-01
            -1.7708E-02
 GRADIENT:   2.9786E+00 -5.7932E-02 -5.7649E-02 -1.3543E-01 -1.0583E-01  9.8771E-01  3.2725E-01  0.0000E+00  8.2551E-02  6.0542E-02
             2.3037E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1492
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.0794E-04  1.0524E-02 -5.3152E-04 -1.4925E-02 -6.9181E-03
 SE:             2.9890E-02  2.0164E-02  2.0146E-04  2.3625E-02  2.5446E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9712E-01  6.0173E-01  8.3328E-03  5.2755E-01  7.8572E-01

 ETASHRINKSD(%)  1.0000E-10  3.2448E+01  9.9325E+01  2.0855E+01  1.4754E+01
 ETASHRINKVR(%)  1.0000E-10  5.4367E+01  9.9995E+01  3.7360E+01  2.7331E+01
 EBVSHRINKSD(%)  2.9949E-01  3.3793E+01  9.9413E+01  2.0316E+01  1.1469E+01
 EBVSHRINKVR(%)  5.9809E-01  5.6167E+01  9.9997E+01  3.6504E+01  2.1623E+01
 RELATIVEINF(%)  9.8796E+01  3.8155E+00  5.5187E-04  5.8908E+00  1.0563E+01
 EPSSHRINKSD(%)  4.4677E+01
 EPSSHRINKVR(%)  6.9394E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1688.1457773856309     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -952.99495082189276     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.09
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1688.146       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  6.93E-01  8.15E-01  1.21E+00  7.51E-01  1.07E+00  1.61E+00  1.00E-02  6.98E-01  1.18E+00  8.89E-01
 


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
+        9.30E+02
 
 TH 2
+       -8.83E+00  4.16E+02
 
 TH 3
+        1.82E+01  1.35E+02  6.74E+02
 
 TH 4
+       -6.33E+00  4.56E+02 -2.74E+02  9.98E+02
 
 TH 5
+       -2.17E+00 -2.71E+02 -7.09E+02  2.48E+02  1.11E+03
 
 TH 6
+       -3.70E-01 -2.53E+00  3.25E+00 -1.19E+00 -9.74E-01  1.71E+02
 
 TH 7
+        1.35E+00  2.88E+01  9.29E-01 -2.54E+00 -7.41E+00  1.97E-01  1.90E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.04E+00 -1.67E+01 -3.15E+01 -1.87E+01  3.10E+01 -5.47E-01  2.43E+01  0.00E+00  1.76E+02
 
 TH10
+       -9.16E-01 -9.44E-01 -7.34E+01 -2.14E+01 -3.03E+01  1.35E-01  6.67E+00  0.00E+00  4.08E+00  7.91E+01
 
 TH11
+       -6.02E+00 -1.52E+01 -4.78E+01 -1.09E+00  1.28E+01  1.25E+00  8.56E-01  0.00E+00  1.43E+01  1.80E+01  2.71E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.696
Stop Time:
Wed Sep 29 11:11:00 CDT 2021
