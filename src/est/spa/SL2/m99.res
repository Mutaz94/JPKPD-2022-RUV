Sat Sep 18 12:35:32 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat99.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1681.20668012724        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.9366E+01 -7.8404E+01 -3.8254E+01 -4.1250E+01  6.2595E+01  3.3742E+01 -1.7568E+01 -2.2853E+00  6.0549E+00  1.1102E+01
             5.6240E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1688.42755065593        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9523E-01  1.1979E+00  1.0245E+00  9.1365E-01  1.0273E+00  8.6556E-01  1.2712E+00  1.1347E+00  9.2772E-01  7.9684E-01
             1.0008E+00
 PARAMETER:  9.5216E-02  2.8055E-01  1.2423E-01  9.6908E-03  1.2698E-01 -4.4383E-02  3.3997E-01  2.2636E-01  2.4975E-02 -1.2710E-01
             1.0079E-01
 GRADIENT:   9.4306E+00  2.3392E+01  1.6206E+01 -1.1164E+01  5.6621E+00 -2.2447E+01  1.5697E+01 -9.6459E+00 -4.7992E+00 -6.2087E+00
            -6.6772E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1689.48119733074        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.9419E-01  1.1256E+00  1.0576E+00  9.5475E-01  9.9865E-01  8.7344E-01  1.2743E+00  1.4115E+00  8.8550E-01  6.9607E-01
             9.9509E-01
 PARAMETER:  9.4171E-02  2.1828E-01  1.5601E-01  5.3694E-02  9.8645E-02 -3.5321E-02  3.4239E-01  4.4463E-01 -2.1600E-02 -2.6231E-01
             9.5081E-02
 GRADIENT:   8.3304E+00  7.2711E+00  5.8726E+00 -9.2447E+00  7.4185E+00 -1.8420E+01  8.2418E+00 -2.4141E+00 -4.2155E+00 -1.0082E+01
            -3.1758E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1690.79519823657        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      260
 NPARAMETR:  9.9073E-01  1.2692E+00  1.1695E+00  8.8400E-01  1.0997E+00  9.0012E-01  1.0681E+00  1.8712E+00  1.0104E+00  8.2487E-01
             9.9339E-01
 PARAMETER:  9.0686E-02  3.3836E-01  2.5654E-01 -2.3298E-02  1.9507E-01 -5.2279E-03  1.6586E-01  7.2657E-01  1.1036E-01 -9.2526E-02
             9.3364E-02
 GRADIENT:  -4.2744E+01  3.9470E-01 -3.5921E+00  6.7825E+00  8.0129E+00 -1.0870E+01  1.3981E+00  3.9164E+00  1.3056E+00 -3.9294E+00
            -1.0072E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1691.42573548345        NO. OF FUNC. EVALS.: 144
 CUMULATIVE NO. OF FUNC. EVALS.:      404
 NPARAMETR:  1.0020E+00  1.2615E+00  1.1777E+00  8.7889E-01  1.0963E+00  9.2079E-01  1.0541E+00  1.7648E+00  1.0033E+00  8.6047E-01
             9.9549E-01
 PARAMETER:  1.0205E-01  3.3231E-01  2.6359E-01 -2.9101E-02  1.9196E-01  1.7476E-02  1.5268E-01  6.6804E-01  1.0329E-01 -5.0276E-02
             9.5484E-02
 GRADIENT:   3.3086E+01  1.1356E+01  3.6416E+00 -4.6421E+00 -3.6684E+00  3.6229E+00 -5.0231E-02  6.4691E-01 -8.2109E-01 -1.2783E-01
             4.9803E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1691.55544187234        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      581
 NPARAMETR:  1.0055E+00  1.2674E+00  1.1424E+00  8.8202E-01  1.0935E+00  9.2311E-01  1.0578E+00  1.6838E+00  1.0136E+00  8.6976E-01
             9.9446E-01
 PARAMETER:  1.0547E-01  3.3700E-01  2.3313E-01 -2.5539E-02  1.8942E-01  1.9988E-02  1.5615E-01  6.2106E-01  1.1348E-01 -3.9541E-02
             9.4443E-02
 GRADIENT:  -2.6017E+00  4.9162E-01  8.5994E-01  1.1442E+00  9.7031E-01 -9.0130E-02 -2.7970E-02  3.1456E-03  1.3625E-01  5.4991E-01
             1.2715E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1691.56480764782        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      745
 NPARAMETR:  1.0066E+00  1.2674E+00  1.1379E+00  8.8039E-01  1.0908E+00  9.2311E-01  1.0592E+00  1.6779E+00  1.0120E+00  8.6091E-01
             9.9422E-01
 PARAMETER:  1.0657E-01  3.3698E-01  2.2919E-01 -2.7389E-02  1.8687E-01  1.9992E-02  1.5755E-01  6.1752E-01  1.1194E-01 -4.9765E-02
             9.4203E-02
 GRADIENT:   2.5175E-01 -4.0592E-01  1.7012E+00 -9.8015E-01  1.1076E-01 -8.4518E-02 -5.7273E-02 -2.4143E-01 -7.5722E-02 -1.9904E-01
            -1.3787E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1691.56906465574        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      921            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0065E+00  1.2672E+00  1.1341E+00  8.8117E-01  1.0907E+00  9.2331E-01  1.0596E+00  1.6801E+00  1.0121E+00  8.6248E-01
             9.9443E-01
 PARAMETER:  1.0648E-01  3.3681E-01  2.2588E-01 -2.6507E-02  1.8685E-01  2.0212E-02  1.5786E-01  6.1885E-01  1.1204E-01 -4.7946E-02
             9.4414E-02
 GRADIENT:   4.6093E+01  1.7681E+01  8.8586E-01  5.9451E+00  2.5866E+00  4.5919E+00  7.9564E-01  3.7746E-01  6.6671E-01  1.1387E-01
             1.2964E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1691.60443576938        NO. OF FUNC. EVALS.:  80
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  1.0063E+00  1.2719E+00  1.0749E+00  8.7563E-01  1.0768E+00  9.2322E-01  1.0626E+00  1.6180E+00  1.0096E+00  8.5151E-01
             9.9375E-01
 PARAMETER:  1.0630E-01  3.4049E-01  1.7220E-01 -3.2817E-02  1.7401E-01  2.0117E-02  1.6076E-01  5.8119E-01  1.0952E-01 -6.0738E-02
             9.3734E-02
 GRADIENT:   4.5046E+01  1.5538E+01 -2.0161E+00  7.1179E+00  5.5323E+00  4.4895E+00  7.5381E-01  1.2668E+00  6.0252E-01  4.4519E-01
             1.3407E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1691.60836929455        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1168
 NPARAMETR:  1.0063E+00  1.2721E+00  1.0717E+00  8.7393E-01  1.0759E+00  9.2347E-01  1.0628E+00  1.6134E+00  1.0096E+00  8.4783E-01
             9.9371E-01
 PARAMETER:  1.0632E-01  3.4070E-01  1.6928E-01 -3.4752E-02  1.7317E-01  2.0388E-02  1.6086E-01  5.7837E-01  1.0951E-01 -6.5078E-02
             9.3694E-02
 GRADIENT:  -9.3542E-01 -4.2918E+00 -1.9316E+00 -3.0254E-02  4.5517E+00  9.9451E-04 -1.1184E-01  9.5324E-01 -2.4503E-02  4.6982E-02
            -2.3149E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1691.61112734351        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  1.0064E+00  1.2724E+00  1.0723E+00  8.7395E-01  1.0757E+00  9.2347E-01  1.0628E+00  1.6135E+00  1.0096E+00  8.4739E-01
             9.9371E-01
 PARAMETER:  1.0634E-01  3.4093E-01  1.6979E-01 -3.4727E-02  1.7298E-01  2.0385E-02  1.6092E-01  5.7838E-01  1.0953E-01 -6.5497E-02
             9.3699E-02
 GRADIENT:  -2.5609E+05  1.5975E+05  3.2081E+05  5.3786E-02 -3.1487E+05 -8.2456E-04  3.3847E+05 -9.4268E+04  4.9728E+05  2.8531E-02
             5.4467E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         4.6         3.3         3.3         3.3         2.2
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1309
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.7034E-04 -1.0669E-02 -3.7409E-02  7.4921E-03 -4.3032E-02
 SE:             2.9840E-02  2.2049E-02  1.5787E-02  2.2121E-02  1.9498E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7673E-01  6.2847E-01  1.7811E-02  7.3485E-01  2.7310E-02

 ETASHRINKSD(%)  3.0800E-02  2.6132E+01  4.7110E+01  2.5892E+01  3.4681E+01
 ETASHRINKVR(%)  6.1590E-02  4.5436E+01  7.2027E+01  4.5080E+01  5.7334E+01
 EBVSHRINKSD(%)  4.9347E-01  2.6524E+01  4.8952E+01  2.7050E+01  3.2303E+01
 EBVSHRINKVR(%)  9.8450E-01  4.6013E+01  7.3941E+01  4.6782E+01  5.4172E+01
 RELATIVEINF(%)  9.8440E+01  1.6543E+00  1.9551E+00  1.6588E+00  1.0994E+01
 EPSSHRINKSD(%)  4.5206E+01
 EPSSHRINKVR(%)  6.9976E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1691.6111273435147     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -956.46030077977650     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.24
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1691.611       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.27E+00  1.07E+00  8.74E-01  1.08E+00  9.23E-01  1.06E+00  1.61E+00  1.01E+00  8.47E-01  9.94E-01
 


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
+        1.19E+09
 
 TH 2
+        1.02E+03  7.23E+07
 
 TH 3
+        2.47E+03 -1.20E+04  4.11E+08
 
 TH 4
+        5.23E+04 -1.25E+04 -3.08E+04  1.78E+09
 
 TH 5
+       -1.54E+04  3.64E+03  8.85E+03  3.02E+04  3.93E+08
 
 TH 6
+       -8.05E+03  1.99E+03  4.74E+03  1.69E+09 -4.63E+03  1.60E+09
 
 TH 7
+        9.86E+03 -2.41E+03 -5.79E+03  9.11E+08  5.67E+03  8.62E+08  4.66E+08
 
 TH 8
+       -4.80E+02  2.35E+03 -1.04E+04  6.01E+03 -1.77E+03 -9.24E+02  1.13E+03  1.57E+07
 
 TH 9
+        4.24E+04 -1.05E+04 -2.49E+04 -5.06E+04  2.44E+04  7.80E+03 -9.52E+03  4.87E+03  1.11E+09
 
 TH10
+       -1.61E+03  3.97E+02  9.41E+02 -1.29E+01 -1.02E+03  1.74E+09  1.01E+03 -1.75E+02  1.56E+03  7.04E+01
 
 TH11
+       -1.28E+09  3.37E+03  8.05E+03 -5.63E+04  1.65E+04  8.67E+03 -1.06E+04 -1.57E+03 -4.56E+04  1.74E+03  1.38E+09
 
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
 #CPUT: Total CPU Time in Seconds,       23.086
Stop Time:
Sat Sep 18 12:35:57 CDT 2021
