Wed Sep 29 22:02:38 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat35.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m35.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1292.60757044675        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.3315E+02  4.7665E+01  5.3357E+01  7.9234E+01  1.3835E+02  5.8124E+01 -4.0753E+01 -6.7215E+01 -1.1835E+01 -8.4127E+01
            -1.4750E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1776.26856716904        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0726E+00  8.9884E-01  1.1055E+00  1.0875E+00  9.2213E-01  9.8231E-01  1.0286E+00  6.6709E-01  1.0821E+00  9.5970E-01
             2.1108E+00
 PARAMETER:  1.7012E-01 -6.6520E-03  2.0028E-01  1.8392E-01  1.8929E-02  8.2147E-02  1.2816E-01 -3.0483E-01  1.7887E-01  5.8868E-02
             8.4708E-01
 GRADIENT:   1.9776E+02  2.0490E+01  2.5028E+01  2.8532E+01 -1.5270E+01  2.7779E+01  2.0077E+00  2.7840E+00  1.8528E+01 -1.0027E+01
            -5.6245E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1783.56675793456        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0681E+00  1.0127E+00  4.7526E-01  9.6736E-01  6.4620E-01  8.9228E-01  1.1750E+00  1.5535E-01  9.5776E-01  7.3207E-01
             2.1691E+00
 PARAMETER:  1.6588E-01  1.1266E-01 -6.4390E-01  6.6819E-02 -3.3665E-01 -1.3970E-02  2.6123E-01 -1.7621E+00  5.6838E-02 -2.1188E-01
             8.7431E-01
 GRADIENT:   1.4809E+02  2.9189E+01 -1.2024E+01  2.4455E+01  3.4625E+01 -1.2648E+01  1.1888E+01  4.7057E-01 -7.1557E-01  1.1370E+01
            -2.6051E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1784.90568530781        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.0326E+00  8.5897E-01  3.4642E-01  1.0167E+00  4.9184E-01  9.0992E-01  1.3099E+00  1.1811E-01  9.1509E-01  4.9478E-01
             2.1699E+00
 PARAMETER:  1.3211E-01 -5.2025E-02 -9.6011E-01  1.1658E-01 -6.0960E-01  5.6047E-03  3.6993E-01 -2.0361E+00  1.1271E-02 -6.0365E-01
             8.7466E-01
 GRADIENT:   1.9099E+01  4.2894E+01 -3.8153E+01  5.5689E+01  6.3530E+01 -1.0372E+01  2.1684E+01  3.4539E-01 -3.7245E+00  5.6297E+00
             1.8564E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1789.82940060285        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  1.0095E+00  7.2525E-01  3.6481E-01  1.0811E+00  4.5562E-01  9.4255E-01  1.3979E+00  1.2158E-01  8.9464E-01  4.0331E-01
             2.1938E+00
 PARAMETER:  1.0948E-01 -2.2124E-01 -9.0837E-01  1.7795E-01 -6.8609E-01  4.0835E-02  4.3495E-01 -2.0072E+00 -1.1338E-02 -8.0806E-01
             8.8563E-01
 GRADIENT:  -1.4174E+02  3.7906E+01 -1.8292E+01  3.4069E+01  2.4404E+00 -9.0737E+00  8.0065E+00  2.8228E-01 -4.4856E+00 -5.1474E-01
             4.0853E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1796.54692882139        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      540
 NPARAMETR:  1.0616E+00  5.8076E-01  4.0470E-01  1.1392E+00  4.4086E-01  9.1967E-01  1.7357E+00  9.4550E-02  8.6870E-01  3.3450E-01
             2.2327E+00
 PARAMETER:  1.5974E-01 -4.4341E-01 -8.0461E-01  2.3029E-01 -7.1903E-01  1.6260E-02  6.5141E-01 -2.2586E+00 -4.0756E-02 -9.9511E-01
             9.0321E-01
 GRADIENT:  -1.4121E+01  2.1579E+01  1.5795E+01 -3.2511E+00 -2.0520E+01 -8.7308E+00  1.0335E+01  5.9577E-02  4.2620E-01 -5.2903E+00
             9.1660E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1799.86542291907        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0706E+00  4.5885E-01  4.1919E-01  1.1938E+00  4.2090E-01  9.4201E-01  1.8256E+00  3.6580E-02  8.5841E-01  5.0499E-01
             2.1372E+00
 PARAMETER:  1.6820E-01 -6.7904E-01 -7.6943E-01  2.7718E-01 -7.6537E-01  4.0264E-02  7.0193E-01 -3.2082E+00 -5.2669E-02 -5.8322E-01
             8.5952E-01
 GRADIENT:   1.5110E+01  4.3051E+00  7.2020E+00 -6.3563E-01 -1.1953E+01  5.3132E-01 -7.3190E-01  2.3703E-02 -1.3592E+00 -3.7567E-01
            -1.5037E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1800.32729109688        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  1.0616E+00  3.6792E-01  4.7076E-01  1.2564E+00  4.3426E-01  9.3614E-01  2.2011E+00  2.1546E-02  8.3526E-01  5.5833E-01
             2.1765E+00
 PARAMETER:  1.5980E-01 -8.9988E-01 -6.5340E-01  3.2827E-01 -7.3412E-01  3.4010E-02  8.8895E-01 -3.7376E+00 -8.0006E-02 -4.8280E-01
             8.7772E-01
 GRADIENT:   2.5835E+00  2.6334E+00  4.4236E+00  2.7910E-01 -8.1355E+00  1.6471E-01  1.2281E+00  8.1622E-03 -5.3750E-01  1.6629E-01
            -1.5899E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1800.37212976386        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  1.0577E+00  3.1583E-01  5.0821E-01  1.2955E+00  4.4804E-01  9.3296E-01  2.4157E+00  1.4754E-02  8.2430E-01  5.9523E-01
             2.1853E+00
 PARAMETER:  1.5608E-01 -1.0525E+00 -5.7686E-01  3.5891E-01 -7.0288E-01  3.0606E-02  9.8198E-01 -4.1162E+00 -9.3222E-02 -4.1880E-01
             8.8177E-01
 GRADIENT:  -7.6530E-01  4.6440E-01  8.0672E-01  1.0067E+00 -1.1302E+00 -4.7014E-02  2.4812E-02  3.6752E-03 -3.3020E-02 -5.0263E-02
            -5.8247E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1800.37396286746        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  1.0583E+00  3.1476E-01  5.0792E-01  1.2945E+00  4.4819E-01  9.3311E-01  2.4177E+00  1.0000E-02  8.2442E-01  5.9563E-01
             2.1854E+00
 PARAMETER:  1.5665E-01 -1.0559E+00 -5.7743E-01  3.5812E-01 -7.0254E-01  3.0769E-02  9.8284E-01 -5.1233E+00 -9.3070E-02 -4.1814E-01
             8.8180E-01
 GRADIENT:   8.1818E-01 -2.9382E-01 -8.8844E-01 -1.5813E+00  1.8880E+00  3.7594E-02 -6.6684E-02  0.0000E+00  6.9525E-02  5.3063E-02
             1.3688E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1162
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6529E-03  3.8310E-02 -2.6411E-04 -2.3234E-02  9.4016E-03
 SE:             2.9458E-02  1.8006E-02  2.0571E-04  2.6133E-02  1.8082E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5525E-01  3.3372E-02  1.9919E-01  3.7397E-01  6.0309E-01

 ETASHRINKSD(%)  1.3125E+00  3.9676E+01  9.9311E+01  1.2451E+01  3.9424E+01
 ETASHRINKVR(%)  2.6078E+00  6.3610E+01  9.9995E+01  2.3352E+01  6.3306E+01
 EBVSHRINKSD(%)  1.6223E+00  4.7713E+01  9.9210E+01  1.1034E+01  3.5775E+01
 EBVSHRINKVR(%)  3.2183E+00  7.2661E+01  9.9994E+01  2.0850E+01  5.8752E+01
 RELATIVEINF(%)  9.5628E+01  7.1167E+00  3.5834E-04  2.7008E+01  2.2954E+00
 EPSSHRINKSD(%)  2.7922E+01
 EPSSHRINKVR(%)  4.8048E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1800.3739628674596     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -881.43542966278687     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.82
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1800.374       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.06E+00  3.15E-01  5.08E-01  1.29E+00  4.48E-01  9.33E-01  2.42E+00  1.00E-02  8.24E-01  5.96E-01  2.19E+00
 


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
+        1.10E+03
 
 TH 2
+       -4.98E+01  5.64E+02
 
 TH 3
+        1.64E+01  5.80E+02  2.76E+03
 
 TH 4
+       -1.99E+01  2.96E+02 -4.06E+02  7.79E+02
 
 TH 5
+        2.19E+01 -1.19E+03 -4.03E+03  1.22E+02  6.59E+03
 
 TH 6
+       -2.57E-01 -7.34E+00  4.56E+00 -4.64E+00  5.45E+00  2.14E+02
 
 TH 7
+        1.10E+00  3.82E+01 -6.21E-01 -2.52E+00 -1.08E+01  7.82E-02  8.39E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.49E+00 -1.85E+01  8.20E-01 -3.34E+00  4.30E+01 -2.71E+00  4.13E+00  0.00E+00  1.87E+02
 
 TH10
+       -5.57E-01  2.53E+01 -1.23E+02 -1.53E+01  7.34E+01  5.02E-01  3.94E+00  0.00E+00 -7.29E+00  1.13E+02
 
 TH11
+       -1.27E+01  3.11E+00 -3.16E+01 -2.03E+01  3.38E+01  3.65E+00  2.08E+00  0.00E+00  1.21E+01  2.26E+01  9.31E+01
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       24.995
Stop Time:
Wed Sep 29 22:03:05 CDT 2021
