Sat Sep 25 02:22:48 CDT 2021
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
$DATA ../../../../data/int/SL3/dat52.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 NO. OF DATA RECS IN DATA SET:      981
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      881
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -81.3707302167559        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.6875E+01 -1.1145E+01  7.9848E+01 -9.6786E+01  3.0215E+01  8.9290E-01 -1.0841E+02 -1.6301E+02 -4.3683E+01 -6.1415E+00
            -7.1385E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2660.49700285866        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0281E+00  1.3551E+00  1.0255E+00  9.3394E-01  1.2197E+00  9.2592E-01  1.2735E+00  9.2319E-01  7.3225E-01  8.4981E-01
             3.0501E+00
 PARAMETER:  1.2776E-01  4.0391E-01  1.2516E-01  3.1662E-02  2.9857E-01  2.3030E-02  3.4176E-01  2.0083E-02 -2.1163E-01 -6.2740E-02
             1.2152E+00
 GRADIENT:   4.6517E+00  4.4182E+01 -1.9847E+01  4.9966E+01  2.5354E+01 -1.9598E+01  1.3776E+01  3.6627E+00  8.1047E-01 -1.3097E+01
             1.4939E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2663.77666859379        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0321E+00  1.3965E+00  1.2454E+00  9.1778E-01  1.3358E+00  9.4482E-01  1.1985E+00  6.2196E-01  6.4839E-01  1.1676E+00
             3.0469E+00
 PARAMETER:  1.3162E-01  4.3396E-01  3.1947E-01  1.4202E-02  3.8951E-01  4.3235E-02  2.8108E-01 -3.7488E-01 -3.3327E-01  2.5499E-01
             1.2141E+00
 GRADIENT:   1.4044E+01  4.8226E+01 -8.8380E+00  6.9924E+01  9.8635E+00 -1.1680E+01  3.3107E+00  1.2633E+00 -2.1425E+00  1.2506E+01
             1.5403E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2672.20992024635        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0193E+00  1.3614E+00  1.2029E+00  8.8727E-01  1.3189E+00  9.7159E-01  1.1723E+00  3.3157E-01  7.2407E-01  1.0916E+00
             2.8192E+00
 PARAMETER:  1.1916E-01  4.0850E-01  2.8473E-01 -1.9604E-02  3.7676E-01  7.1182E-02  2.5900E-01 -1.0039E+00 -2.2286E-01  1.8764E-01
             1.1364E+00
 GRADIENT:  -2.9494E+00 -7.4263E-01 -9.7955E-01  1.9371E+00  2.4072E+00 -2.1711E+00  1.4744E+00  2.8739E-01  1.6486E-01 -7.4550E-01
            -7.0721E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2672.89744241867        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  1.0218E+00  1.4672E+00  1.2652E+00  8.2858E-01  1.4234E+00  9.7954E-01  1.1034E+00  1.2549E-01  7.1819E-01  1.1855E+00
             2.8298E+00
 PARAMETER:  1.2153E-01  4.8335E-01  3.3522E-01 -8.8044E-02  4.5307E-01  7.9332E-02  1.9843E-01 -1.9755E+00 -2.3101E-01  2.7019E-01
             1.1402E+00
 GRADIENT:   1.9737E+00  6.0244E+00 -1.3871E+00  8.8866E+00  2.5508E+00  8.0967E-01 -1.3181E+00  3.2912E-02 -2.5325E-01  1.9225E+00
             3.6440E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2673.44614637220        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      436
 NPARAMETR:  1.0276E+00  1.6280E+00  1.3963E+00  7.3700E-01  1.6017E+00  9.7827E-01  1.0426E+00  3.8992E-02  7.0990E-01  1.2846E+00
             2.8225E+00
 PARAMETER:  1.2719E-01  5.8736E-01  4.3381E-01 -2.0517E-01  5.7106E-01  7.8030E-02  1.4169E-01 -3.1444E+00 -2.4263E-01  3.5044E-01
             1.1376E+00
 GRADIENT:   7.9228E+00  1.0891E+01 -2.5615E+00  1.5756E+01  7.2812E+00 -5.0055E-01  1.8954E-02  1.7345E-03 -3.4286E-01  4.4845E-02
            -2.2254E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2673.94993702504        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      611
 NPARAMETR:  1.0243E+00  1.8432E+00  1.5753E+00  5.8941E-01  1.8240E+00  9.7975E-01  9.5344E-01  1.0000E-02  7.0188E-01  1.4202E+00
             2.8180E+00
 PARAMETER:  1.2404E-01  7.1153E-01  5.5447E-01 -4.2863E-01  7.0105E-01  7.9538E-02  5.2325E-02 -7.0054E+00 -2.5399E-01  4.5078E-01
             1.1360E+00
 GRADIENT:   1.5514E+00 -6.1547E-01  1.9185E-02 -4.2004E-01 -6.5052E-02 -8.5112E-02  2.3637E-01  0.0000E+00  5.3686E-02 -6.6857E-02
             3.4550E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2673.95078719079        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      786
 NPARAMETR:  1.0236E+00  1.8452E+00  1.5783E+00  5.8851E-01  1.8262E+00  9.8001E-01  9.5237E-01  1.0000E-02  6.9900E-01  1.4219E+00
             2.8180E+00
 PARAMETER:  1.2333E-01  7.1257E-01  5.5637E-01 -4.3015E-01  7.0224E-01  7.9811E-02  5.1203E-02 -7.0271E+00 -2.5811E-01  4.5199E-01
             1.1360E+00
 GRADIENT:  -5.5100E-02 -3.8314E-02 -1.3680E-02 -2.0657E-02  1.0770E-02  3.3070E-03 -2.6031E-03  0.0000E+00 -3.8982E-05  7.0775E-04
             8.7481E-03

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2673.95081641477        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      949
 NPARAMETR:  1.0236E+00  1.8444E+00  1.5847E+00  5.8909E-01  1.8269E+00  9.8000E-01  9.5283E-01  1.0000E-02  6.9817E-01  1.4223E+00
             2.8180E+00
 PARAMETER:  1.2336E-01  7.1215E-01  5.6042E-01 -4.2918E-01  7.0262E-01  7.9801E-02  5.1682E-02 -7.0083E+00 -2.5930E-01  4.5230E-01
             1.1360E+00
 GRADIENT:   6.5134E-03  4.3646E-03 -4.1204E-04  2.4873E-03  2.2661E-03 -3.4617E-04  1.9686E-03  0.0000E+00  7.2898E-04 -1.6643E-04
            -8.3437E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      949
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3013E-03 -3.3677E-03 -5.6029E-05 -1.1315E-02 -1.3928E-02
 SE:             2.9275E-02  2.6749E-02  3.7968E-05  1.0700E-02  2.5477E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6454E-01  8.9981E-01  1.4003E-01  2.9031E-01  5.8460E-01

 ETASHRINKSD(%)  1.9244E+00  1.0387E+01  9.9873E+01  6.4152E+01  1.4649E+01
 ETASHRINKVR(%)  3.8118E+00  1.9696E+01  1.0000E+02  8.7149E+01  2.7151E+01
 EBVSHRINKSD(%)  1.8882E+00  1.0439E+01  9.9876E+01  6.4725E+01  1.3566E+01
 EBVSHRINKVR(%)  3.7406E+00  1.9789E+01  1.0000E+02  8.7557E+01  2.5291E+01
 RELATIVEINF(%)  9.6174E+01  6.0186E+00  1.9288E-05  7.7042E-01  1.3460E+01
 EPSSHRINKSD(%)  1.5299E+01
 EPSSHRINKVR(%)  2.8257E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          881
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1619.1696955066332     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2673.9508164147742     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1054.7811209081410     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    11.30
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2673.951       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.84E+00  1.58E+00  5.89E-01  1.83E+00  9.80E-01  9.53E-01  1.00E-02  6.98E-01  1.42E+00  2.82E+00
 


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
+        1.06E+03
 
 TH 2
+       -1.25E+01  2.90E+02
 
 TH 3
+        1.70E+00  2.76E+00  7.22E+00
 
 TH 4
+       -2.55E+01  4.01E+02 -3.18E+01  8.15E+02
 
 TH 5
+       -4.49E+00 -2.76E+01 -2.61E+01  1.04E+02  1.23E+02
 
 TH 6
+        3.42E+00 -3.99E+00  5.77E-01 -8.12E+00 -2.05E+00  1.88E+02
 
 TH 7
+        3.87E+00  1.17E+01 -1.26E+00 -6.38E+01  1.21E+00 -2.34E-01  1.43E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.91E-01 -5.71E+00  9.03E-01 -1.58E+01  9.19E-01 -8.77E-03  1.88E+01  0.00E+00  1.08E+01
 
 TH10
+       -4.44E-01 -1.97E+00 -4.51E+00  2.14E+01 -1.03E+01 -3.35E-01  4.83E-01  0.00E+00  2.54E+00  5.35E+01
 
 TH11
+       -1.63E+01 -1.18E+01 -6.52E-01 -1.86E+01  1.77E+00  3.23E+00  6.21E+00  0.00E+00  3.33E+00  7.39E+00  1.47E+02
 
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
 #CPUT: Total CPU Time in Seconds,       30.890
Stop Time:
Sat Sep 25 02:23:21 CDT 2021
