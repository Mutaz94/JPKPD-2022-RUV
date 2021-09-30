Wed Sep 29 11:40:38 CDT 2021
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
$DATA ../../../../data/spa/B/dat97.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m97.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1645.14302947729        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.5573E+02 -3.2303E+01 -5.1506E+01  3.8111E+01  6.4280E+01  3.2783E+01  1.9479E+00  1.3476E+01  2.5499E+01  1.7980E+01
             2.1937E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1653.37432808857        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:      173
 NPARAMETR:  8.7113E-01  1.0813E+00  1.0961E+00  9.8788E-01  1.0160E+00  1.1067E+00  1.0012E+00  9.5232E-01  9.3212E-01  9.2907E-01
             9.4312E-01
 PARAMETER: -3.7962E-02  1.7820E-01  1.9179E-01  8.7810E-02  1.1583E-01  2.0142E-01  1.0123E-01  5.1145E-02  2.9707E-02  2.6430E-02
             4.1434E-02
 GRADIENT:  -1.1626E+02  8.5225E+00  4.3582E-01  1.1390E+01 -6.4691E+00  1.7673E+01 -3.3956E+00  4.4397E+00  2.0715E+00  1.8750E+00
            -7.5010E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1656.68210624428        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      351
 NPARAMETR:  8.7928E-01  9.4004E-01  8.9232E-01  1.0748E+00  8.8022E-01  1.0606E+00  1.3471E+00  4.4304E-01  8.0670E-01  8.4910E-01
             9.2781E-01
 PARAMETER: -2.8651E-02  3.8164E-02 -1.3930E-02  1.7215E-01 -2.7580E-02  1.5888E-01  3.9798E-01 -7.1409E-01 -1.1480E-01 -6.3573E-02
             2.5073E-02
 GRADIENT:  -1.1005E+02  1.8606E+01 -1.9991E+01  3.8414E+01  2.4952E+01  1.9848E+00  8.9082E+00  1.9277E+00  6.1587E-01  6.4520E+00
            -1.0636E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1662.19915144499        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      528
 NPARAMETR:  9.3512E-01  7.8920E-01  8.1286E-01  1.1434E+00  7.6778E-01  1.0447E+00  1.4687E+00  2.6198E-01  7.5668E-01  7.2946E-01
             9.5369E-01
 PARAMETER:  3.2924E-02 -1.3673E-01 -1.0720E-01  2.3402E-01 -1.6425E-01  1.4376E-01  4.8436E-01 -1.2395E+00 -1.7882E-01 -2.1545E-01
             5.2586E-02
 GRADIENT:   1.0492E+01  9.2275E+00 -1.1096E+00  1.6501E+01  3.8409E-01  1.9078E+00 -9.8935E-02  5.9035E-01 -1.2774E+00  2.2283E-01
             1.2248E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1662.68914445381        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      704
 NPARAMETR:  9.2910E-01  6.6737E-01  7.8421E-01  1.1998E+00  7.1106E-01  1.0364E+00  1.6688E+00  1.3040E-01  7.2975E-01  7.0552E-01
             9.4818E-01
 PARAMETER:  2.6464E-02 -3.0441E-01 -1.4308E-01  2.8212E-01 -2.4100E-01  1.3579E-01  6.1213E-01 -1.9372E+00 -2.1505E-01 -2.4882E-01
             4.6789E-02
 GRADIENT:  -1.4834E+00  8.3502E-01  7.6795E-01  2.5389E-01 -1.7950E+00 -1.2288E+00 -6.9628E-02  1.3173E-01  6.2535E-02 -1.4942E-01
             3.5398E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1662.73189160596        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      861
 NPARAMETR:  9.3025E-01  6.5061E-01  7.8882E-01  1.2078E+00  7.0907E-01  1.0427E+00  1.6998E+00  5.2737E-02  7.2625E-01  7.1155E-01
             9.4705E-01
 PARAMETER:  2.7700E-02 -3.2984E-01 -1.3721E-01  2.8884E-01 -2.4381E-01  1.4183E-01  6.3052E-01 -2.8424E+00 -2.1986E-01 -2.4031E-01
             4.5591E-02
 GRADIENT:   4.1855E+02  5.1941E+01  8.1318E+00  3.1602E+02  3.1114E+01  1.0135E+02  2.5399E+01  4.7109E-02  1.3956E+01  1.0275E+00
             6.1709E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1662.73837484801        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1039
 NPARAMETR:  9.2931E-01  6.5195E-01  7.8769E-01  1.2088E+00  7.0890E-01  1.0390E+00  1.6960E+00  2.9842E-02  7.2640E-01  7.1585E-01
             9.4687E-01
 PARAMETER:  2.6685E-02 -3.2778E-01 -1.3865E-01  2.8962E-01 -2.4404E-01  1.3827E-01  6.2828E-01 -3.4119E+00 -2.1965E-01 -2.3428E-01
             4.5404E-02
 GRADIENT:  -6.7983E-01  1.0632E-01  1.8166E-02  4.6721E-01 -4.2651E-02 -1.9837E-01 -2.0652E-02  6.1231E-03 -7.7910E-03  1.5076E-01
            -8.2725E-03

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1662.74678209104        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1213
 NPARAMETR:  9.3033E-01  6.5288E-01  7.8601E-01  1.2067E+00  7.0858E-01  1.0411E+00  1.6969E+00  1.0000E-02  7.2682E-01  7.1407E-01
             9.4702E-01
 PARAMETER:  2.7786E-02 -3.2636E-01 -1.4078E-01  2.8788E-01 -2.4449E-01  1.4029E-01  6.2880E-01 -4.7028E+00 -2.1908E-01 -2.3678E-01
             4.5563E-02
 GRADIENT:   1.5338E+00 -6.7489E-01  1.5175E-01 -2.3486E+00  1.5416E-01  6.0994E-01  2.1223E-01  0.0000E+00  1.1555E-01  1.5351E-01
             7.4963E-02

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1662.74678209104        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  9.3033E-01  6.5288E-01  7.8601E-01  1.2067E+00  7.0858E-01  1.0411E+00  1.6969E+00  1.0000E-02  7.2682E-01  7.1407E-01
             9.4702E-01
 PARAMETER:  2.7786E-02 -3.2636E-01 -1.4078E-01  2.8788E-01 -2.4449E-01  1.4029E-01  6.2880E-01 -4.7028E+00 -2.1908E-01 -2.3678E-01
             4.5563E-02
 GRADIENT:   1.5338E+00 -6.7489E-01  1.5175E-01 -2.3486E+00  1.5416E-01  6.0994E-01  2.1223E-01  0.0000E+00  1.1555E-01  1.5351E-01
             7.4963E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1235
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5506E-04  1.8857E-02 -5.0301E-04 -1.8501E-02  2.0765E-03
 SE:             2.9860E-02  2.1546E-02  2.2781E-04  2.4327E-02  2.2574E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9318E-01  3.8146E-01  2.7240E-02  4.4696E-01  9.2671E-01

 ETASHRINKSD(%)  1.0000E-10  2.7818E+01  9.9237E+01  1.8501E+01  2.4374E+01
 ETASHRINKVR(%)  1.0000E-10  4.7897E+01  9.9994E+01  3.3579E+01  4.2807E+01
 EBVSHRINKSD(%)  3.6122E-01  2.8648E+01  9.9263E+01  1.7818E+01  2.2755E+01
 EBVSHRINKVR(%)  7.2113E-01  4.9088E+01  9.9995E+01  3.2462E+01  4.0332E+01
 RELATIVEINF(%)  9.8880E+01  5.9375E+00  4.6709E-04  9.5396E+00  4.2254E+00
 EPSSHRINKSD(%)  4.3529E+01
 EPSSHRINKVR(%)  6.8110E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1662.7467820910413     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -927.59595552730309     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.52
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1662.747       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.30E-01  6.53E-01  7.86E-01  1.21E+00  7.09E-01  1.04E+00  1.70E+00  1.00E-02  7.27E-01  7.14E-01  9.47E-01
 


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
+        1.18E+03
 
 TH 2
+       -9.13E+00  4.84E+02
 
 TH 3
+        1.97E+01  2.92E+02  1.05E+03
 
 TH 4
+       -7.19E+00  3.65E+02 -3.89E+02  9.90E+02
 
 TH 5
+       -5.63E+00 -5.81E+02 -1.51E+03  4.86E+02  2.56E+03
 
 TH 6
+        1.51E-01 -1.74E+00  4.30E+00 -2.63E+00 -2.40E+00  1.81E+02
 
 TH 7
+        1.29E+00  3.79E+01 -3.20E+00 -1.00E+01 -2.36E+00  6.72E-02  2.42E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.12E+00 -2.33E+01 -3.75E+01  7.63E+00  2.64E+01 -1.97E-02  1.33E+01  0.00E+00  1.92E+02
 
 TH10
+       -2.07E+00 -5.54E-01 -1.02E+02 -4.66E+01 -2.39E+01 -4.87E-02  1.13E+01  0.00E+00  1.24E+01  1.42E+02
 
 TH11
+       -7.09E+00 -1.23E+01 -4.57E+01 -4.58E+00  1.88E+01  1.57E+00  2.39E+00  0.00E+00  1.17E+01  3.13E+01  2.45E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.623
Stop Time:
Wed Sep 29 11:41:02 CDT 2021
