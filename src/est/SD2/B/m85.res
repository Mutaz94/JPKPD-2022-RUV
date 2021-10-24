Sat Oct 23 17:08:27 CDT 2021
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
$DATA ../../../../data/SD2/B/dat85.csv ignore=@
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

NM-TRAN MESSAGES
  
 WARNINGS AND ERRORS (IF ANY) FOR PROBLEM    1
             
 (WARNING  2) NM-TRAN INFERS THAT THE DATA ARE POPULATION.

License Registered to: University of Minnesota
Expiration Date:    14 APR 2022
Current Date:       23 OCT 2021
Days until program expires : 176
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
 NO. OF DATA RECS IN DATA SET:      800
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

 TOT. NO. OF OBS RECS:      700
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

 #PARA: PARAFILE=../../../../rmpi.pnm, PROTOCOL=MPI, NODES= 10

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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2947.24846226942        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6231E+02 -5.4895E+01 -1.1565E+00  4.1317E+01  2.7664E+01  7.2189E+01 -1.6745E+01  1.3185E+01 -6.9105E+00  3.0219E+01
            -4.5492E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2958.43672341735        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:      181
 NPARAMETR:  9.8643E-01  1.1336E+00  1.0002E+00  9.9301E-01  1.0344E+00  9.5285E-01  1.0388E+00  9.7280E-01  1.0283E+00  9.4495E-01
             1.0397E+00
 PARAMETER:  8.6340E-02  2.2539E-01  1.0024E-01  9.2985E-02  1.3387E-01  5.1699E-02  1.3809E-01  7.2424E-02  1.2791E-01  4.3381E-02
             1.3892E-01
 GRADIENT:   1.5465E+01  6.2809E+00  3.8730E+00  4.8709E+01 -2.7885E+01  8.8055E+00 -9.7571E+00  1.0375E+01 -8.8350E+00  1.6845E+01
             1.3526E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2958.91579202502        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  9.7476E-01  1.1387E+00  9.9577E-01  9.9747E-01  1.0436E+00  9.3173E-01  1.0494E+00  9.5767E-01  1.0355E+00  9.3647E-01
             1.0359E+00
 PARAMETER:  7.4440E-02  2.2991E-01  9.5766E-02  9.7463E-02  1.4264E-01  2.9290E-02  1.4822E-01  5.6751E-02  1.3487E-01  3.4366E-02
             1.3529E-01
 GRADIENT:  -1.4182E+01  8.1049E+00 -2.1276E+00  6.3189E+01 -1.2318E+01  9.8417E-02 -8.4858E+00  9.8401E+00 -7.4047E+00  1.5235E+01
             7.9961E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2963.74676713522        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  9.7552E-01  1.1371E+00  9.6973E-01  9.4614E-01  1.0574E+00  9.2804E-01  1.1476E+00  6.7963E-01  1.0743E+00  8.2479E-01
             1.0317E+00
 PARAMETER:  7.5210E-02  2.2844E-01  6.9259E-02  4.4639E-02  1.5584E-01  2.5317E-02  2.3770E-01 -2.8621E-01  1.7165E-01 -9.2623E-02
             1.3118E-01
 GRADIENT:  -1.1823E+01 -2.7869E+01  1.1006E+01  1.5887E+00  1.2579E+01 -1.4747E+00 -5.2253E-01  1.0255E+00  8.3232E-01 -1.0453E+00
            -3.9342E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2965.01568932739        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      681
 NPARAMETR:  9.8012E-01  1.1370E+00  8.6679E-01  9.4273E-01  1.0042E+00  9.3172E-01  1.1637E+00  4.3446E-01  1.0773E+00  7.8589E-01
             1.0375E+00
 PARAMETER:  7.9915E-02  2.2841E-01 -4.2955E-02  4.1022E-02  1.0416E-01  2.9279E-02  2.5161E-01 -7.3365E-01  1.7442E-01 -1.4093E-01
             1.3683E-01
 GRADIENT:  -2.6141E-01 -1.4813E+01 -7.4160E-01  2.1334E-01  1.2393E-01  6.5820E-02  6.1077E-01  7.9357E-02  3.1431E-01  2.4694E-01
             8.5886E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2965.18345825422        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      838
 NPARAMETR:  9.8071E-01  1.1535E+00  8.6460E-01  9.3958E-01  1.0091E+00  9.3192E-01  1.1553E+00  4.3685E-01  1.0787E+00  7.8892E-01
             1.0376E+00
 PARAMETER:  8.0521E-02  2.4277E-01 -4.5493E-02  3.7680E-02  1.0905E-01  2.9494E-02  2.4434E-01 -7.2816E-01  1.7572E-01 -1.3709E-01
             1.3695E-01
 GRADIENT:   1.2049E+00 -6.1093E+00 -1.6135E+00  6.0565E+00 -2.1627E+00  1.3001E-01  7.8066E-01  1.3312E-01 -4.1476E-01 -2.7202E-01
             1.7181E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2965.74037372346        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  9.8016E-01  1.2526E+00  8.8450E-01  8.8295E-01  1.0814E+00  9.3155E-01  1.0722E+00  4.8829E-01  1.1299E+00  8.5811E-01
             1.0364E+00
 PARAMETER:  7.9962E-02  3.2519E-01 -2.2733E-02 -2.4486E-02  1.7822E-01  2.9095E-02  1.6970E-01 -6.1685E-01  2.2215E-01 -5.3023E-02
             1.3571E-01
 GRADIENT:  -2.2026E-01  2.9070E-01 -3.1281E-01  3.1871E-01  2.3908E-01 -7.6528E-02  4.8883E-03  1.5763E-02 -1.3256E-02 -1.0454E-01
            -8.1624E-02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -2965.74037372346        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1037
 NPARAMETR:  9.8016E-01  1.2526E+00  8.8450E-01  8.8295E-01  1.0814E+00  9.3155E-01  1.0722E+00  4.8829E-01  1.1299E+00  8.5811E-01
             1.0364E+00
 PARAMETER:  7.9962E-02  3.2519E-01 -2.2733E-02 -2.4486E-02  1.7822E-01  2.9095E-02  1.6970E-01 -6.1685E-01  2.2215E-01 -5.3023E-02
             1.3571E-01
 GRADIENT:  -2.2026E-01  2.9070E-01 -3.1281E-01  3.1871E-01  2.3908E-01 -7.6528E-02  4.8883E-03  1.5763E-02 -1.3256E-02 -1.0454E-01
            -8.1624E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1037
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0249E-03 -1.6635E-02 -1.2001E-02  1.1935E-02 -2.0662E-02
 SE:             2.9897E-02  2.3592E-02  8.1992E-03  2.5465E-02  2.2945E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7265E-01  4.8072E-01  1.4330E-01  6.3930E-01  3.6784E-01

 ETASHRINKSD(%)  1.0000E-10  2.0965E+01  7.2532E+01  1.4689E+01  2.3133E+01
 ETASHRINKVR(%)  1.0000E-10  3.7535E+01  9.2455E+01  2.7220E+01  4.0915E+01
 EBVSHRINKSD(%)  3.3558E-01  2.1183E+01  7.4128E+01  1.5159E+01  2.2918E+01
 EBVSHRINKVR(%)  6.7004E-01  3.7878E+01  9.3307E+01  2.8020E+01  4.0584E+01
 RELATIVEINF(%)  9.9327E+01  1.4537E+01  3.9279E+00  2.0336E+01  1.3147E+01
 EPSSHRINKSD(%)  2.4243E+01
 EPSSHRINKVR(%)  4.2609E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          700
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1286.5139464865417     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2965.7403737234617     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1679.2264272369200     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     8.82
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2965.740       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.25E+00  8.84E-01  8.83E-01  1.08E+00  9.32E-01  1.07E+00  4.88E-01  1.13E+00  8.58E-01  1.04E+00
 


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
 
 Elapsed finaloutput time in seconds:     0.00
 #CPUT: Total CPU Time in Seconds,       64.287
Stop Time:
Sat Oct 23 17:08:39 CDT 2021
