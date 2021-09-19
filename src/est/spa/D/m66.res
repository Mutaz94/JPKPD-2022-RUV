Sat Sep 18 15:32:20 CDT 2021
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
$DATA ../../../../data/spa/D/dat66.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m66.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   16857.1861389260        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.0180E+02  3.3059E+02 -2.2244E+01  1.0926E+02  1.8294E+02 -2.5183E+03 -9.3647E+02 -9.5576E+01 -1.5750E+03 -8.5398E+02
            -3.0967E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -516.770062460791        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2718E+00  9.9455E-01  8.0345E-01  1.6337E+00  1.6512E+00  2.0267E+00  1.1856E+00  9.7496E-01  9.6709E-01  1.1057E+00
             1.4743E+01
 PARAMETER:  3.4046E-01  9.4531E-02 -1.1884E-01  5.9083E-01  6.0152E-01  8.0642E-01  2.7028E-01  7.4643E-02  6.6532E-02  2.0047E-01
             2.7907E+00
 GRADIENT:  -1.7326E+00  1.8499E+01 -3.9383E+00  2.5955E+01 -5.5501E+00  4.1808E+01  8.0385E-02  5.5066E+00  5.8071E+00  2.1005E-01
             4.2927E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -533.630238721589        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.2147E+00  9.7825E-01  6.3004E-01  1.4732E+00  3.8605E+00  1.6181E+00  1.7268E+00  1.8086E-01  5.4423E-01  7.9775E+00
             1.4140E+01
 PARAMETER:  2.9447E-01  7.8014E-02 -3.6197E-01  4.8743E-01  1.4508E+00  5.8124E-01  6.4628E-01 -1.6101E+00 -5.0838E-01  2.1766E+00
             2.7490E+00
 GRADIENT:   1.1694E+01  2.2664E+01 -2.2813E+01  3.7262E+01 -1.8322E+01 -4.4108E+00  4.0407E+00  1.7649E-01  4.5670E+00  8.9658E+00
             1.2569E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -558.104928600431        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  1.0659E+00  4.7697E-01  6.3493E-01  1.3932E+00  1.1907E+01  1.3913E+00  1.7797E+00  2.2558E-01  3.8413E-01  1.1029E+01
             1.3848E+01
 PARAMETER:  1.6383E-01 -6.4030E-01 -3.5423E-01  4.3164E-01  2.5771E+00  4.3026E-01  6.7642E-01 -1.3891E+00 -8.5678E-01  2.5005E+00
             2.7281E+00
 GRADIENT:  -2.8576E+01  1.0046E+01  6.4628E+00 -2.6681E-01 -1.0259E+00 -1.6394E+01  2.8105E+00  8.0353E-02  5.6982E+00 -2.7556E-01
             4.6877E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -585.701184180291        NO. OF FUNC. EVALS.: 100
 CUMULATIVE NO. OF FUNC. EVALS.:      325
 NPARAMETR:  1.0232E+00  1.8577E-01  4.9319E-01  1.5764E+00  2.4793E+01  1.4681E+00  2.2688E+00  6.2179E-02  1.9439E-01  1.4026E+01
             1.3535E+01
 PARAMETER:  1.2292E-01 -1.5832E+00 -6.0685E-01  5.5512E-01  3.3106E+00  4.8398E-01  9.1924E-01 -2.6777E+00 -1.5379E+00  2.7409E+00
             2.7053E+00
 GRADIENT:  -4.5618E+01  6.2626E+00  1.9422E+01  8.2512E+01  6.5953E+00 -3.7062E+01  4.3180E-01  1.1121E-02  1.3947E+00 -1.8535E+01
             9.5920E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -589.510539985243        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:      426
 NPARAMETR:  1.0700E+00  1.8581E-01  4.9299E-01  1.5761E+00  2.4749E+01  1.6047E+00  2.1879E+00  2.4773E-02  1.0000E-02  1.4041E+01
             1.3526E+01
 PARAMETER:  1.6767E-01 -1.5830E+00 -6.0727E-01  5.5498E-01  3.3088E+00  5.7294E-01  8.8294E-01 -3.5980E+00 -4.7595E+00  2.7420E+00
             2.7046E+00
 GRADIENT:   1.0500E+00  4.1668E+00  2.7591E+01  4.3238E+01  8.2735E+00 -1.6223E+00  4.7225E-01  1.6218E-03  0.0000E+00 -6.9875E+00
             5.0582E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -589.972687413208        NO. OF FUNC. EVALS.:  81
 CUMULATIVE NO. OF FUNC. EVALS.:      507
 NPARAMETR:  1.0431E+00  1.8541E-01  4.9111E-01  1.5708E+00  2.4391E+01  1.6671E+00  1.3278E+00  1.0000E-02  1.0000E-02  1.4139E+01
             1.3483E+01
 PARAMETER:  1.4221E-01 -1.5852E+00 -6.1108E-01  5.5159E-01  3.2942E+00  6.1111E-01  3.8354E-01 -1.2563E+01 -2.5413E+01  2.7489E+00
             2.7014E+00
 GRADIENT:  -7.7476E+00  3.9478E+00  1.0571E+01  1.9778E+01  6.1619E+00  4.4658E+00  1.7537E-01  0.0000E+00  0.0000E+00  5.2597E+00
             2.2275E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -590.013403151172        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      582
 NPARAMETR:  1.0533E+00  1.8546E-01  4.9137E-01  1.5715E+00  2.4443E+01  1.6424E+00  1.4353E+00  1.0000E-02  1.0000E-02  1.4119E+01
             1.3484E+01
 PARAMETER:  1.5196E-01 -1.5849E+00 -6.1056E-01  5.5203E-01  3.2963E+00  5.9614E-01  4.6137E-01 -1.1303E+01 -2.2602E+01  2.7475E+00
             2.7015E+00
 GRADIENT:  -2.8084E+00  3.9154E+00  1.3645E+01  2.4567E+01  6.5389E+00  1.3156E+00  2.0534E-01  0.0000E+00  0.0000E+00  3.4034E+00
             1.7030E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -590.020769314289        NO. OF FUNC. EVALS.: 112
 CUMULATIVE NO. OF FUNC. EVALS.:      694
 NPARAMETR:  1.0514E+00  1.8515E-01  4.9159E-01  1.5704E+00  2.4511E+01  1.6368E+00  1.4370E+00  1.0000E-02  1.0000E-02  1.4073E+01
             1.3505E+01
 PARAMETER:  1.5030E-01 -1.5851E+00 -6.1069E-01  5.5188E-01  3.2961E+00  5.9220E-01  4.6209E-01 -1.1470E+01 -2.3027E+01  2.7468E+00
             2.7006E+00
 GRADIENT:   3.3803E+03  1.6458E+02 -8.2439E+02  4.8323E+02 -1.4829E+02 -7.8287E+02 -1.1016E+03  0.0000E+00  0.0000E+00  1.7532E+02
            -3.2714E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      694
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5993E-02 -5.3229E-03 -6.5371E-05 -7.1836E-04 -3.6181E-02
 SE:             2.7404E-02  1.1038E-03  8.2664E-05  2.5303E-04  8.1159E-03
 N:                     100         100         100         100         100

 P VAL.:         5.5949E-01  1.4197E-06  4.2906E-01  4.5248E-03  8.2788E-06

 ETASHRINKSD(%)  8.1932E+00  9.6302E+01  9.9723E+01  9.9152E+01  7.2811E+01
 ETASHRINKVR(%)  1.5715E+01  9.9863E+01  9.9999E+01  9.9993E+01  9.2607E+01
 EBVSHRINKSD(%)  6.6075E+00  9.5251E+01  9.9694E+01  9.9119E+01  7.8353E+01
 EBVSHRINKVR(%)  1.2778E+01  9.9774E+01  9.9999E+01  9.9992E+01  9.5314E+01
 RELATIVEINF(%)  3.4240E+01  2.5033E-02  2.3210E-04  8.0511E-04  4.0724E+00
 EPSSHRINKSD(%)  3.3009E+00
 EPSSHRINKVR(%)  6.4929E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -590.02076931428871     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       145.13005724944946     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.84
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.40
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -590.021       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  1.85E-01  4.91E-01  1.57E+00  2.44E+01  1.64E+00  1.44E+00  1.00E-02  1.00E-02  1.41E+01  1.35E+01
 


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
+        5.09E+06
 
 TH 2
+        2.24E+06  1.47E+06
 
 TH 3
+       -2.18E+06  1.17E+06  1.43E+06
 
 TH 4
+       -1.35E+03 -6.60E+03  4.03E+05  1.67E+05
 
 TH 5
+        1.38E+01 -6.18E+00 -4.31E+03  1.52E+03  1.94E+01
 
 TH 6
+        1.23E+03 -3.65E+05 -1.25E+03 -1.25E+05 -1.35E+03  1.35E+05
 
 TH 7
+        6.45E+01 -5.34E+05  5.21E+05 -2.09E+00  4.92E-02 -5.78E-01  2.89E+05
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -2.93E+01 -5.76E+00 -9.77E+01 -4.78E+01  3.43E+01  1.14E+01 -9.72E-02  0.00E+00  0.00E+00  8.25E+01
 
 TH11
+        1.90E+01 -1.54E+01  9.91E+03  4.03E+01 -2.66E+00  3.04E+03  4.45E+03  0.00E+00  0.00E+00  7.85E+01  1.77E+02
 
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
 #CPUT: Total CPU Time in Seconds,       24.337
Stop Time:
Sat Sep 18 15:32:45 CDT 2021
