Sat Sep 25 07:35:42 CDT 2021
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
$DATA ../../../../data/spa/B/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1698.24229520817        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -1.9777E+01 -5.6765E+01 -5.0352E+01 -1.0144E+01  7.8186E+01  9.1591E+00  6.7487E+00  8.3360E+00  1.9774E+01  1.6474E+00
            -1.7826E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1704.47782952364        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  1.0290E+00  1.0473E+00  1.0614E+00  9.9162E-01  9.8454E-01  9.7077E-01  9.5531E-01  9.6744E-01  9.2485E-01  9.6612E-01
             1.0604E+00
 PARAMETER:  1.2860E-01  1.4625E-01  1.5959E-01  9.1587E-02  8.4420E-02  7.0332E-02  5.4284E-02  6.6897E-02  2.1881E-02  6.5532E-02
             1.5860E-01
 GRADIENT:  -3.2607E+00 -8.0967E+00 -1.7176E+00 -1.0833E+01  3.6386E+00 -5.3330E+00  3.5708E+00  2.2355E+00 -4.0407E-01 -1.1176E+00
             3.5856E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1704.76732226584        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      311
 NPARAMETR:  1.0289E+00  1.1607E+00  1.0271E+00  9.2768E-01  1.0220E+00  9.7015E-01  8.1482E-01  8.2192E-01  1.0181E+00  1.0144E+00
             1.0383E+00
 PARAMETER:  1.2852E-01  2.4903E-01  1.2671E-01  2.4936E-02  1.2180E-01  6.9693E-02 -1.0479E-01 -9.6107E-02  1.1790E-01  1.1432E-01
             1.3758E-01
 GRADIENT:  -4.6753E+00  3.9048E+00  5.9578E+00  4.4187E+00  2.1378E+00 -5.9854E+00 -6.7601E-02 -2.2067E+00  1.3567E+00 -2.3400E+00
            -6.6431E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1705.31620001214        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  1.0340E+00  1.2550E+00  8.3969E-01  8.6164E-01  9.7760E-01  9.9247E-01  8.0644E-01  6.9867E-01  1.0357E+00  9.5520E-01
             1.0558E+00
 PARAMETER:  1.3345E-01  3.2710E-01 -7.4721E-02 -4.8918E-02  7.7341E-02  9.2444E-02 -1.1513E-01 -2.5858E-01  1.3511E-01  5.4163E-02
             1.5427E-01
 GRADIENT:   2.9075E+00  5.1052E+00  7.1790E-01  6.2788E+00 -3.4066E+00  2.2892E+00 -7.8309E-01  1.9097E-01 -1.4318E+00  2.6276E-01
             1.8762E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1705.95810287750        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      662
 NPARAMETR:  1.0338E+00  1.5846E+00  6.0850E-01  6.3684E-01  1.0475E+00  9.9194E-01  7.0876E-01  4.5778E-01  1.2867E+00  9.7567E-01
             1.0534E+00
 PARAMETER:  1.3321E-01  5.6031E-01 -3.9676E-01 -3.5124E-01  1.4642E-01  9.1907E-02 -2.4424E-01 -6.8136E-01  3.5211E-01  7.5369E-02
             1.5199E-01
 GRADIENT:  -2.0471E-01 -4.3339E+00 -1.2290E+00  1.7999E-01  1.0018E-01  1.1904E+00  1.8615E-01  3.0617E-01  1.6388E+00  1.8485E+00
             1.4925E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1706.04191747032        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      837
 NPARAMETR:  1.0338E+00  1.6787E+00  5.3825E-01  5.7338E-01  1.0686E+00  9.8951E-01  6.9280E-01  3.1927E-01  1.3580E+00  9.7096E-01
             1.0516E+00
 PARAMETER:  1.3322E-01  6.1801E-01 -5.1942E-01 -4.5620E-01  1.6639E-01  8.9457E-02 -2.6701E-01 -1.0417E+00  4.0601E-01  7.0531E-02
             1.5033E-01
 GRADIENT:  -5.7051E-01 -3.4793E+00 -1.3799E+00 -5.0377E-01  2.1556E+00  4.9483E-02  2.6608E-01  1.9398E-01  9.5919E-02  4.4791E-01
             3.5293E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1706.11250085133        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  1.0338E+00  1.6498E+00  5.3948E-01  5.9284E-01  1.0498E+00  9.8933E-01  7.0060E-01  9.7641E-02  1.3238E+00  9.5534E-01
             1.0509E+00
 PARAMETER:  1.3326E-01  6.0067E-01 -5.1716E-01 -4.2282E-01  1.4856E-01  8.9275E-02 -2.5582E-01 -2.2265E+00  3.8050E-01  5.4313E-02
             1.4961E-01
 GRADIENT:  -6.3511E-01 -4.5438E-01 -7.5061E-01  1.1536E+00  2.3192E+00 -4.6219E-02 -1.5715E-01  1.0948E-02 -3.0850E-01 -1.7376E-01
             2.8254E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1706.12237480265        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1188
 NPARAMETR:  1.0341E+00  1.6661E+00  5.3159E-01  5.8155E-01  1.0545E+00  9.8945E-01  6.9690E-01  2.3702E-02  1.3432E+00  9.5856E-01
             1.0508E+00
 PARAMETER:  1.3351E-01  6.1047E-01 -5.3189E-01 -4.4205E-01  1.5303E-01  8.9392E-02 -2.6112E-01 -3.6422E+00  3.9505E-01  5.7679E-02
             1.4959E-01
 GRADIENT:  -3.8061E-02 -1.0074E-01 -2.9542E-02  1.4217E-02  6.1676E-02 -9.5238E-03  1.4524E-02  7.0698E-04  1.1926E-02  8.9657E-03
             2.6967E-03

0ITERATION NO.:   39    OBJECTIVE VALUE:  -1706.12264288885        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1315
 NPARAMETR:  1.0341E+00  1.6673E+00  5.3096E-01  5.8076E-01  1.0549E+00  9.8947E-01  6.9660E-01  1.0000E-02  1.3443E+00  9.5872E-01
             1.0509E+00
 PARAMETER:  1.3352E-01  6.1120E-01 -5.3306E-01 -4.4341E-01  1.5340E-01  8.9417E-02 -2.6154E-01 -4.6586E+00  3.9590E-01  5.7845E-02
             1.4961E-01
 GRADIENT:  -7.7728E-03 -2.0538E-02 -1.5582E-02  8.0255E-03  1.5442E-02 -9.7324E-04  1.0193E-02  0.0000E+00  3.5998E-03  3.8932E-03
             4.7720E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1315
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7208E-04 -3.1849E-02 -2.7382E-04  2.3527E-02 -3.3469E-02
 SE:             2.9825E-02  2.2181E-02  1.1053E-04  2.3505E-02  2.2903E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9540E-01  1.5105E-01  1.3238E-02  3.1685E-01  1.4393E-01

 ETASHRINKSD(%)  8.0945E-02  2.5691E+01  9.9630E+01  2.1257E+01  2.3271E+01
 ETASHRINKVR(%)  1.6182E-01  4.4781E+01  9.9999E+01  3.7995E+01  4.1126E+01
 EBVSHRINKSD(%)  4.8768E-01  2.5057E+01  9.9683E+01  2.2785E+01  2.1793E+01
 EBVSHRINKVR(%)  9.7299E-01  4.3836E+01  9.9999E+01  4.0379E+01  3.8837E+01
 RELATIVEINF(%)  9.8976E+01  3.2995E+00  1.1771E-04  3.8010E+00  1.1147E+01
 EPSSHRINKSD(%)  4.3947E+01
 EPSSHRINKVR(%)  6.8580E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1706.1226428888492     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -970.97181632511104     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.21
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.94
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1706.123       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.67E+00  5.31E-01  5.81E-01  1.05E+00  9.89E-01  6.97E-01  1.00E-02  1.34E+00  9.59E-01  1.05E+00
 


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
+        1.05E+03
 
 TH 2
+       -7.18E+00  4.87E+02
 
 TH 3
+        8.63E+00  1.93E+02  4.37E+02
 
 TH 4
+       -1.90E+01  4.12E+02 -3.03E+02  1.10E+03
 
 TH 5
+       -4.01E+00 -2.58E+02 -4.04E+02  2.70E+02  6.63E+02
 
 TH 6
+        4.64E-01 -1.30E+00  1.62E+00 -4.98E+00 -3.09E+00  2.01E+02
 
 TH 7
+        1.68E+00  3.39E+00 -7.22E+00 -1.36E+01 -1.67E+01 -1.13E+00  1.47E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.19E+00 -2.15E+01 -3.30E+01  5.67E+01  9.28E-01 -1.94E-01  2.39E+01  0.00E+00  4.76E+01
 
 TH10
+       -3.47E+00 -1.88E+01 -4.26E+01 -2.73E-01 -6.35E+01  1.06E+00  1.21E+01  0.00E+00  6.81E+00  8.90E+01
 
 TH11
+       -7.39E+00 -1.99E+01 -2.53E+01 -1.50E+00  3.13E-02  2.04E+00  1.16E+01  0.00E+00  5.56E+00  1.99E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       22.231
Stop Time:
Sat Sep 25 07:36:05 CDT 2021
