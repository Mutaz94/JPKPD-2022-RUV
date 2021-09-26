Sat Sep 25 12:13:41 CDT 2021
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
$DATA ../../../../data/spa/S2/dat31.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m31.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1653.58265425753        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0094E+02 -4.5733E+01 -4.6174E+01 -3.8630E+00  8.8258E+01  4.8717E-01  6.6299E+00  7.2386E+00  5.7787E+00  2.2244E+00
            -2.3300E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1658.53861473314        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.6282E-01  9.4334E-01  9.8506E-01  1.0344E+00  8.9360E-01  1.0084E+00  8.4346E-01  9.0554E-01  1.0078E+00  8.9350E-01
             1.0289E+00
 PARAMETER:  6.2116E-02  4.1670E-02  8.4943E-02  1.3386E-01 -1.2499E-02  1.0833E-01 -7.0246E-02  7.7817E-04  1.0773E-01 -1.2612E-02
             1.2854E-01
 GRADIENT:   1.5204E+01 -9.4898E+00  1.8814E+00 -1.2545E+01 -6.7764E+00  5.8434E+00  5.4314E-01  3.1118E+00  4.3463E+00 -1.2800E+00
             8.6039E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1659.05912289188        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:      248
 NPARAMETR:  9.6974E-01  9.4915E-01  9.6978E-01  1.0438E+00  8.9088E-01  9.9962E-01  7.9005E-01  7.7463E-01  1.0125E+00  9.5118E-01
             1.0288E+00
 PARAMETER:  6.9271E-02  4.7808E-02  6.9314E-02  1.4286E-01 -1.5546E-02  9.9617E-02 -1.3566E-01 -1.5536E-01  1.1245E-01  4.9948E-02
             1.2836E-01
 GRADIENT:  -9.3884E+00  1.7513E+00  1.6392E+00  1.2314E+00 -1.1679E+01 -1.9403E+00 -4.5766E-01  1.6509E+00  2.1801E+00  3.5435E+00
             8.2539E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1659.74956186091        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  9.7781E-01  1.0504E+00  8.1249E-01  9.7278E-01  8.6907E-01  1.0071E+00  9.3790E-01  5.2760E-01  1.0036E+00  8.7207E-01
             1.0030E+00
 PARAMETER:  7.7562E-02  1.4912E-01 -1.0765E-01  7.2400E-02 -4.0329E-02  1.0712E-01  3.5886E-02 -5.3942E-01  1.0360E-01 -3.6887E-02
             1.0302E-01
 GRADIENT:   5.8891E+00 -2.3761E+00 -3.5221E+00  3.3578E+00  7.4504E+00  3.7921E-01  8.6805E-01  7.5850E-01  4.1341E-01 -1.5240E+00
            -2.5021E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1660.18490783962        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      599
 NPARAMETR:  9.7587E-01  1.2383E+00  7.0191E-01  8.4543E-01  9.0366E-01  1.0072E+00  8.2774E-01  2.2060E-01  1.1106E+00  8.9281E-01
             1.0090E+00
 PARAMETER:  7.5572E-02  3.1370E-01 -2.5395E-01 -6.7913E-02 -1.2999E-03  1.0721E-01 -8.9058E-02 -1.4114E+00  2.0489E-01 -1.3381E-02
             1.0896E-01
 GRADIENT:  -3.6272E-01 -5.9403E+00 -2.6323E+00 -1.5392E+00  5.6460E+00  1.8031E-01 -2.6716E-01  1.3417E-01  4.6549E-02 -1.3301E-01
            -2.4348E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1660.22668473884        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  9.7611E-01  1.3063E+00  6.7399E-01  8.0224E-01  9.2204E-01  1.0068E+00  7.9908E-01  1.3057E-01  1.1560E+00  8.9946E-01
             1.0101E+00
 PARAMETER:  7.5816E-02  3.6717E-01 -2.9454E-01 -1.2035E-01  1.8830E-02  1.0675E-01 -1.2430E-01 -1.9358E+00  2.4499E-01 -5.9650E-03
             1.1003E-01
 GRADIENT:  -9.9804E-02 -2.0708E+00 -7.3655E-01 -1.1013E+00  1.8391E+00 -2.3750E-02 -1.3711E-01  4.6714E-02  4.0467E-02 -8.8054E-02
            -1.0864E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1660.25926978706        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      951
 NPARAMETR:  9.7612E-01  1.2645E+00  6.8360E-01  8.2976E-01  9.0436E-01  1.0068E+00  8.2141E-01  3.3014E-02  1.1247E+00  8.9294E-01
             1.0100E+00
 PARAMETER:  7.5833E-02  3.3472E-01 -2.8038E-01 -8.6621E-02 -5.3135E-04  1.0681E-01 -9.6734E-02 -3.3108E+00  2.1748E-01 -1.3241E-02
             1.0998E-01
 GRADIENT:  -5.5325E-02 -2.9917E-01 -3.0858E-01  1.0896E-02  1.2388E-01 -1.1218E-02  3.0146E-03  2.9208E-03  6.7777E-02  1.8465E-01
             8.4150E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1660.26077074150        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1126
 NPARAMETR:  9.7614E-01  1.2631E+00  6.8494E-01  8.3081E-01  9.0428E-01  1.0069E+00  8.2246E-01  1.0000E-02  1.1235E+00  8.9252E-01
             1.0100E+00
 PARAMETER:  7.5853E-02  3.3356E-01 -2.7842E-01 -8.5358E-02 -6.1921E-04  1.0683E-01 -9.5458E-02 -5.1503E+00  2.1646E-01 -1.3708E-02
             1.0995E-01
 GRADIENT:   7.3180E-03 -1.1267E-02  9.9374E-05 -2.4188E-03  3.5967E-03 -2.6601E-04  4.6162E-04  0.0000E+00  4.8450E-03 -2.2954E-03
            -3.5426E-04

0ITERATION NO.:   36    OBJECTIVE VALUE:  -1660.26077074150        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1155
 NPARAMETR:  9.7616E-01  1.2631E+00  6.8494E-01  8.3081E-01  9.0427E-01  1.0069E+00  8.2240E-01  1.0000E-02  1.1234E+00  8.9261E-01
             1.0100E+00
 PARAMETER:  7.5853E-02  3.3356E-01 -2.7842E-01 -8.5358E-02 -6.1921E-04  1.0683E-01 -9.5458E-02 -5.1503E+00  2.1646E-01 -1.3708E-02
             1.0995E-01
 GRADIENT:  -7.6646E-03 -5.0937E-03  2.3515E-04 -2.0195E-03  3.9786E-03 -1.7643E-03  1.0211E-03  0.0000E+00  4.2428E-03 -2.6132E-03
            -2.5936E-04

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1155
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.1550E-05 -2.1874E-02 -3.3920E-04  1.0951E-02 -2.4001E-02
 SE:             2.9829E-02  2.0098E-02  1.5111E-04  2.5091E-02  2.3563E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9782E-01  2.7644E-01  2.4786E-02  6.6252E-01  3.0840E-01

 ETASHRINKSD(%)  7.0370E-02  3.2669E+01  9.9494E+01  1.5940E+01  2.1061E+01
 ETASHRINKVR(%)  1.4069E-01  5.4665E+01  9.9997E+01  2.9340E+01  3.7687E+01
 EBVSHRINKSD(%)  4.4153E-01  3.2237E+01  9.9545E+01  1.6286E+01  2.0323E+01
 EBVSHRINKVR(%)  8.8111E-01  5.4082E+01  9.9998E+01  2.9919E+01  3.6515E+01
 RELATIVEINF(%)  9.8997E+01  1.9166E+00  2.1826E-04  3.8638E+00  6.9519E+00
 EPSSHRINKSD(%)  4.4101E+01
 EPSSHRINKVR(%)  6.8753E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1660.2607707414995     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -925.10994417776135     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.48
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1660.261       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  1.26E+00  6.85E-01  8.31E-01  9.04E-01  1.01E+00  8.22E-01  1.00E-02  1.12E+00  8.93E-01  1.01E+00
 


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
+        1.14E+03
 
 TH 2
+       -6.71E+00  4.97E+02
 
 TH 3
+        1.15E+01  2.56E+02  5.26E+02
 
 TH 4
+       -1.40E+01  3.82E+02 -2.35E+02  8.99E+02
 
 TH 5
+        4.12E-01 -4.08E+02 -6.21E+02  2.58E+02  1.07E+03
 
 TH 6
+       -1.49E+00 -1.36E+00  3.98E+00 -3.81E+00 -1.30E+00  1.92E+02
 
 TH 7
+        3.12E+00  1.13E+01  3.62E-01 -9.36E+00 -1.99E+01  5.20E+00  6.15E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.86E+00 -2.76E+01 -2.76E+01  3.82E+01 -6.05E-01 -9.44E-02  2.51E+01  0.00E+00  8.12E+01
 
 TH10
+       -1.45E+00 -1.09E+01 -5.69E+01 -1.13E+01 -5.89E+01 -4.43E+00  2.28E+01  0.00E+00  6.38E+00  1.06E+02
 
 TH11
+       -6.88E+00 -1.99E+01 -3.40E+01 -1.51E-01  5.33E-01  3.57E+00  6.74E+00  0.00E+00  7.90E+00  2.03E+01  2.07E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.159
Stop Time:
Sat Sep 25 12:14:01 CDT 2021
