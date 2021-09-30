Wed Sep 29 12:47:11 CDT 2021
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
$DATA ../../../../data/spa/A2/dat38.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m38.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1074.56414610110        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8794E+02  9.4913E+01  1.4885E+01  1.4001E+02  1.2921E+02  3.5190E+01 -4.6508E+01 -6.6624E+00 -5.6091E+01 -7.0647E+01
            -9.5934E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1392.27655524519        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0638E+00  9.4867E-01  1.0270E+00  1.0200E+00  9.0073E-01  1.0771E+00  1.1281E+00  9.3885E-01  1.1988E+00  9.9013E-01
             2.4503E+00
 PARAMETER:  1.6186E-01  4.7303E-02  1.2666E-01  1.1982E-01 -4.5450E-03  1.7424E-01  2.2055E-01  3.6900E-02  2.8132E-01  9.0077E-02
             9.9620E-01
 GRADIENT:   2.3697E+02  2.1304E+01  7.1255E+00  2.2462E+01 -2.2883E+01  3.1789E+01  2.0823E+00  5.6555E+00  1.3529E+01  7.1722E+00
             3.0594E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1399.26330483861        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0417E+00  7.2180E-01  6.9562E-01  1.1638E+00  6.8716E-01  1.0794E+00  1.6544E+00  2.9417E-01  1.0418E+00  8.3078E-01
             2.2597E+00
 PARAMETER:  1.4089E-01 -2.2600E-01 -2.6296E-01  2.5167E-01 -2.7519E-01  1.7644E-01  6.0346E-01 -1.1236E+00  1.4095E-01 -8.5388E-02
             9.1522E-01
 GRADIENT:   1.9037E+02  2.8478E+01 -3.0218E+01  1.1393E+02  4.7897E+01  3.5341E+01  1.1444E+01  1.3018E+00  1.1242E+01  4.3784E+00
             9.0443E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1403.32010856551        NO. OF FUNC. EVALS.: 120
 CUMULATIVE NO. OF FUNC. EVALS.:      279
 NPARAMETR:  1.0133E+00  6.9422E-01  6.7736E-01  1.1463E+00  6.5620E-01  1.0394E+00  1.4622E+00  1.9255E-01  1.0464E+00  8.4062E-01
             2.2250E+00
 PARAMETER:  1.1320E-01 -2.6496E-01 -2.8956E-01  2.3654E-01 -3.2129E-01  1.3868E-01  4.7994E-01 -1.5474E+00  1.4537E-01 -7.3617E-02
             8.9976E-01
 GRADIENT:   2.1001E+01  1.1374E+01 -1.4213E+01  2.8381E+01  1.4669E+01  1.0761E+01 -4.5610E+00  5.4465E-01  5.1689E+00  6.0369E+00
            -3.3651E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1412.92627305780        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      454
 NPARAMETR:  1.0127E+00  2.8721E-01  3.5588E-01  1.2231E+00  3.4905E-01  1.0220E+00  2.5224E+00  1.0000E-02  9.3313E-01  6.7920E-01
             2.0924E+00
 PARAMETER:  1.1263E-01 -1.1475E+00 -9.3316E-01  3.0137E-01 -9.5254E-01  1.2176E-01  1.0252E+00 -4.5417E+00  3.0786E-02 -2.8685E-01
             8.3831E-01
 GRADIENT:   1.6385E+01  7.1491E+00 -3.8653E+01  3.4844E+01  3.6976E+01  1.1369E+00  6.4500E+00  1.8583E-04 -1.6224E+00  7.6416E-02
             5.8639E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1416.06352134057        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      632
 NPARAMETR:  1.0024E+00  1.7515E-01  3.7792E-01  1.2478E+00  3.3971E-01  1.0180E+00  3.1002E+00  1.0000E-02  9.6182E-01  7.3900E-01
             2.0515E+00
 PARAMETER:  1.0239E-01 -1.6421E+00 -8.7308E-01  3.2140E-01 -9.7967E-01  1.1784E-01  1.2315E+00 -5.7349E+00  6.1069E-02 -2.0246E-01
             8.1857E-01
 GRADIENT:   7.2149E+00  1.2564E+00  1.5981E+01 -1.9243E+01 -1.6439E+01  1.4780E+00 -2.1645E+00  0.0000E+00  7.9454E+00  9.6996E-01
             1.2045E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1417.48696089436        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  9.9159E-01  1.0651E-01  3.8941E-01  1.2931E+00  3.4330E-01  1.0098E+00  4.2636E+00  1.0000E-02  9.2037E-01  7.2079E-01
             2.0510E+00
 PARAMETER:  9.1551E-02 -2.1395E+00 -8.4311E-01  3.5701E-01 -9.6915E-01  1.0979E-01  1.5501E+00 -6.6100E+00  1.7018E-02 -2.2741E-01
             8.1831E-01
 GRADIENT:  -5.4227E+00  5.6354E+00 -1.5165E+00  1.1818E+00  1.9603E+00 -5.1443E-01  6.4274E+00  0.0000E+00 -3.2179E+00 -4.0752E+00
            -1.4506E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1417.71824374547        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  9.9155E-01  8.9611E-02  3.8554E-01  1.2952E+00  3.4053E-01  1.0105E+00  4.7337E+00  1.0000E-02  9.2430E-01  7.0973E-01
             2.0455E+00
 PARAMETER:  9.1519E-02 -2.3123E+00 -8.5310E-01  3.5866E-01 -9.7725E-01  1.1042E-01  1.6547E+00 -7.0032E+00  2.1279E-02 -2.4288E-01
             8.1566E-01
 GRADIENT:   6.7721E-01  2.1416E+01 -1.7656E+01 -1.8492E+01  3.0799E+01  4.6803E-01  3.3969E+01  0.0000E+00 -1.0718E+01 -1.7363E+01
            -6.8852E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1417.92876417003        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  9.9325E-01  7.6245E-02  3.6776E-01  1.2812E+00  3.2995E-01  1.0149E+00  5.1266E+00  1.0000E-02  9.3582E-01  6.8700E-01
             2.0378E+00
 PARAMETER:  9.3223E-02 -2.4738E+00 -9.0033E-01  3.4776E-01 -1.0088E+00  1.1475E-01  1.7344E+00 -7.4782E+00  3.3670E-02 -2.7542E-01
             8.1187E-01
 GRADIENT:   2.4261E+00  4.3207E+00 -2.8239E+00 -4.4238E+00  8.4188E+00  1.0954E+00  7.2162E+00  0.0000E+00  7.4374E-01 -6.0788E+00
            -2.2425E+00

0ITERATION NO.:   41    OBJECTIVE VALUE:  -1417.92876417003        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     1193
 NPARAMETR:  9.9293E-01  7.5958E-02  3.6830E-01  1.2819E+00  3.2944E-01  1.0137E+00  5.1138E+00  1.0000E-02  9.3599E-01  6.8734E-01
             2.0407E+00
 PARAMETER:  9.3223E-02 -2.4738E+00 -9.0033E-01  3.4776E-01 -1.0088E+00  1.1475E-01  1.7344E+00 -7.4782E+00  3.3670E-02 -2.7542E-01
             8.1187E-01
 GRADIENT:   1.5645E+00  4.6209E+01 -1.4460E+02 -3.8599E+02  1.2528E+02  9.0816E-01  6.3574E+01  0.0000E+00 -1.3758E+03 -2.5116E+02
            -1.6942E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1193
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         8.5675E-04  4.1338E-02 -1.5133E-04 -1.8802E-02  1.0259E-02
 SE:             2.9289E-02  1.4618E-02  2.2074E-04  2.6625E-02  2.1795E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7666E-01  4.6845E-03  4.9300E-01  4.8008E-01  6.3785E-01

 ETASHRINKSD(%)  1.8773E+00  5.1029E+01  9.9261E+01  1.0801E+01  2.6984E+01
 ETASHRINKVR(%)  3.7195E+00  7.6019E+01  9.9995E+01  2.0436E+01  4.6686E+01
 EBVSHRINKSD(%)  1.6074E+00  6.5195E+01  9.9177E+01  7.1528E+00  2.3338E+01
 EBVSHRINKVR(%)  3.1890E+00  8.7886E+01  9.9993E+01  1.3794E+01  4.1230E+01
 RELATIVEINF(%)  9.5975E+01  8.0045E+00  3.0941E-04  3.2012E+01  2.7900E+00
 EPSSHRINKSD(%)  3.8134E+01
 EPSSHRINKVR(%)  6.1726E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1417.9287641700298     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -682.77793760629163     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    15.99
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.05
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1417.929       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  7.62E-02  3.68E-01  1.28E+00  3.30E-01  1.01E+00  5.13E+00  1.00E-02  9.36E-01  6.87E-01  2.04E+00
 


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
+        1.82E+05
 
 TH 2
+       -1.39E+02  9.44E+04
 
 TH 3
+        6.18E+01  2.24E+02  3.65E+04
 
 TH 4
+        4.38E+01  3.66E+02 -7.51E+02  1.78E+04
 
 TH 5
+        3.20E+01 -1.38E+03 -7.42E+03  1.62E+02  4.26E+04
 
 TH 6
+        1.54E+05 -1.17E+01  8.35E+00 -1.07E+01  3.66E+00  2.85E+02
 
 TH 7
+       -1.14E+00 -2.62E+01  1.35E+00  3.90E+00 -1.97E+01 -1.23E-01  4.49E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.91E+05  1.22E+02 -1.14E+02 -1.78E+02  1.84E+02 -2.66E+01  3.79E+00  0.00E+00  3.93E+05
 
 TH10
+        2.19E+02  7.82E+01 -3.27E+02 -1.92E+02  2.93E+02 -3.01E+01  2.61E+00  0.00E+00  9.38E+04  9.59E+04
 
 TH11
+        1.11E+01 -6.99E+00 -6.72E+01 -3.35E+01  4.15E+01  9.92E-02  1.43E-01  0.00E+00 -5.72E+01  5.28E+03  1.31E+03
 
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
 #CPUT: Total CPU Time in Seconds,       23.120
Stop Time:
Wed Sep 29 12:47:36 CDT 2021
