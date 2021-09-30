Wed Sep 29 06:30:41 CDT 2021
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
$DATA ../../../../data/int/TD1/dat60.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3266.44836540557        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1383E+02  2.5479E+01  7.9530E+01  1.6113E+00  1.3865E+02  4.7013E+01 -8.2562E+00 -2.6937E+02 -6.1141E+01  1.2425E-01
            -8.9089E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3612.63734831837        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:       88
 NPARAMETR:  8.5757E-01  9.1922E-01  9.8367E-01  1.1188E+00  8.3780E-01  8.5197E-01  9.1385E-01  1.0810E+00  1.0892E+00  8.8785E-01
             1.7003E+00
 PARAMETER: -5.3653E-02  1.5772E-02  8.3531E-02  2.1228E-01 -7.6974E-02 -6.0206E-02  9.9078E-03  1.7789E-01  1.8549E-01 -1.8957E-02
             6.3083E-01
 GRADIENT:  -3.3116E+02  2.7871E+01  3.9785E+01  1.3193E+02 -9.2770E+01 -1.1499E+02 -5.1708E+00  1.6821E+01  4.2737E+01  8.9674E-01
             7.0642E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3636.34669399078        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      184
 NPARAMETR:  8.2309E-01  6.6440E-01  5.8698E-01  1.2443E+00  5.7998E-01  9.0593E-01  8.4525E-01  4.6522E-01  9.6563E-01  9.2687E-01
             1.5438E+00
 PARAMETER: -9.4684E-02 -3.0888E-01 -4.3276E-01  3.1855E-01 -4.4475E-01  1.2047E-03 -6.8127E-02 -6.6524E-01  6.5026E-02  2.4063E-02
             5.3423E-01
 GRADIENT:  -5.4066E+02  4.1748E+00 -7.8634E+01  1.2242E+02 -3.2255E+01 -1.6175E+02 -1.8295E+01  6.8736E+00  4.5475E+00  7.3134E+00
             5.7226E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3734.53555863563        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      365
 NPARAMETR:  8.9786E-01  1.0767E+00  1.0775E+00  1.0658E+00  9.8616E-01  8.4543E-01  5.1105E-01  8.5326E-01  1.1057E+00  1.2376E+00
             1.4370E+00
 PARAMETER: -7.7422E-03  1.7386E-01  1.7465E-01  1.6369E-01  8.6061E-02 -6.7911E-02 -5.7129E-01 -5.8688E-02  2.0047E-01  3.1321E-01
             4.6253E-01
 GRADIENT:  -3.5358E+02  5.3652E+01  3.2918E+01  1.0514E+02 -4.4257E+01 -1.2998E+02 -1.7869E+01 -1.7214E+00 -6.1652E+00  2.5037E+01
             5.0607E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3761.09499069906        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      537
 NPARAMETR:  9.5674E-01  1.0764E+00  1.0357E+00  1.0655E+00  1.0017E+00  9.4426E-01  5.1188E-01  8.7974E-01  1.1052E+00  1.1607E+00
             1.4363E+00
 PARAMETER:  5.5773E-02  1.7365E-01  1.3507E-01  1.6346E-01  1.0169E-01  4.2644E-02 -5.6966E-01 -2.8130E-02  2.0004E-01  2.4902E-01
             4.6205E-01
 GRADIENT:  -1.2708E+02  4.0538E+01  7.2908E+00  1.0573E+02 -9.1295E+00 -3.9941E+01 -1.9795E+01 -1.3211E-01 -6.0029E+00  9.6529E+00
             5.0883E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3766.09878030195        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      713
 NPARAMETR:  1.0106E+00  1.0765E+00  1.0238E+00  1.0655E+00  1.0056E+00  1.0146E+00  5.1195E-01  8.8597E-01  1.1053E+00  1.1169E+00
             1.4358E+00
 PARAMETER:  1.1051E-01  1.7368E-01  1.2348E-01  1.6348E-01  1.0558E-01  1.1450E-01 -5.6952E-01 -2.1074E-02  2.0008E-01  2.1059E-01
             4.6172E-01
 GRADIENT:   7.3450E+00  3.6651E+01 -7.1645E-01  1.0540E+02  1.9442E+00 -2.6150E+00 -2.1472E+01 -1.1397E-01 -5.5688E+00 -5.7302E-01
             5.0826E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3766.84034219256        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      912
 NPARAMETR:  1.0146E+00  1.0762E+00  1.0163E+00  1.0651E+00  1.0026E+00  1.0093E+00  5.1525E-01  9.3765E-01  1.1052E+00  1.1192E+00
             1.4339E+00
 PARAMETER:  1.1446E-01  1.7347E-01  1.1622E-01  1.6308E-01  1.0259E-01  1.0922E-01 -5.6310E-01  3.5624E-02  2.0002E-01  2.1262E-01
             4.6038E-01
 GRADIENT:   1.6046E+01  3.7139E+01 -5.1686E+00  1.0457E+02  3.3364E-01 -4.8592E+00 -2.1038E+01  1.7040E+00 -5.4379E+00  1.0663E+00
             5.0967E+02

0ITERATION NO.:   31    OBJECTIVE VALUE:  -3766.84034219256        NO. OF FUNC. EVALS.:  34
 CUMULATIVE NO. OF FUNC. EVALS.:      946
 NPARAMETR:  1.0140E+00  1.0762E+00  1.0164E+00  1.0651E+00  1.0026E+00  1.0104E+00  5.1519E-01  9.3763E-01  1.1051E+00  1.1189E+00
             1.4340E+00
 PARAMETER:  1.1446E-01  1.7347E-01  1.1622E-01  1.6308E-01  1.0259E-01  1.0922E-01 -5.6310E-01  3.5624E-02  2.0002E-01  2.1262E-01
             4.6038E-01
 GRADIENT:   1.4167E+01  5.3507E+04 -7.9814E+04  5.6980E+04  3.6340E-01 -5.1310E+00  3.2928E+04  9.2720E+04  9.2722E+04  9.1250E-01
            -4.0013E+04
 NUMSIGDIG:         0.9         2.3         2.3         2.3         2.6         0.6         2.3         2.3         2.3         1.5
                    2.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      946
 NO. OF SIG. DIGITS IN FINAL EST.:  0.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.9276E-03 -3.7653E-02 -1.4303E-02 -4.4317E-02 -1.5326E-02
 SE:             3.0208E-02  1.7321E-02  1.4062E-02  2.8174E-02  2.6824E-02
 N:                     100         100         100         100         100

 P VAL.:         8.4443E-01  2.9720E-02  3.0910E-01  1.1573E-01  5.6775E-01

 ETASHRINKSD(%)  1.0000E-10  4.1972E+01  5.2890E+01  5.6122E+00  1.0136E+01
 ETASHRINKVR(%)  1.0000E-10  6.6327E+01  7.7806E+01  1.0909E+01  1.9245E+01
 EBVSHRINKSD(%)  4.9855E-01  5.3324E+01  5.3453E+01  5.2963E+00  9.5743E+00
 EBVSHRINKVR(%)  9.9462E-01  7.8213E+01  7.8334E+01  1.0312E+01  1.8232E+01
 RELATIVEINF(%)  9.9000E+01  6.0344E+00  1.3049E+01  4.5398E+01  2.2740E+01
 EPSSHRINKSD(%)  3.9035E+01
 EPSSHRINKVR(%)  6.2833E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3766.8403421925623     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -2112.7509824241515     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.11
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    12.06
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3766.840       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.08E+00  1.02E+00  1.07E+00  1.00E+00  1.01E+00  5.15E-01  9.38E-01  1.11E+00  1.12E+00  1.43E+00
 


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
+        6.55E+02  1.33E+07
 
 TH 3
+        3.38E+07  9.95E+03  3.32E+07
 
 TH 4
+        7.04E+02 -5.74E+03  1.06E+04  1.54E+07
 
 TH 5
+       -5.73E-01 -6.10E+02 -2.18E+02  2.60E+07  4.39E+07
 
 TH 6
+       -3.62E+07  1.14E+03  6.62E-01  1.23E+03 -4.09E+07  3.82E+07
 
 TH 7
+        4.25E+02 -1.16E+03  6.37E+03 -3.95E+03 -9.77E+01  7.36E+02  5.51E+06
 
 TH 8
+        1.31E+03 -3.12E+04 -4.18E+07  2.85E+07  4.81E+07  2.28E+03 -2.01E+04  5.27E+07
 
 TH 9
+        5.61E+02 -8.28E+03  8.38E+03 -5.21E+03 -1.44E+02  9.70E+02 -9.18E+02 -2.63E+04  9.48E+06
 
 TH10
+       -1.68E+07 -7.21E+01 -1.65E+07 -5.96E+01 -9.73E+00 -1.77E+07  5.19E+00 -1.14E+02 -5.59E+01  8.19E+06
 
 TH11
+       -1.95E+02  1.44E+03 -2.83E+03  1.72E+03  4.79E+01 -3.22E+02  3.40E+02  8.87E+03  2.33E+03  1.89E+01  1.08E+06
 
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
 #CPUT: Total CPU Time in Seconds,       36.284
Stop Time:
Wed Sep 29 06:31:19 CDT 2021
