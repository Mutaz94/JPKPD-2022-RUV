Sat Sep 25 01:16:14 CDT 2021
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
$DATA ../../../../data/int/SL2/dat45.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      995
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

 TOT. NO. OF OBS RECS:      895
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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2132.20084459365        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0667E+02  3.0657E+00  1.3452E+02  2.4579E+01  1.6194E+02 -6.4989E+01 -8.6254E+01 -5.8981E+02 -1.9764E+02 -2.5941E+01
            -2.4855E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3141.82089767067        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.7181E-01  1.1341E+00  1.0304E+00  9.0828E-01  1.0265E+00  1.1668E+00  9.2611E-01  1.1823E+00  1.0949E+00  8.2381E-01
             1.8380E+00
 PARAMETER:  7.1405E-02  2.2584E-01  1.2997E-01  3.8000E-03  1.2613E-01  2.5431E-01  2.3238E-02  2.6747E-01  1.9064E-01 -9.3814E-02
             7.0868E-01
 GRADIENT:   1.9472E+00 -1.3264E+01 -9.5112E-01 -2.1753E+01  6.9686E+00  7.7160E+00 -1.5995E+00 -9.6829E+00 -4.5456E+00 -8.5218E+00
            -1.6284E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3143.17260638110        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.7845E-01  1.2514E+00  1.1035E+00  8.5818E-01  1.1246E+00  1.1722E+00  7.7621E-01  1.3288E+00  1.1534E+00  9.7933E-01
             1.8081E+00
 PARAMETER:  7.8213E-02  3.2428E-01  1.9852E-01 -5.2938E-02  2.1742E-01  2.5890E-01 -1.5334E-01  3.8428E-01  2.4273E-01  7.9113E-02
             6.9228E-01
 GRADIENT:   1.3745E+01  7.7012E+00  1.7684E+00 -1.7332E+00  4.9753E+00  9.7671E+00 -3.8569E+00 -9.8096E+00 -2.0470E+00 -1.5009E+00
            -5.7436E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3145.02226002359        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      340
 NPARAMETR:  9.8156E-01  1.2580E+00  1.1704E+00  8.6476E-01  1.1433E+00  1.1576E+00  8.1367E-01  1.5861E+00  1.1273E+00  9.7169E-01
             1.8414E+00
 PARAMETER:  8.1390E-02  3.2954E-01  2.5732E-01 -4.5307E-02  2.3390E-01  2.4634E-01 -1.0620E-01  5.6131E-01  2.1980E-01  7.1280E-02
             7.1053E-01
 GRADIENT:   3.5272E+00  4.4802E+00  4.8651E-01  4.8730E+00  9.7905E-01 -2.1731E+00 -4.8070E-01 -4.8023E+00 -3.9190E+00  3.5604E-02
            -8.8048E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3145.17109741271        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      506
 NPARAMETR:  9.8136E-01  1.2601E+00  1.1807E+00  8.5744E-01  1.1475E+00  1.1595E+00  8.1513E-01  1.6265E+00  1.1581E+00  9.7019E-01
             1.8422E+00
 PARAMETER:  8.1187E-02  3.3118E-01  2.6612E-01 -5.3806E-02  2.3761E-01  2.4799E-01 -1.0440E-01  5.8645E-01  2.4681E-01  6.9742E-02
             7.1094E-01
 GRADIENT:   1.7492E+01  8.4240E+00  1.5896E+00 -2.5783E-01  2.5665E+00  5.0296E+00  6.3741E-01 -3.7420E+00  2.1235E+00  5.4424E-02
            -3.9659E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3145.18354851477        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      684
 NPARAMETR:  9.8136E-01  1.2601E+00  1.1807E+00  8.5959E-01  1.1475E+00  1.1640E+00  8.1247E-01  1.6265E+00  1.1494E+00  9.7151E-01
             1.8422E+00
 PARAMETER:  8.1187E-02  3.3118E-01  2.6612E-01 -5.1298E-02  2.3761E-01  2.5186E-01 -1.0768E-01  5.8645E-01  2.3926E-01  7.1099E-02
             7.1094E-01
 GRADIENT:   3.2954E+00 -7.2630E-01  9.3854E-01 -2.5667E-02 -4.5915E-02  2.9374E-03  2.3594E-03 -3.9882E+00 -6.4759E-03  4.6279E-03
            -5.7032E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3145.23426854180        NO. OF FUNC. EVALS.: 151
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  9.7993E-01  1.2601E+00  1.1800E+00  8.5951E-01  1.1474E+00  1.1630E+00  8.1222E-01  1.6469E+00  1.1490E+00  9.7131E-01
             1.8429E+00
 PARAMETER:  7.9727E-02  3.3123E-01  2.6549E-01 -5.1393E-02  2.3751E-01  2.5098E-01 -1.0798E-01  5.9889E-01  2.3890E-01  7.0891E-02
             7.1133E-01
 GRADIENT:   1.5037E+01  1.0305E+01  5.1484E-01  2.1441E+00  2.8918E+00  6.2926E+00  2.3372E-01 -3.1483E+00  6.2822E-01  2.1548E-01
            -2.7762E+00

0ITERATION NO.:   33    OBJECTIVE VALUE:  -3145.23427395662        NO. OF FUNC. EVALS.:  83
 CUMULATIVE NO. OF FUNC. EVALS.:      918
 NPARAMETR:  9.7982E-01  1.2601E+00  1.1800E+00  8.5951E-01  1.1474E+00  1.1629E+00  8.1219E-01  1.6469E+00  1.1490E+00  9.7121E-01
             1.8429E+00
 PARAMETER:  7.9612E-02  3.3123E-01  2.6549E-01 -5.1412E-02  2.3749E-01  2.5077E-01 -1.0799E-01  5.9889E-01  2.3889E-01  7.0885E-02
             7.1133E-01
 GRADIENT:  -8.5840E+05 -1.7571E+01  2.9347E+01 -1.3613E-01 -1.8073E+05 -3.9805E-01  1.5585E-02 -7.1695E+04  7.4948E-04  9.6233E-02
            -2.3649E+02
 NUMSIGDIG:         3.3         7.8         7.6         2.6         3.3         2.1         2.5         3.3         4.8         1.9
                    6.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      918
 NO. OF SIG. DIGITS IN FINAL EST.:  1.9

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.6588E-04 -3.0864E-02 -2.1808E-02  1.7402E-02 -3.0137E-02
 SE:             2.9782E-02  2.0143E-02  1.7345E-02  2.5822E-02  2.3163E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7948E-01  1.2547E-01  2.0865E-01  5.0037E-01  1.9324E-01

 ETASHRINKSD(%)  2.2565E-01  3.2517E+01  4.1892E+01  1.3493E+01  2.2399E+01
 ETASHRINKVR(%)  4.5078E-01  5.4460E+01  6.6235E+01  2.5165E+01  3.9782E+01
 EBVSHRINKSD(%)  6.5823E-01  3.3095E+01  4.6442E+01  1.4359E+01  2.1156E+01
 EBVSHRINKVR(%)  1.3121E+00  5.5237E+01  7.1316E+01  2.6657E+01  3.7837E+01
 RELATIVEINF(%)  9.8679E+01  1.1316E+01  1.8312E+01  2.5383E+01  1.7628E+01
 EPSSHRINKSD(%)  1.8737E+01
 EPSSHRINKVR(%)  3.3963E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          895
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1644.8999744363641     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3145.2342739566202     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1500.3342995202561     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.76
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.77
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3145.234       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.80E-01  1.26E+00  1.18E+00  8.59E-01  1.15E+00  1.16E+00  8.12E-01  1.65E+00  1.15E+00  9.71E-01  1.84E+00
 


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
+        4.47E+09
 
 TH 2
+        5.25E+08  2.46E+08
 
 TH 3
+        3.98E-02  1.64E+08  2.19E+08
 
 TH 4
+        4.54E+04  3.49E+02 -5.12E+01  2.90E+09
 
 TH 5
+        4.80E+03 -2.78E+02 -1.24E+02 -9.16E+08  2.89E+08
 
 TH 6
+       -5.46E+03 -1.92E+00  6.72E-01  8.56E+08 -1.96E+03  1.44E+02
 
 TH 7
+       -8.42E+03  2.62E+01 -4.63E+00 -2.85E+09 -3.03E+03  3.75E-01  6.94E+01
 
 TH 8
+       -2.22E+08  5.21E+07  4.88E+04  6.78E+00  7.98E+07 -8.46E-01  8.39E-01  4.41E+07
 
 TH 9
+        2.51E+03 -1.58E+01  2.18E+00  3.50E+01  9.01E+02 -7.32E-02  1.31E+01  2.54E+02  7.99E+01
 
 TH10
+        5.10E+02 -2.06E+01  2.54E+00 -2.57E+09  8.11E+08 -7.58E+08  2.52E+09  2.24E+08  1.60E+00  2.27E+09
 
 TH11
+       -1.67E+08 -3.92E+07  1.01E+05  1.91E+08 -5.52E+00  5.62E+07  1.87E+08 -1.66E+07  1.23E+01  7.31E+00  2.50E+07
 
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
 #CPUT: Total CPU Time in Seconds,       38.631
Stop Time:
Sat Sep 25 01:16:54 CDT 2021
