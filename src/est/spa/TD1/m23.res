Sat Sep 18 13:56:00 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1709.19822584279        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:  -4.7812E+00  6.9706E+00 -1.6404E+01  3.2801E+01 -1.3705E+01  2.1611E+01  3.6841E+00  7.3128E+00  1.7470E+01  9.7536E+00
             1.4234E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1712.40156885295        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  1.0229E+00  1.0429E+00  1.0675E+00  9.7057E-01  1.0519E+00  9.3946E-01  9.9066E-01  9.7111E-01  9.3577E-01  9.7078E-01
             1.0030E+00
 PARAMETER:  1.2265E-01  1.4197E-01  1.6529E-01  7.0128E-02  1.5063E-01  3.7550E-02  9.0614E-02  7.0680E-02  3.3615E-02  7.0344E-02
             1.0297E-01
 GRADIENT:  -9.9873E-01  7.0260E+00  1.2971E-01  1.1176E+01 -7.5584E+00 -6.0406E+00 -1.3898E+00  1.3781E-01  6.5818E-01 -2.5997E+00
            -2.1521E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1712.69995851705        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  1.0215E+00  1.1193E+00  1.1397E+00  9.2934E-01  1.1249E+00  9.4900E-01  9.6211E-01  1.1738E+00  9.5659E-01  1.0267E+00
             1.0069E+00
 PARAMETER:  1.2125E-01  2.1268E-01  2.3080E-01  2.6720E-02  2.1768E-01  4.7656E-02  6.1371E-02  2.6026E-01  5.5618E-02  1.2633E-01
             1.0686E-01
 GRADIENT:  -4.5731E+00  1.0499E+01 -3.7032E+00  1.5476E+01 -6.6719E-02 -1.9693E+00  8.6129E-01  1.6869E+00 -1.3529E+00 -5.5375E-01
            -3.3695E-01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1713.57740345956        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      495
 NPARAMETR:  1.0231E+00  1.1618E+00  1.4178E+00  9.0254E-01  1.2567E+00  9.5612E-01  6.9305E-01  1.5542E+00  1.1307E+00  1.1527E+00
             1.0063E+00
 PARAMETER:  1.2286E-01  2.5001E-01  4.4912E-01 -2.5368E-03  3.2846E-01  5.5131E-02 -2.6665E-01  5.4096E-01  2.2280E-01  2.4213E-01
             1.0628E-01
 GRADIENT:   8.9886E-01  2.0291E+00 -5.3689E+00  1.2229E+01  1.3256E+01  1.2827E+00  3.2626E+00 -6.4337E-02  5.2842E+00  1.0662E+00
            -1.4080E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1713.78492137478        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:      694             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0231E+00  1.1626E+00  1.4196E+00  9.0178E-01  1.2572E+00  9.5443E-01  6.3452E-01  1.5586E+00  1.1191E+00  1.1493E+00
             1.0065E+00
 PARAMETER:  1.2288E-01  2.5065E-01  4.5036E-01 -3.3883E-03  3.2886E-01  5.3362E-02 -3.5488E-01  5.4377E-01  2.1255E-01  2.3911E-01
             1.0648E-01
 GRADIENT:   5.8757E+01  2.0838E+01 -5.3630E+00  1.9311E+01  1.5581E+01  5.7417E+00  1.4668E+00 -3.2483E-01 -1.6943E-01  1.2231E-01
            -6.0633E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1713.81362701737        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:      823             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0232E+00  1.1625E+00  1.4194E+00  9.0180E-01  1.2573E+00  9.5277E-01  6.1140E-01  1.5584E+00  1.1252E+00  1.1500E+00
             1.0075E+00
 PARAMETER:  1.2291E-01  2.5059E-01  4.5025E-01 -3.3633E-03  3.2894E-01  5.1617E-02 -3.9200E-01  5.4363E-01  2.1799E-01  2.3973E-01
             1.0745E-01
 GRADIENT:   5.8698E+01  2.1829E+01 -5.6007E+00  2.0295E+01  1.5699E+01  5.0502E+00  1.0405E+00 -4.3006E-01 -5.7446E-01  1.1587E-03
            -2.4911E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1713.99403393895        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1003
 NPARAMETR:  1.0228E+00  1.1555E+00  1.4610E+00  8.9275E-01  1.2405E+00  9.5295E-01  6.0238E-01  1.5717E+00  1.1342E+00  1.1512E+00
             1.0083E+00
 PARAMETER:  1.2258E-01  2.4454E-01  4.7914E-01 -1.3453E-02  3.1549E-01  5.1808E-02 -4.0686E-01  5.5218E-01  2.2589E-01  2.4082E-01
             1.0824E-01
 GRADIENT:   7.2240E-01 -2.5993E+00  2.5424E+00 -7.0542E+00 -6.1655E+00  1.1803E-01  1.0898E+00 -1.7598E+00 -5.6046E-01  1.7720E+00
             2.0595E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1714.11065049617        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1185
 NPARAMETR:  1.0252E+00  1.1578E+00  1.4654E+00  8.9478E-01  1.2459E+00  9.5414E-01  5.5201E-01  1.6392E+00  1.1502E+00  1.1364E+00
             1.0083E+00
 PARAMETER:  1.2490E-01  2.4656E-01  4.8211E-01 -1.1177E-02  3.1983E-01  5.3054E-02 -4.9418E-01  5.9419E-01  2.3993E-01  2.2784E-01
             1.0822E-01
 GRADIENT:   6.0699E+00  8.7542E-01 -1.2691E+00  1.0001E+00 -1.4905E+00  5.4549E-01  1.8104E-01  2.5876E-01 -7.5407E-01 -8.3788E-01
            -1.0033E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1714.11613029135        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1384            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0252E+00  1.1577E+00  1.4653E+00  8.9470E-01  1.2457E+00  9.5280E-01  5.4473E-01  1.6331E+00  1.1557E+00  1.1365E+00
             1.0085E+00
 PARAMETER:  1.2486E-01  2.4643E-01  4.8209E-01 -1.1265E-02  3.1973E-01  5.1650E-02 -5.0747E-01  5.9047E-01  2.4467E-01  2.2796E-01
             1.0846E-01
 GRADIENT:   6.4782E+01  1.5341E+01 -6.9373E-01  7.2579E+00  6.1101E-02  5.0915E+00  9.4823E-01  1.1020E-01  1.6267E+00 -7.6120E-01
             9.1053E-02

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1714.11613029135        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:     1455
 NPARAMETR:  1.0252E+00  1.1577E+00  1.4653E+00  8.9470E-01  1.2457E+00  9.5280E-01  5.4447E-01  1.6331E+00  1.1557E+00  1.1368E+00
             1.0085E+00
 PARAMETER:  1.2486E-01  2.4643E-01  4.8209E-01 -1.1265E-02  3.1973E-01  5.1650E-02 -5.0747E-01  5.9047E-01  2.4467E-01  2.2796E-01
             1.0846E-01
 GRADIENT:  -8.4574E+05  2.1426E+05  2.1906E+05  5.2798E+05  3.3029E+05 -8.4055E-03  1.7746E-01 -1.1926E-02 -2.6633E-02 -9.2295E-01
             3.5597E-03
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         3.7         1.3         3.2         3.1         1.3
                    4.1

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1455
 NO. OF SIG. DIGITS IN FINAL EST.:  1.3
 ADDITIONAL PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.0847E-03 -3.2063E-02 -3.4713E-02  7.4354E-03 -4.2144E-02
 SE:             2.9879E-02  1.1981E-02  1.4572E-02  2.7007E-02  2.2110E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4438E-01  7.4501E-03  1.7211E-02  7.8308E-01  5.6638E-02

 ETASHRINKSD(%)  1.0000E-10  5.9861E+01  5.1183E+01  9.5231E+00  2.5929E+01
 ETASHRINKVR(%)  1.0000E-10  8.3888E+01  7.6169E+01  1.8139E+01  4.5135E+01
 EBVSHRINKSD(%)  4.6198E-01  6.0996E+01  5.5558E+01  9.6803E+00  2.2943E+01
 EBVSHRINKVR(%)  9.2184E-01  8.4787E+01  8.0249E+01  1.8424E+01  4.0623E+01
 RELATIVEINF(%)  9.8755E+01  8.9204E-01  5.9231E+00  5.3615E+00  2.1182E+01
 EPSSHRINKSD(%)  4.4408E+01
 EPSSHRINKVR(%)  6.9096E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1714.1161302913538     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -978.96530372761561     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.81
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.48
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1714.116       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.16E+00  1.47E+00  8.95E-01  1.25E+00  9.53E-01  5.45E-01  1.63E+00  1.16E+00  1.14E+00  1.01E+00
 


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
+        1.61E+09
 
 TH 2
+        1.85E+03  3.24E+08
 
 TH 3
+        7.57E+02  2.54E+04  5.29E+07
 
 TH 4
+       -2.31E+09 -2.31E+04  4.18E+08  3.30E+09
 
 TH 5
+       -5.18E+08  2.49E+04  1.00E+04  7.99E+04  1.66E+08
 
 TH 6
+       -6.26E+03  2.81E+03  1.13E+03  8.96E+03  2.01E+03  2.16E+02
 
 TH 7
+        1.83E+05 -8.20E+04 -3.31E+04 -2.61E+05 -5.87E+04  8.48E-01  2.42E+01
 
 TH 8
+        2.67E+05 -1.20E+05 -4.83E+04 -3.81E+05 -8.57E+04  9.72E-02  1.73E+00  1.30E+01
 
 TH 9
+        4.47E+04 -2.01E+04 -8.10E+03 -6.39E+04 -1.44E+04  6.78E-01  3.27E+01  3.29E-01  1.00E+02
 
 TH10
+        7.78E+04 -3.49E+04 -1.41E+04 -1.11E+05  2.56E+08 -8.15E-01 -3.69E+08  1.06E+08  1.05E-01  3.93E+08
 
 TH11
+        1.40E+05 -6.30E+04 -2.54E+04 -2.01E+05 -4.51E+04  7.64E-01  2.40E+00  3.10E+00  7.51E+00 -9.32E+08  2.08E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.344
Stop Time:
Sat Sep 18 13:56:29 CDT 2021
