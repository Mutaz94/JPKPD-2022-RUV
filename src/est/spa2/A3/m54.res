Thu Sep 30 06:32:05 CDT 2021
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
$DATA ../../../../data/spa2/A3/dat54.csv ignore=@
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
Current Date:       30 SEP 2021
Days until program expires : 199
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m54.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   743.447182152694        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8874E+02  2.2140E+02  4.2670E+02  6.4257E+01  3.2885E+02  4.0259E+01 -2.0767E+02 -4.1987E+02 -1.6690E+02 -2.0263E+02
            -5.3693E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1596.21955361400        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0976E+00  9.9459E-01  8.3742E-01  1.1100E+00  8.8882E-01  8.8879E-01  9.9635E-01  1.0346E+00  9.8152E-01  9.0927E-01
             5.2951E+00
 PARAMETER:  1.9313E-01  9.4579E-02 -7.7425E-02  2.0432E-01 -1.7856E-02 -1.7892E-02  9.6339E-02  1.3400E-01  8.1350E-02  4.8817E-03
             1.7668E+00
 GRADIENT:   5.2480E+01 -6.8939E+00 -1.7039E+01  2.3392E+01  9.3566E+00 -3.0774E+00  1.2180E+01  1.0376E+01  1.7526E+01  2.7924E+01
             3.5740E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1621.50296373512        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0760E+00  6.2544E-01  3.2184E-01  1.2717E+00  4.1200E-01  1.0261E+00  9.7619E-01  7.2648E-01  1.1496E+00  2.6956E-01
             4.8849E+00
 PARAMETER:  1.7326E-01 -3.6930E-01 -1.0337E+00  3.4034E-01 -7.8674E-01  1.2572E-01  7.5901E-02 -2.1954E-01  2.3944E-01 -1.2110E+00
             1.6862E+00
 GRADIENT:   4.8195E+00  8.9832E+01  4.5868E+00  1.2711E+02 -4.3139E+01  2.4392E+01 -4.2193E+00  4.9641E+00  1.3875E+01  7.9844E-01
             3.1485E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1685.95952286077        NO. OF FUNC. EVALS.: 139
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0308E+00  5.1240E-01  2.8752E-01  1.2688E+00  3.5988E-01  9.0234E-01  1.2711E+00  3.4655E-02  1.1995E+00  3.4583E-01
             3.7898E+00
 PARAMETER:  1.3031E-01 -5.6865E-01 -1.1464E+00  3.3803E-01 -9.2199E-01 -2.7636E-03  3.3989E-01 -3.2623E+00  2.8194E-01 -9.6180E-01
             1.4323E+00
 GRADIENT:  -7.5578E+01  6.1964E+01 -9.1238E+00  1.4115E+02 -2.0598E+01 -2.8313E+01 -6.7043E+00 -6.4740E-03 -3.1485E+00 -4.3109E+00
             1.2100E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1726.55379843487        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      475
 NPARAMETR:  1.0524E+00  3.1916E-01  1.8613E-01  1.1061E+00  2.3990E-01  1.0061E+00  8.5475E-01  1.0000E-02  1.4050E+00  7.1982E-01
             3.0590E+00
 PARAMETER:  1.5108E-01 -1.0421E+00 -1.5813E+00  2.0084E-01 -1.3275E+00  1.0613E-01 -5.6949E-02 -1.1788E+01  4.4007E-01 -2.2875E-01
             1.2181E+00
 GRADIENT:   5.6174E+00  3.0029E+01  2.8876E+01  5.1628E+01 -6.0199E+01  4.3191E+00  1.4643E+00  0.0000E+00  3.0981E+00  9.2967E-01
            -4.1701E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1731.62431419100        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      650
 NPARAMETR:  1.0509E+00  2.7599E-01  1.5027E-01  9.8253E-01  2.1248E-01  9.9563E-01  5.7474E-01  1.0000E-02  1.4744E+00  7.4722E-01
             3.1562E+00
 PARAMETER:  1.4963E-01 -1.1874E+00 -1.7953E+00  8.2377E-02 -1.4489E+00  9.5623E-02 -4.5385E-01 -1.6133E+01  4.8822E-01 -1.9139E-01
             1.2494E+00
 GRADIENT:   2.9151E-01 -1.4149E+00 -2.6441E+00 -3.6681E+00  6.2600E+00  6.5271E-01  1.0598E+00  0.0000E+00  4.6318E-01 -1.1681E+00
             1.4829E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1732.09450862651        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      827
 NPARAMETR:  1.0508E+00  2.7630E-01  1.5293E-01  9.8947E-01  2.1400E-01  9.9363E-01  1.9818E-01  1.0000E-02  1.4499E+00  7.7514E-01
             3.1784E+00
 PARAMETER:  1.4957E-01 -1.1863E+00 -1.7778E+00  8.9418E-02 -1.4418E+00  9.3605E-02 -1.5186E+00 -1.5894E+01  4.7147E-01 -1.5472E-01
             1.2564E+00
 GRADIENT:  -1.0892E+00 -2.5960E-01 -4.1484E-01 -2.5252E-01  5.1925E-01  4.2252E-01  1.2587E-01  0.0000E+00 -8.8532E-01  1.0649E+00
             2.6782E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1732.15518512851        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1002
 NPARAMETR:  1.0514E+00  2.7605E-01  1.5303E-01  9.8983E-01  2.1402E-01  9.9236E-01  2.2176E-02  1.0000E-02  1.4535E+00  7.7431E-01
             3.1817E+00
 PARAMETER:  1.5013E-01 -1.1872E+00 -1.7771E+00  8.9782E-02 -1.4417E+00  9.2328E-02 -3.7088E+00 -1.5693E+01  4.7401E-01 -1.5579E-01
             1.2574E+00
 GRADIENT:   6.3798E-02 -7.6715E-02 -2.0928E-02 -1.7365E-02  7.9370E-02  1.0000E-03  1.5041E-03  0.0000E+00 -1.4997E-02 -7.5097E-02
            -2.1039E-02

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1732.15580164065        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     1094
 NPARAMETR:  1.0514E+00  2.7613E-01  1.5305E-01  9.8988E-01  2.1405E-01  9.9238E-01  1.0000E-02  1.0000E-02  1.4535E+00  7.7452E-01
             3.1818E+00
 PARAMETER:  1.5011E-01 -1.1869E+00 -1.7770E+00  8.9827E-02 -1.4416E+00  9.2354E-02 -4.6763E+00 -1.5592E+01  4.7397E-01 -1.5551E-01
             1.2574E+00
 GRADIENT:   1.6037E-02  1.9730E-02 -1.2144E-02 -2.7471E-02  1.5001E-02  1.0221E-02  0.0000E+00  0.0000E+00 -6.8587E-03 -1.8704E-02
             3.1638E-05

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1094
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.8219E-04 -2.7009E-05  2.6437E-04 -1.0680E-02  3.3526E-03
 SE:             2.8963E-02  1.7773E-04  2.1701E-04  2.6442E-02  2.6731E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8396E-01  8.7921E-01  2.2315E-01  6.8627E-01  9.0019E-01

 ETASHRINKSD(%)  2.9703E+00  9.9405E+01  9.9273E+01  1.1417E+01  1.0448E+01
 ETASHRINKVR(%)  5.8525E+00  9.9996E+01  9.9995E+01  2.1530E+01  1.9805E+01
 EBVSHRINKSD(%)  2.8369E+00  9.9400E+01  9.9380E+01  8.5563E+00  1.0698E+01
 EBVSHRINKVR(%)  5.5934E+00  9.9996E+01  9.9996E+01  1.6380E+01  2.0252E+01
 RELATIVEINF(%)  9.4160E+01  6.9117E-04  4.6027E-04  3.8989E+01  5.5678E+00
 EPSSHRINKSD(%)  2.2465E+01
 EPSSHRINKVR(%)  3.9883E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1732.1558016406466     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -629.42956179503949     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.65
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.56
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1732.156       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.05E+00  2.76E-01  1.53E-01  9.90E-01  2.14E-01  9.92E-01  1.00E-02  1.00E-02  1.45E+00  7.75E-01  3.18E+00
 


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
+        9.66E+02
 
 TH 2
+        4.98E+01  5.01E+03
 
 TH 3
+       -5.03E+01  2.72E+03  1.61E+04
 
 TH 4
+       -1.65E+01 -3.36E+01 -6.76E+02  4.48E+02
 
 TH 5
+        6.56E+01 -8.93E+03 -1.77E+04 -2.53E+02  3.24E+04
 
 TH 6
+        2.21E+00 -1.61E+01  5.71E+01 -7.94E+00 -5.08E+00  1.80E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.21E+00 -1.73E+01  2.17E+02 -6.33E+00  4.45E+01 -1.54E+00  0.00E+00  0.00E+00  5.93E+01
 
 TH10
+       -1.73E+00 -5.42E+01 -4.55E+00  1.09E+01  4.03E+01  1.06E+00  0.00E+00  0.00E+00  5.51E-01  2.08E+02
 
 TH11
+       -1.80E+01 -8.32E+00 -3.46E+01 -5.97E+00  6.78E+01  3.11E+00  0.00E+00  0.00E+00  4.97E+00  1.23E+01  6.28E+01
 
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
 #CPUT: Total CPU Time in Seconds,       32.235
Stop Time:
Thu Sep 30 06:32:39 CDT 2021
