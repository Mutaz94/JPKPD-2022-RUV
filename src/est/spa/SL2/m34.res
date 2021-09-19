Sat Sep 18 12:12:58 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat34.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m34.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1575.98013112698        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.5735E+02 -1.1715E+02 -3.7161E+01 -7.5649E+01  1.1947E+02 -3.2402E+01 -9.2553E+00 -1.3123E+01  3.4282E+01 -3.6180E+01
            -4.5248E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1598.75692674751        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.4803E-01  1.1877E+00  9.3394E-01  9.3373E-01  9.4507E-01  1.0552E+00  1.0684E+00  1.1785E+00  6.6389E-01  1.1336E+00
             1.0909E+00
 PARAMETER:  4.6628E-02  2.7202E-01  3.1661E-02  3.1431E-02  4.3502E-02  1.5373E-01  1.6620E-01  2.6421E-01 -3.0964E-01  2.2538E-01
             1.8699E-01
 GRADIENT:   2.8973E+01 -5.8864E+00  1.5868E+01 -3.6717E+01 -3.8747E+01 -6.6593E-01 -8.4890E+00 -5.3761E+00 -6.8084E+00  6.6937E+00
            -3.9170E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1601.54322507794        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.5104E-01  1.1351E+00  1.4891E+00  9.8184E-01  1.1562E+00  1.0772E+00  1.2224E+00  1.8843E+00  6.1512E-01  1.2046E+00
             1.0922E+00
 PARAMETER:  4.9804E-02  2.2671E-01  4.9816E-01  8.1675E-02  2.4515E-01  1.7439E-01  3.0080E-01  7.3357E-01 -3.8594E-01  2.8614E-01
             1.8816E-01
 GRADIENT:   4.0095E+01 -1.1013E+01  5.7823E-01 -2.7912E+01  1.3669E+01  8.6781E+00  6.8756E+00 -2.9739E+00 -5.1851E-01 -1.1196E+01
            -4.3148E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1602.64175319672        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      255
 NPARAMETR:  9.3908E-01  1.0376E+00  1.6779E+00  1.0596E+00  1.1656E+00  1.0571E+00  1.2698E+00  1.9543E+00  6.0926E-01  1.2855E+00
             1.0983E+00
 PARAMETER:  3.7142E-02  1.3690E-01  6.1754E-01  1.5792E-01  2.5324E-01  1.5550E-01  3.3885E-01  7.7003E-01 -3.9551E-01  3.5114E-01
             1.9375E-01
 GRADIENT:  -1.4962E+01 -5.1528E+00 -2.6263E+00 -9.5641E+00  5.7885E+00 -4.3096E+00 -1.2405E-02 -1.2149E+00 -1.3150E+00 -2.3318E+00
            -1.5079E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1602.84105978558        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      419
 NPARAMETR:  9.4784E-01  1.0316E+00  1.7354E+00  1.0693E+00  1.1662E+00  1.0689E+00  1.2512E+00  1.9539E+00  6.3385E-01  1.3043E+00
             1.1011E+00
 PARAMETER:  4.6433E-02  1.3113E-01  6.5125E-01  1.6699E-01  2.5377E-01  1.6663E-01  3.2409E-01  7.6984E-01 -3.5594E-01  3.6564E-01
             1.9627E-01
 GRADIENT:   3.3595E+00 -4.2781E-01  4.6515E-01 -4.2023E+00 -8.7527E-01  2.3589E-01 -9.2393E-01 -2.9128E+00 -1.5018E-01 -7.6875E-02
            -3.9355E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.01982868181        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      600
 NPARAMETR:  9.4257E-01  9.3671E-01  1.8591E+00  1.1324E+00  1.1614E+00  1.0631E+00  1.3703E+00  1.9539E+00  5.9699E-01  1.3124E+00
             1.1000E+00
 PARAMETER:  4.0856E-02  3.4615E-02  7.2008E-01  2.2431E-01  2.4960E-01  1.6122E-01  4.1505E-01  7.6984E-01 -4.1585E-01  3.7182E-01
             1.9530E-01
 GRADIENT:  -6.2204E+00  2.1750E+00  1.5674E+00 -1.3866E+00 -5.9565E-01 -1.6899E+00 -6.6066E-02 -4.8155E+00 -1.2216E+00  5.1923E-01
            -7.6678E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.04674729205        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:      790
 NPARAMETR:  9.4280E-01  9.2472E-01  1.8698E+00  1.1405E+00  1.1597E+00  1.0632E+00  1.3838E+00  1.9539E+00  5.9565E-01  1.3112E+00
             1.1000E+00
 PARAMETER:  4.1098E-02  2.1731E-02  7.2583E-01  2.3143E-01  2.4814E-01  1.6133E-01  4.2480E-01  7.6984E-01 -4.1810E-01  3.7094E-01
             1.9533E-01
 GRADIENT:  -5.5952E+00  2.4014E+00  1.4896E+00 -6.4157E-01 -4.9616E-01 -1.6198E+00 -2.4281E-02 -4.7933E+00 -1.1735E+00  4.9384E-01
            -7.2820E-01

0ITERATION NO.:   34    OBJECTIVE VALUE:  -1603.04929012490        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:      939
 NPARAMETR:  9.4282E-01  9.2471E-01  1.8700E+00  1.1405E+00  1.1600E+00  1.0676E+00  1.3837E+00  1.9535E+00  5.9562E-01  1.3112E+00
             1.0999E+00
 PARAMETER:  4.1123E-02  2.1731E-02  7.2583E-01  2.3143E-01  2.4814E-01  1.6531E-01  4.2480E-01  7.6964E-01 -4.1810E-01  3.7094E-01
             1.9529E-01
 GRADIENT:  -4.2467E+00  2.4717E+05 -6.8090E+04 -2.3885E+01 -4.5562E-01 -5.2168E-02  1.1637E+05 -1.8087E+02  1.1823E+05  4.5223E-01
             2.5314E+05
 NUMSIGDIG:         8.7         3.3         3.3         7.6         2.3         2.8         3.3         6.2         3.3         9.1
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      939
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.2427E-03  4.1207E-03 -4.8792E-02 -1.8375E-02 -4.4424E-02
 SE:             2.9808E-02  2.2614E-02  1.6834E-02  1.8645E-02  2.1247E-02
 N:                     100         100         100         100         100

 P VAL.:         9.1337E-01  8.5541E-01  3.7516E-03  3.2438E-01  3.6546E-02

 ETASHRINKSD(%)  1.3867E-01  2.4241E+01  4.3603E+01  3.7537E+01  2.8820E+01
 ETASHRINKVR(%)  2.7715E-01  4.2605E+01  6.8193E+01  6.0983E+01  4.9333E+01
 EBVSHRINKSD(%)  4.8687E-01  2.4468E+01  5.3667E+01  3.8547E+01  2.2759E+01
 EBVSHRINKVR(%)  9.7136E-01  4.2949E+01  7.8533E+01  6.2235E+01  4.0338E+01
 RELATIVEINF(%)  9.8762E+01  2.9370E+00  5.5261E+00  1.8839E+00  2.1783E+01
 EPSSHRINKSD(%)  4.4296E+01
 EPSSHRINKVR(%)  6.8971E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.0492901249022     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -867.89846356116402     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.45
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.049       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.43E-01  9.25E-01  1.87E+00  1.14E+00  1.16E+00  1.07E+00  1.38E+00  1.95E+00  5.96E-01  1.31E+00  1.10E+00
 


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
+        1.09E+03
 
 TH 2
+       -1.16E+01  1.45E+09
 
 TH 3
+       -9.66E+07  9.85E+07  1.34E+07
 
 TH 4
+       -5.39E+04  5.06E+08 -3.37E+01  3.55E+08
 
 TH 5
+       -8.19E+03 -6.31E+01 -3.16E+07  2.96E+03  1.49E+08
 
 TH 6
+        2.16E+00  4.23E+03 -2.87E+02 -2.00E+00 -1.94E+00  1.71E+02
 
 TH 7
+        2.23E+08  2.27E+08 -4.24E+02  7.97E+07  7.31E+07  6.64E+02  4.14E+01
 
 TH 8
+       -3.43E-01  8.90E+07 -2.17E+01  1.48E+01  5.08E+02 -1.77E-01 -1.40E+07  5.45E+06
 
 TH 9
+       -1.58E+00 -1.23E+01 -2.70E+00 -1.88E+08  6.98E+00  2.81E+08  3.48E+01  9.95E+00  3.99E+08
 
 TH10
+       -4.11E-01  2.75E+08 -1.87E+07 -9.63E+07 -8.83E+07  2.17E-01  1.08E+01 -2.92E+00  1.02E+08  5.22E+07
 
 TH11
+       -1.90E+04  6.22E+08  1.32E+03 -2.18E+08  2.00E+08  1.82E+03  9.79E+07 -3.83E+07  2.31E+08 -1.18E+08  1.76E+02
 
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
 #CPUT: Total CPU Time in Seconds,       19.136
Stop Time:
Sat Sep 18 12:13:19 CDT 2021
