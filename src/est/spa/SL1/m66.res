Sat Sep 25 10:39:18 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat66.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1624.15112092105        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4616E+02 -7.5944E+01 -5.1299E+00 -8.1309E+01  1.5764E+01  1.0569E+01 -9.2794E+00  3.9512E+00  2.2572E+01 -1.4267E+01
            -1.9064E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1633.91112456174        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7070E-01  1.1580E+00  9.8826E-01  9.5491E-01  1.0854E+00  9.3256E-01  1.1254E+00  9.4741E-01  7.5002E-01  1.1587E+00
             1.0876E+00
 PARAMETER:  7.0265E-02  2.4670E-01  8.8190E-02  5.3863E-02  1.8197E-01  3.0179E-02  2.1814E-01  4.5978E-02 -1.8765E-01  2.4732E-01
             1.8397E-01
 GRADIENT:   7.7142E+01 -5.9217E+00 -2.5910E+00 -1.3695E+01  9.8691E+00 -1.2370E+01 -3.9691E+00  1.2586E+00 -6.3569E+00  1.6461E+00
             1.2434E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1635.46142242829        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.6915E-01  9.2139E-01  9.8699E-01  1.1129E+00  9.5273E-01  9.3133E-01  1.4881E+00  7.0817E-01  6.4439E-01  1.0481E+00
             1.0868E+00
 PARAMETER:  6.8663E-02  1.8133E-02  8.6906E-02  2.0699E-01  5.1578E-02  2.8855E-02  4.9752E-01 -2.4506E-01 -3.3946E-01  1.4702E-01
             1.8319E-01
 GRADIENT:   7.5402E+01  1.6894E+01  4.7459E+00  1.5787E+01  9.0321E-01 -1.2521E+01  6.5380E+00 -6.2999E-01 -6.1933E+00 -2.0720E+00
             1.0832E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1636.81805020043        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.4139E-01  9.5299E-01  8.0977E-01  1.0738E+00  8.6856E-01  9.5348E-01  1.3366E+00  4.6788E-01  7.3159E-01  9.6976E-01
             1.0526E+00
 PARAMETER:  3.9605E-02  5.1848E-02 -1.1100E-01  1.7119E-01 -4.0918E-02  5.2365E-02  3.9016E-01 -6.5955E-01 -2.1253E-01  6.9290E-02
             1.5130E-01
 GRADIENT:   8.3522E-01  1.3052E+00 -6.9309E+00  7.5662E+00  6.9023E+00 -2.9591E+00  2.7998E-01  1.1043E+00 -7.9669E-02  1.8195E+00
             2.7169E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1637.71857521947        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.4264E-01  8.1327E-01  7.2055E-01  1.1388E+00  7.4958E-01  9.5161E-01  1.5185E+00  1.9448E-01  6.9500E-01  8.4907E-01
             1.0512E+00
 PARAMETER:  4.0933E-02 -1.0670E-01 -2.2774E-01  2.2994E-01 -1.8824E-01  5.0402E-02  5.1773E-01 -1.5374E+00 -2.6384E-01 -6.3616E-02
             1.4996E-01
 GRADIENT:   2.8099E+00 -1.0535E-01 -2.6733E+00  4.0614E+00  4.4150E+00 -3.7847E+00  4.3024E-01  2.0657E-01 -2.8422E-01 -6.5777E-01
             1.3084E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1637.87423237853        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  9.4161E-01  7.2481E-01  6.6127E-01  1.1741E+00  6.7580E-01  9.6298E-01  1.6468E+00  7.7781E-02  6.7846E-01  7.8027E-01
             1.0472E+00
 PARAMETER:  3.9839E-02 -2.2185E-01 -3.1360E-01  2.6054E-01 -2.9186E-01  6.2276E-02  5.9883E-01 -2.4539E+00 -2.8792E-01 -1.4812E-01
             1.4608E-01
 GRADIENT:   3.5540E-01 -7.6799E-01 -1.1495E+00 -1.1978E+00  2.3942E-01  8.4056E-01 -3.5435E-02  6.0752E-02  2.8028E-01  8.1295E-01
             5.7292E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1638.86918222467        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.5564E-01  6.3133E-01  7.2877E-01  1.2450E+00  6.8740E-01  9.6781E-01  1.8896E+00  3.8531E-02  6.5522E-01  8.3231E-01
             1.0534E+00
 PARAMETER:  5.4623E-02 -3.5993E-01 -2.1640E-01  3.1911E-01 -2.7483E-01  6.7281E-02  7.3635E-01 -3.1563E+00 -3.2279E-01 -8.3553E-02
             1.5201E-01
 GRADIENT:   7.3814E+00  3.1520E+00  1.2541E+00  2.6948E+00 -2.5260E+00  6.3604E-01  1.1450E-01  4.3665E-04 -2.0214E-01 -9.4875E-01
             3.2635E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1638.99560123142        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      693
 NPARAMETR:  9.5157E-01  5.6150E-01  7.4746E-01  1.2823E+00  6.7961E-01  9.6491E-01  2.0651E+00  1.3902E-02  6.4122E-01  8.5113E-01
             1.0536E+00
 PARAMETER:  5.0360E-02 -4.7714E-01 -1.9108E-01  3.4869E-01 -2.8623E-01  6.4284E-02  8.2518E-01 -4.1757E+00 -3.4439E-01 -6.1196E-02
             1.5220E-01
 GRADIENT:  -4.1076E-03  9.0324E-02  1.1160E-01  8.2989E-02 -1.5248E-01 -5.5208E-03 -2.5707E-03 -2.7556E-04 -2.9191E-02 -1.5990E-02
             1.7900E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1638.99563560418        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:      750
 NPARAMETR:  9.5156E-01  5.6029E-01  7.4754E-01  1.2829E+00  6.7937E-01  9.6492E-01  2.0683E+00  1.3620E-02  6.4111E-01  8.5125E-01
             1.0535E+00
 PARAMETER:  5.0344E-02 -4.7930E-01 -1.9096E-01  3.4914E-01 -2.8658E-01  6.4289E-02  8.2673E-01 -4.1962E+00 -3.4456E-01 -6.1052E-02
             1.5213E-01
 GRADIENT:  -1.1875E-03  5.2964E-03 -2.5588E-03  1.8298E-02 -1.2107E-03  3.7446E-03 -2.2570E-03 -2.6337E-04  2.4116E-03  3.4382E-03
            -5.7720E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      750
 NO. OF SIG. DIGITS IN FINAL EST.:  3.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.2613E-04  2.6516E-02 -6.9200E-04 -2.7922E-02  4.8345E-03
 SE:             2.9821E-02  2.1393E-02  3.1694E-04  2.3131E-02  2.3421E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8325E-01  2.1517E-01  2.9010E-02  2.2739E-01  8.3646E-01

 ETASHRINKSD(%)  9.6393E-02  2.8332E+01  9.8938E+01  2.2508E+01  2.1538E+01
 ETASHRINKVR(%)  1.9269E-01  4.8637E+01  9.9989E+01  3.9950E+01  3.8437E+01
 EBVSHRINKSD(%)  5.1279E-01  2.9120E+01  9.9020E+01  2.1409E+01  1.9415E+01
 EBVSHRINKVR(%)  1.0230E+00  4.9760E+01  9.9990E+01  3.8235E+01  3.5061E+01
 RELATIVEINF(%)  9.8376E+01  8.0790E+00  8.3211E-04  1.0902E+01  4.9955E+00
 EPSSHRINKSD(%)  4.2938E+01
 EPSSHRINKVR(%)  6.7440E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1638.9956356041760     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -903.84480904043778     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:     7.26
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.87
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1638.996       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  5.60E-01  7.48E-01  1.28E+00  6.79E-01  9.65E-01  2.07E+00  1.36E-02  6.41E-01  8.51E-01  1.05E+00
 


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
+        1.30E+03
 
 TH 2
+       -1.55E+01  4.43E+02
 
 TH 3
+        2.64E+01  2.11E+02  1.12E+03
 
 TH 4
+       -8.96E+00  3.76E+02 -4.06E+02  1.05E+03
 
 TH 5
+       -6.72E+00 -4.18E+02 -1.43E+03  4.14E+02  2.17E+03
 
 TH 6
+       -2.39E+00 -5.04E+00  5.44E+00 -2.79E+00 -5.21E+00  2.05E+02
 
 TH 7
+        1.87E+00  3.79E+01 -5.45E+00 -1.52E+01  3.90E+00  3.08E-01  1.80E+01
 
 TH 8
+        2.56E+00  1.96E+00 -5.14E+00  1.47E-01  2.55E+00 -7.65E+00 -5.65E-02 -1.83E+01
 
 TH 9
+        1.20E+00 -2.37E+01 -4.10E+01 -1.69E+01  4.93E+01 -1.85E+00  1.30E+01 -2.43E+00  1.96E+02
 
 TH10
+        1.33E-01 -5.55E+00 -8.76E+01 -3.15E+01 -3.91E+01 -5.90E-01  3.48E+00 -2.57E+00  1.88E+01  1.15E+02
 
 TH11
+       -8.21E+00 -8.05E+00 -4.72E+01 -9.70E+00  1.80E+01  5.78E+00  1.26E+00  6.76E-01  1.81E+01  2.53E+01  1.96E+02
 
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
 #CPUT: Total CPU Time in Seconds,       13.200
Stop Time:
Sat Sep 25 10:39:33 CDT 2021
