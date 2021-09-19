Sat Sep 18 08:32:31 CDT 2021
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
$DATA ../../../../data/spa/B/dat49.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1587.58044096147        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4999E+02 -4.7723E+01 -3.3303E+01 -7.4505E+00  4.1701E+01 -7.1637E+00 -4.1846E+01 -2.1918E+00 -2.3116E+01  1.0296E+00
            -2.4073E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1596.97998173042        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7311E-01  1.1606E+00  1.1097E+00  8.8686E-01  1.1051E+00  9.9943E-01  1.5720E+00  1.0544E+00  1.1984E+00  9.2642E-01
             1.1318E+00
 PARAMETER:  7.2738E-02  2.4891E-01  2.0410E-01 -2.0073E-02  1.9995E-01  9.9433E-02  5.5237E-01  1.5296E-01  2.8101E-01  2.3575E-02
             2.2377E-01
 GRADIENT:   8.1164E+01 -3.9203E+00  5.3718E+00 -1.9393E+01  2.0133E+01 -3.9569E+00  3.0843E+01 -5.5178E+00  2.2184E+01 -2.3452E+00
             2.7575E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1599.31619336999        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      162
 NPARAMETR:  9.7509E-01  1.0755E+00  1.2556E+00  9.6793E-01  1.0993E+00  1.0100E+00  1.6055E+00  1.4506E+00  1.1089E+00  8.6670E-01
             1.1229E+00
 PARAMETER:  7.4770E-02  1.7281E-01  3.2761E-01  6.7409E-02  1.9470E-01  1.0990E-01  5.7343E-01  4.7199E-01  2.0340E-01 -4.3061E-02
             2.1590E-01
 GRADIENT:   8.6623E+01  2.1234E+00 -3.5499E+00 -1.2834E+00  2.0833E+01  6.3531E-01  2.4701E+01  2.3496E+00  1.9599E+01 -8.4514E+00
             2.3220E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1603.27775378996        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  9.4690E-01  1.0178E+00  1.1703E+00  1.0137E+00  1.0352E+00  1.0019E+00  1.5384E+00  1.2262E+00  9.8348E-01  8.8280E-01
             1.0702E+00
 PARAMETER:  4.5438E-02  1.1760E-01  2.5728E-01  1.1365E-01  1.3461E-01  1.0194E-01  5.3077E-01  3.0390E-01  8.3342E-02 -2.4656E-02
             1.6784E-01
 GRADIENT:   2.6638E+01  8.8562E+00  2.2745E+00  1.3875E+01  9.6865E-01 -9.8005E-01  6.8148E+00 -1.2069E+00  4.3750E+00 -3.6016E+00
             5.7268E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1603.36123784683        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      344
 NPARAMETR:  9.4396E-01  1.0276E+00  1.1636E+00  1.0048E+00  1.0383E+00  1.0024E+00  1.5098E+00  1.2307E+00  9.8223E-01  8.9353E-01
             1.0656E+00
 PARAMETER:  4.2324E-02  1.2727E-01  2.5149E-01  1.0478E-01  1.3762E-01  1.0236E-01  5.1199E-01  3.0762E-01  8.2072E-02 -1.2579E-02
             1.6350E-01
 GRADIENT:  -1.7405E+01  2.6660E+00  1.4403E+00  3.1957E+00 -2.3474E-01 -6.8039E+00  1.5328E+00 -1.0042E+00  2.4048E+00 -2.7179E+00
             4.0739E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1603.54581829565        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      523
 NPARAMETR:  9.4736E-01  1.0302E+00  1.1984E+00  1.0038E+00  1.0537E+00  1.0102E+00  1.4963E+00  1.2851E+00  9.7789E-01  9.2060E-01
             1.0602E+00
 PARAMETER:  4.5921E-02  1.2974E-01  2.8098E-01  1.0382E-01  1.5228E-01  1.1019E-01  5.0301E-01  3.5087E-01  7.7638E-02  1.7275E-02
             1.5844E-01
 GRADIENT:  -9.1022E+00  1.3618E+00  8.3979E-01  1.6563E+00 -2.6922E-01 -3.5271E+00  9.3439E-01 -5.7401E-01  1.2791E+00 -1.2646E+00
             2.2017E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1603.56408959678        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:      692
 NPARAMETR:  9.4746E-01  1.0302E+00  1.1983E+00  1.0038E+00  1.0537E+00  1.0193E+00  1.4962E+00  1.2854E+00  9.7762E-01  9.2092E-01
             1.0600E+00
 PARAMETER:  4.6029E-02  1.2971E-01  2.8090E-01  1.0379E-01  1.5232E-01  1.1914E-01  5.0294E-01  3.5110E-01  7.7371E-02  1.7616E-02
             1.5824E-01
 GRADIENT:  -8.7082E+00  1.2721E+00  8.1464E-01  1.5893E+00 -2.5426E-01  6.1699E-03  9.1223E-01 -5.5692E-01  1.2390E+00 -1.2237E+00
             2.1640E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1603.60328637807        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:      856            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5132E-01  1.0276E+00  1.1949E+00  1.0029E+00  1.0552E+00  1.0193E+00  1.4906E+00  1.2953E+00  9.6799E-01  9.3111E-01
             1.0542E+00
 PARAMETER:  5.0095E-02  1.2720E-01  2.7806E-01  1.0286E-01  1.5378E-01  1.1909E-01  4.9915E-01  3.5873E-01  6.7465E-02  2.8621E-02
             1.5282E-01
 GRADIENT:   3.8229E+01  2.5339E+00 -8.6614E-01  7.3313E+00  2.6841E+00  7.0519E+00  3.3121E+00  2.5369E-01  7.9719E-01 -1.3469E-01
             2.1094E-01

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1603.60409013357        NO. OF FUNC. EVALS.: 106
 CUMULATIVE NO. OF FUNC. EVALS.:      962
 NPARAMETR:  9.5132E-01  1.0279E+00  1.1954E+00  1.0029E+00  1.0551E+00  1.0193E+00  1.4907E+00  1.2949E+00  9.6801E-01  9.3125E-01
             1.0542E+00
 PARAMETER:  5.0091E-02  1.2739E-01  2.7831E-01  1.0287E-01  1.5369E-01  1.1909E-01  4.9921E-01  3.5853E-01  6.7478E-02  2.8703E-02
             1.5282E-01
 GRADIENT:  -7.8429E-02 -1.5296E+00 -1.0425E+00 -6.2453E-02  1.4584E+00  1.3473E-02 -6.0460E-02  1.3037E-01 -3.1244E-02 -1.6816E-01
             7.0401E-02
 NUMSIGDIG:         3.4         1.4         1.6         3.0         1.8         3.5         2.8         1.9         2.3         1.5
                    2.7

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      962
 NO. OF SIG. DIGITS IN FINAL EST.:  1.4

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.7803E-04  1.6674E-03 -3.6988E-02 -6.1583E-03 -3.5327E-02
 SE:             2.9835E-02  2.2430E-02  1.4703E-02  2.1745E-02  2.0111E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8722E-01  9.4074E-01  1.1879E-02  7.7702E-01  7.8990E-02

 ETASHRINKSD(%)  4.9330E-02  2.4858E+01  5.0744E+01  2.7152E+01  3.2625E+01
 ETASHRINKVR(%)  9.8636E-02  4.3536E+01  7.5738E+01  4.6932E+01  5.4606E+01
 EBVSHRINKSD(%)  4.6937E-01  2.5419E+01  5.4435E+01  2.7803E+01  2.9846E+01
 EBVSHRINKVR(%)  9.3653E-01  4.4376E+01  7.9238E+01  4.7876E+01  5.0784E+01
 RELATIVEINF(%)  9.8573E+01  3.4489E+00  2.9161E+00  3.2548E+00  1.1369E+01
 EPSSHRINKSD(%)  4.4652E+01
 EPSSHRINKVR(%)  6.9366E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1603.6040901335709     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -868.45326356983276     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.47
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.45
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1603.604       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.51E-01  1.03E+00  1.20E+00  1.00E+00  1.06E+00  1.02E+00  1.49E+00  1.30E+00  9.68E-01  9.31E-01  1.05E+00
 


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
+        1.17E+03
 
 TH 2
+       -7.93E+00  3.37E+09
 
 TH 3
+        5.62E+00  1.33E+09  5.23E+08
 
 TH 4
+       -7.81E+00  4.28E+09  1.69E+09  5.44E+09
 
 TH 5
+       -2.16E+00 -2.72E+09 -1.07E+09 -3.46E+09  2.20E+09
 
 TH 6
+       -1.33E-01 -1.73E+00  9.66E-01 -2.02E+00 -1.20E+00  1.87E+02
 
 TH 7
+        1.11E+00  2.27E+01  9.10E-01 -1.11E+01 -8.95E-01  1.32E-01  3.40E+01
 
 TH 8
+       -2.04E-01 -9.52E+08 -3.75E+08  7.45E+00  7.68E+08 -3.90E-01  1.81E+00  1.69E+01
 
 TH 9
+       -1.28E+00 -1.17E+01 -1.80E+09  6.36E+00  4.62E+00  9.35E-01  1.65E+01  3.70E+00  6.68E+01
 
 TH10
+        2.60E-01 -1.29E+00 -1.87E+09 -1.25E+01 -3.83E+09  3.24E-02  4.37E+00  7.66E+00  3.94E+00  6.73E+01
 
 TH11
+       -6.48E+00 -6.33E+00  1.08E+09 -2.85E+00 -1.07E+01  1.26E+00  1.62E+00  3.48E+00  6.31E+00  1.48E+01  1.93E+02
 
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
 #CPUT: Total CPU Time in Seconds,       18.980
Stop Time:
Sat Sep 18 08:32:52 CDT 2021
