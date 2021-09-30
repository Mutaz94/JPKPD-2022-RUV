Thu Sep 30 00:04:32 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat22.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m22.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   170.325421704882        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2831E+02  9.9713E+01  1.6469E+02  6.7414E+01  2.3693E+02  7.7372E+01 -6.9548E+01 -1.3188E+02 -7.1741E+01 -1.9737E+02
            -3.9832E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1464.41353939080        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9311E-01  9.9116E-01  9.2787E-01  1.0954E+00  9.2149E-01  7.2341E-01  9.8329E-01  9.8927E-01  9.8062E-01  1.0609E+00
             4.5894E+00
 PARAMETER:  9.3090E-02  9.1125E-02  2.5131E-02  1.9114E-01  1.8238E-02 -2.2378E-01  8.3153E-02  8.9217E-02  8.0430E-02  1.5909E-01
             1.6238E+00
 GRADIENT:   4.8374E+01 -5.3269E+00 -1.2157E+01  5.6568E+00 -1.2229E+01 -9.4141E+00  1.2762E+01  7.3675E+00  2.2646E+01  2.7969E+01
             2.4196E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1490.95577654417        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  9.7138E-01  6.1506E-01  2.8429E-01  1.2806E+00  3.7485E-01  7.5601E-01  5.7975E-01  8.3808E-02  1.0727E+00  4.7275E-01
             4.0364E+00
 PARAMETER:  7.0966E-02 -3.8604E-01 -1.1577E+00  3.4733E-01 -8.8124E-01 -1.7970E-01 -4.4516E-01 -2.3792E+00  1.7014E-01 -6.4918E-01
             1.4954E+00
 GRADIENT:  -3.1006E+01  3.0370E+01 -5.0477E+01  1.9616E+02  5.5370E+01 -1.6852E+01  2.1632E+00  1.4669E-01 -2.0175E+01  1.0681E+01
             1.7633E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1525.20194730611        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      233
 NPARAMETR:  9.1991E-01  5.1340E-01  2.1561E-01  1.1609E+00  2.9580E-01  8.1825E-01  9.7253E-01  1.0000E-02  1.1596E+00  2.6295E-01
             3.0526E+00
 PARAMETER:  1.6519E-02 -5.6669E-01 -1.4343E+00  2.4920E-01 -1.1181E+00 -1.0059E-01  7.2142E-02 -5.3870E+00  2.4806E-01 -1.2358E+00
             1.2160E+00
 GRADIENT:  -1.0844E+02  2.5127E+01 -2.7224E+01  1.1852E+02  8.3360E+01 -2.4339E+00 -6.9286E+00  0.0000E+00 -3.4855E+01 -7.2310E+00
            -9.3818E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1533.61621668823        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      389
 NPARAMETR:  9.7696E-01  4.7465E-01  2.2159E-01  1.1018E+00  2.8515E-01  7.8508E-01  1.3304E+00  1.0000E-02  1.2274E+00  2.5689E-01
             2.9832E+00
 PARAMETER:  7.6691E-02 -6.4518E-01 -1.4069E+00  1.9698E-01 -1.1547E+00 -1.4196E-01  3.8551E-01 -7.7357E+00  3.0491E-01 -1.2591E+00
             1.1930E+00
 GRADIENT:   5.5423E+01  2.3776E+01  1.5898E+01  1.1110E+01 -2.7688E+01 -1.0138E+01  5.8997E+00  0.0000E+00 -7.9214E+00 -6.2223E+00
            -3.3217E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1541.22483779063        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      564
 NPARAMETR:  9.5989E-01  4.0061E-01  1.9501E-01  1.0765E+00  2.5173E-01  8.1128E-01  7.7920E-01  1.0000E-02  1.3374E+00  5.4072E-01
             2.9790E+00
 PARAMETER:  5.9068E-02 -8.1477E-01 -1.5347E+00  1.7372E-01 -1.2794E+00 -1.0914E-01 -1.4949E-01 -1.0185E+01  3.9071E-01 -5.1485E-01
             1.1916E+00
 GRADIENT:  -3.2353E+00  1.8832E+00  3.1717E+00 -1.0153E+00 -6.1444E+00  1.6162E-01  1.9963E+00  0.0000E+00 -2.1192E-01  4.2452E-01
            -8.2760E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1541.93334033544        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      739
 NPARAMETR:  9.6255E-01  4.0233E-01  1.8819E-01  1.0696E+00  2.4806E-01  8.1269E-01  2.5069E-01  1.0000E-02  1.3662E+00  5.8363E-01
             2.9762E+00
 PARAMETER:  6.1833E-02 -8.1048E-01 -1.5703E+00  1.6732E-01 -1.2941E+00 -1.0740E-01 -1.2836E+00 -9.4864E+00  4.1206E-01 -4.3849E-01
             1.1906E+00
 GRADIENT:   5.4882E+00  2.6488E+00  3.7340E-01  2.2900E+00 -2.1634E+00  1.2093E-01  1.0329E-01  0.0000E+00 -2.5001E-01 -6.2498E-01
            -3.1461E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1542.01522892047        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      914
 NPARAMETR:  9.6095E-01  3.9421E-01  1.8651E-01  1.0662E+00  2.4554E-01  8.1197E-01  9.6245E-02  1.0000E-02  1.3685E+00  5.9074E-01
             2.9854E+00
 PARAMETER:  6.0166E-02 -8.3087E-01 -1.5792E+00  1.6409E-01 -1.3043E+00 -1.0830E-01 -2.2409E+00 -8.8266E+00  4.1375E-01 -4.2637E-01
             1.1937E+00
 GRADIENT:   5.6723E-02 -6.9926E-01 -4.9938E-01  9.4826E-02  1.2849E+00  4.0736E-02  1.8541E-02  0.0000E+00 -1.1273E-01  1.7792E-01
             4.7835E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1542.02534148812        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  9.6090E-01  3.9468E-01  1.8616E-01  1.0655E+00  2.4533E-01  8.1207E-01  1.1046E-02  1.0000E-02  1.3706E+00  5.9122E-01
             2.9837E+00
 PARAMETER:  6.0116E-02 -8.2967E-01 -1.5811E+00  1.6344E-01 -1.3052E+00 -1.0816E-01 -4.4057E+00 -6.9463E+00  4.1527E-01 -4.2556E-01
             1.1932E+00
 GRADIENT:  -1.5701E-02 -2.3665E-02 -3.9642E-03  3.9880E-02  2.0336E-02  3.0907E-02  2.5141E-04  0.0000E+00  1.3404E-03 -2.4969E-03
            -1.6604E-02

0ITERATION NO.:   43    OBJECTIVE VALUE:  -1542.02537511465        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1183
 NPARAMETR:  9.6091E-01  3.9480E-01  1.8614E-01  1.0653E+00  2.4533E-01  8.1188E-01  1.0000E-02  1.0000E-02  1.3707E+00  5.9131E-01
             2.9838E+00
 PARAMETER:  6.0130E-02 -8.2937E-01 -1.5813E+00  1.6328E-01 -1.3052E+00 -1.0840E-01 -5.2327E+00 -6.2239E+00  4.1531E-01 -4.2541E-01
             1.1932E+00
 GRADIENT:   2.0715E-02  3.6480E-02  7.2152E-03 -6.4302E-02 -2.9967E-02 -5.0295E-02  0.0000E+00  0.0000E+00 -4.0685E-03  7.5403E-03
             2.9854E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1183
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6760E-03 -1.1259E-04  1.9033E-04 -1.0583E-02  1.9276E-03
 SE:             2.8554E-02  1.4583E-04  2.2069E-04  2.6859E-02  2.2206E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5319E-01  4.4007E-01  3.8845E-01  6.9357E-01  9.3083E-01

 ETASHRINKSD(%)  4.3398E+00  9.9511E+01  9.9261E+01  1.0017E+01  2.5607E+01
 ETASHRINKVR(%)  8.4912E+00  9.9998E+01  9.9995E+01  1.9031E+01  4.4658E+01
 EBVSHRINKSD(%)  4.0736E+00  9.9512E+01  9.9290E+01  7.9899E+00  2.5754E+01
 EBVSHRINKVR(%)  7.9813E+00  9.9998E+01  9.9995E+01  1.5341E+01  4.4875E+01
 RELATIVEINF(%)  9.1711E+01  2.7788E-04  3.6553E-04  5.4280E+01  1.7632E+00
 EPSSHRINKSD(%)  2.4675E+01
 EPSSHRINKVR(%)  4.3262E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1542.0253751146524     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -623.08684190997974     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    17.16
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1542.025       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.61E-01  3.95E-01  1.86E-01  1.07E+00  2.45E-01  8.12E-01  1.00E-02  1.00E-02  1.37E+00  5.91E-01  2.98E+00
 


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
+        1.69E+03
 
 TH 2
+       -7.43E+01  1.69E+03
 
 TH 3
+       -1.62E+02  3.43E+03  1.57E+04
 
 TH 4
+       -2.07E+01  1.27E+02 -4.07E+02  4.45E+02
 
 TH 5
+        3.09E+02 -6.06E+03 -1.92E+04 -3.02E+02  2.93E+04
 
 TH 6
+        5.61E+00 -1.79E+01  7.87E+01 -1.33E+01 -7.78E+00  2.56E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.64E+01 -3.85E+01  1.48E+02 -8.21E+00  5.10E+01 -1.11E+00  0.00E+00  0.00E+00  7.04E+01
 
 TH10
+       -1.01E+01 -6.16E+01 -2.82E+01  6.95E+00  1.68E+02  4.73E+00  0.00E+00  0.00E+00  2.37E+00  1.67E+02
 
 TH11
+       -2.69E+01 -6.80E+00 -5.80E+01 -5.50E+00  4.19E+01  4.28E+00  0.00E+00  0.00E+00  5.92E+00  2.49E+01  5.36E+01
 
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
 
 Elapsed finaloutput time in seconds:     1.45
 #CPUT: Total CPU Time in Seconds,       24.875
Stop Time:
Thu Sep 30 00:05:03 CDT 2021
