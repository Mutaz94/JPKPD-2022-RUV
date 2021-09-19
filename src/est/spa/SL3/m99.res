Sat Sep 18 13:10:50 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat99.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m99.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1624.94099245751        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.4964E+01 -9.2719E+01 -4.8227E+01 -3.7851E+01  1.0670E+02  3.3430E+01 -2.3557E+01 -4.9938E+00  3.7685E+00 -2.2094E+01
            -6.4390E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1637.91461822853        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.9848E-01  1.1746E+00  1.0402E+00  9.6030E-01  9.8340E-01  8.7138E-01  1.1964E+00  1.0666E+00  9.5045E-01  1.0557E+00
             1.2025E+00
 PARAMETER:  9.8481E-02  2.6089E-01  1.3940E-01  5.9488E-02  8.3266E-02 -3.7679E-02  2.7933E-01  1.6448E-01  4.9176E-02  1.5422E-01
             2.8444E-01
 GRADIENT:  -3.4871E+00  3.6318E+01  1.1555E+01  2.2169E+01 -2.8919E+01 -1.9462E+01  4.3281E+00 -4.7558E+00  1.1053E+00  3.6773E+00
             2.1014E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1639.52337005694        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.9824E-01  1.0852E+00  1.3090E+00  1.0255E+00  1.0509E+00  8.9539E-01  1.2798E+00  1.6249E+00  8.7203E-01  1.0492E+00
             1.2003E+00
 PARAMETER:  9.8238E-02  1.8181E-01  3.6925E-01  1.2517E-01  1.4969E-01 -1.0501E-02  3.4668E-01  5.8546E-01 -3.6932E-02  1.4804E-01
             2.8255E-01
 GRADIENT:   1.1364E+00  2.3415E+01 -3.9031E+00  3.0783E+01 -4.9547E+00 -7.6680E+00  5.4147E+00  3.2564E+00 -2.7227E+00 -1.3753E+00
             1.9994E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1640.79219638141        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      239
 NPARAMETR:  9.9603E-01  1.0784E+00  1.4002E+00  1.0141E+00  1.0842E+00  9.0920E-01  1.2029E+00  1.6457E+00  9.1474E-01  1.1053E+00
             1.1523E+00
 PARAMETER:  9.6026E-02  1.7552E-01  4.3658E-01  1.1402E-01  1.8085E-01  4.8103E-03  2.8472E-01  5.9819E-01  1.0886E-02  2.0013E-01
             2.4179E-01
 GRADIENT:   7.2842E-01  4.8213E+00 -6.0650E-01  5.8514E+00  8.2928E-02 -1.4855E+00  1.5008E+00  1.0153E-01 -6.5516E-01  8.6577E-02
             3.7418E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1640.84676042995        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  9.9695E-01  1.0786E+00  1.4103E+00  1.0140E+00  1.0872E+00  9.1064E-01  1.1976E+00  1.6585E+00  9.1777E-01  1.1081E+00
             1.1505E+00
 PARAMETER:  9.6950E-02  1.7564E-01  4.4383E-01  1.1389E-01  1.8360E-01  6.3935E-03  2.8029E-01  6.0591E-01  1.4191E-02  2.0266E-01
             2.4022E-01
 GRADIENT:  -2.8836E+01 -7.7469E-01 -9.2126E-01 -7.4094E-01 -5.2927E-01 -4.2707E+00  3.6564E-01 -3.6455E-02 -8.9128E-01 -6.2610E-02
             2.8586E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1640.88738461008        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  9.9735E-01  1.0789E+00  1.4130E+00  1.0139E+00  1.0880E+00  9.2021E-01  1.1965E+00  1.6619E+00  9.2383E-01  1.1088E+00
             1.1504E+00
 PARAMETER:  9.7347E-02  1.7597E-01  4.4575E-01  1.1379E-01  1.8434E-01  1.6850E-02  2.7941E-01  6.0795E-01  2.0777E-02  2.0329E-01
             2.4007E-01
 GRADIENT:  -2.7193E+01 -9.5328E-01 -9.3167E-01 -6.9916E-01 -4.9638E-01 -1.5667E-01  6.9943E-01  6.5491E-03 -4.7948E-02 -3.5448E-02
             2.9431E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1640.92459929555        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      732
 NPARAMETR:  9.9866E-01  1.0791E+00  1.4138E+00  1.0139E+00  1.0881E+00  9.2058E-01  1.1954E+00  1.6619E+00  9.2416E-01  1.1088E+00
             1.1494E+00
 PARAMETER:  9.8661E-02  1.7613E-01  4.4632E-01  1.1385E-01  1.8440E-01  1.7251E-02  2.7851E-01  6.0799E-01  2.1129E-02  2.0332E-01
             2.3921E-01
 GRADIENT:  -2.3754E+01 -8.0909E-01 -8.1230E-01 -5.4719E-01 -6.2705E-01  5.7920E-02  6.3819E-01 -4.8518E-02 -6.8956E-02 -7.1436E-02
             2.5391E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1640.94763214657        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:      917
 NPARAMETR:  9.9957E-01  1.0792E+00  1.4145E+00  1.0139E+00  1.0882E+00  9.2043E-01  1.1946E+00  1.6618E+00  9.2440E-01  1.1091E+00
             1.1489E+00
 PARAMETER:  9.9592E-02  1.7621E-01  4.4667E-01  1.1387E-01  1.8453E-01  1.7185E-02  2.7790E-01  6.0799E-01  2.1372E-02  2.0348E-01
             2.3872E-01
 GRADIENT:   5.2414E+05  2.9745E+05 -5.8665E+04  4.6032E+05 -1.4203E+05  5.7506E-02  9.4306E+04  8.6111E+04 -2.6208E+05 -1.2880E+05
            -2.1956E+05
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         2.6         3.3         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:      917
 NO. OF SIG. DIGITS IN FINAL EST.:  2.6

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.1412E-03 -8.5564E-03 -3.9070E-02  6.9464E-04 -4.1453E-02
 SE:             2.9761E-02  2.0074E-02  1.5328E-02  2.2465E-02  2.1096E-02
 N:                     100         100         100         100         100

 P VAL.:         7.5872E-01  6.6993E-01  1.0803E-02  9.7533E-01  4.9415E-02

 ETASHRINKSD(%)  2.9719E-01  3.2750E+01  4.8650E+01  2.4738E+01  2.9327E+01
 ETASHRINKVR(%)  5.9350E-01  5.4774E+01  7.3632E+01  4.3356E+01  5.0054E+01
 EBVSHRINKSD(%)  6.5099E-01  3.2769E+01  5.3549E+01  2.6060E+01  2.5724E+01
 EBVSHRINKVR(%)  1.2977E+00  5.4800E+01  7.8423E+01  4.5328E+01  4.4831E+01
 RELATIVEINF(%)  9.7692E+01  1.1188E+00  2.6773E+00  1.3965E+00  1.4685E+01
 EPSSHRINKSD(%)  4.4750E+01
 EPSSHRINKVR(%)  6.9474E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1640.9476321465736     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -905.79680558283542     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    12.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.54
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1640.948       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.08E+00  1.41E+00  1.01E+00  1.09E+00  9.21E-01  1.19E+00  1.66E+00  9.24E-01  1.11E+00  1.15E+00
 


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
+        1.31E+09
 
 TH 2
+        2.59E+03  3.62E+08
 
 TH 3
+       -7.80E+02  3.83E+04  3.28E+07
 
 TH 4
+        4.28E+03 -2.59E+04  7.86E+03  9.83E+08
 
 TH 5
+       -6.53E+08  7.87E+03 -2.50E+03  1.32E+04  3.25E+08
 
 TH 6
+        9.20E+03  4.82E+03 -1.45E+03  7.94E+03 -4.56E+03  2.31E+02
 
 TH 7
+        1.49E+03 -5.04E+03  1.53E+03 -8.34E+03  4.56E+03  2.76E+03  1.19E+08
 
 TH 8
+        1.30E+08 -2.42E+03 -2.05E+07  1.12E+08 -6.47E+07  9.06E+02  3.91E+07  1.28E+07
 
 TH 9
+       -4.00E+04 -2.10E+04  6.32E+03 -3.46E+04  7.06E+08 -9.91E+03 -1.20E+04 -1.40E+08  1.53E+09
 
 TH10
+       -5.81E+08 -1.16E+03  3.44E+02 -1.91E+03  2.89E+08 -4.06E+03 -6.59E+02 -5.75E+07  6.28E+08  2.57E+08
 
 TH11
+       -1.82E+03  8.88E+03 -2.65E+04  1.82E+04 -5.54E+03 -3.34E+03  3.51E+03  3.84E+03  1.46E+04  8.12E+02  1.74E+08
 
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
 #CPUT: Total CPU Time in Seconds,       18.619
Stop Time:
Sat Sep 18 13:11:11 CDT 2021
