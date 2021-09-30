Wed Sep 29 22:54:43 CDT 2021
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
$DATA ../../../../data/spa1/A1/dat91.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m91.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1939.49139317951        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.3833E+02 -3.7730E+01 -5.7842E+01  2.8483E+01  8.8289E+01  7.5479E+00 -3.5207E+00  1.7084E+01  2.6271E+00  2.3103E+01
            -2.6314E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1963.36982923894        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  9.9611E-01  1.0434E+00  1.3077E+00  1.0335E+00  1.0335E+00  1.1511E+00  9.9310E-01  7.7006E-01  1.0098E+00  6.1452E-01
             1.5876E+00
 PARAMETER:  9.6104E-02  1.4250E-01  3.6824E-01  1.3295E-01  1.3299E-01  2.4075E-01  9.3079E-02 -1.6129E-01  1.0976E-01 -3.8692E-01
             5.6221E-01
 GRADIENT:   2.5208E+02  3.8855E+01  3.3691E+01  3.2207E+01 -2.0641E+01  6.5955E+01 -5.8747E-02 -2.3595E+00  1.3972E+00 -1.1222E+01
             1.1979E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1972.06729715001        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:      182
 NPARAMETR:  9.7469E-01  9.1261E-01  5.8356E-01  1.1004E+00  6.5297E-01  1.0451E+00  1.1740E+00  3.6634E-01  1.0299E+00  2.8699E-01
             1.5022E+00
 PARAMETER:  7.4367E-02  8.5503E-03 -4.3861E-01  1.9571E-01 -3.2622E-01  1.4412E-01  2.6042E-01 -9.0419E-01  1.2948E-01 -1.1483E+00
             5.0690E-01
 GRADIENT:   6.6598E+01  3.5928E+01 -4.0325E+00  8.1924E+01 -3.9875E+01 -2.8363E+00 -1.1959E+01  2.4364E+00  2.6336E+01 -1.4344E+00
             1.0333E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1990.65684996826        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      358
 NPARAMETR:  9.2635E-01  8.1259E-01  7.7968E-01  1.1332E+00  7.4131E-01  1.0404E+00  1.3627E+00  2.3943E-01  9.0048E-01  5.7055E-01
             1.2874E+00
 PARAMETER:  2.3493E-02 -1.0753E-01 -1.4887E-01  2.2501E-01 -1.9934E-01  1.3957E-01  4.0946E-01 -1.3295E+00 -4.8313E-03 -4.6115E-01
             3.5259E-01
 GRADIENT:  -2.6137E+01  5.8908E+00 -3.6146E+00  1.0779E+01  7.4088E-01 -2.2824E+00 -1.1324E+00  8.5505E-01  1.9351E+00 -1.3259E+00
             1.6674E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1991.62243952367        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      533
 NPARAMETR:  9.3499E-01  6.4796E-01  8.4133E-01  1.2298E+00  7.1902E-01  1.0367E+00  1.6346E+00  1.0909E-01  8.4202E-01  6.2969E-01
             1.2815E+00
 PARAMETER:  3.2784E-02 -3.3393E-01 -7.2770E-02  3.0686E-01 -2.2987E-01  1.3602E-01  5.9139E-01 -2.1156E+00 -7.1950E-02 -3.6252E-01
             3.4801E-01
 GRADIENT:  -2.5164E+00  2.5773E+00  7.9875E-01  5.7335E+00 -3.1628E+00 -2.0863E+00 -1.9557E-01  1.1029E-01  3.2717E-01 -6.4143E-02
            -1.4312E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1991.68382243465        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      692
 NPARAMETR:  9.3627E-01  6.2143E-01  8.5297E-01  1.2412E+00  7.1849E-01  1.0424E+00  1.6931E+00  1.1353E-02  8.3320E-01  6.3728E-01
             1.2831E+00
 PARAMETER:  3.4152E-02 -3.7573E-01 -5.9027E-02  3.1607E-01 -2.3060E-01  1.4155E-01  6.2655E-01 -4.3782E+00 -8.2482E-02 -3.5054E-01
             3.4930E-01
 GRADIENT:   2.3177E+02  3.2219E+01  3.3749E+00  2.1194E+02  1.5080E+01  3.5339E+01  1.5153E+01  2.0936E-03  5.1505E+00  4.6110E-01
             2.8691E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1991.68873232834        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  9.3513E-01  6.2657E-01  8.5086E-01  1.2384E+00  7.1923E-01  1.0413E+00  1.6801E+00  1.0000E-02  8.3449E-01  6.3821E-01
             1.2837E+00
 PARAMETER:  3.2933E-02 -3.6750E-01 -6.1511E-02  3.1379E-01 -2.2957E-01  1.4046E-01  6.1886E-01 -1.4363E+01 -8.0934E-02 -3.4908E-01
             3.4976E-01
 GRADIENT:  -1.4742E+00 -6.6662E-01 -1.1173E-01 -1.6461E+00  7.7802E-01 -1.5429E-01  1.2018E-01  0.0000E+00 -6.2434E-02 -7.9317E-02
             2.0367E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1991.69188407781        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      987
 NPARAMETR:  9.3570E-01  6.2940E-01  8.4937E-01  1.2377E+00  7.1898E-01  1.0416E+00  1.6732E+00  1.0000E-02  8.3549E-01  6.3803E-01
             1.2837E+00
 PARAMETER:  3.3542E-02 -3.6299E-01 -6.3256E-02  3.1324E-01 -2.2992E-01  1.4080E-01  6.1471E-01 -1.2580E+01 -7.9733E-02 -3.4937E-01
             3.4975E-01
 GRADIENT:  -3.4199E-01 -1.2195E-02  6.2226E-02  1.2293E-01 -4.3696E-02 -4.1404E-02  3.6742E-02  0.0000E+00 -5.2043E-03  1.3187E-02
             3.4245E-02

0ITERATION NO.:   37    OBJECTIVE VALUE:  -1991.69913999437        NO. OF FUNC. EVALS.:  67
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  9.3630E-01  6.4090E-01  8.4341E-01  1.2290E+00  7.1988E-01  1.0424E+00  1.6545E+00  1.0000E-02  8.3953E-01  6.3463E-01
             1.2836E+00
 PARAMETER:  3.4404E-02 -3.4447E-01 -6.9866E-02  3.0725E-01 -2.2903E-01  1.4145E-01  6.0016E-01 -2.2873E+01 -7.5908E-02 -3.5624E-01
             3.4974E-01
 GRADIENT:   1.9298E-01  6.4545E-02  2.2473E-01  1.0337E+00 -4.0465E-01 -1.4815E-02 -1.6886E-01  0.0000E+00 -9.7000E-02 -6.2001E-02
             3.0038E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1054
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.7029E-04  1.8328E-02 -3.6882E-04 -1.5407E-02 -7.4557E-04
 SE:             2.9767E-02  2.0288E-02  2.1090E-04  2.5612E-02  2.0429E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7936E-01  3.6633E-01  8.0329E-02  5.4747E-01  9.7089E-01

 ETASHRINKSD(%)  2.7574E-01  3.2032E+01  9.9293E+01  1.4197E+01  3.1559E+01
 ETASHRINKVR(%)  5.5072E-01  5.3803E+01  9.9995E+01  2.6379E+01  5.3159E+01
 EBVSHRINKSD(%)  5.6457E-01  3.3500E+01  9.9273E+01  1.3588E+01  3.0542E+01
 EBVSHRINKVR(%)  1.1260E+00  5.5777E+01  9.9995E+01  2.5330E+01  5.1756E+01
 RELATIVEINF(%)  9.8447E+01  5.1384E+00  4.9420E-04  1.1986E+01  4.1302E+00
 EPSSHRINKSD(%)  3.0921E+01
 EPSSHRINKVR(%)  5.2281E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1991.6991399943709     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1072.7606067896982     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    14.88
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.04
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1991.699       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.37E-01  6.41E-01  8.44E-01  1.23E+00  7.20E-01  1.04E+00  1.65E+00  1.00E-02  8.39E-01  6.34E-01  1.28E+00
 


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
+        1.15E+03
 
 TH 2
+       -1.38E+01  4.65E+02
 
 TH 3
+        7.32E+00  3.02E+02  8.67E+02
 
 TH 4
+       -2.25E+00  3.04E+02 -2.71E+02  7.91E+02
 
 TH 5
+        7.18E-01 -6.37E+02 -1.38E+03  3.30E+02  2.57E+03
 
 TH 6
+        1.98E+00 -3.01E+00  3.52E+00 -2.81E+00 -5.07E+00  1.78E+02
 
 TH 7
+        1.36E+00  3.84E+01  1.72E-01 -8.20E+00 -5.52E+00  3.64E-02  2.23E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.15E+00 -2.33E+01 -1.47E+01  1.31E+01  1.27E+01 -4.44E-01  1.11E+01  0.00E+00  1.65E+02
 
 TH10
+       -1.20E-02  3.38E+00 -7.76E+01 -2.67E+01 -1.47E+01  1.01E+00  9.07E+00  0.00E+00  8.23E+00  1.21E+02
 
 TH11
+       -1.02E+01 -9.14E+00 -3.42E+01 -1.43E+01  1.47E+01  2.52E+00  2.95E+00  0.00E+00  9.98E+00  3.62E+01  2.58E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.995
Stop Time:
Wed Sep 29 22:55:07 CDT 2021
