Sat Sep 18 12:48:50 CDT 2021
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
$DATA ../../../../data/spa/SL3/dat36.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m36.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1633.17979280834        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.5375E+01 -7.1732E+01 -4.6110E+01 -1.9884E+01  9.1856E+01 -3.6217E+01 -8.2351E+00  4.8059E+00  2.1328E+01 -1.2354E+01
            -5.1840E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1644.38326580210        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9236E-01  1.1311E+00  1.0333E+00  9.5319E-01  1.0035E+00  1.0717E+00  1.0989E+00  9.5304E-01  8.2225E-01  9.9769E-01
             1.1339E+00
 PARAMETER:  9.2332E-02  2.2315E-01  1.3276E-01  5.2062E-02  1.0346E-01  1.6927E-01  1.9428E-01  5.1899E-02 -9.5707E-02  9.7691E-02
             2.2566E-01
 GRADIENT:   9.5513E-01  1.1290E+01  8.6847E+00  1.6212E+00  2.9163E+00 -2.0288E+00 -5.0388E+00 -5.4899E-01 -4.8758E+00 -3.6327E+00
             8.6604E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1645.86433554144        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  1.0061E+00  1.0048E+00  9.0575E-01  1.0252E+00  8.8990E-01  1.0869E+00  1.3489E+00  7.7536E-01  6.9902E-01  8.8737E-01
             1.1155E+00
 PARAMETER:  1.0608E-01  1.0477E-01  1.0113E-03  1.2486E-01 -1.6645E-02  1.8331E-01  3.9928E-01 -1.5443E-01 -2.5808E-01 -1.9495E-02
             2.0932E-01
 GRADIENT:   2.9304E+01  1.6001E+01  5.3478E+00  1.1196E+01  2.5597E+00  4.5064E+00  3.7944E+00  5.3004E-01 -8.1392E+00 -4.4875E+00
            -4.6867E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1647.28130503086        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.9287E-01  9.9995E-01  6.8422E-01  1.0005E+00  7.7808E-01  1.0824E+00  1.2772E+00  3.3623E-01  7.7066E-01  8.2140E-01
             1.1141E+00
 PARAMETER:  9.2842E-02  9.9948E-02 -2.7948E-01  1.0054E-01 -1.5092E-01  1.7917E-01  3.4470E-01 -9.8997E-01 -1.6051E-01 -9.6741E-02
             2.0803E-01
 GRADIENT:  -2.7741E+00 -1.9080E+00 -1.0083E+01  9.9222E+00  1.2741E+01  8.0955E-01  5.6579E-01  1.0672E+00  1.0671E+00  3.5576E+00
             2.2615E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1647.39306911871        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  9.9395E-01  9.7364E-01  6.0669E-01  1.0007E+00  7.1485E-01  1.0813E+00  1.3005E+00  2.3200E-01  7.5800E-01  7.3027E-01
             1.1095E+00
 PARAMETER:  9.3934E-02  7.3290E-02 -3.9974E-01  1.0068E-01 -2.3568E-01  1.7812E-01  3.6276E-01 -1.3610E+00 -1.7708E-01 -2.1434E-01
             2.0389E-01
 GRADIENT:  -2.2895E+00 -1.0574E+00 -3.3748E+00  2.8449E+00  2.6809E+00 -2.4688E-01 -1.0834E-01  6.0214E-01  8.8910E-01  1.8052E+00
             8.2418E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1648.37098597962        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      487
 NPARAMETR:  1.0139E+00  8.2979E-01  7.0501E-01  1.1001E+00  7.1852E-01  1.1084E+00  1.5113E+00  1.3963E-01  7.1287E-01  8.0648E-01
             1.1062E+00
 PARAMETER:  1.1382E-01 -8.6588E-02 -2.4955E-01  1.9536E-01 -2.3057E-01  2.0291E-01  5.1299E-01 -1.8688E+00 -2.3846E-01 -1.1507E-01
             2.0094E-01
 GRADIENT:   3.5274E+00  1.0542E+00 -1.5409E+00  3.7497E+00  2.9676E-01  5.6740E-01 -7.6314E-02  1.2090E-01 -3.7202E-01  5.5232E-01
            -5.5270E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1648.40582412636        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      663
 NPARAMETR:  1.0127E+00  8.1816E-01  7.1113E-01  1.1041E+00  7.1936E-01  1.1077E+00  1.5287E+00  8.9903E-02  7.1052E-01  8.1055E-01
             1.1073E+00
 PARAMETER:  1.1266E-01 -1.0070E-01 -2.4089E-01  1.9905E-01 -2.2939E-01  2.0229E-01  5.2443E-01 -2.3090E+00 -2.4176E-01 -1.1005E-01
             2.0190E-01
 GRADIENT:   1.7597E+00 -1.7975E+00 -7.7839E-01 -2.2642E+00  1.4279E+00  4.3779E-01  1.2122E-02  4.2188E-02 -3.6468E-01  5.4447E-02
            -1.9814E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1648.43647996602        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      839
 NPARAMETR:  1.0119E+00  8.4530E-01  7.0318E-01  1.0897E+00  7.2364E-01  1.1069E+00  1.4871E+00  1.3719E-02  7.2078E-01  8.0711E-01
             1.1081E+00
 PARAMETER:  1.1187E-01 -6.8063E-02 -2.5215E-01  1.8592E-01 -2.2346E-01  2.0152E-01  4.9685E-01 -4.1889E+00 -2.2742E-01 -1.1430E-01
             2.0268E-01
 GRADIENT:  -1.0573E-01  1.6340E-01  2.0938E-01  6.9895E-02 -2.9895E-01 -2.4235E-02 -1.8280E-02  9.9031E-04 -9.2917E-03 -4.3273E-02
            -6.7255E-03

0ITERATION NO.:   38    OBJECTIVE VALUE:  -1648.43675830995        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      931
 NPARAMETR:  1.0120E+00  8.4477E-01  7.0312E-01  1.0899E+00  7.2358E-01  1.1069E+00  1.4880E+00  1.0000E-02  7.2062E-01  8.0738E-01
             1.1081E+00
 PARAMETER:  1.1193E-01 -6.8689E-02 -2.5223E-01  1.8608E-01 -2.2354E-01  2.0159E-01  4.9741E-01 -4.8132E+00 -2.2765E-01 -1.1396E-01
             2.0264E-01
 GRADIENT:   8.7122E-03 -3.4339E-02 -5.3414E-02 -1.3563E-02  6.4908E-02  5.3363E-03  6.0962E-03  0.0000E+00  4.9059E-03  1.4461E-02
             6.9898E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:      931
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0438E-04  1.0315E-02 -4.4006E-04 -1.4392E-02 -3.3474E-03
 SE:             2.9842E-02  2.2641E-02  1.9343E-04  2.2983E-02  2.2611E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9721E-01  6.4867E-01  2.2907E-02  5.3119E-01  8.8231E-01

 ETASHRINKSD(%)  2.5878E-02  2.4151E+01  9.9352E+01  2.3005E+01  2.4250E+01
 ETASHRINKVR(%)  5.1749E-02  4.2470E+01  9.9996E+01  4.0717E+01  4.2620E+01
 EBVSHRINKSD(%)  4.2875E-01  2.3871E+01  9.9391E+01  2.3232E+01  2.3144E+01
 EBVSHRINKVR(%)  8.5567E-01  4.2044E+01  9.9996E+01  4.1066E+01  4.0932E+01
 RELATIVEINF(%)  9.8790E+01  5.8561E+00  3.2020E-04  6.2093E+00  4.4145E+00
 EPSSHRINKSD(%)  4.2935E+01
 EPSSHRINKVR(%)  6.7436E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1648.4367583099477     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -913.28593174620949     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    10.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.08
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1648.437       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  8.45E-01  7.03E-01  1.09E+00  7.24E-01  1.11E+00  1.49E+00  1.00E-02  7.21E-01  8.07E-01  1.11E+00
 


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
+        8.78E+02
 
 TH 2
+       -7.15E+00  4.25E+02
 
 TH 3
+        1.80E+01  2.13E+02  9.12E+02
 
 TH 4
+       -1.31E+01  3.63E+02 -4.59E+02  1.08E+03
 
 TH 5
+       -4.67E+00 -4.23E+02 -1.20E+03  5.24E+02  1.91E+03
 
 TH 6
+       -2.38E-01 -2.90E+00  4.74E+00 -3.00E+00 -1.91E+00  1.60E+02
 
 TH 7
+        1.37E+00  3.46E+01 -1.17E+01 -1.62E+01  4.93E+00  2.55E-01  3.55E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.69E+00 -2.35E+01 -4.22E+01  1.71E+01  2.05E+01 -8.32E-01  1.70E+01  0.00E+00  1.46E+02
 
 TH10
+       -3.15E+00 -5.47E+00 -7.52E+01 -3.32E+01 -4.41E+01 -9.76E-02  9.43E+00  0.00E+00  1.94E+01  1.06E+02
 
 TH11
+       -5.91E+00 -1.23E+01 -3.38E+01 -4.57E+00  6.38E+00  1.66E+00  3.23E+00  0.00E+00  1.65E+01  2.51E+01  1.77E+02
 
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
 #CPUT: Total CPU Time in Seconds,       16.803
Stop Time:
Sat Sep 18 12:49:09 CDT 2021
