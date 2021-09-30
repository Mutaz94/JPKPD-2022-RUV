Thu Sep 30 04:16:03 CDT 2021
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
$DATA ../../../../data/spa2/B/dat49.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1831.96985892623        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9951E+02 -1.8428E+01  8.0962E+01  1.8753E+01  6.2183E+01  2.1455E+01 -5.0150E+01 -2.3935E+02 -5.3716E+01  2.7126E+00
            -8.4581E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2227.00484885749        NO. OF FUNC. EVALS.: 121
 CUMULATIVE NO. OF FUNC. EVALS.:      134
 NPARAMETR:  8.6004E-01  1.0019E+00  9.4249E-01  9.4386E-01  1.0358E+00  8.0102E-01  1.2467E+00  1.0939E+00  1.0232E+00  9.6735E-01
             1.8408E+00
 PARAMETER: -5.0773E-02  1.0187E-01  4.0766E-02  4.2224E-02  1.3520E-01 -1.2187E-01  3.2051E-01  1.8978E-01  1.2296E-01  6.6804E-02
             7.1019E-01
 GRADIENT:  -3.6449E+02 -1.1095E+02 -4.5087E+01 -4.3571E+01  7.2876E+01 -1.5047E+02 -3.1791E+01  1.1793E+01  1.2621E+01  1.1529E+01
             4.0771E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2265.06074282046        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  8.5064E-01  6.0991E-01  9.1804E-01  1.2295E+00  7.5882E-01  8.1938E-01  2.5605E+00  4.2756E-01  8.0474E-01  7.0343E-01
             1.6586E+00
 PARAMETER: -6.1768E-02 -3.9445E-01  1.4486E-02  3.0663E-01 -1.7599E-01 -9.9203E-02  1.0402E+00 -7.4966E-01 -1.1724E-01 -2.5179E-01
             6.0597E-01
 GRADIENT:  -3.7586E+02 -1.6703E+00  5.7778E-01  7.4051E+01 -1.0716E+01 -1.4451E+02  2.3960E+01 -8.5891E-01 -1.5809E+00  5.9218E+00
             3.2789E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2392.80508921223        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      491
 NPARAMETR:  9.0197E-01  8.3823E-01  7.9049E-01  1.0081E+00  8.3953E-01  9.0431E-01  1.7870E+00  5.8013E-01  9.8042E-01  7.2857E-01
             1.0639E+00
 PARAMETER: -3.1743E-03 -7.6465E-02 -1.3511E-01  1.0803E-01 -7.4909E-02 -5.8388E-04  6.8056E-01 -4.4450E-01  8.0229E-02 -2.1667E-01
             1.6196E-01
 GRADIENT:  -1.4465E+02 -4.6899E+01 -3.0559E+01 -5.8496E+01  3.5032E+01 -6.2398E+01  1.3296E+01 -8.8435E-01  5.7104E+00  4.5695E+00
             2.9863E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2399.02479935119        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:      680             RESET HESSIAN, TYPE I
 NPARAMETR:  9.1694E-01  8.9414E-01  8.2770E-01  1.0015E+00  8.7477E-01  9.2376E-01  1.6965E+00  6.7142E-01  9.8198E-01  7.5579E-01
             1.0503E+00
 PARAMETER:  1.3286E-02 -1.1893E-02 -8.9102E-02  1.0152E-01 -3.3791E-02  2.0693E-02  6.2856E-01 -2.9836E-01  8.1811E-02 -1.7999E-01
             1.4909E-01
 GRADIENT:   2.8442E+02  9.7194E+00 -1.6226E+01  4.6026E+01  5.2171E+01  1.6945E+01  7.9135E+01  4.9956E-01  1.5370E+01  4.3070E+00
             1.9787E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2400.73795354045        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      861
 NPARAMETR:  9.3000E-01  9.1362E-01  8.2770E-01  1.0015E+00  8.7477E-01  9.2375E-01  1.6965E+00  6.7143E-01  9.8198E-01  7.5579E-01
             1.0503E+00
 PARAMETER:  2.7424E-02  9.6588E-03 -8.9106E-02  1.0153E-01 -3.3795E-02  2.0689E-02  6.2859E-01 -2.9835E-01  8.1811E-02 -1.7999E-01
             1.4909E-01
 GRADIENT:  -5.9615E+01 -2.1329E+01 -1.8161E+01 -3.4836E+01  1.6341E+01 -4.5659E+01  1.1951E+01 -2.6641E-01  5.6413E+00  1.5028E+00
             1.7225E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2403.53240925492        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     1019
 NPARAMETR:  9.3734E-01  9.2474E-01  8.2770E-01  1.0015E+00  8.7448E-01  1.0075E+00  1.7001E+00  6.7077E-01  9.8230E-01  7.5579E-01
             1.0508E+00
 PARAMETER:  3.5287E-02  2.1756E-02 -8.9108E-02  1.0153E-01 -3.4130E-02  1.0746E-01  6.3069E-01 -2.9933E-01  8.2144E-02 -1.7999E-01
             1.4958E-01
 GRADIENT:   3.4755E+02  2.9586E+01 -9.8074E+00  5.4473E+01  3.4123E+01  5.9516E+01  8.1551E+01  3.4540E-01  1.5124E+01  3.8731E+00
             1.9437E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2404.95256777751        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1200
 NPARAMETR:  9.5198E-01  9.4796E-01  8.4427E-01  1.0204E+00  8.7012E-01  1.0250E+00  1.5985E+00  6.7849E-01  9.5756E-01  7.4586E-01
             1.0337E+00
 PARAMETER:  5.0790E-02  4.6559E-02 -6.9281E-02  1.2022E-01 -3.9126E-02  1.2470E-01  5.6907E-01 -2.8788E-01  5.6636E-02 -1.9322E-01
             1.3314E-01
 GRADIENT:   1.0527E+00  9.3102E+00  4.2375E+00 -3.9530E+00 -1.6829E+01  1.1131E+00  2.3380E-01 -2.3851E+00 -6.3670E-01 -4.3002E+00
            -1.6200E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2405.10482299569        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1396             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5200E-01  9.4608E-01  8.4443E-01  1.0207E+00  8.7152E-01  1.0238E+00  1.6007E+00  6.8900E-01  9.5763E-01  7.7447E-01
             1.0342E+00
 PARAMETER:  5.0811E-02  4.4576E-02 -6.9090E-02  1.2051E-01 -3.7515E-02  1.2350E-01  5.7042E-01 -2.7252E-01  5.6711E-02 -1.5557E-01
             1.3362E-01
 GRADIENT:   3.9406E+02  5.0975E+01  5.5849E+00  1.0156E+02  1.0343E+01  7.6166E+01  6.0306E+01 -8.1541E-01  1.0357E+01  2.3097E+00
             1.7098E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2405.12074244311        NO. OF FUNC. EVALS.: 173
 CUMULATIVE NO. OF FUNC. EVALS.:     1569
 NPARAMETR:  9.5222E-01  9.4608E-01  8.4443E-01  1.0209E+00  8.7192E-01  1.0240E+00  1.5907E+00  6.9041E-01  9.5759E-01  7.7448E-01
             1.0342E+00
 PARAMETER:  5.1040E-02  4.4576E-02 -6.9090E-02  1.2064E-01 -3.7054E-02  1.2376E-01  5.6418E-01 -2.7047E-01  5.6664E-02 -1.5557E-01
             1.3362E-01
 GRADIENT:   6.0094E+00 -1.2387E+00  3.0905E+05 -7.5466E+01  3.0898E+05 -2.4967E+05  3.2145E+00  1.1409E+05  3.0892E+05  9.6881E-01
            -7.8013E+02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1569
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         4.2139E-04 -6.3635E-03 -2.1877E-02  2.7067E-03 -9.9505E-03
 SE:             2.9856E-02  2.4746E-02  1.2922E-02  2.5811E-02  2.0965E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8874E-01  7.9707E-01  9.0466E-02  9.1648E-01  6.3506E-01

 ETASHRINKSD(%)  1.0000E-10  1.7096E+01  5.6708E+01  1.3530E+01  2.9764E+01
 ETASHRINKVR(%)  1.0000E-10  3.1270E+01  8.1258E+01  2.5229E+01  5.0669E+01
 EBVSHRINKSD(%)  3.4025E-01  1.5904E+01  5.9625E+01  1.4630E+01  2.9943E+01
 EBVSHRINKVR(%)  6.7934E-01  2.9278E+01  8.3699E+01  2.7119E+01  5.0920E+01
 RELATIVEINF(%)  9.9313E+01  2.8871E+01  7.1166E+00  3.4536E+01  1.1381E+01
 EPSSHRINKSD(%)  3.0101E+01
 EPSSHRINKVR(%)  5.1142E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2405.1207424431063     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1302.3945025974992     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.81
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2405.121       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.52E-01  9.46E-01  8.44E-01  1.02E+00  8.72E-01  1.02E+00  1.59E+00  6.90E-01  9.58E-01  7.74E-01  1.03E+00
 


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
+        8.52E+07
 
 TH 2
+        8.57E+07  8.63E+07
 
 TH 3
+        6.58E+00 -9.67E+07  2.17E+08
 
 TH 4
+       -1.86E+04  1.76E+02 -2.10E+02  5.09E+07
 
 TH 5
+       -9.30E+07 -1.38E+03  1.05E+08 -7.19E+07  1.02E+08
 
 TH 6
+        6.40E+07 -6.44E+07  3.30E+00  2.75E-01  2.53E+03  4.81E+07
 
 TH 7
+       -9.04E+06  9.10E+06 -1.50E+01 -6.99E+06 -9.87E+06  6.79E+06  9.59E+05
 
 TH 8
+        4.35E+07 -9.96E+02  4.91E+07  3.36E+07  5.37E+02  1.18E+03  4.62E+06  2.21E+07
 
 TH 9
+        1.68E+03  8.53E+07 -9.55E+07  6.55E+07  8.10E+00  7.92E+00  8.99E+06 -4.33E+07  1.68E+08
 
 TH10
+       -6.73E+07 -1.72E+01  3.03E+04  1.46E+04 -5.34E+01  3.91E-01 -2.73E+01  2.10E+02  1.12E+01  1.06E+08
 
 TH11
+       -1.16E+03  1.33E+03 -2.64E+04 -4.53E+07 -6.40E+07 -4.40E+07 -3.32E+02  2.99E+07  5.83E+07 -2.59E+02  4.06E+07
 
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
 #CPUT: Total CPU Time in Seconds,       38.748
Stop Time:
Thu Sep 30 04:16:43 CDT 2021
