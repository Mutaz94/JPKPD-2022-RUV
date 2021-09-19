Sat Sep 18 00:55:14 CDT 2021
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
$DATA ../../../../data/int/A2/dat60.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m60.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2276.53492817609        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7669E+01  1.2443E+02  1.6818E+02 -6.6350E+01  2.0603E+02  8.0446E+00 -1.2544E+02 -2.3615E+02 -1.3553E+01 -9.7127E+01
            -2.7231E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3082.80516161180        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0510E+00  8.1366E-01  8.1540E-01  1.1649E+00  7.4903E-01  9.2156E-01  1.1341E+00  9.7410E-01  7.8617E-01  8.8061E-01
             1.9436E+00
 PARAMETER:  1.4974E-01 -1.0621E-01 -1.0408E-01  2.5260E-01 -1.8897E-01  1.8312E-02  2.2588E-01  7.3763E-02 -1.4058E-01 -2.7145E-02
             7.6455E-01
 GRADIENT:   1.2098E+02  4.1344E+01  2.6867E+01  3.8294E+01 -9.2254E+00 -2.4608E+01 -1.0735E+01  6.1362E+00 -3.1850E+01  3.1093E+00
            -7.6324E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3097.13003299687        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0457E+00  5.9981E-01  4.9051E-01  1.3140E+00  4.8463E-01  1.0134E+00  1.4209E+00  7.7328E-01  7.9584E-01  5.8904E-01
             1.8771E+00
 PARAMETER:  1.4465E-01 -4.1115E-01 -6.1231E-01  3.7310E-01 -6.2437E-01  1.1332E-01  4.5132E-01 -1.5712E-01 -1.2836E-01 -4.2926E-01
             7.2975E-01
 GRADIENT:   9.1701E+01  1.3135E+02 -1.1610E+01  2.8082E+02  1.2522E+01  1.0963E+01 -9.3956E+00 -2.6883E+00 -6.9660E+01 -2.5011E+01
            -2.3369E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3121.29303822522        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      232
 NPARAMETR:  9.5686E-01  3.3041E-01  2.2895E-01  1.3359E+00  2.6147E-01  1.2379E+00  1.5986E+00  1.3951E+00  1.0124E+00  2.4029E-01
             1.8025E+00
 PARAMETER:  5.5904E-02 -1.0074E+00 -1.3742E+00  3.8959E-01 -1.2414E+00  3.1338E-01  5.6913E-01  4.3296E-01  1.1235E-01 -1.3259E+00
             6.8915E-01
 GRADIENT:  -6.5747E+01  1.0996E+02  1.4977E+01  3.4482E+02 -1.1251E+02  6.9211E+01  1.1449E+00 -1.4854E+01 -9.1088E+01 -2.6414E+00
             1.7411E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3173.86282439578        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.7175E-01  2.4962E-01  1.7977E-01  1.0702E+00  2.2560E-01  1.0213E+00  1.5557E+00  1.4220E+00  1.1928E+00  2.0496E-01
             1.6393E+00
 PARAMETER:  7.1342E-02 -1.2878E+00 -1.6161E+00  1.6787E-01 -1.3890E+00  1.2104E-01  5.4195E-01  4.5207E-01  2.7630E-01 -1.4850E+00
             5.9430E-01
 GRADIENT:  -6.2049E+01  4.3707E+00 -7.7368E+00  6.3950E+00  1.6772E+01  1.1253E+01 -2.1746E+00 -1.3127E+01 -8.3339E+00 -1.4699E+01
             2.4715E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3182.05107301503        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:      460
 NPARAMETR:  1.0120E+00  2.6827E-01  2.0228E-01  1.1052E+00  2.4057E-01  9.8129E-01  1.5219E+00  1.4546E+00  1.1770E+00  3.4060E-01
             1.6066E+00
 PARAMETER:  1.1193E-01 -1.2158E+00 -1.4981E+00  2.0003E-01 -1.3248E+00  8.1108E-02  5.1998E-01  4.7471E-01  2.6295E-01 -9.7706E-01
             5.7411E-01
 GRADIENT:   1.3525E+01  1.7369E+00  2.5021E+00  7.8159E+00  1.7174E+01 -2.2200E+00 -8.3570E-01  1.0599E+00  3.5233E+00 -2.3385E+00
            -1.8573E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3188.83599947310        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  1.0057E+00  2.0626E-01  1.3600E-01  9.7271E-01  1.8448E-01  9.9864E-01  1.3592E+00  1.3136E+00  1.2168E+00  6.9881E-01
             1.6136E+00
 PARAMETER:  1.0571E-01 -1.4786E+00 -1.8951E+00  7.2328E-02 -1.5902E+00  9.8636E-02  4.0689E-01  3.7274E-01  2.9621E-01 -2.5837E-01
             5.7849E-01
 GRADIENT:  -2.8200E+00 -7.6613E+00  2.7784E+00  9.6658E+00  5.0025E+00  5.7567E+00 -2.2147E+00  1.5790E-01 -4.5497E-01 -1.3226E-01
            -1.6060E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3190.50570414300        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:      826             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0079E+00  1.9452E-01  1.2091E-01  9.1094E-01  1.7079E-01  9.8360E-01  1.3530E+00  1.2421E+00  1.2467E+00  7.7486E-01
             1.6202E+00
 PARAMETER:  1.0789E-01 -1.5372E+00 -2.0127E+00  6.7193E-03 -1.6673E+00  8.3464E-02  4.0234E-01  3.1679E-01  3.2051E-01 -1.5507E-01
             5.8255E-01
 GRADIENT:   1.2810E+01  1.8177E+01  1.7208E+01  1.1502E+01  9.9041E+01  1.0116E-01 -8.0057E-01  1.5969E+01  1.9700E+00  2.5295E-01
            -5.1345E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3190.53364843343        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  1.0050E+00  1.9451E-01  1.2089E-01  9.0608E-01  1.7074E-01  9.8340E-01  1.3551E+00  1.2421E+00  1.2467E+00  7.7441E-01
             1.6224E+00
 PARAMETER:  1.0496E-01 -1.5373E+00 -2.0128E+00  1.3764E-03 -1.6676E+00  8.3259E-02  4.0384E-01  3.1679E-01  3.2051E-01 -1.5565E-01
             5.8392E-01
 GRADIENT:   5.7448E+00  1.8881E+01  1.8423E+01  5.2666E+00  1.0040E+02  3.0562E-02 -3.8105E-01  1.6099E+01  2.0303E+00  9.9137E-02
            -2.3661E+00

0ITERATION NO.:   44    OBJECTIVE VALUE:  -3190.54352195149        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1035
 NPARAMETR:  1.0058E+00  1.9452E-01  1.2090E-01  9.0557E-01  1.7074E-01  9.8411E-01  1.3566E+00  1.2421E+00  1.2467E+00  7.7496E-01
             1.6231E+00
 PARAMETER:  1.0581E-01 -1.5373E+00 -2.0128E+00  8.0251E-04 -1.6676E+00  8.3886E-02  4.0502E-01  3.1679E-01  3.2051E-01 -1.5503E-01
             5.8435E-01
 GRADIENT:  -4.6310E+03 -6.3226E+02 -4.8103E+02 -9.7982E+03  5.7336E+02 -6.7699E-01  2.4175E+03 -3.0874E+03 -7.2554E-01 -3.4500E-01
            -8.3747E+02
 NUMSIGDIG:         3.3         3.3         3.3         3.3         3.3         1.5         3.3         3.3         7.2         1.7
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1035
 NO. OF SIG. DIGITS IN FINAL EST.:  1.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.5816E-03  1.3568E-02  1.8484E-02  2.7970E-03  2.6726E-02
 SE:             2.9775E-02  2.5806E-02  1.9826E-02  2.7044E-02  2.3712E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3091E-01  5.9905E-01  3.5119E-01  9.1763E-01  2.5969E-01

 ETASHRINKSD(%)  2.4887E-01  1.3546E+01  3.3580E+01  9.4003E+00  2.0563E+01
 ETASHRINKVR(%)  4.9712E-01  2.5258E+01  5.5883E+01  1.7917E+01  3.6897E+01
 EBVSHRINKSD(%)  7.1034E-01  1.2609E+01  2.7103E+01  8.6906E+00  2.1250E+01
 EBVSHRINKVR(%)  1.4156E+00  2.3628E+01  4.6860E+01  1.6626E+01  3.7984E+01
 RELATIVEINF(%)  9.8584E+01  2.5859E+01  9.8267E+00  3.2649E+01  1.0770E+01
 EPSSHRINKSD(%)  2.2763E+01
 EPSSHRINKVR(%)  4.0344E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3190.5435219514911     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1536.4541621830804     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.35
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.57
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3190.544       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.95E-01  1.21E-01  9.06E-01  1.71E-01  9.84E-01  1.36E+00  1.24E+00  1.25E+00  7.75E-01  1.62E+00
 


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
+        2.16E+07
 
 TH 2
+       -6.83E+02  2.74E+06
 
 TH 3
+       -8.33E+02 -3.41E+06  4.09E+06
 
 TH 4
+        4.61E+04  1.61E+04  1.91E+04  2.99E+07
 
 TH 5
+        7.24E+02  6.04E+03  4.17E+04 -1.85E+04  3.09E+06
 
 TH 6
+       -3.55E+03 -1.27E+03 -1.55E+03 -4.18E+03  1.32E+03  2.10E+02
 
 TH 7
+        3.75E+02 -8.71E+02 -1.12E+03 -8.93E+03  9.61E+02  6.88E+02  8.11E+05
 
 TH 8
+       -5.21E+02  1.65E+04  2.03E+04  1.25E+04 -1.74E+04 -9.61E+02  1.13E+06  1.59E+06
 
 TH 9
+        5.76E+06  1.10E+03  1.51E+03  6.76E+06 -1.02E+03  3.98E-01 -1.11E+06  1.55E+06  3.07E+06
 
 TH10
+       -2.15E+04 -7.61E+03 -9.31E+03 -2.53E+04  8.34E+03 -3.27E-01  4.18E+03 -5.82E+03  1.91E+01  1.63E+02
 
 TH11
+       -2.30E+02 -5.49E+03 -6.75E+03  5.18E+03  5.80E+03 -3.97E+02  4.71E+05  6.56E+05  6.48E+05 -2.41E+03  2.71E+05
 
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
 #CPUT: Total CPU Time in Seconds,       38.038
Stop Time:
Sat Sep 18 00:55:53 CDT 2021
