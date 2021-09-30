Wed Sep 29 07:18:11 CDT 2021
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
$DATA ../../../../data/int/TD2/dat34.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -3322.01797063723        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.0717E+02 -9.8865E+01  2.3396E+01  1.1295E+02  1.6480E+02 -7.1378E+00 -5.2001E+01 -4.7859E+01 -4.1231E+01 -2.6260E+01
            -8.9408E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3474.10056183065        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      127
 NPARAMETR:  9.5759E-01  1.1970E+00  9.6002E-01  9.6606E-01  9.6898E-01  1.0647E+00  1.0896E+00  1.0434E+00  1.0426E+00  1.0849E+00
             1.4586E+00
 PARAMETER:  5.6663E-02  2.7986E-01  5.9201E-02  6.5469E-02  6.8488E-02  1.6265E-01  1.8580E-01  1.4245E-01  1.4175E-01  1.8151E-01
             4.7751E-01
 GRADIENT:   1.2358E+01  2.2742E+01 -3.7580E+00  4.2929E+01 -5.5800E+01 -4.3369E+00  4.6519E+00 -1.5248E+00  2.0716E+00 -1.6840E+01
             1.1235E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3497.17018916833        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.5746E-01  1.7133E+00  1.5603E+00  6.9330E-01  1.4241E+00  1.0555E+00  7.5418E-01  2.1603E+00  1.0958E+00  1.6078E+00
             1.4156E+00
 PARAMETER:  5.6529E-02  6.3841E-01  5.4487E-01 -2.6630E-01  4.5356E-01  1.5400E-01 -1.8213E-01  8.7025E-01  1.9146E-01  5.7487E-01
             4.4755E-01
 GRADIENT:   1.1997E+01  1.1662E+02  2.6186E+01  4.9480E+01 -3.5497E+01 -8.4191E+00 -1.8050E+01 -1.8426E+01  4.3207E+00  8.1861E+00
             7.1903E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3498.73078403838        NO. OF FUNC. EVALS.: 205
 CUMULATIVE NO. OF FUNC. EVALS.:      510             RESET HESSIAN, TYPE I
 NPARAMETR:  9.5738E-01  1.7183E+00  1.6042E+00  6.9058E-01  1.4387E+00  1.0658E+00  7.7704E-01  2.3108E+00  1.0248E+00  1.6226E+00
             1.4148E+00
 PARAMETER:  5.6443E-02  6.4131E-01  5.7262E-01 -2.7022E-01  4.6377E-01  1.6376E-01 -1.5227E-01  9.3761E-01  1.2454E-01  5.8402E-01
             4.4697E-01
 GRADIENT:   2.0798E+02  5.6710E+02  3.0293E+01  1.0090E+02  5.7902E+01  2.7249E+01 -1.3803E+01 -1.3210E+01  5.1962E+00  3.5419E+01
             8.2935E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3499.83686919035        NO. OF FUNC. EVALS.: 172
 CUMULATIVE NO. OF FUNC. EVALS.:      682            RESET HESSIAN, TYPE II
 NPARAMETR:  9.5726E-01  1.7154E+00  1.5788E+00  6.8813E-01  1.4441E+00  1.0766E+00  7.8327E-01  2.3048E+00  1.0229E+00  1.5648E+00
             1.4129E+00
 PARAMETER:  5.6324E-02  6.3964E-01  5.5668E-01 -2.7377E-01  4.6748E-01  1.7377E-01 -1.4428E-01  9.3498E-01  1.2266E-01  5.4777E-01
             4.4566E-01
 GRADIENT:   2.0808E+02  5.5785E+02  2.6546E+01  9.7651E+01  6.1471E+01  3.3366E+01 -1.2004E+01 -8.4901E+00  4.4408E+00  2.3807E+01
             7.6855E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3500.53559962302        NO. OF FUNC. EVALS.: 150
 CUMULATIVE NO. OF FUNC. EVALS.:      832
 NPARAMETR:  9.5765E-01  1.7041E+00  1.5770E+00  6.8528E-01  1.4515E+00  1.0740E+00  7.8447E-01  2.2833E+00  1.0232E+00  1.5595E+00
             1.4179E+00
 PARAMETER:  5.6726E-02  6.3301E-01  5.5552E-01 -2.7792E-01  4.7263E-01  1.7136E-01 -1.4275E-01  9.2561E-01  1.2293E-01  5.4435E-01
             4.4915E-01
 GRADIENT:   2.0767E+02  5.2972E+02  2.8322E+01  8.5868E+01  6.8321E+01  3.1607E+01 -1.1583E+01 -1.2677E+01  5.7629E+00  2.3249E+01
             8.7265E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3502.24792560989        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  9.5769E-01  1.6459E+00  1.5725E+00  6.8345E-01  1.4595E+00  1.1337E+00  7.8448E-01  2.2983E+00  1.0233E+00  1.5749E+00
             1.4076E+00
 PARAMETER:  5.6772E-02  5.9827E-01  5.5266E-01 -2.8060E-01  4.7806E-01  2.2545E-01 -1.4273E-01  9.3218E-01  1.2303E-01  5.5419E-01
             4.4191E-01
 GRADIENT:   1.1853E+01  1.8998E+00  2.1339E+01 -7.4385E+00 -9.5923E+00  1.8776E+01 -1.5834E+01 -1.4670E+01  4.5864E+00  3.9617E+00
             6.8355E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3504.46525572884        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1150
 NPARAMETR:  9.5612E-01  1.6299E+00  1.5154E+00  6.8391E-01  1.4865E+00  1.0764E+00  7.9490E-01  2.3944E+00  1.0140E+00  1.5497E+00
             1.3998E+00
 PARAMETER:  5.5130E-02  5.8852E-01  5.1570E-01 -2.7993E-01  4.9643E-01  1.7366E-01 -1.2953E-01  9.7313E-01  1.1386E-01  5.3805E-01
             4.3634E-01
 GRADIENT:   2.1061E+02  3.9031E+02  1.9532E+01  4.5092E+01  1.1232E+02  3.4096E+01 -1.0099E+01 -6.1509E+00  8.6065E+00  2.4317E+01
             7.0615E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3504.72817466493        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1328
 NPARAMETR:  9.5602E-01  1.6298E+00  1.5128E+00  6.8395E-01  1.4868E+00  1.0763E+00  8.0041E-01  2.4174E+00  1.0126E+00  1.5498E+00
             1.3994E+00
 PARAMETER:  5.5025E-02  5.8843E-01  5.1397E-01 -2.7986E-01  4.9661E-01  1.7357E-01 -1.2263E-01  9.8268E-01  1.1256E-01  5.3814E-01
             4.3601E-01
 GRADIENT:   9.6224E+00 -2.4969E+01  1.3575E+01 -1.1886E+01  1.3054E+01 -5.5938E-01 -1.4564E+01 -1.0902E+01  7.3135E+00 -1.4274E+00
             6.3385E+01

0ITERATION NO.:   41    OBJECTIVE VALUE:  -3504.72817466493        NO. OF FUNC. EVALS.:  27
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  9.5602E-01  1.6298E+00  1.5128E+00  6.8395E-01  1.4868E+00  1.0763E+00  8.0041E-01  2.4174E+00  1.0126E+00  1.5498E+00
             1.3994E+00
 PARAMETER:  5.5025E-02  5.8843E-01  5.1397E-01 -2.7986E-01  4.9661E-01  1.7357E-01 -1.2263E-01  9.8268E-01  1.1256E-01  5.3814E-01
             4.3601E-01
 GRADIENT:  -3.7264E+04  6.3015E+03  3.5103E+01 -2.4976E+01  1.1371E+01  2.1474E+04  3.0380E+04  3.7221E+03  2.2863E+00 -3.6782E+00
             4.5032E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1355
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.6455E-03 -1.9933E-02 -4.3707E-02  3.4771E-02 -3.5528E-02
 SE:             2.9841E-02  2.4436E-02  1.8590E-02  2.1385E-02  2.6286E-02
 N:                     100         100         100         100         100

 P VAL.:         9.0277E-01  4.1467E-01  1.8721E-02  1.0397E-01  1.7650E-01

 ETASHRINKSD(%)  2.7616E-02  1.8135E+01  3.7720E+01  2.8357E+01  1.1940E+01
 ETASHRINKVR(%)  5.5225E-02  3.2982E+01  6.1212E+01  4.8673E+01  2.2454E+01
 EBVSHRINKSD(%)  5.1171E-01  2.1962E+01  4.6942E+01  3.0950E+01  8.7837E+00
 EBVSHRINKVR(%)  1.0208E+00  3.9100E+01  7.1848E+01  5.2321E+01  1.6796E+01
 RELATIVEINF(%)  9.8975E+01  1.4242E+01  1.9464E+01  1.0786E+01  5.0154E+01
 EPSSHRINKSD(%)  2.2156E+01
 EPSSHRINKVR(%)  3.9403E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3504.7281746649287     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1850.6388148965179     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    43.67
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    15.23
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3504.728       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.56E-01  1.63E+00  1.51E+00  6.84E-01  1.49E+00  1.08E+00  8.00E-01  2.42E+00  1.01E+00  1.55E+00  1.40E+00
 


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
+        2.04E+07
 
 TH 2
+        8.30E+00  1.01E+05
 
 TH 3
+        6.18E-02  1.20E+01  1.54E+05
 
 TH 4
+       -5.26E+00 -5.08E+05 -4.26E+01  5.09E+06
 
 TH 5
+       -7.68E-01 -5.55E+01 -2.36E+01  6.60E+05  3.42E+05
 
 TH 6
+       -1.41E+02  1.14E+01  3.48E-01 -2.44E+00 -5.67E-01  2.67E+06
 
 TH 7
+       -4.96E+01  2.89E+01 -3.69E-01  2.60E+01 -3.32E+00  5.08E+06  9.67E+06
 
 TH 8
+        4.17E+05  4.02E+04 -1.16E+01  7.97E+00  5.40E+04  5.31E+00  1.71E+00  1.60E+04
 
 TH 9
+        8.55E+06  8.54E+05  1.05E+06  2.61E+01 -1.11E+06  7.63E-02  4.75E+01 -3.50E+05  1.43E+07
 
 TH10
+        3.67E-01 -7.79E+00  1.44E+05 -5.84E+05  1.51E+05  1.90E-01  1.35E+00 -4.78E+04  3.37E+00  2.68E+05
 
 TH11
+       -9.95E+00 -2.27E+01  4.12E+02  1.71E+01 -7.95E+00  3.70E+00  1.53E+01  1.47E+01  1.34E+06  6.46E+00  5.01E+05
 
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
 #CPUT: Total CPU Time in Seconds,       59.000
Stop Time:
Wed Sep 29 07:19:12 CDT 2021
