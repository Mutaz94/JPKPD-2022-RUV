Thu Sep 30 00:07:54 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat29.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m29.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   395.371384477235        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7997E+02  8.0669E+01  9.2194E+01  4.4887E+01  2.5330E+02  2.6831E+01 -8.2981E+01 -3.4983E+01 -5.8088E+01 -1.4670E+02
            -4.6015E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1414.03134384769        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0113E+00  9.9069E-01  1.0226E+00  1.0961E+00  9.1599E-01  8.0277E-01  1.0425E+00  9.3203E-01  1.0687E+00  9.0856E-01
             5.2561E+00
 PARAMETER:  1.1119E-01  9.0645E-02  1.2238E-01  1.9173E-01  1.2248E-02 -1.1969E-01  1.4160E-01  2.9611E-02  1.6640E-01  4.1007E-03
             1.7594E+00
 GRADIENT:  -4.5212E+01 -1.2859E+01 -9.5594E+00 -1.2943E+01 -1.1713E+01 -4.1460E+01  6.7682E+00  5.5773E+00  1.7302E+01  1.9900E+01
             3.3101E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1452.30099777299        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8433E-01  5.7514E-01  6.7685E-01  1.4375E+00  6.0481E-01  9.3256E-01  1.6525E+00  1.2685E-01  1.1494E+00  3.9326E-01
             4.5884E+00
 PARAMETER:  8.4207E-02 -4.5314E-01 -2.9030E-01  4.6289E-01 -4.0284E-01  3.0182E-02  6.0231E-01 -1.9647E+00  2.3925E-01 -8.3327E-01
             1.6235E+00
 GRADIENT:  -7.0697E+01  2.9671E+01 -3.5987E+01  1.4032E+02  2.3138E+01 -7.8147E+00  9.9796E+00  3.3360E-01  3.3164E+01  5.8618E+00
             2.5370E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1501.22055088942        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.8324E-01  6.2646E-01  7.5299E-01  1.2421E+00  6.4152E-01  9.6242E-01  1.2938E+00  3.3914E-02  1.0184E+00  2.6684E-01
             3.4708E+00
 PARAMETER:  8.3100E-02 -3.6767E-01 -1.8371E-01  3.1678E-01 -3.4391E-01  6.1699E-02  3.5759E-01 -3.2839E+00  1.1821E-01 -1.2211E+00
             1.3444E+00
 GRADIENT:   9.8045E-01  1.3284E+01  1.4465E+01  1.2551E+01 -2.3405E+01  1.7450E+00 -7.7747E-01  1.1064E-02  4.4181E+00  4.8406E-01
             3.5631E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1503.17706097986        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:      370
 NPARAMETR:  9.9507E-01  4.9098E-01  7.3137E-01  1.3367E+00  6.0092E-01  9.6266E-01  1.5777E+00  5.5304E-02  9.6766E-01  1.6991E-01
             3.5169E+00
 PARAMETER:  9.5059E-02 -6.1136E-01 -2.1284E-01  3.9021E-01 -4.0929E-01  6.1948E-02  5.5597E-01 -2.7949E+00  6.7129E-02 -1.6725E+00
             1.3576E+00
 GRADIENT:  -3.0930E+00  9.8278E+00 -3.1281E-01  1.3482E+01 -1.0110E+01 -5.0636E-01  2.4610E-02  2.6479E-02  1.3606E+00  2.0152E-02
            -1.1653E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1505.47945338084        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      545
 NPARAMETR:  9.8714E-01  3.0346E-01  9.5150E-01  1.4728E+00  6.7314E-01  9.5575E-01  2.1669E+00  7.3433E-02  8.7233E-01  8.1048E-02
             3.5221E+00
 PARAMETER:  8.7056E-02 -1.0925E+00  5.0283E-02  4.8714E-01 -2.9581E-01  5.4745E-02  8.7331E-01 -2.5114E+00 -3.6587E-02 -2.4127E+00
             1.3591E+00
 GRADIENT:  -7.9787E+00  4.2577E+00  3.7417E+00  1.2006E+01 -6.4656E+00 -1.3629E+00 -3.2968E-02  5.5516E-02 -4.3329E+00  3.9805E-02
            -5.7382E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1507.08192126860        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      722
 NPARAMETR:  9.8171E-01  8.0217E-02  6.6591E-01  1.5209E+00  4.8842E-01  9.6106E-01  4.1954E+00  1.5823E-01  8.9774E-01  1.0000E-02
             3.4978E+00
 PARAMETER:  8.1539E-02 -2.4230E+00 -3.0659E-01  5.1928E-01 -6.1657E-01  6.0278E-02  1.5340E+00 -1.7437E+00 -7.8701E-03 -7.1875E+00
             1.3521E+00
 GRADIENT:  -4.3002E+00  6.5798E-01 -1.8455E+00  1.4220E+01 -7.1274E-01  4.7723E-01 -8.1158E-01  1.3862E-01 -7.0425E-01  0.0000E+00
            -5.9713E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1507.51716897134        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      898
 NPARAMETR:  9.8424E-01  3.8176E-02  6.3862E-01  1.5264E+00  4.6680E-01  9.5718E-01  6.4759E+00  6.4512E-02  8.9962E-01  1.0000E-02
             3.5227E+00
 PARAMETER:  8.4116E-02 -3.1655E+00 -3.4844E-01  5.2289E-01 -6.6184E-01  5.6240E-02  1.9681E+00 -2.6409E+00 -5.7877E-03 -9.8381E+00
             1.3592E+00
 GRADIENT:   5.0164E+00  4.9836E-01  1.4043E+00  7.5165E+00 -4.1005E+00 -5.4165E-01  1.3894E-01  1.7571E-02 -1.3062E-01  0.0000E+00
             1.5483E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1507.96589269404        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  9.7923E-01  1.0000E-02  6.5668E-01  1.5425E+00  4.7248E-01  9.5807E-01  1.3003E+01  1.0000E-02  8.8985E-01  1.0000E-02
             3.5147E+00
 PARAMETER:  7.9014E-02 -4.5384E+00 -3.2055E-01  5.3343E-01 -6.4976E-01  5.7168E-02  2.6652E+00 -4.7144E+00 -1.6700E-02 -1.4577E+01
             1.3570E+00
 GRADIENT:  -4.2432E+00 -1.8972E-02  6.7824E+00  1.3141E+01 -7.3430E+00  5.7138E-01  6.3464E-01  0.0000E+00 -2.4718E-01  0.0000E+00
             6.6237E+00

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1507.98579660309        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  9.8023E-01  1.0000E-02  6.5492E-01  1.5430E+00  4.7210E-01  9.5783E-01  1.3104E+01  1.0000E-02  8.9124E-01  1.0000E-02
             3.5214E+00
 PARAMETER:  7.9030E-02 -4.5626E+00 -3.2376E-01  5.3288E-01 -6.5203E-01  5.7172E-02  2.6768E+00 -4.7560E+00 -1.6074E-02 -1.4667E+01
             1.3569E+00
 GRADIENT:  -3.5188E+00  0.0000E+00 -3.3108E+02 -3.5744E+02 -4.2611E+00  1.6152E-01  6.9762E+01  0.0000E+00 -6.9567E-01  0.0000E+00
            -1.3917E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1167
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.6408E-03  4.0596E-03  8.7226E-05 -1.6729E-02  3.8982E-05
 SE:             2.8657E-02  4.2839E-03  1.9307E-04  2.7048E-02  2.9966E-04
 N:                     100         100         100         100         100

 P VAL.:         9.2657E-01  3.4332E-01  6.5142E-01  5.3625E-01  8.9650E-01

 ETASHRINKSD(%)  3.9969E+00  8.5648E+01  9.9353E+01  9.3867E+00  9.8996E+01
 ETASHRINKVR(%)  7.8341E+00  9.7940E+01  9.9996E+01  1.7892E+01  9.9990E+01
 EBVSHRINKSD(%)  3.6328E+00  8.8094E+01  9.9326E+01  9.3316E+00  9.9055E+01
 EBVSHRINKVR(%)  7.1336E+00  9.8582E+01  9.9995E+01  1.7792E+01  9.9991E+01
 RELATIVEINF(%)  9.2142E+01  1.0194E+00  2.5564E-04  3.8647E+01  4.9267E-04
 EPSSHRINKSD(%)  1.8486E+01
 EPSSHRINKVR(%)  3.3555E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1507.9857966030947     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -589.04726339842205     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    18.83
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.88
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1507.986       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  6.55E-01  1.54E+00  4.71E-01  9.58E-01  1.32E+01  1.00E-02  8.90E-01  1.00E-02  3.51E+00
 


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
+        1.19E+03
 
 TH 2
+        0.00E+00  2.68E+06
 
 TH 3
+        1.59E+01  0.00E+00  6.17E+04
 
 TH 4
+       -1.20E+01  0.00E+00 -3.05E+02  7.67E+03
 
 TH 5
+        1.01E+02  0.00E+00 -4.37E+04  1.31E+02  4.20E+03
 
 TH 6
+       -1.23E+01  0.00E+00  1.35E+05 -5.52E+01 -3.34E+00  2.25E+02
 
 TH 7
+       -1.38E+00  0.00E+00  1.19E+01  1.08E+01 -1.66E+01  1.42E+00  4.24E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.11E+01  0.00E+00  1.45E+05 -9.24E+00  6.05E+01  1.49E+02 -3.24E-02  0.00E+00  3.46E+05
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.52E+01  0.00E+00 -5.27E+01  5.25E+02  8.60E+01 -3.47E+00  1.22E-01  0.00E+00  1.02E+01  0.00E+00  2.61E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.778
Stop Time:
Thu Sep 30 00:08:23 CDT 2021
