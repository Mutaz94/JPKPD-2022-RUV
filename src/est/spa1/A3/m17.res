Thu Sep 30 00:02:14 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat17.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m17.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -118.302403210791        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7410E+02  3.9191E+01  2.3672E+02 -6.6144E+01  2.1477E+02  4.6049E+01 -7.5532E+01 -2.8282E+02 -1.1881E+02 -1.5001E+02
            -3.3513E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1376.18565422262        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.1184E+00  1.0228E+00  8.7758E-01  1.2226E+00  8.8497E-01  8.0270E-01  1.0650E+00  9.8776E-01  1.0896E+00  1.0749E+00
             6.6201E+00
 PARAMETER:  2.1193E-01  1.2255E-01 -3.0592E-02  3.0095E-01 -2.2203E-02 -1.1977E-01  1.6295E-01  8.7682E-02  1.8584E-01  1.7224E-01
             1.9901E+00
 GRADIENT:   2.3637E+01 -2.1200E+00 -1.2938E+01  1.6351E+01 -1.4466E+01 -1.3667E+01  9.9821E+00  7.4205E+00  3.0441E+01  2.2110E+01
             4.3975E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1531.39648941967        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  1.0285E+00  6.8608E-01  2.1461E-01  1.2112E+00  3.4964E-01  9.5838E-01  9.2171E-01  1.0000E-02  1.5941E+00  2.7722E-01
             3.8591E+00
 PARAMETER:  1.2809E-01 -2.7676E-01 -1.4389E+00  2.9159E-01 -9.5084E-01  5.7489E-02  1.8473E-02 -6.8219E+00  5.6633E-01 -1.1830E+00
             1.4504E+00
 GRADIENT:  -2.6932E+01 -2.6261E+01 -1.0437E+02  9.5547E+01  1.7132E+02  4.3832E+00 -2.9109E+00  0.0000E+00  5.4304E+01  1.2388E-01
             1.9210E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1578.01125842552        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  9.8237E-01  4.6034E-01  1.5794E-01  1.2495E+00  2.3381E-01  1.0260E+00  7.8956E-01  1.0000E-02  1.5199E+00  7.0827E-01
             2.7638E+00
 PARAMETER:  8.2211E-02 -6.7580E-01 -1.7455E+00  3.2278E-01 -1.3532E+00  1.2569E-01 -1.3628E-01 -7.3152E+00  5.1863E-01 -2.4493E-01
             1.1166E+00
 GRADIENT:  -1.1874E+02  9.6967E+01 -5.1694E+01  1.4833E+02 -1.6703E+01  1.1821E+01 -1.9765E+00  0.0000E+00 -1.1158E+01  5.3615E+00
             8.6198E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1605.26476271113        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      472
 NPARAMETR:  1.0319E+00  3.9842E-01  1.8819E-01  1.1103E+00  2.4116E-01  9.6005E-01  8.6538E-01  1.0000E-02  1.3103E+00  6.8683E-01
             2.6617E+00
 PARAMETER:  1.3140E-01 -8.2025E-01 -1.5703E+00  2.0463E-01 -1.3223E+00  5.9234E-02 -4.4582E-02 -9.0198E+00  3.7026E-01 -2.7567E-01
             1.0790E+00
 GRADIENT:  -2.2643E+01  4.1337E+01  2.7717E+01  1.9690E+01 -3.2608E+01  1.2929E+00 -3.5803E+00  0.0000E+00 -1.1959E+00 -3.6582E+00
            -2.5486E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1614.51677657039        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      648
 NPARAMETR:  1.0380E+00  2.7019E-01  1.5579E-01  1.0618E+00  2.0250E-01  9.5204E-01  1.7730E+00  1.0000E-02  1.3474E+00  6.6439E-01
             2.5636E+00
 PARAMETER:  1.3729E-01 -1.2086E+00 -1.7593E+00  1.5998E-01 -1.4970E+00  5.0853E-02  6.7269E-01 -1.2149E+01  3.9815E-01 -3.0889E-01
             1.0414E+00
 GRADIENT:  -7.0696E-01  7.8163E+00  6.0302E+00  1.9686E+01 -5.7689E+00  5.6285E-01  4.9031E+00  0.0000E+00 -1.1267E+01 -1.1109E+00
            -2.0982E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1619.66875115956        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      824
 NPARAMETR:  1.0331E+00  2.0611E-01  1.0054E-01  8.8269E-01  1.5795E-01  9.3691E-01  1.6457E+00  1.0000E-02  1.8235E+00  7.7062E-01
             2.4966E+00
 PARAMETER:  1.3256E-01 -1.4793E+00 -2.1972E+00 -2.4785E-02 -1.7455E+00  3.4836E-02  5.9814E-01 -1.7109E+01  7.0075E-01 -1.6056E-01
             1.0149E+00
 GRADIENT:  -1.0673E-02 -8.8033E+00 -1.9127E+01 -3.3938E+00  2.5427E+01  7.1733E-01  4.0840E+00  0.0000E+00 -1.7794E+00  4.5435E-01
             4.9865E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1621.33102183846        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  1.0340E+00  2.4049E-01  1.0610E-01  9.1546E-01  1.6383E-01  9.4005E-01  1.1868E+00  1.0000E-02  1.7949E+00  7.7012E-01
             2.5010E+00
 PARAMETER:  1.3341E-01 -1.3251E+00 -2.1433E+00  1.1674E-02 -1.7089E+00  3.8175E-02  2.7122E-01 -1.5461E+01  6.8497E-01 -1.6121E-01
             1.0167E+00
 GRADIENT:  -4.8491E-01  2.9064E-02 -2.9662E-01  8.4753E-01  3.7158E-01 -2.6190E-01  2.4095E-01  0.0000E+00  3.0224E-01  4.2384E-01
            -5.8909E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1621.34400623651        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1178
 NPARAMETR:  1.0345E+00  2.4206E-01  1.0640E-01  9.1490E-01  1.6416E-01  9.4082E-01  1.1481E+00  1.0000E-02  1.7909E+00  7.6921E-01
             2.5065E+00
 PARAMETER:  1.3392E-01 -1.3186E+00 -2.1406E+00  1.1065E-02 -1.7069E+00  3.8993E-02  2.3813E-01 -1.5448E+01  6.8274E-01 -1.6239E-01
             1.0189E+00
 GRADIENT:   3.5744E-01  2.8091E-01  8.0776E-01 -2.8497E-01 -4.0083E-02  3.2486E-03  3.5610E-02  0.0000E+00  2.7311E-01 -9.0585E-02
            -1.0763E-01

0ITERATION NO.:   42    OBJECTIVE VALUE:  -1621.34478516415        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.0344E+00  2.4193E-01  1.0635E-01  9.1531E-01  1.6418E-01  9.4082E-01  1.1427E+00  1.0000E-02  1.7893E+00  7.6959E-01
             2.5069E+00
 PARAMETER:  1.3384E-01 -1.3191E+00 -2.1410E+00  1.1504E-02 -1.7068E+00  3.9000E-02  2.3336E-01 -1.5448E+01  6.8185E-01 -1.6189E-01
             1.0190E+00
 GRADIENT:   1.4669E-01 -9.2923E-02  5.7072E-02  1.1589E-01  1.0624E+00  1.2072E-02 -1.3679E-02  0.0000E+00  3.2288E-02 -6.8978E-02
            -1.7401E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1235
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.7120E-03  1.7800E-02  2.5377E-04 -1.1195E-02  1.1423E-02
 SE:             2.9105E-02  1.2294E-02  1.9465E-04  2.5825E-02  2.6110E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5310E-01  1.4768E-01  1.9232E-01  6.6466E-01  6.6176E-01

 ETASHRINKSD(%)  2.4945E+00  5.8812E+01  9.9348E+01  1.3482E+01  1.2528E+01
 ETASHRINKVR(%)  4.9268E+00  8.3036E+01  9.9996E+01  2.5147E+01  2.3486E+01
 EBVSHRINKSD(%)  2.5197E+00  5.9922E+01  9.9376E+01  8.4385E+00  1.3818E+01
 EBVSHRINKVR(%)  4.9759E+00  8.3938E+01  9.9996E+01  1.6165E+01  2.5727E+01
 RELATIVEINF(%)  9.4810E+01  7.4276E+00  7.4444E-04  4.7497E+01  1.0554E+01
 EPSSHRINKSD(%)  2.9951E+01
 EPSSHRINKVR(%)  5.0932E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1621.3447851641481     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -702.40625195947541     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.87
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1621.345       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  2.42E-01  1.06E-01  9.15E-01  1.64E-01  9.41E-01  1.14E+00  1.00E-02  1.79E+00  7.70E-01  2.51E+00
 


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
+        1.12E+03
 
 TH 2
+       -1.84E+01  2.38E+03
 
 TH 3
+       -1.36E+02  2.28E+03  4.22E+04
 
 TH 4
+        5.08E+00 -5.98E+01 -1.02E+03  3.54E+02
 
 TH 5
+        1.64E+02 -6.02E+03 -3.48E+04 -3.67E+02  5.29E+04
 
 TH 6
+        4.74E+00 -2.92E+01  3.82E+01 -7.00E+00 -2.24E+01  2.01E+02
 
 TH 7
+        4.25E-01  5.01E+01 -4.54E+01 -1.71E+00 -2.05E+01 -2.07E-01  4.38E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.06E+01 -1.38E+01  2.64E+02 -1.02E+01  7.19E+01 -4.91E-01 -9.22E-01  0.00E+00  3.83E+01
 
 TH10
+       -3.91E+00  5.34E+01  1.53E+02  1.00E+01  1.10E+02  5.25E+00  9.19E+00  0.00E+00  4.07E+00  1.93E+02
 
 TH11
+       -2.02E+01  2.72E+01 -6.41E+01 -1.01E+00  7.52E+01  2.88E+00  7.35E+00  0.00E+00  4.47E+00  9.31E+00  6.62E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.569
Stop Time:
Thu Sep 30 00:02:45 CDT 2021
