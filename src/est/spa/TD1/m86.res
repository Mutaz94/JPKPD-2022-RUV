Wed Sep 29 18:33:45 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat86.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m86.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1697.25877960219        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9266E+02 -2.0789E+01 -5.0280E+01  6.0756E+01  8.8990E+01  5.6214E+01  9.8361E+00  8.0840E+00  2.4949E+01  2.8360E+00
            -1.7386E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1703.70191498332        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      195
 NPARAMETR:  1.0289E+00  1.0532E+00  1.0639E+00  1.0017E+00  9.8033E-01  9.8376E-01  9.4817E-01  9.6774E-01  9.2918E-01  9.6362E-01
             1.0518E+00
 PARAMETER:  1.2852E-01  1.5187E-01  1.6195E-01  1.0168E-01  8.0138E-02  8.3624E-02  4.6778E-02  6.7204E-02  2.6551E-02  6.2938E-02
             1.5050E-01
 GRADIENT:   3.4898E-02  9.3561E+00 -2.5096E-01  1.2459E+01  1.7996E-01 -2.0765E-01  3.5272E+00  1.4839E+00 -4.6980E-01 -7.7493E-01
            -3.6879E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1703.91782623429        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0319E+00  1.0981E+00  1.0422E+00  9.6912E-01  9.9299E-01  9.8745E-01  8.1534E-01  8.7228E-01  9.8534E-01  1.0212E+00
             1.0547E+00
 PARAMETER:  1.3140E-01  1.9357E-01  1.4137E-01  6.8636E-02  9.2962E-02  8.7369E-02 -1.0415E-01 -3.6643E-02  8.5227E-02  1.2098E-01
             1.5321E-01
 GRADIENT:   5.9234E+00  4.9976E+00  1.6827E+00  8.0563E+00 -4.1896E+00  1.2523E+00 -9.3854E-01 -5.5817E-01 -6.2178E-02  4.0023E+00
             1.2656E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1704.49778540587        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0297E+00  1.3095E+00  8.0132E-01  8.2577E-01  9.8135E-01  9.8402E-01  7.8383E-01  7.2024E-01  1.0763E+00  9.3139E-01
             1.0503E+00
 PARAMETER:  1.2924E-01  3.6967E-01 -1.2149E-01 -9.1444E-02  8.1170E-02  8.3888E-02 -1.4356E-01 -2.2817E-01  1.7355E-01  2.8922E-02
             1.4906E-01
 GRADIENT:  -3.4754E+00  7.5804E+00  1.4640E+00  7.8412E+00 -3.4907E+00 -1.3250E+00 -6.2261E-01  1.7609E-01 -6.3061E-01 -1.5859E+00
            -9.4843E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1704.89828340407        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0335E+00  1.5859E+00  6.0260E-01  6.4346E-01  1.0338E+00  9.9202E-01  7.1687E-01  5.1149E-01  1.2650E+00  9.5400E-01
             1.0526E+00
 PARAMETER:  1.3296E-01  5.6116E-01 -4.0650E-01 -3.4089E-01  1.3324E-01  9.1990E-02 -2.3286E-01 -5.7042E-01  3.3510E-01  5.2913E-02
             1.5124E-01
 GRADIENT:   2.8375E+00  1.1668E+01  5.8289E-01  7.4427E+00 -4.8188E+00  1.1643E+00  8.9729E-01  2.5025E-01 -4.5547E-01  1.3107E+00
             4.2155E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1705.00543031717        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      908
 NPARAMETR:  1.0331E+00  1.6564E+00  5.5975E-01  5.9079E-01  1.0593E+00  9.8935E-01  6.9553E-01  3.8837E-01  1.3428E+00  9.6084E-01
             1.0525E+00
 PARAMETER:  1.3256E-01  6.0466E-01 -4.8027E-01 -4.2630E-01  1.5766E-01  8.9294E-02 -2.6308E-01 -8.4580E-01  3.9478E-01  6.0056E-02
             1.5114E-01
 GRADIENT:   1.9573E+00 -2.6556E-02  1.1793E+00  1.3087E+00 -6.0428E-01  5.0136E-02  1.2513E-01  3.7853E-02 -2.1800E-01 -7.7060E-02
            -5.3721E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1705.02189119509        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1085
 NPARAMETR:  1.0290E+00  1.6880E+00  5.2290E-01  5.6805E-01  1.0602E+00  9.8874E-01  6.8845E-01  2.1068E-01  1.3776E+00  9.5877E-01
             1.0524E+00
 PARAMETER:  1.2854E-01  6.2354E-01 -5.4836E-01 -4.6555E-01  1.5843E-01  8.8673E-02 -2.7332E-01 -1.4574E+00  4.2034E-01  5.7893E-02
             1.5105E-01
 GRADIENT:  -7.4992E+00 -3.7799E+00 -1.8165E-01  1.8782E+00  9.5041E-01 -3.2025E-01 -8.4289E-02  2.1344E-02  6.6223E-01  6.5199E-01
             2.4403E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1705.02253292919        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  1.0286E+00  1.7007E+00  5.1075E-01  5.5934E-01  1.0610E+00  9.8877E-01  6.8621E-01  1.6597E-01  1.3892E+00  9.5684E-01
             1.0524E+00
 PARAMETER:  1.2820E-01  6.3103E-01 -5.7188E-01 -4.8100E-01  1.5923E-01  8.8704E-02 -2.7657E-01 -1.6959E+00  4.2874E-01  5.5881E-02
             1.5109E-01
 GRADIENT:  -8.3444E+00 -4.0963E+00 -5.3809E-01  1.9724E+00  7.6705E-01 -3.4207E-01 -3.8551E-02  2.2783E-02  7.6875E-01  7.1690E-01
             3.1400E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1705.06086568151        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  1.0340E+00  1.7006E+00  5.1052E-01  5.5809E-01  1.0602E+00  9.8960E-01  6.8730E-01  8.4159E-02  1.3831E+00  9.5221E-01
             1.0521E+00
 PARAMETER:  1.3343E-01  6.3096E-01 -5.7232E-01 -4.8324E-01  1.5847E-01  8.9546E-02 -2.7498E-01 -2.3750E+00  4.2434E-01  5.1025E-02
             1.5077E-01
 GRADIENT:   3.7165E+00 -4.2428E+00  1.0570E+00 -4.7249E-01 -4.7789E-01  4.4054E-02 -1.6514E-01  2.4014E-03 -3.1447E-01 -1.8724E-01
            -1.8213E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1705.06691349883        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1606             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0337E+00  1.7001E+00  5.0780E-01  5.5744E-01  1.0590E+00  9.8987E-01  6.8787E-01  1.6751E-02  1.3845E+00  9.5161E-01
             1.0521E+00
 PARAMETER:  1.3311E-01  6.3070E-01 -5.7766E-01 -4.8440E-01  1.5736E-01  8.9814E-02 -2.7416E-01 -3.9893E+00  4.2535E-01  5.0399E-02
             1.5074E-01
 GRADIENT:   5.6076E+02  6.1863E+02  4.6264E+00  1.0782E+02  1.0650E+01  4.5031E+01  1.2097E+01  8.3867E-04  1.6467E+01  7.1206E-01
             1.1648E+00

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1705.06768371824        NO. OF FUNC. EVALS.: 132
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  1.0337E+00  1.6997E+00  5.0650E-01  5.5783E-01  1.0590E+00  9.8986E-01  6.8801E-01  1.0000E-02  1.3850E+00  9.5121E-01
             1.0522E+00
 PARAMETER:  1.3310E-01  6.3047E-01 -5.8024E-01 -4.8369E-01  1.5730E-01  8.9807E-02 -2.7395E-01 -5.3379E+00  4.2568E-01  4.9981E-02
             1.5084E-01
 GRADIENT:   2.9382E+00 -7.4035E+00 -1.8810E-01 -6.8279E-02  1.6012E+00  1.3524E-01  6.1767E-02  0.0000E+00  2.6431E-01  1.4164E-01
             7.8417E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1738
 NO. OF SIG. DIGITS IN FINAL EST.:  2.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0732E-04 -3.2697E-02 -2.6554E-04  2.4753E-02 -3.5012E-02
 SE:             2.9828E-02  2.2368E-02  1.0836E-04  2.3441E-02  2.2774E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9445E-01  1.4380E-01  1.4264E-02  2.9097E-01  1.2420E-01

 ETASHRINKSD(%)  7.3685E-02  2.5065E+01  9.9637E+01  2.1471E+01  2.3704E+01
 ETASHRINKVR(%)  1.4732E-01  4.3847E+01  9.9999E+01  3.8332E+01  4.1789E+01
 EBVSHRINKSD(%)  4.8813E-01  2.4414E+01  9.9693E+01  2.2981E+01  2.2185E+01
 EBVSHRINKVR(%)  9.7387E-01  4.2868E+01  9.9999E+01  4.0681E+01  3.9448E+01
 RELATIVEINF(%)  9.8985E+01  3.5530E+00  1.1122E-04  3.9620E+00  1.1428E+01
 EPSSHRINKSD(%)  4.4034E+01
 EPSSHRINKVR(%)  6.8678E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1705.0676837182350     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.91685715449682     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.03
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.12
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1705.068       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.70E+00  5.06E-01  5.58E-01  1.06E+00  9.90E-01  6.88E-01  1.00E-02  1.38E+00  9.51E-01  1.05E+00
 


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
 ********************                                     T MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.21E+03
 
 TH 2
+        1.37E+02  5.23E+02
 
 TH 3
+        2.27E+02  2.51E+02  4.91E+02
 
 TH 4
+       -1.46E+02  4.07E+02 -3.46E+02  1.12E+03
 
 TH 5
+       -3.07E+02 -3.52E+02 -4.90E+02  2.55E+02  7.42E+02
 
 TH 6
+        5.04E+01 -1.90E+01  3.31E+00 -2.96E+01  3.46E+01  1.72E+02
 
 TH 7
+        2.17E+01 -5.33E+00 -1.04E+01 -8.32E+00 -4.54E+01 -1.44E+01  2.41E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -7.04E+00 -3.90E+00 -4.99E+01  5.68E+01  5.25E+00 -1.23E+01  1.17E+01  0.00E+00  1.32E+01
 
 TH10
+        1.75E+01 -9.01E+00 -4.32E+01  2.18E+01 -6.45E+01 -1.22E+01  2.75E+01  0.00E+00  2.30E+01  5.24E+01
 
 TH11
+        6.87E+01  7.48E+00  1.92E+01 -2.84E+01 -4.08E+01  1.29E+01  5.08E+01  0.00E+00  5.12E+00  1.41E+01  2.22E+02
 
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
+        1.05E+03
 
 TH 2
+       -7.55E+00  4.89E+02
 
 TH 3
+        7.96E+00  1.96E+02  4.51E+02
 
 TH 4
+       -1.85E+01  4.07E+02 -3.18E+02  1.12E+03
 
 TH 5
+       -3.40E+00 -2.57E+02 -4.08E+02  2.78E+02  6.63E+02
 
 TH 6
+        2.99E-01 -1.64E+00  1.91E+00 -4.21E+00 -1.13E+00  1.99E+02
 
 TH 7
+        1.28E+00  2.23E+00 -8.51E+00 -1.30E+01 -1.77E+01 -3.99E-01  1.56E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.01E+00 -2.11E+01 -3.35E+01  5.78E+01  4.71E-01 -3.24E-01  2.32E+01  0.00E+00  4.45E+01
 
 TH10
+       -6.18E-01 -1.91E+01 -4.34E+01 -1.36E+00 -6.35E+01  8.35E-01  1.38E+01  0.00E+00  6.63E+00  8.98E+01
 
 TH11
+       -6.98E+00 -1.98E+01 -2.54E+01 -6.60E-01 -1.35E+00  2.71E+00  1.18E+01  0.00E+00  5.31E+00  1.83E+01  1.92E+02
 
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
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                                     S MATRIX                                   ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        1.05E+03
 
 TH 2
+       -1.23E+02  4.75E+02
 
 TH 3
+       -1.81E+02  1.70E+02  4.39E+02
 
 TH 4
+        9.96E+01  4.06E+02 -3.53E+02  1.12E+03
 
 TH 5
+        2.49E+02 -2.14E+02 -3.93E+02  3.24E+02  6.78E+02
 
 TH 6
+       -6.76E+01  2.54E+01 -9.88E-01 -7.48E+00 -5.96E+01  2.37E+02
 
 TH 7
+       -3.34E+01 -4.84E+01  1.23E+00 -7.71E+01 -1.88E+01  1.70E+00  1.18E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.71E+01 -4.38E+01 -5.15E+01  3.05E+01 -9.45E-02 -3.20E+00  2.69E+01  0.00E+00  4.15E+01
 
 TH10
+       -6.10E+01 -1.99E+01 -2.85E+01 -2.72E+01 -8.95E+01  3.64E+01  2.67E+01  0.00E+00  8.16E+00  1.03E+02
 
 TH11
+       -4.61E+01 -2.18E+01 -4.66E+01  3.14E+01  2.14E+01  1.69E-01  2.07E+01  0.00E+00  1.09E+01  1.94E+01  1.61E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.02
 #CPUT: Total CPU Time in Seconds,       29.210
Stop Time:
Wed Sep 29 18:34:15 CDT 2021
