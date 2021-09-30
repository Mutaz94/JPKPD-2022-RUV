Wed Sep 29 12:02:05 CDT 2021
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
$DATA ../../../../data/spa/A1/dat28.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m28.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1436.66352029042        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2540E+02 -7.0457E+01 -1.9312E+01 -9.0148E+01  9.7901E+01  3.8620E+01 -1.2735E+00  6.4679E+00 -2.9405E+01 -7.3945E+00
            -3.7419E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1495.85442908762        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9616E-01  1.0073E+00  1.0394E+00  1.1441E+00  9.1452E-01  9.0649E-01  8.6902E-01  8.7973E-01  1.0343E+00  7.7084E-01
             2.5697E+00
 PARAMETER:  9.6153E-02  1.0728E-01  1.3864E-01  2.3459E-01  1.0642E-02  1.8230E-03 -4.0394E-02 -2.8141E-02  1.3376E-01 -1.6028E-01
             1.0438E+00
 GRADIENT:   6.2406E+01  3.4715E+01  6.0658E+00  5.5672E+01 -3.2162E+01 -1.3906E+01  7.7017E+00  6.6307E+00  1.1793E+01  1.6627E+01
             1.6050E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1510.12544692890        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      158
 NPARAMETR:  9.8448E-01  6.8771E-01  5.4123E-01  1.3127E+00  5.4399E-01  9.3321E-01  8.2111E-01  3.1585E-01  9.3801E-01  4.7455E-01
             2.2890E+00
 PARAMETER:  8.4353E-02 -2.7439E-01 -5.1390E-01  3.7205E-01 -5.0882E-01  3.0872E-02 -9.7094E-02 -1.0525E+00  3.6003E-02 -6.4539E-01
             9.2812E-01
 GRADIENT:   3.7606E+01  3.0937E+01 -4.3427E+01  1.7202E+02  4.6290E+01 -6.7397E+00 -4.6554E+00  1.5718E+00  4.8844E+00  3.0681E+00
             1.1680E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1535.79667215284        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:      273
 NPARAMETR:  9.6735E-01  7.2292E-01  5.5617E-01  1.2135E+00  5.6222E-01  9.6104E-01  1.0846E+00  1.1353E-01  9.3011E-01  5.5930E-01
             1.7378E+00
 PARAMETER:  6.6810E-02 -2.2445E-01 -4.8669E-01  2.9353E-01 -4.7585E-01  6.0259E-02  1.8121E-01 -2.0757E+00  2.7552E-02 -4.8107E-01
             6.5265E-01
 GRADIENT:  -4.9237E+01  2.2559E+00 -9.3199E+00  1.2829E+00  7.3243E+00 -5.0752E+00 -2.3713E+00  9.4129E-02  4.6026E+00 -2.3017E+00
             2.7206E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1539.99055700612        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      450
 NPARAMETR:  9.8687E-01  5.3431E-01  6.3112E-01  1.3414E+00  5.4533E-01  9.7003E-01  1.2586E+00  6.4128E-02  8.7340E-01  7.0729E-01
             1.6142E+00
 PARAMETER:  8.6779E-02 -5.2678E-01 -3.6026E-01  3.9370E-01 -5.0636E-01  6.9572E-02  3.3000E-01 -2.6469E+00 -3.5360E-02 -2.4632E-01
             5.7881E-01
 GRADIENT:   3.9138E+00  1.6455E+01  1.1957E+01  1.0443E+01 -2.3953E+01 -1.2830E+00 -4.7075E-01  4.5867E-02 -1.9249E+00  5.6262E+00
             1.3670E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1544.66906250216        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      627
 NPARAMETR:  9.8010E-01  2.1004E-01  6.4605E-01  1.5095E+00  4.8965E-01  9.7441E-01  2.2751E+00  1.0000E-02  8.1611E-01  7.3491E-01
             1.6081E+00
 PARAMETER:  7.9897E-02 -1.4604E+00 -3.3687E-01  5.1181E-01 -6.1407E-01  7.4079E-02  9.2200E-01 -4.7751E+00 -1.0321E-01 -2.0800E-01
             5.7507E-01
 GRADIENT:   5.0223E+00  5.4492E+00  1.2904E+01  1.1950E+01 -1.5748E+01  1.7012E+00  1.1932E+00  0.0000E+00  3.2000E+00  4.9773E+00
             3.1015E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1546.79841202265        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      802
 NPARAMETR:  9.7497E-01  6.7253E-02  5.3314E-01  1.5310E+00  4.0587E-01  9.7135E-01  3.9698E+00  1.0000E-02  7.9295E-01  6.7875E-01
             1.5883E+00
 PARAMETER:  7.4654E-02 -2.5993E+00 -5.2897E-01  5.2591E-01 -8.0173E-01  7.0931E-02  1.4787E+00 -8.2099E+00 -1.3199E-01 -2.8750E-01
             5.6264E-01
 GRADIENT:   1.5663E+00  1.2600E+00  3.9392E+00  1.0454E+01 -1.2814E+01  7.4368E-01 -7.0404E-02  0.0000E+00  4.6452E-01 -5.1870E-02
            -2.1813E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1547.39501642175        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      977
 NPARAMETR:  9.7102E-01  1.2773E-02  5.8357E-01  1.5766E+00  4.2901E-01  9.6778E-01  7.6789E+00  1.0000E-02  7.7694E-01  7.0290E-01
             1.5853E+00
 PARAMETER:  7.0589E-02 -4.2605E+00 -4.3859E-01  5.5528E-01 -7.4628E-01  6.7254E-02  2.1385E+00 -1.2401E+01 -1.5239E-01 -2.5254E-01
             5.6075E-01
 GRADIENT:  -1.0982E+00  1.6974E-01  4.0888E-01  8.3906E+00 -8.2829E-01 -1.6268E-02 -4.1256E-02  0.0000E+00 -8.0863E-01 -3.2010E-01
            -2.5086E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1547.44437357467        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:     1124
 NPARAMETR:  9.6862E-01  1.0000E-02  5.8258E-01  1.5707E+00  4.2803E-01  9.6695E-01  9.0556E+00  1.0000E-02  7.7672E-01  7.0171E-01
             1.5863E+00
 PARAMETER:  6.8115E-02 -4.5176E+00 -4.4029E-01  5.5153E-01 -7.4857E-01  6.6396E-02  2.3034E+00 -1.3093E+01 -1.5268E-01 -2.5424E-01
             5.6142E-01
 GRADIENT:  -6.3900E+00  1.2899E-01  2.4375E+00 -5.6678E+00 -6.7497E-01 -2.7387E-01 -2.8642E-02  0.0000E+00 -5.8901E-01 -4.7608E-01
            -2.0171E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1547.46465382912        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1312             RESET HESSIAN, TYPE I
 NPARAMETR:  9.7135E-01  1.0000E-02  5.8160E-01  1.5705E+00  4.2755E-01  9.6762E-01  1.0603E+01  1.0000E-02  7.7788E-01  7.0189E-01
             1.5930E+00
 PARAMETER:  7.0932E-02 -4.5852E+00 -4.4198E-01  5.5140E-01 -7.4968E-01  6.7088E-02  2.4611E+00 -1.3093E+01 -1.5118E-01 -2.5398E-01
             5.6563E-01
 GRADIENT:   1.2918E+02  0.0000E+00  1.5554E+01  2.8476E+02  5.8141E+01  1.0940E+01  2.8382E-01  0.0000E+00  6.8235E+00  1.1433E+00
             4.6035E+00

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1547.46537154524        NO. OF FUNC. EVALS.:  93
 CUMULATIVE NO. OF FUNC. EVALS.:     1405
 NPARAMETR:  9.7135E-01  1.0000E-02  5.8152E-01  1.5709E+00  4.2755E-01  9.6761E-01  1.0345E+01  1.0000E-02  7.7785E-01  7.0187E-01
             1.5930E+00
 PARAMETER:  7.0930E-02 -4.5852E+00 -4.4211E-01  5.5163E-01 -7.4968E-01  6.7074E-02  2.4365E+00 -1.3093E+01 -1.5122E-01 -2.5401E-01
             5.6561E-01
 GRADIENT:   6.6800E-03  0.0000E+00  1.3918E+00 -4.7832E+00  1.9121E-01  2.9126E-02  1.1855E-02  0.0000E+00  3.9229E-02  3.9091E-02
             6.5296E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1405
 NO. OF SIG. DIGITS IN FINAL EST.:  2.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         3.7690E-04  6.1369E-04 -6.6010E-05 -6.5080E-03 -9.6984E-03
 SE:             2.9606E-02  1.7856E-03  2.6889E-04  2.8375E-02  2.4412E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8984E-01  7.3108E-01  8.0608E-01  8.1859E-01  6.9116E-01

 ETASHRINKSD(%)  8.1535E-01  9.4018E+01  9.9099E+01  4.9397E+00  1.8217E+01
 ETASHRINKVR(%)  1.6240E+00  9.9642E+01  9.9992E+01  9.6354E+00  3.3115E+01
 EBVSHRINKSD(%)  1.0379E+00  9.4468E+01  9.9095E+01  4.8869E+00  1.7531E+01
 EBVSHRINKVR(%)  2.0650E+00  9.9694E+01  9.9992E+01  9.5349E+00  3.1989E+01
 RELATIVEINF(%)  8.8347E+01  1.7707E-02  3.5006E-04  1.0088E+01  2.0372E+00
 EPSSHRINKSD(%)  3.8826E+01
 EPSSHRINKVR(%)  6.2577E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1547.4653715452439     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -812.31454498150572     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.38
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1547.465       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.71E-01  1.00E-02  5.82E-01  1.57E+00  4.28E-01  9.68E-01  1.03E+01  1.00E-02  7.78E-01  7.02E-01  1.59E+00
 


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
+        1.23E+03
 
 TH 2
+        0.00E+00  7.64E+02
 
 TH 3
+       -6.26E+00  0.00E+00  2.70E+03
 
 TH 4
+       -1.87E+01  0.00E+00 -2.13E+02  6.74E+02
 
 TH 5
+        4.53E+01  0.00E+00 -4.55E+03 -1.31E+02  8.40E+03
 
 TH 6
+        7.40E-01  0.00E+00  5.15E+00 -5.21E+00 -1.94E-01  2.04E+02
 
 TH 7
+        1.89E-03  0.00E+00 -6.89E-03 -1.63E-02  3.42E-02  7.78E-04  4.17E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        5.80E+00  0.00E+00  2.98E+01 -8.80E+00 -1.77E+00  7.04E-01  3.50E-02  0.00E+00  2.72E+02
 
 TH10
+       -3.06E+00  0.00E+00 -3.27E+01  8.06E-01 -7.59E+01 -9.37E-01  2.88E-03  0.00E+00  1.48E+00  1.82E+02
 
 TH11
+       -1.05E+01  0.00E+00 -2.82E+01 -7.43E+00  1.15E+01  2.67E+00 -1.74E-04  0.00E+00  1.01E+01  3.73E+01  9.44E+01
 
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
 #CPUT: Total CPU Time in Seconds,       22.063
Stop Time:
Wed Sep 29 12:02:28 CDT 2021
