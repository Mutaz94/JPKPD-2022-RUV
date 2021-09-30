Wed Sep 29 11:22:47 CDT 2021
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
$DATA ../../../../data/spa/B/dat60.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1663.92620476356        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4154E+02 -4.3286E+01 -3.5142E+01  1.7917E+01  4.1523E+01  8.2250E+01 -1.4588E+01  6.7787E+00  1.2585E+01  1.6317E+01
            -2.1663E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.84717992145        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      175
 NPARAMETR:  1.0043E+00  1.1778E+00  1.1540E+00  9.4155E-01  1.1003E+00  1.0860E+00  1.1550E+00  9.5689E-01  9.4673E-01  8.5284E-01
             1.1805E+00
 PARAMETER:  1.0425E-01  2.6361E-01  2.4323E-01  3.9773E-02  1.9561E-01  1.8251E-01  2.4409E-01  5.5937E-02  4.5263E-02 -5.9185E-02
             2.6595E-01
 GRADIENT:   7.8315E-01 -1.1898E-01  1.8142E+01 -1.1600E+01  4.1869E+00  2.1711E+01 -3.2345E+00 -6.6436E+00 -1.1057E+00 -6.0886E+00
             3.1851E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1670.60975864455        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      357
 NPARAMETR:  1.0060E+00  1.1822E+00  1.0613E+00  9.3219E-01  1.0438E+00  1.0424E+00  1.2994E+00  1.0245E+00  8.9663E-01  7.2035E-01
             1.1596E+00
 PARAMETER:  1.0603E-01  2.6736E-01  1.5948E-01  2.9779E-02  1.4289E-01  1.4148E-01  3.6194E-01  1.2424E-01 -9.1097E-03 -2.2802E-01
             2.4807E-01
 GRADIENT:   4.0531E+00  1.2657E+01  2.3894E+01 -2.5890E+01 -1.1711E+01  5.7276E+00  1.2367E+01 -5.7595E+00 -4.5709E+00 -1.0447E+01
             2.5254E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1674.21944697133        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      536
 NPARAMETR:  1.0057E+00  1.0868E+00  9.3696E-01  9.9213E-01  9.6660E-01  1.0324E+00  1.3024E+00  7.8197E-01  9.1138E-01  7.9497E-01
             1.0668E+00
 PARAMETER:  1.0570E-01  1.8321E-01  3.4883E-02  9.2102E-02  6.6032E-02  1.3189E-01  3.6420E-01 -1.4594E-01  7.2040E-03 -1.2945E-01
             1.6469E-01
 GRADIENT:   3.5336E+00  3.6084E+00  2.2321E+00  3.8334E+00 -2.6869E-01  1.4218E+00  3.1061E+00 -6.2809E-01  3.6089E+00  5.3469E-01
             1.5676E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1674.57494346788        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      715
 NPARAMETR:  1.0059E+00  1.2891E+00  7.1111E-01  8.5591E-01  9.4253E-01  1.0334E+00  1.1219E+00  5.5417E-01  9.6467E-01  7.3743E-01
             1.0638E+00
 PARAMETER:  1.0590E-01  3.5396E-01 -2.4093E-01 -5.5588E-02  4.0810E-02  1.3283E-01  2.1506E-01 -4.9029E-01  6.4035E-02 -2.0458E-01
             1.6189E-01
 GRADIENT:  -5.4744E-01  8.3567E+00  2.1285E+00  5.1263E+00 -7.7750E+00  7.6010E-01 -5.4232E-01  3.4727E-01 -1.9100E+00  1.4050E-02
             7.2563E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1674.83773453548        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      890
 NPARAMETR:  1.0065E+00  1.4887E+00  5.7937E-01  7.2513E-01  9.7683E-01  1.0350E+00  1.0006E+00  3.1943E-01  1.0716E+00  7.4218E-01
             1.0622E+00
 PARAMETER:  1.0646E-01  4.9794E-01 -4.4581E-01 -2.2141E-01  7.6553E-02  1.3437E-01  1.0057E-01 -1.0412E+00  1.6918E-01 -1.9817E-01
             1.6034E-01
 GRADIENT:  -9.2449E-01  1.0091E+01  1.4427E+00  3.7729E+00 -7.8469E+00  9.3546E-01 -2.3566E-01  2.2442E-01 -2.7057E+00  1.9889E-01
             3.1945E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1674.94707596873        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1066
 NPARAMETR:  1.0070E+00  1.6250E+00  5.3766E-01  6.4233E-01  1.0382E+00  1.0322E+00  9.2755E-01  2.0442E-01  1.2039E+00  7.8170E-01
             1.0607E+00
 PARAMETER:  1.0701E-01  5.8549E-01 -5.2052E-01 -3.4265E-01  1.3750E-01  1.3172E-01  2.4786E-02 -1.4876E+00  2.8555E-01 -1.4628E-01
             1.5893E-01
 GRADIENT:   1.3207E-01  1.1029E+01  8.2232E-01  7.6791E+00 -3.5550E+00 -1.9208E-01 -1.2302E+00  8.1008E-02  5.7439E-02 -4.6061E-01
            -6.5887E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1674.94925651500        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1245
 NPARAMETR:  1.0071E+00  1.6556E+00  5.2596E-01  6.2252E-01  1.0516E+00  1.0321E+00  9.1382E-01  1.7928E-01  1.2338E+00  7.9017E-01
             1.0609E+00
 PARAMETER:  1.0704E-01  6.0416E-01 -5.4253E-01 -3.7398E-01  1.5034E-01  1.3156E-01  9.8785E-03 -1.6188E+00  3.1012E-01 -1.3551E-01
             1.5912E-01
 GRADIENT:   2.0366E-01  1.0073E+01  7.3261E-01  7.2696E+00 -2.9326E+00 -2.5552E-01 -1.1744E+00  6.2978E-02  2.7486E-01 -4.6647E-01
            -6.5906E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1675.05227845201        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1428
 NPARAMETR:  1.0072E+00  1.6496E+00  5.2457E-01  6.1747E-01  1.0538E+00  1.0336E+00  9.1832E-01  8.8009E-02  1.2339E+00  7.9342E-01
             1.0617E+00
 PARAMETER:  1.0718E-01  6.0052E-01 -5.4518E-01 -3.8212E-01  1.5238E-01  1.3309E-01  1.4796E-02 -2.3303E+00  3.1018E-01 -1.3140E-01
             1.5991E-01
 GRADIENT:   6.4539E-01 -3.2201E+00  2.9701E-01  7.4754E-03 -2.7509E-01  4.1157E-01  1.3003E-02  1.6609E-02 -1.3529E-01 -3.2853E-02
             3.6542E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1675.06338692031        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1616
 NPARAMETR:  1.0079E+00  1.6477E+00  5.2397E-01  6.1709E-01  1.0536E+00  1.0355E+00  9.1836E-01  1.1702E-02  1.2355E+00  7.9358E-01
             1.0615E+00
 PARAMETER:  1.0783E-01  5.9938E-01 -5.4633E-01 -3.8274E-01  1.5217E-01  1.3485E-01  1.4833E-02 -4.3480E+00  3.1149E-01 -1.3120E-01
             1.5967E-01
 GRADIENT:   2.0411E+00 -5.8674E+00  1.7732E-02 -8.9205E-01  4.2302E-01  1.1349E+00 -1.1148E-02  3.1971E-04  7.4662E-02  5.1595E-02
             3.7818E-02

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1675.06376718327        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1717
 NPARAMETR:  1.0079E+00  1.6424E+00  5.2298E-01  6.1713E-01  1.0518E+00  1.0356E+00  9.1947E-01  1.0000E-02  1.2319E+00  7.9267E-01
             1.0623E+00
 PARAMETER:  1.0781E-01  5.9937E-01 -5.4654E-01 -3.8193E-01  1.5165E-01  1.3483E-01  1.5045E-02 -6.0910E+00  3.1117E-01 -1.3214E-01
             1.5949E-01
 GRADIENT:  -2.3017E-02  8.5102E-01  6.6887E-02  8.2858E-02  2.6004E-01 -9.1380E-03 -2.9787E-02  0.0000E+00  5.6344E-02  2.9646E-03
            -5.7209E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1717
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2143E-04 -1.8859E-02 -2.9834E-04  1.8065E-02 -2.9688E-02
 SE:             2.9848E-02  2.5348E-02  1.1598E-04  2.2375E-02  2.0469E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9675E-01  4.5689E-01  1.0102E-02  4.1947E-01  1.4695E-01

 ETASHRINKSD(%)  5.7294E-03  1.5080E+01  9.9611E+01  2.5040E+01  3.1427E+01
 ETASHRINKVR(%)  1.1458E-02  2.7885E+01  9.9998E+01  4.3810E+01  5.2977E+01
 EBVSHRINKSD(%)  4.4416E-01  1.5193E+01  9.9654E+01  2.6109E+01  3.0853E+01
 EBVSHRINKVR(%)  8.8635E-01  2.8077E+01  9.9999E+01  4.5401E+01  5.2187E+01
 RELATIVEINF(%)  9.9067E+01  5.4134E+00  1.0727E-04  3.5455E+00  7.6860E+00
 EPSSHRINKSD(%)  4.3548E+01
 EPSSHRINKVR(%)  6.8132E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1675.0637671832737     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.91294061953556     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.80
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1675.064       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.01E+00  1.65E+00  5.24E-01  6.18E-01  1.05E+00  1.04E+00  9.19E-01  1.00E-02  1.24E+00  7.93E-01  1.06E+00
 


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
+        1.01E+03
 
 TH 2
+       -4.94E+00  3.60E+02
 
 TH 3
+        9.28E+00  1.72E+02  5.35E+02
 
 TH 4
+       -1.53E+01  2.68E+02 -4.44E+02  1.07E+03
 
 TH 5
+       -5.11E+00 -2.21E+02 -5.05E+02  4.22E+02  7.71E+02
 
 TH 6
+       -5.16E-02 -9.10E-01  2.40E+00 -3.60E+00 -9.36E-01  1.83E+02
 
 TH 7
+        4.08E-01  1.75E+01 -3.00E+01 -1.71E+01  7.35E+00 -9.61E-02  1.29E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        9.52E-01 -1.68E+01 -4.07E+01  6.37E+01 -7.28E+00 -3.66E-01  1.13E+01  0.00E+00  5.13E+01
 
 TH10
+       -3.68E-01 -1.35E+01 -3.49E+01 -1.39E+01 -7.69E+01  3.07E-01  1.75E+01  0.00E+00  1.07E+01  8.54E+01
 
 TH11
+       -6.20E+00 -1.50E+01 -2.35E+01 -4.64E-01 -1.18E+01  2.22E+00 -6.47E+06  0.00E+00  8.27E+00  2.14E+01  1.90E+02
 
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
 #CPUT: Total CPU Time in Seconds,       28.872
Stop Time:
Wed Sep 29 11:23:17 CDT 2021
