Thu Sep 30 00:46:34 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat100.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m100.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   392.625252349541        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.2122E+02  4.9339E+01  2.7623E+02 -7.2038E+00  2.4631E+02  2.3534E+01 -6.1017E+01 -3.3840E+02 -7.9137E+00 -1.5916E+02
            -4.2735E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1429.25780069412        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.0148E+00  1.0638E+00  9.1700E-01  1.1139E+00  9.9299E-01  8.5106E-01  9.1853E-01  1.0021E+00  8.0927E-01  9.1766E-01
             4.8488E+00
 PARAMETER:  1.1469E-01  1.6185E-01  1.3355E-02  2.0782E-01  9.2964E-02 -6.1278E-02  1.5021E-02  1.0213E-01 -1.1162E-01  1.4073E-02
             1.6787E+00
 GRADIENT:   1.0838E+02 -4.2227E+00 -2.1509E+01  2.3089E+01  4.9990E+00 -2.0565E+01  7.5323E+00  6.6174E+00  1.2626E+01  1.6512E+01
             2.6167E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1458.13921985766        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.7017E-01  8.7757E-01  3.0553E-01  1.1150E+00  4.3639E-01  9.2120E-01  9.1485E-01  2.3474E-01  8.6640E-01  2.6409E-01
             4.2288E+00
 PARAMETER:  6.9721E-02 -3.0603E-02 -1.0857E+00  2.0884E-01 -7.2921E-01  1.7918E-02  1.1008E-02 -1.3493E+00 -4.3408E-02 -1.2314E+00
             1.5419E+00
 GRADIENT:  -8.1847E+00  8.9469E+01  3.5831E+01  6.4904E+01 -9.8384E+01 -8.5720E+00 -3.5533E+00 -2.0805E-01 -2.0988E+00 -3.1393E-01
             1.8073E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1481.38575905379        NO. OF FUNC. EVALS.:  97
 CUMULATIVE NO. OF FUNC. EVALS.:      254
 NPARAMETR:  9.2902E-01  7.3792E-01  2.5738E-01  1.1142E+00  3.8400E-01  9.4014E-01  1.0538E+00  8.2555E-02  9.8768E-01  2.7227E-01
             3.3683E+00
 PARAMETER:  2.6378E-02 -2.0392E-01 -1.2572E+00  2.0814E-01 -8.5712E-01  3.8276E-02  1.5241E-01 -2.3943E+00  8.7603E-02 -1.2010E+00
             1.3144E+00
 GRADIENT:  -7.8403E+01  1.7090E+01 -1.7949E+01  4.2399E+01 -6.2700E+00 -1.1527E+01 -4.1574E+00 -1.7930E-01 -3.9364E+00 -3.8782E+00
             1.1753E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1486.85247784380        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      429
 NPARAMETR:  9.6263E-01  6.8060E-01  2.8428E-01  1.1397E+00  3.8216E-01  9.6530E-01  1.0358E+00  3.0409E-02  1.0068E+00  5.7587E-01
             3.1854E+00
 PARAMETER:  6.1915E-02 -2.8478E-01 -1.1578E+00  2.3077E-01 -8.6192E-01  6.4679E-02  1.3521E-01 -3.3930E+00  1.0675E-01 -4.5188E-01
             1.2586E+00
 GRADIENT:   1.0393E+01  1.4567E+01 -9.0231E-02  2.3047E+01  3.9561E+00  1.1130E+00  2.5546E+00 -8.9242E-03  1.2734E+00  2.8852E-01
            -1.9894E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1492.73980678807        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      605
 NPARAMETR:  9.5809E-01  4.6296E-01  2.1008E-01  1.1281E+00  2.7206E-01  9.7304E-01  6.1721E-01  1.0000E-02  1.0746E+00  6.8189E-01
             3.1427E+00
 PARAMETER:  5.7182E-02 -6.7012E-01 -1.4603E+00  2.2055E-01 -1.2017E+00  7.2667E-02 -3.8254E-01 -5.4362E+00  1.7192E-01 -2.8288E-01
             1.2451E+00
 GRADIENT:  -1.2066E+00  8.5553E+00  6.6614E+00  1.5332E+00 -7.7023E+00  1.4327E+00  7.8145E-01  0.0000E+00 -4.3452E+00 -6.6825E+00
             4.3647E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1495.46300272334        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      780
 NPARAMETR:  9.5724E-01  3.5766E-01  1.5555E-01  1.0563E+00  2.1461E-01  9.7047E-01  3.2769E-01  1.0000E-02  1.2108E+00  7.7521E-01
             3.0272E+00
 PARAMETER:  5.6301E-02 -9.2818E-01 -1.7608E+00  1.5481E-01 -1.4389E+00  7.0024E-02 -1.0157E+00 -7.5774E+00  2.9124E-01 -1.5462E-01
             1.2076E+00
 GRADIENT:   1.1225E+00 -3.5034E+00 -2.4329E+00 -3.7734E+00  5.2346E+00  5.5084E-02  1.0157E+00  0.0000E+00  7.7452E-01  1.4676E-01
             2.6893E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1495.98149105027        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      957
 NPARAMETR:  9.5668E-01  3.6807E-01  1.5880E-01  1.0662E+00  2.1863E-01  9.6988E-01  6.3129E-02  1.0000E-02  1.2025E+00  7.7981E-01
             3.0183E+00
 PARAMETER:  5.5716E-02 -8.9948E-01 -1.7401E+00  1.6405E-01 -1.4204E+00  6.9413E-02 -2.6626E+00 -5.0419E+00  2.8439E-01 -1.4870E-01
             1.2047E+00
 GRADIENT:  -1.6091E-02 -1.0608E+00 -1.1048E+00  1.5916E+00  2.2555E+00 -4.9181E-01  3.6541E-02  0.0000E+00  1.3179E-01  3.9063E-01
            -1.6298E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1498.32284578892        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1135
 NPARAMETR:  9.5607E-01  3.6410E-01  1.5357E-01  1.0563E+00  2.1384E-01  9.7207E-01  1.0000E-02  1.0216E+00  1.2164E+00  7.4921E-01
             2.9840E+00
 PARAMETER:  5.5074E-02 -9.1032E-01 -1.7736E+00  1.5480E-01 -1.4425E+00  7.1668E-02 -6.5832E+00  1.2134E-01  2.9586E-01 -1.8873E-01
             1.1932E+00
 GRADIENT:  -2.3814E+00  4.2519E-01  7.0680E+00 -5.5888E+00 -2.0298E+00  5.5428E-02  0.0000E+00  4.8231E+00  2.8893E-02  1.8803E+01
             2.3331E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1501.71422852103        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1312
 NPARAMETR:  9.5296E-01  3.5177E-01  1.4595E-01  1.0646E+00  2.0874E-01  9.7427E-01  1.0000E-02  1.2777E+00  1.3147E+00  5.5931E-01
             2.8421E+00
 PARAMETER:  5.1818E-02 -9.4477E-01 -1.8245E+00  1.6263E-01 -1.4666E+00  7.3930E-02 -6.6342E+00  3.4503E-01  3.7361E-01 -4.8105E-01
             1.1445E+00
 GRADIENT:   1.9416E-01 -1.7911E-01 -4.6616E-01 -1.1022E-02  6.9024E-01 -1.2474E-01  0.0000E+00  5.2228E-01  4.8424E-01  7.0806E-01
            -2.9422E-01

0ITERATION NO.:   48    OBJECTIVE VALUE:  -1501.71823478028        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     1410
 NPARAMETR:  9.5302E-01  3.5385E-01  1.4754E-01  1.0673E+00  2.1008E-01  9.7478E-01  1.0000E-02  1.2799E+00  1.3086E+00  5.4691E-01
             2.8449E+00
 PARAMETER:  5.1804E-02 -9.3928E-01 -1.8142E+00  1.6571E-01 -1.4598E+00  7.4353E-02 -6.5531E+00  3.4792E-01  3.6898E-01 -4.9849E-01
             1.1455E+00
 GRADIENT:  -1.2123E-01 -1.4348E-01 -2.3827E-01  4.3252E-01  8.4936E-01 -2.5797E-02  0.0000E+00  9.0159E-02 -1.0762E-03  9.9751E-02
            -3.6666E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1410
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.4509E-03 -2.0750E-04  1.9128E-02 -7.3870E-03  2.0180E-02
 SE:             2.8966E-02  1.3907E-04  2.0696E-02  2.6146E-02  1.8166E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3257E-01  1.3567E-01  3.5537E-01  7.7754E-01  2.6664E-01

 ETASHRINKSD(%)  2.9597E+00  9.9534E+01  3.0664E+01  1.2409E+01  3.9141E+01
 ETASHRINKVR(%)  5.8317E+00  9.9998E+01  5.1925E+01  2.3278E+01  6.2962E+01
 EBVSHRINKSD(%)  2.8424E+00  9.9525E+01  3.0775E+01  1.0948E+01  4.0120E+01
 EBVSHRINKVR(%)  5.6041E+00  9.9998E+01  5.2080E+01  2.0697E+01  6.4144E+01
 RELATIVEINF(%)  9.4289E+01  4.4455E-04  7.9827E+00  5.0129E+01  2.6202E+00
 EPSSHRINKSD(%)  2.8418E+01
 EPSSHRINKVR(%)  4.8760E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1501.7182347802841     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -582.77970157561140     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.13
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.74
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1501.718       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.53E-01  3.54E-01  1.47E-01  1.07E+00  2.10E-01  9.75E-01  1.00E-02  1.28E+00  1.31E+00  5.50E-01  2.84E+00
 


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
+        1.22E+03
 
 TH 2
+       -2.59E+01  2.15E+03
 
 TH 3
+       -1.14E+02  2.36E+03  1.47E+04
 
 TH 4
+       -2.00E+00  1.19E+02 -2.85E+02  4.42E+02
 
 TH 5
+        1.57E+02 -6.45E+03 -1.66E+04 -7.12E+02  3.06E+04
 
 TH 6
+        4.61E+00 -1.88E+01  1.20E+01 -6.71E+00  2.02E+01  1.85E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.91E+00  1.74E+01 -1.96E+01 -3.02E+00 -5.05E+01  3.21E+00  0.00E+00  3.54E+01
 
 TH 9
+        1.39E+01 -4.47E+01  1.49E+02 -1.64E+01  1.86E+02  5.41E-01  0.00E+00 -2.54E+00  6.43E+01
 
 TH10
+       -2.47E+00 -1.06E+02  4.38E+01  1.13E+01  3.77E+02  2.46E+00  0.00E+00  2.61E+01  1.56E+01  9.63E+01
 
 TH11
+       -2.18E+01 -2.22E+01 -2.92E+01 -9.90E-01  6.26E+01  1.33E+00  0.00E+00  4.94E+00  8.56E+00  8.18E+00  5.39E+01
 
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
 #CPUT: Total CPU Time in Seconds,       30.948
Stop Time:
Thu Sep 30 00:47:06 CDT 2021
