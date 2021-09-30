Thu Sep 30 03:32:45 CDT 2021
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
$DATA ../../../../data/spa1/D/dat72.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m72.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15869.2354075598        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.6404E+02  3.8189E+02 -2.4126E+01 -4.2007E+00  2.1475E+02 -2.2053E+03 -7.9919E+02 -9.4399E+01 -1.7667E+03 -6.8636E+02
            -3.0021E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -701.974263341438        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.4157E+00  9.4627E-01  8.2933E-01  2.2016E+00  1.2896E+00  2.8376E+00  1.2623E+00  9.9243E-01  2.1028E+00  1.2504E+00
             1.2661E+01
 PARAMETER:  4.4761E-01  4.4772E-02 -8.7138E-02  8.8917E-01  3.5432E-01  1.1430E+00  3.3290E-01  9.2405E-02  8.4329E-01  3.2347E-01
             2.6385E+00
 GRADIENT:   1.1975E+00  4.6123E+01 -3.0101E+01  1.0638E+02 -8.6018E+00  6.2012E+01 -2.2103E+00  7.8398E+00 -2.7001E+01  3.2092E+00
             1.4849E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -751.468736629235        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.3652E+00  9.9774E-01  2.2264E+00  2.3442E+00  3.2319E+00  3.2304E+00  2.6241E+00  5.5414E-01  3.7390E+00  5.5220E+00
             1.0259E+01
 PARAMETER:  4.1127E-01  9.7738E-02  9.0038E-01  9.5195E-01  1.2731E+00  1.2726E+00  1.0647E+00 -4.9034E-01  1.4188E+00  1.8087E+00
             2.4282E+00
 GRADIENT:   1.7701E+01  9.3358E+00 -1.1490E-01  6.9674E+01 -1.8909E+01  1.0969E+02  1.1826E+01 -1.3527E-01  6.6802E+01  1.0075E+01
             1.3671E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -792.330458367096        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1748E+00  6.0959E-01  8.4409E-01  1.3871E+00  9.2341E+00  2.3590E+00  1.5074E+00  1.4591E-02  2.0465E+00  1.4009E+01
             1.0020E+01
 PARAMETER:  2.6110E-01 -3.9497E-01 -6.9500E-02  4.2719E-01  2.3229E+00  9.5823E-01  5.1041E-01 -4.1273E+00  8.1615E-01  2.7397E+00
             2.4046E+00
 GRADIENT:  -8.9864E+00  7.2875E+00  2.6917E+01 -4.2408E+01 -9.2824E+00 -8.0742E+00  5.5514E+00 -2.3666E-03 -2.0692E+01  3.9960E+01
             1.3436E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -889.072634283059        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      299
 NPARAMETR:  8.1147E-01  6.0315E-02  1.0611E-01  8.9198E-01  7.8481E+00  1.7005E+00  5.1147E-02  3.4777E-02  1.1732E+00  2.1563E+00
             7.8897E+00
 PARAMETER: -1.0890E-01 -2.7082E+00 -2.1433E+00 -1.4307E-02  2.1603E+00  6.3094E-01 -2.8731E+00 -3.2588E+00  2.5976E-01  8.6842E-01
             2.1656E+00
 GRADIENT:   3.3227E+01  5.0142E+00 -4.9358E+01  2.1279E+02 -1.2713E+01 -1.0919E+02  8.3208E-03 -2.4910E-02 -5.6173E-01  4.3127E+00
            -1.0403E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -936.705582315271        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      440
 NPARAMETR:  5.9520E-01  1.6787E-02  4.9669E-02  4.5463E-01  6.8247E+00  2.0944E+00  1.0000E-02  7.7372E-01  9.4533E-01  9.8993E-01
             6.9837E+00
 PARAMETER: -4.1886E-01 -3.9872E+00 -2.9024E+00 -6.8827E-01  2.0205E+00  8.3928E-01 -7.2162E+00 -1.5654E-01  4.3778E-02  8.9879E-02
             2.0436E+00
 GRADIENT:   5.8928E+00  1.5547E+00  3.4786E+00  2.4726E+00 -2.2410E+01 -3.5631E-01  0.0000E+00  2.4805E+00  3.5472E+01  1.1233E+00
            -1.9603E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -963.012414400438        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      616
 NPARAMETR:  5.5133E-01  1.1166E-02  3.9721E-02  3.8584E-01  7.9496E+00  2.1360E+00  1.0000E-02  8.1565E-01  6.4400E-01  7.7199E-01
             8.1413E+00
 PARAMETER: -4.9543E-01 -4.3948E+00 -3.1259E+00 -8.5234E-01  2.1731E+00  8.5894E-01 -8.6721E+00 -1.0377E-01 -3.4005E-01 -1.5879E-01
             2.1970E+00
 GRADIENT:  -4.9072E+00  1.7858E-01  1.4483E+01 -1.4823E+01 -4.1605E+00  2.0316E+01  0.0000E+00 -3.1968E+00  8.3769E+00  9.8701E-02
            -1.2526E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -967.125732939855        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      791
 NPARAMETR:  5.3480E-01  1.0000E-02  3.6133E-02  3.6002E-01  1.1259E+01  1.9701E+00  1.0000E-02  1.0301E+00  4.2905E-01  7.0726E-01
             8.2649E+00
 PARAMETER: -5.2587E-01 -4.7939E+00 -3.2206E+00 -9.2159E-01  2.5211E+00  7.7810E-01 -8.4460E+00  1.2964E-01 -7.4618E-01 -2.4635E-01
             2.2120E+00
 GRADIENT:  -4.4535E+00  0.0000E+00  1.1662E+00 -4.0725E+00 -9.1013E-02 -3.2427E+00  0.0000E+00  2.3404E+00  4.3229E+00  1.6213E-04
             2.4564E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -968.644683975434        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      968
 NPARAMETR:  5.7091E-01  1.0000E-02  4.2913E-02  4.0669E-01  1.1062E+02  1.9901E+00  1.0000E-02  1.1564E+00  9.7104E-02  1.3889E+00
             8.2946E+00
 PARAMETER: -4.6053E-01 -6.0198E+00 -3.0486E+00 -7.9971E-01  4.8061E+00  7.8818E-01 -6.3115E+00  2.4531E-01 -2.2320E+00  4.2852E-01
             2.2156E+00
 GRADIENT:   5.2962E-01  0.0000E+00 -2.8666E+00  5.5510E+00  6.4766E-03 -1.3751E+00  0.0000E+00 -4.9167E+00  1.2202E-01  9.1126E-07
            -7.1106E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -968.779610150452        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1143
 NPARAMETR:  5.7152E-01  1.0000E-02  4.3434E-02  4.0899E-01  7.1121E+02  1.9963E+00  1.0000E-02  1.1992E+00  2.6750E-02  2.1323E+00
             8.2811E+00
 PARAMETER: -4.5945E-01 -7.1976E+00 -3.0365E+00 -7.9406E-01  6.6670E+00  7.9131E-01 -5.0146E+00  2.8170E-01 -3.5212E+00  8.5720E-01
             2.2140E+00
 GRADIENT:  -8.8958E-02  0.0000E+00  1.3405E-01 -8.4895E-02  5.3481E-04 -2.5957E-02  0.0000E+00 -1.4807E-01  1.3950E-02  1.0825E-07
            -4.8273E-02

0ITERATION NO.:   49    OBJECTIVE VALUE:  -968.784720544839        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1272
 NPARAMETR:  5.7166E-01  1.0000E-02  4.3442E-02  4.0907E-01  4.4509E+03  1.9964E+00  2.1084E-02  1.2009E+00  1.0000E-02  3.2288E+00
             8.2815E+00
 PARAMETER: -4.5922E-01 -8.3614E+00 -3.0363E+00 -7.9387E-01  8.5009E+00  7.9133E-01 -3.7592E+00  2.8307E-01 -4.7948E+00  1.2721E+00
             2.2140E+00
 GRADIENT:   1.9497E-02  0.0000E+00 -3.2030E-02  3.6355E-02  8.6165E-05 -1.0112E-02  3.6511E-07  1.2359E-02  0.0000E+00  6.5149E-09
             7.7124E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1272
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.5684E-03 -1.6421E-06  6.5370E-03 -3.8440E-04 -1.4402E-06
 SE:             2.9295E-02  1.6278E-06  2.4149E-02  2.8652E-04  8.0288E-07
 N:                     100         100         100         100         100

 P VAL.:         8.4924E-01  3.1307E-01  7.8662E-01  1.7972E-01  7.2848E-02

 ETASHRINKSD(%)  1.8589E+00  9.9995E+01  1.9099E+01  9.9040E+01  9.9997E+01
 ETASHRINKVR(%)  3.6833E+00  1.0000E+02  3.4550E+01  9.9991E+01  1.0000E+02
 EBVSHRINKSD(%)  1.7851E+00  9.9993E+01  1.8989E+01  9.9016E+01  9.9997E+01
 EBVSHRINKVR(%)  3.5383E+00  1.0000E+02  3.4371E+01  9.9990E+01  1.0000E+02
 RELATIVEINF(%)  5.4133E+00  3.0474E-08  6.1893E-01  6.0101E-05  4.7477E-09
 EPSSHRINKSD(%)  1.1496E+01
 EPSSHRINKVR(%)  2.1671E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -968.78472054483939     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -49.846187340166694     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.12
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     9.25
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -968.785       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         5.72E-01  1.00E-02  4.34E-02  4.09E-01  4.45E+03  2.00E+00  2.11E-02  1.20E+00  1.00E-02  3.23E+00  8.28E+00
 


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
+        8.16E+02
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -9.46E+02  0.00E+00  2.61E+05
 
 TH 4
+       -4.09E+02  0.00E+00 -3.82E+04  6.21E+03
 
 TH 5
+        1.05E-06  0.00E+00 -1.96E-05  2.65E-06 -2.05E-12
 
 TH 6
+        6.03E+00  0.00E+00  1.46E+02 -3.68E+01  5.68E-08  4.48E+01
 
 TH 7
+        7.88E-03  0.00E+00  1.99E-02  3.78E-03 -6.13E-08  1.31E-03  3.41E-02
 
 TH 8
+        1.86E+00  0.00E+00 -1.82E+02 -4.88E+01 -2.45E-07  1.96E+00  4.49E-03  6.16E+01
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -7.19E-05  0.00E+00  9.78E-04 -6.88E-05 -1.29E-09  5.83E-05 -1.17E-03  2.96E-04  0.00E+00  6.46E-06
 
 TH11
+       -1.05E+01  0.00E+00  1.75E+02 -2.72E+01 -3.66E-08  1.32E+00  1.27E-05  3.04E+00  0.00E+00  9.22E-06  7.38E+00
 
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
 #CPUT: Total CPU Time in Seconds,       30.444
Stop Time:
Thu Sep 30 03:33:17 CDT 2021
