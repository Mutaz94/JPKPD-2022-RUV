Thu Sep 30 01:49:15 CDT 2021
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
$DATA ../../../../data/spa1/TD2/dat2.csv ignore=@
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
 (E4.0,E3.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m2.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2140.18576235256        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.4684E+02 -2.3086E+01 -1.8401E+01  8.3882E+00 -2.2037E+01  3.6300E+01  1.0963E+01  1.8932E+01  4.4151E+01  9.3626E+00
             7.8887E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2151.55824301586        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.7154E-01  9.9378E-01  1.2065E+00  1.0455E+00  1.1472E+00  9.2386E-01  9.1697E-01  7.9565E-01  7.5818E-01  1.0861E+00
             9.8982E-01
 PARAMETER:  7.1122E-02  9.3759E-02  2.8771E-01  1.4453E-01  2.3735E-01  2.0809E-02  1.3320E-02 -1.2859E-01 -1.7684E-01  1.8264E-01
             8.9771E-02
 GRADIENT:   3.9196E+02  2.0457E-01 -1.9000E+01  8.1970E+01  3.9267E+01  7.7052E+00 -7.9331E+00  4.7585E+00 -9.4390E+00 -7.4432E+00
            -5.0321E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2153.00715927601        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      202
 NPARAMETR:  9.6535E-01  1.0420E+00  1.2698E+00  1.0076E+00  1.2242E+00  9.4041E-01  9.3109E-01  6.0613E-01  8.1749E-01  1.2062E+00
             9.7313E-01
 PARAMETER:  6.4735E-02  1.4111E-01  3.3884E-01  1.0755E-01  3.0225E-01  3.8557E-02  2.8599E-02 -4.0066E-01 -1.0152E-01  2.8750E-01
             7.2765E-02
 GRADIENT:  -3.5808E+01 -4.1941E+01 -7.4644E+00 -4.7025E+01  2.8509E+01 -2.5159E+01 -2.4237E+00 -1.4406E+00 -5.9135E+00 -3.2524E+00
            -2.0608E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2155.62443301184        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.8010E-01  1.0841E+00  1.3553E+00  1.0105E+00  1.2385E+00  9.9757E-01  8.5925E-01  8.0830E-01  8.7339E-01  1.2298E+00
             9.9777E-01
 PARAMETER:  7.9900E-02  1.8071E-01  4.0405E-01  1.1044E-01  3.1394E-01  9.7567E-02 -5.1695E-02 -1.1282E-01 -3.5376E-02  3.0684E-01
             9.7771E-02
 GRADIENT:   1.0723E+00 -5.6591E-01  1.1820E-01  1.4830E+00 -3.1883E-01  6.5611E-02  3.2112E-01 -3.0285E-01  8.9248E-01  1.4578E-01
             6.4895E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2155.69617063055        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      552
 NPARAMETR:  9.8008E-01  1.1852E+00  1.3639E+00  9.4588E-01  1.2906E+00  9.9751E-01  7.3418E-01  9.0016E-01  9.4627E-01  1.2654E+00
             9.9976E-01
 PARAMETER:  7.9877E-02  2.6988E-01  4.1038E-01  4.4366E-02  3.5510E-01  9.7510E-02 -2.0900E-01 -5.1864E-03  4.4769E-02  3.3541E-01
             9.9758E-02
 GRADIENT:  -4.3325E-01  8.2086E-01  4.8802E-01  2.0748E+00 -1.2726E+00 -1.3040E-01 -9.8217E-01 -1.1088E-01 -3.4648E-01  2.1979E-01
             1.7824E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2155.84195498114        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      730
 NPARAMETR:  9.7975E-01  1.3814E+00  1.0932E+00  8.2149E-01  1.2825E+00  9.9550E-01  7.4215E-01  6.9285E-01  1.0160E+00  1.2317E+00
             1.0010E+00
 PARAMETER:  7.9547E-02  4.2308E-01  1.8912E-01 -9.6631E-02  3.4878E-01  9.5495E-02 -1.9820E-01 -2.6694E-01  1.1590E-01  3.0842E-01
             1.0096E-01
 GRADIENT:  -4.9694E+00  1.2205E+01 -9.0833E-01  1.3107E+01 -3.6793E+00 -1.6152E+00 -3.0708E-01  5.3864E-01  5.0841E-01  7.8778E-01
             3.9801E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2156.05126700780        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      906
 NPARAMETR:  9.8450E-01  1.5723E+00  9.1088E-01  6.9462E-01  1.3201E+00  1.0020E+00  7.1378E-01  4.9603E-01  1.1151E+00  1.2289E+00
             9.9016E-01
 PARAMETER:  8.4379E-02  5.5251E-01  6.6548E-03 -2.6438E-01  3.7771E-01  1.0197E-01 -2.3717E-01 -6.0113E-01  2.0899E-01  3.0609E-01
             9.0113E-02
 GRADIENT:   3.3551E+00  1.5286E+01  9.3927E-01  1.0295E+01 -1.9708E+00  4.6447E-01  3.0230E-01  1.8263E-01 -1.1424E+00 -1.0060E+00
            -5.1949E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2156.24678629761        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1082
 NPARAMETR:  9.8669E-01  1.7809E+00  7.3709E-01  5.5196E-01  1.3876E+00  1.0067E+00  6.5752E-01  2.5952E-01  1.3004E+00  1.2619E+00
             9.9163E-01
 PARAMETER:  8.6602E-02  6.7714E-01 -2.0505E-01 -4.9429E-01  4.2756E-01  1.0666E-01 -3.1928E-01 -1.2489E+00  3.6264E-01  3.3265E-01
             9.1598E-02
 GRADIENT:   6.2194E+00  5.7443E+00 -5.1297E-01  4.9949E+00 -1.2551E-03  1.9894E+00 -1.5968E+00  1.2226E-01 -1.8261E+00 -1.5324E-01
            -3.7153E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2156.35287112385        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1264
 NPARAMETR:  9.8423E-01  1.7825E+00  7.3229E-01  5.4667E-01  1.3897E+00  1.0029E+00  6.5823E-01  1.0237E-01  1.3206E+00  1.2630E+00
             9.9500E-01
 PARAMETER:  8.4103E-02  6.7804E-01 -2.1158E-01 -5.0390E-01  4.2911E-01  1.0293E-01 -3.1821E-01 -2.1792E+00  3.7805E-01  3.3352E-01
             9.4983E-02
 GRADIENT:   7.5069E-01 -3.7941E+00 -2.0697E-01  1.2600E+00  7.1290E-01  5.3526E-01 -5.2259E-01  1.7586E-02 -4.5931E-01 -2.0821E-02
            -6.5493E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2156.36382769840        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1445
 NPARAMETR:  9.8444E-01  1.7841E+00  7.3070E-01  5.4410E-01  1.3896E+00  1.0016E+00  6.6066E-01  3.4530E-02  1.3276E+00  1.2632E+00
             9.9582E-01
 PARAMETER:  8.4321E-02  6.7890E-01 -2.1375E-01 -5.0863E-01  4.2905E-01  1.0158E-01 -3.1452E-01 -3.2659E+00  3.8334E-01  3.3366E-01
             9.5807E-02
 GRADIENT:   1.2266E+00 -6.8582E+00  2.1005E-01 -4.7448E-01 -2.5020E-01  7.5132E-03  2.9233E-01  2.0225E-03  1.1668E-01  9.3758E-02
             1.1267E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2156.36523919074        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1612
 NPARAMETR:  9.8459E-01  1.7838E+00  7.3019E-01  5.4431E-01  1.3891E+00  1.0019E+00  6.5972E-01  1.0000E-02  1.3274E+00  1.2627E+00
             9.9569E-01
 PARAMETER:  8.4469E-02  6.7876E-01 -2.1445E-01 -5.0823E-01  4.2865E-01  1.0185E-01 -3.1595E-01 -4.9588E+00  3.8325E-01  3.3329E-01
             9.5684E-02
 GRADIENT:   1.5443E+00 -6.7649E+00  1.2071E-01 -3.4230E-01 -1.9928E-01  1.1203E-01  1.0245E-01  0.0000E+00  6.7136E-02  9.2052E-02
             8.6000E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1612
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0195E-04 -3.3461E-02 -2.6252E-04  2.7248E-02 -3.9730E-02
 SE:             2.9867E-02  2.1566E-02  9.5620E-05  2.2811E-02  2.3513E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8659E-01  1.2076E-01  6.0419E-03  2.3228E-01  9.1082E-02

 ETASHRINKSD(%)  1.0000E-10  2.7751E+01  9.9680E+01  2.3580E+01  2.1229E+01
 ETASHRINKVR(%)  1.0000E-10  4.7801E+01  9.9999E+01  4.1600E+01  3.7951E+01
 EBVSHRINKSD(%)  3.2945E-01  2.5475E+01  9.9725E+01  2.7039E+01  1.7645E+01
 EBVSHRINKVR(%)  6.5782E-01  4.4460E+01  9.9999E+01  4.6766E+01  3.2177E+01
 RELATIVEINF(%)  9.9216E+01  2.9126E+00  1.1941E-04  2.8290E+00  2.0746E+01
 EPSSHRINKSD(%)  3.2515E+01
 EPSSHRINKVR(%)  5.4457E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2156.3652391907431     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1237.4267059860704     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    24.49
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.43
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2156.365       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  1.78E+00  7.30E-01  5.44E-01  1.39E+00  1.00E+00  6.60E-01  1.00E-02  1.33E+00  1.26E+00  9.96E-01
 


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
+        1.13E+03
 
 TH 2
+       -6.00E+00  4.72E+02
 
 TH 3
+        4.07E+00  9.76E+01  1.91E+02
 
 TH 4
+       -6.19E+00  5.18E+02 -1.71E+02  1.14E+03
 
 TH 5
+        1.71E+00 -1.01E+02 -1.30E+02  1.38E+02  2.46E+02
 
 TH 6
+        6.69E-01 -1.25E+00  9.95E-01 -2.77E+00 -4.34E-01  1.96E+02
 
 TH 7
+        1.30E+00 -8.55E+00  9.74E+00 -1.37E+01 -1.74E+01 -4.48E-01  1.46E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.36E-01 -1.77E+01 -2.49E+01  5.72E+01  4.96E+00 -2.99E-01  3.68E+01  0.00E+00  4.05E+01
 
 TH10
+        1.34E+00 -1.07E+01 -2.17E+01  8.05E+00 -3.72E+01  1.11E-01  4.29E+00  0.00E+00  4.39E+00  6.18E+01
 
 TH11
+       -8.19E+00 -2.50E+01 -2.78E+01  2.50E+00  2.91E+00  1.86E+00  1.01E+01  0.00E+00  5.18E+00  1.34E+01  4.25E+02
 
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
 #CPUT: Total CPU Time in Seconds,       31.982
Stop Time:
Thu Sep 30 01:49:48 CDT 2021
