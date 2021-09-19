Sat Sep 18 11:45:18 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat55.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m55.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1668.05694245072        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2376E+01 -9.4626E+01 -1.8608E+01 -9.6841E+01  7.5148E+01  5.5689E+01 -9.5537E+00 -3.2228E+00 -1.8604E+00 -1.7422E+01
             5.7267E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1677.23567922605        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  1.0075E+00  1.0528E+00  9.2520E-01  1.0269E+00  9.1440E-01  7.8299E-01  1.0486E+00  1.0394E+00  9.6134E-01  1.0702E+00
             9.7396E-01
 PARAMETER:  1.0746E-01  1.5144E-01  2.2257E-02  1.2653E-01  1.0515E-02 -1.4464E-01  1.4746E-01  1.3866E-01  6.0577E-02  1.6789E-01
             7.3619E-02
 GRADIENT:   7.1647E+01 -2.8201E+00  4.7827E+00 -1.5172E+01 -2.0655E+01 -3.3617E+01 -1.0646E+00  3.6115E+00 -1.5821E+00  1.0754E+01
            -3.6958E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1678.24364119171        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.0016E+00  1.0163E+00  6.6267E-01  1.0421E+00  7.6200E-01  8.0166E-01  1.2031E+00  7.4572E-01  8.6911E-01  8.1810E-01
             9.8635E-01
 PARAMETER:  1.0155E-01  1.1613E-01 -3.1147E-01  1.4122E-01 -1.7181E-01 -1.2108E-01  2.8488E-01 -1.9340E-01 -4.0282E-02 -1.0077E-01
             8.6254E-02
 GRADIENT:   3.7082E+01  7.7312E+00 -1.5444E+01  2.3589E+01  1.1753E+01 -2.4599E+01  3.7880E+00  4.8933E+00 -8.1292E+00  3.9489E+00
             1.8112E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1680.24911113884        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      276
 NPARAMETR:  9.9788E-01  1.0221E+00  7.4073E-01  1.0408E+00  8.1185E-01  8.5764E-01  1.1655E+00  6.5170E-01  9.3234E-01  8.9228E-01
             9.8355E-01
 PARAMETER:  9.7874E-02  1.2190E-01 -2.0012E-01  1.4000E-01 -1.0844E-01 -5.3571E-02  2.5317E-01 -3.2817E-01  2.9945E-02 -1.3979E-02
             8.3411E-02
 GRADIENT:  -1.1045E+01  5.7706E-01  4.5626E-01 -4.9108E-01 -3.5828E-01  6.7371E-01  1.2805E-01  4.5945E-01  2.4313E-01  5.8162E-01
             6.5570E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1680.42644102784        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      451
 NPARAMETR:  1.0022E+00  9.4377E-01  6.6355E-01  1.0751E+00  7.3308E-01  8.5689E-01  1.2500E+00  4.4251E-01  8.9609E-01  8.2178E-01
             9.8182E-01
 PARAMETER:  1.0221E-01  4.2130E-02 -3.1015E-01  1.7242E-01 -2.1050E-01 -5.4442E-02  3.2312E-01 -7.1528E-01 -9.7154E-03 -9.6283E-02
             8.1648E-02
 GRADIENT:   6.9062E-02 -2.4130E-01 -3.8155E-02 -2.4892E-01  8.4176E-02  6.5718E-02  1.4897E-02  2.3175E-02  1.9598E-02  2.5116E-02
            -4.1410E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1680.43306471758        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:      633
 NPARAMETR:  1.0022E+00  1.0072E+00  6.3605E-01  1.0352E+00  7.4294E-01  8.5668E-01  1.1832E+00  3.8711E-01  9.2181E-01  8.2443E-01
             9.8203E-01
 PARAMETER:  1.0222E-01  1.0721E-01 -3.5249E-01  1.3457E-01 -1.9714E-01 -5.4691E-02  2.6819E-01 -8.4906E-01  1.8582E-02 -9.3065E-02
             8.1867E-02
 GRADIENT:  -8.6848E-01  7.2664E-01  6.0025E-01 -1.4157E-02 -1.5197E+00 -1.7122E-01 -1.7888E-03  5.2401E-02  4.8371E-02  3.5097E-01
             1.5341E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1680.45431159557        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      810
 NPARAMETR:  1.0029E+00  1.1028E+00  5.8153E-01  9.7091E-01  7.5364E-01  8.5725E-01  1.0994E+00  2.3274E-01  9.6390E-01  8.1689E-01
             9.8141E-01
 PARAMETER:  1.0292E-01  1.9787E-01 -4.4209E-01  7.0474E-02 -1.8285E-01 -5.4027E-02  1.9476E-01 -1.3578E+00  6.3228E-02 -1.0226E-01
             8.1234E-02
 GRADIENT:   2.2337E-01 -9.8317E-01 -1.1782E+00 -1.6322E-01  9.8343E-01 -1.0299E-01  3.5316E-01  6.2068E-02  2.2133E-01  3.1040E-01
             2.1279E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1680.47206916395        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      986
 NPARAMETR:  1.0029E+00  1.2132E+00  5.5539E-01  9.0380E-01  7.9076E-01  8.5743E-01  1.0160E+00  1.0615E-01  1.0187E+00  8.4007E-01
             9.8159E-01
 PARAMETER:  1.0286E-01  2.9322E-01 -4.8809E-01 -1.1500E-03 -1.3476E-01 -5.3811E-02  1.1583E-01 -2.1429E+00  1.1857E-01 -7.4271E-02
             8.1423E-02
 GRADIENT:   2.5561E-02  2.4071E-01  3.9503E-02  1.6754E-03 -1.7369E-01 -1.1917E-02  7.8314E-02  1.1104E-02 -2.0588E-02  7.0725E-02
            -4.7468E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1680.47340940904        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1164
 NPARAMETR:  1.0029E+00  1.2242E+00  5.5353E-01  8.9679E-01  7.9542E-01  8.5746E-01  1.0087E+00  8.2053E-02  1.0257E+00  8.4317E-01
             9.8154E-01
 PARAMETER:  1.0289E-01  3.0227E-01 -4.9144E-01 -8.9284E-03 -1.2888E-01 -5.3775E-02  1.0869E-01 -2.4004E+00  1.2540E-01 -7.0592E-02
             8.1366E-02
 GRADIENT:   1.4327E-01 -9.6858E-02  2.3848E-01 -4.1614E-01 -1.2046E-01  8.7401E-03  8.6680E-02  5.8025E-03  9.9297E-02 -1.5969E-02
            -9.8679E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1680.47738476396        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1340
 NPARAMETR:  1.0029E+00  1.2012E+00  5.5286E-01  9.1018E-01  7.8317E-01  8.5750E-01  1.0236E+00  2.2697E-02  1.0123E+00  8.3314E-01
             9.8165E-01
 PARAMETER:  1.0287E-01  2.8332E-01 -4.9265E-01  5.8827E-03 -1.4441E-01 -5.3735E-02  1.2332E-01 -3.6855E+00  1.1219E-01 -8.2557E-02
             8.1478E-02
 GRADIENT:  -8.7880E-03 -8.1465E-02 -1.4357E-01  1.2039E-01  1.9553E-01 -1.0139E-03 -1.7178E-03  4.7198E-04  6.7100E-03 -1.7488E-03
             3.0939E-04

0ITERATION NO.:   49    OBJECTIVE VALUE:  -1680.47758577701        NO. OF FUNC. EVALS.: 127
 CUMULATIVE NO. OF FUNC. EVALS.:     1467
 NPARAMETR:  1.0029E+00  1.2034E+00  5.5275E-01  9.0889E-01  7.8414E-01  8.5750E-01  1.0222E+00  1.0000E-02  1.0134E+00  8.3401E-01
             9.8167E-01
 PARAMETER:  1.0287E-01  2.8512E-01 -4.9285E-01  4.4699E-03 -1.4317E-01 -5.3738E-02  1.2191E-01 -4.6189E+00  1.1334E-01 -8.1505E-02
             8.1502E-02
 GRADIENT:  -3.9910E-03  1.1732E-02  3.3790E-03  5.8465E-03 -4.8960E-03  1.2317E-04 -1.9814E-03  0.0000E+00 -5.2959E-03 -1.3257E-03
             2.0640E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1467
 NO. OF SIG. DIGITS IN FINAL EST.:  3.0
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4430E-04 -1.2384E-02 -4.1771E-04  7.0781E-03 -1.7656E-02
 SE:             2.9799E-02  2.2691E-02  1.7518E-04  2.5023E-02  2.2297E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9614E-01  5.8522E-01  1.7107E-02  7.7728E-01  4.2845E-01

 ETASHRINKSD(%)  1.6993E-01  2.3983E+01  9.9413E+01  1.6168E+01  2.5301E+01
 ETASHRINKVR(%)  3.3957E-01  4.2214E+01  9.9997E+01  2.9722E+01  4.4200E+01
 EBVSHRINKSD(%)  5.5443E-01  2.3742E+01  9.9497E+01  1.6448E+01  2.5136E+01
 EBVSHRINKVR(%)  1.1058E+00  4.1847E+01  9.9997E+01  3.0190E+01  4.3953E+01
 RELATIVEINF(%)  9.8859E+01  3.5170E+00  2.7485E-04  5.1619E+00  4.8138E+00
 EPSSHRINKSD(%)  4.5152E+01
 EPSSHRINKVR(%)  6.9917E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1680.4775857770119     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -945.32675921327370     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    16.39
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.35
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1680.478       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.20E+00  5.53E-01  9.09E-01  7.84E-01  8.57E-01  1.02E+00  1.00E-02  1.01E+00  8.34E-01  9.82E-01
 


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
+        1.49E+03
 
 TH 2
+       -6.82E+00  4.43E+02
 
 TH 3
+        9.30E+00  3.14E+02  1.02E+03
 
 TH 4
+       -1.77E+01  2.86E+02 -4.50E+02  9.17E+02
 
 TH 5
+       -3.06E+00 -4.64E+02 -1.03E+03  4.32E+02  1.45E+03
 
 TH 6
+        2.40E+00 -1.61E+00  2.70E+00 -2.38E+00 -2.58E+00  2.65E+02
 
 TH 7
+        1.59E+00  2.36E+01 -4.14E+01 -4.28E+00  4.68E+00 -4.78E-01  7.03E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.92E+00 -2.51E+01 -3.43E+01  4.08E+01 -1.09E+01  9.64E-01  1.56E+01  0.00E+00  1.04E+02
 
 TH10
+       -4.70E+00 -1.53E+01 -7.85E+01 -1.81E+01 -5.99E+01  1.03E-01  2.43E+01  0.00E+00  1.01E+01  9.68E+01
 
 TH11
+       -1.17E+01 -1.23E+01 -3.17E+01 -2.96E+00  3.28E+00  1.79E+00  5.41E+00  0.00E+00  7.43E+00  1.77E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       21.805
Stop Time:
Sat Sep 18 11:45:41 CDT 2021
