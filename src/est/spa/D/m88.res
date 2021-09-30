Wed Sep 29 20:23:58 CDT 2021
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
$DATA ../../../../data/spa/D/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22542.9955871431        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   7.8001E+02  5.2351E+02  7.9985E-01  5.4813E+02  2.4602E+02 -2.6943E+03 -1.1019E+03 -8.4407E+01 -1.5472E+03 -7.5490E+02
            -4.1639E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -473.904654226163        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.2374E+00  1.0078E+00  8.1148E-01  1.4795E+00  1.6624E+00  1.9743E+00  1.0280E+00  9.5412E-01  7.8992E-01  9.4971E-01
             1.5135E+01
 PARAMETER:  3.1305E-01  1.0777E-01 -1.0890E-01  4.9169E-01  6.0828E-01  7.8021E-01  1.2766E-01  5.3037E-02 -1.3582E-01  4.8399E-02
             2.8170E+00
 GRADIENT:   5.8970E+00  1.4581E+01 -3.6961E+00  2.1041E+01 -3.9407E+00  5.8730E+01 -4.2202E-01  4.0129E+00  3.5110E+00  2.1118E-01
             1.5636E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -489.244836448482        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.2279E+00  9.0292E-01  7.2350E-01  1.4452E+00  3.4695E+00  1.4713E+00  1.2583E+00  2.4782E-01  3.5451E-01  1.4672E+00
             1.6239E+01
 PARAMETER:  3.0531E-01 -2.1224E-03 -2.2366E-01  4.6826E-01  1.3440E+00  4.8617E-01  3.2972E-01 -1.2951E+00 -9.3703E-01  4.8337E-01
             2.8874E+00
 GRADIENT:  -7.4538E-01  2.3941E+01 -2.5527E+00  3.4626E+01 -8.9894E+00 -1.0100E+01  8.6642E-01  2.7641E-01  1.6310E+00  6.3364E-01
             4.1524E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -512.277211472209        NO. OF FUNC. EVALS.:  79
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  1.0194E+00  3.0096E-01  6.9445E-01  1.4959E+00  7.5115E+00  1.4708E+00  1.6851E+00  1.5194E-01  3.4814E-01  1.0220E+01
             1.4005E+01
 PARAMETER:  1.1925E-01 -1.1008E+00 -2.6464E-01  5.0270E-01  2.1164E+00  4.8583E-01  6.2183E-01 -1.7842E+00 -9.5515E-01  2.4244E+00
             2.7394E+00
 GRADIENT:  -5.6044E+01  9.3206E+00  1.2892E+01  1.6705E+01 -4.3202E-01 -6.7421E+00  6.8471E-01  1.9100E-02  4.2751E+00  1.3101E+01
             1.6701E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -520.552596538312        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      307
 NPARAMETR:  1.0686E+00  7.4554E-02  6.3601E-01  1.5783E+00  5.3875E+00  1.4504E+00  2.9241E+00  6.0607E-02  4.7396E-02  8.3882E+00
             1.4306E+01
 PARAMETER:  1.6634E-01 -2.4962E+00 -3.5254E-01  5.5636E-01  1.7841E+00  4.7183E-01  1.1730E+00 -2.7034E+00 -2.9492E+00  2.2268E+00
             2.7607E+00
 GRADIENT:  -2.8326E+00  2.2814E+00  3.4174E+00  2.1876E+01 -5.7334E+00  2.6770E+00  1.3138E-01  5.6466E-04  9.9193E-02 -1.4355E+00
            -7.7356E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -523.052160800816        NO. OF FUNC. EVALS.: 117
 CUMULATIVE NO. OF FUNC. EVALS.:      424
 NPARAMETR:  1.0732E+00  1.9265E-02  5.9377E-01  1.6235E+00  6.9655E+00  1.4945E+00  3.5884E+00  3.1599E-02  1.0000E-02  9.2529E+00
             1.4165E+01
 PARAMETER:  1.7062E-01 -3.8494E+00 -4.2126E-01  5.8457E-01  2.0410E+00  5.0176E-01  1.3777E+00 -3.3546E+00 -4.5156E+00  2.3249E+00
             2.7508E+00
 GRADIENT:  -1.1641E+01  4.3620E-01 -4.8790E+00  1.2570E+01  7.9906E+00  6.9210E+00  1.0532E-02 -1.5540E-03  1.8105E-03 -2.3166E+01
            -1.4401E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -524.187373946655        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      570
 NPARAMETR:  1.0780E+00  1.0000E-02  5.9678E-01  1.6205E+00  6.9572E+00  1.5123E+00  2.5213E+00  6.8855E-02  1.0000E-02  9.3986E+00
             1.4468E+01
 PARAMETER:  1.7513E-01 -5.0610E+00 -4.1620E-01  5.8272E-01  2.0398E+00  5.1364E-01  1.0248E+00 -2.5758E+00 -4.5516E+00  2.3406E+00
             2.7720E+00
 GRADIENT:  -3.5242E+00  0.0000E+00 -3.6625E-01 -8.8561E+00  9.0300E+00  1.9508E+01  1.6864E-03 -7.6314E-03  0.0000E+00 -1.4914E+01
             3.6710E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -538.186679527713        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  1.0792E+00  1.0000E-02  5.9778E-01  1.5860E+00  6.8809E+00  1.5161E+00  2.4498E+00  3.2637E+00  1.0000E-02  9.5640E+00
             1.3358E+01
 PARAMETER:  1.7618E-01 -5.5168E+00 -4.1453E-01  5.6121E-01  2.0287E+00  5.1612E-01  9.9602E-01  1.2829E+00 -4.5508E+00  2.3580E+00
             2.6921E+00
 GRADIENT:   4.6203E+00  0.0000E+00  2.2054E+01  4.7979E+01  7.9620E+00 -3.6114E+01 -2.3753E-04  3.9262E-02  0.0000E+00  7.2772E+00
             1.6734E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -540.914307826382        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      716
 NPARAMETR:  1.0790E+00  1.0000E-02  4.9913E-01  1.4795E+00  6.8075E+00  1.5493E+00  2.8511E+00  2.8419E+00  1.0000E-02  9.2847E+00
             1.2999E+01
 PARAMETER:  1.7600E-01 -6.1483E+00 -5.9488E-01  4.9169E-01  2.0180E+00  5.3784E-01  1.1477E+00  1.1445E+00 -4.5620E+00  2.3284E+00
             2.6649E+00
 GRADIENT:   3.5289E+01  0.0000E+00  1.0330E+01  2.0825E+01 -7.7220E+00 -2.8925E+01  1.6909E-03 -1.5471E+01  0.0000E+00  1.8663E+01
            -3.9308E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -541.854814041620        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      788
 NPARAMETR:  1.0642E+00  1.0000E-02  3.8629E-01  1.3497E+00  6.5556E+00  1.6130E+00  2.5332E+00  3.0332E+00  1.0000E-02  8.7514E+00
             1.2669E+01
 PARAMETER:  1.6226E-01 -6.8727E+00 -8.5116E-01  3.9988E-01  1.9803E+00  5.7812E-01  1.0295E+00  1.2096E+00 -4.5728E+00  2.2692E+00
             2.6391E+00
 GRADIENT:   6.6376E+01  0.0000E+00 -6.3311E+00  4.3628E+00 -1.7863E+01 -1.7047E+01  1.5895E-03  5.2864E-01  0.0000E+00  3.8507E+01
            -4.0982E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -543.655418049929        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      864
 NPARAMETR:  9.9101E-01  1.0000E-02  3.0734E-01  1.2651E+00  6.0349E+00  1.7056E+00  9.0836E-01  3.1926E+00  1.0000E-02  7.5685E+00
             1.2711E+01
 PARAMETER:  9.0974E-02 -6.7798E+00 -1.0798E+00  3.3515E-01  1.8976E+00  6.3389E-01  3.8848E-03  1.2608E+00 -4.5584E+00  2.1240E+00
             2.6425E+00
 GRADIENT:   2.7947E+01  0.0000E+00 -1.8578E+01  3.8222E+01 -2.2008E+01  9.3094E+00  1.0731E-04  1.7056E+01  0.0000E+00  2.8486E+01
            -6.0066E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -545.439646982430        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:      941
 NPARAMETR:  8.6791E-01  1.0000E-02  2.2985E-01  1.1432E+00  5.4986E+00  1.8604E+00  9.1295E-02  2.8339E+00  1.0000E-02  5.3583E+00
             1.1829E+01
 PARAMETER: -4.1671E-02 -6.4422E+00 -1.3703E+00  2.3380E-01  1.8045E+00  7.2078E-01 -2.2937E+00  1.1417E+00 -4.5504E+00  1.7786E+00
             2.5706E+00
 GRADIENT:  -1.1018E+01  0.0000E+00 -2.6530E+01  7.8636E+01 -3.6356E+01  3.2203E+01  1.1147E-06  1.4525E+01  1.2086E-05  3.1642E+01
            -4.1367E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -546.216325685715        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1013
 NPARAMETR:  8.3110E-01  1.0000E-02  2.0944E-01  1.0932E+00  5.3980E+00  1.8937E+00  2.6988E-02  2.3548E+00  1.0000E-02  4.4559E+00
             1.1257E+01
 PARAMETER: -8.5009E-02 -6.3209E+00 -1.4633E+00  1.8911E-01  1.7860E+00  7.3852E-01 -3.5123E+00  9.5645E-01 -4.5555E+00  1.5942E+00
             2.5210E+00
 GRADIENT:  -1.5485E+01  0.0000E+00 -2.0385E+01  8.7576E+01 -3.7067E+01  3.2800E+01  1.7028E-07 -2.5851E+00  0.0000E+00  3.0662E+01
            -8.2207E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -547.635147016963        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1088
 NPARAMETR:  8.0465E-01  1.0000E-02  1.4417E-01  9.3449E-01  5.0373E+00  1.8612E+00  1.0000E-02  1.5252E+00  1.0000E-02  3.0086E+00
             1.0981E+01
 PARAMETER: -1.1735E-01 -6.6029E+00 -1.8368E+00  3.2249E-02  1.7169E+00  7.2121E-01 -6.2037E+00  5.2215E-01 -4.6531E+00  1.2015E+00
             2.4961E+00
 GRADIENT:   1.5160E+01  0.0000E+00 -2.4010E+01  9.0136E+01 -5.2204E+01  2.3846E+01  0.0000E+00 -2.4723E+01  0.0000E+00  3.4022E+01
            -1.2332E+02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -557.630537544796        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  7.9782E-01  1.0000E-02  9.9830E-02  7.8877E-01  4.3406E+00  1.5900E+00  1.0000E-02  1.2686E+00  1.0000E-02  1.6796E+00
             1.2509E+01
 PARAMETER: -1.2587E-01 -5.8723E+00 -2.2043E+00 -1.3728E-01  1.5680E+00  5.6377E-01 -1.0972E+01  3.3791E-01 -4.9201E+00  6.1858E-01
             2.6264E+00
 GRADIENT:   4.1954E+01  0.0000E+00 -1.7357E+01  3.3302E+01 -7.2739E+01 -2.5702E+01  0.0000E+00 -5.5123E+00  0.0000E+00  1.5936E+01
            -2.5616E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -563.217542114752        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1232
 NPARAMETR:  7.2196E-01  1.0000E-02  8.6289E-02  7.3454E-01  4.9281E+00  1.5646E+00  1.0000E-02  1.5559E+00  1.0000E-02  1.1717E+00
             1.1638E+01
 PARAMETER: -2.2579E-01 -4.9462E+00 -2.3500E+00 -2.0852E-01  1.6950E+00  5.4764E-01 -1.4277E+01  5.4207E-01 -4.9751E+00  2.5844E-01
             2.5543E+00
 GRADIENT:   1.2366E+01  0.0000E+00 -2.9915E+01  8.3260E+01 -3.2738E+01 -3.2523E+01  0.0000E+00  5.2454E+00  0.0000E+00  3.5855E+00
            -5.8178E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -563.229949923057        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  7.0623E-01  1.0000E-02  7.8599E-02  6.9750E-01  5.0176E+00  1.5281E+00  1.0000E-02  1.5180E+00  1.0000E-02  9.7919E-01
             1.1609E+01
 PARAMETER: -2.4781E-01 -4.6960E+00 -2.4434E+00 -2.6025E-01  1.7130E+00  5.2403E-01 -1.5948E+01  5.1742E-01 -5.0360E+00  7.8970E-02
             2.5518E+00
 GRADIENT:   1.7039E+01  0.0000E+00 -3.4403E+01  8.9179E+01 -2.8970E+01 -4.1195E+01  0.0000E+00  5.0914E+00  0.0000E+00  3.0986E+00
            -6.1508E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -563.242180316904        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:     1400
 NPARAMETR:  6.9263E-01  1.0042E-02  7.3404E-02  6.6998E-01  5.0650E+00  1.5088E+00  1.0000E-02  1.4863E+00  1.0000E-02  8.5289E-01
             1.1585E+01
 PARAMETER: -2.6725E-01 -4.5010E+00 -2.5118E+00 -3.0050E-01  1.7223E+00  5.1132E-01 -1.7266E+01  4.9626E-01 -5.0816E+00 -5.9122E-02
             2.5497E+00
 GRADIENT:   6.0896E+00  3.5960E-01 -5.4174E+01  8.7408E+01 -3.1240E+01 -5.1754E+01  0.0000E+00  4.3740E+00  0.0000E+00  2.8384E+00
            -8.5814E+01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -577.797147321974        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1577
 NPARAMETR:  6.8892E-01  1.2302E-02  7.5701E-02  6.3715E-01  8.1732E+00  1.6507E+00  1.0000E-02  1.2294E+00  1.0000E-02  8.6270E-01
             1.2755E+01
 PARAMETER: -2.7263E-01 -4.2980E+00 -2.4810E+00 -3.5075E-01  2.2009E+00  6.0120E-01 -2.2338E+01  3.0653E-01 -5.3204E+00 -4.7688E-02
             2.6460E+00
 GRADIENT:  -2.2399E+00 -3.0057E-02  2.2023E+00  9.0797E+00 -3.3862E-01 -1.0705E+01  0.0000E+00  7.1051E-01  0.0000E+00  1.1114E-02
             1.0307E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -581.171925363958        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1763             RESET HESSIAN, TYPE I
 NPARAMETR:  5.2298E-01  2.2032E-02  3.3312E-02  3.6282E-01  1.4653E+01  1.6401E+00  1.0000E-02  6.3549E-01  1.0000E-02  1.5801E-01
             1.2943E+01
 PARAMETER: -5.4821E-01 -3.7152E+00 -3.3018E+00 -9.1386E-01  2.7847E+00  5.9474E-01 -4.3533E+01 -3.5335E-01 -6.0794E+00 -1.7451E+00
             2.6605E+00
 GRADIENT:   2.9246E+01  2.2704E+00  3.7275E+01  5.1400E+00 -1.3492E-01  9.3609E+00  0.0000E+00  8.4008E-01  0.0000E+00  2.1931E-04
             1.6563E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -581.587987876218        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:     1836
 NPARAMETR:  5.0133E-01  2.4535E-02  3.3034E-02  3.6364E-01  1.3188E+01  1.5738E+00  1.0000E-02  2.4007E-01  1.0000E-02  1.1840E-01
             1.3117E+01
 PARAMETER: -5.9050E-01 -3.6076E+00 -3.3102E+00 -9.1159E-01  2.6793E+00  5.5351E-01 -4.3533E+01 -1.3268E+00 -6.0794E+00 -2.0337E+00
             2.6739E+00
 GRADIENT:   4.8776E+00  5.6383E+00  4.9125E+01  4.8935E-01 -4.6278E-01  8.0282E-01  0.0000E+00  9.7977E-01  0.0000E+00  4.3879E-04
             1.5905E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -581.964147636227        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1906
 NPARAMETR:  4.6499E-01  2.0164E-02  2.8173E-02  3.2608E-01  9.4985E+00  1.5559E+00  1.0000E-02  1.0000E-02  1.0000E-02  4.1642E-02
             1.3018E+01
 PARAMETER: -6.6573E-01 -3.8039E+00 -3.4694E+00 -1.0206E+00  2.3511E+00  5.4208E-01 -4.3533E+01 -5.4085E+00 -6.0794E+00 -3.0787E+00
             2.6663E+00
 GRADIENT:  -1.7165E+00  4.9725E+00  4.6777E+01  1.6115E+01 -4.8469E-01 -2.1739E+00  0.0000E+00  0.0000E+00  0.0000E+00  2.7451E-04
             1.1782E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -583.282483764644        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2088
 NPARAMETR:  4.8210E-01  1.6454E-02  2.6713E-02  3.1301E-01  1.0208E+01  1.5900E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.9045E-02
             1.3588E+01
 PARAMETER: -6.2959E-01 -4.0072E+00 -3.5226E+00 -1.0615E+00  2.4232E+00  5.6374E-01 -4.3533E+01 -5.7938E+00 -6.0794E+00 -3.1430E+00
             2.7092E+00
 GRADIENT:   3.9856E+00  4.0260E-01 -2.5623E+00  7.8336E-01  2.6073E-01  1.7656E+00  0.0000E+00  0.0000E+00  0.0000E+00  3.5481E-05
             6.0449E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -583.654044444079        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2268
 NPARAMETR:  4.6558E-01  1.3625E-02  2.4765E-02  2.9489E-01  7.1756E+00  1.5834E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.3253E-02
             1.3476E+01
 PARAMETER: -6.6446E-01 -4.1959E+00 -3.5983E+00 -1.1212E+00  2.0707E+00  5.5961E-01 -4.3533E+01 -7.7600E+00 -6.0794E+00 -4.2235E+00
             2.7009E+00
 GRADIENT:  -8.8112E-01  9.5922E-01  9.5937E+00 -1.4098E+01 -8.1825E-02  7.2962E-01  0.0000E+00  0.0000E+00  0.0000E+00  1.6272E-04
             4.9539E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -583.873307445969        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2445
 NPARAMETR:  4.6164E-01  1.0000E-02  2.3873E-02  2.9279E-01  6.5041E+00  1.5874E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.0000E-02
             1.3353E+01
 PARAMETER: -6.7312E-01 -4.5047E+00 -3.6271E+00 -1.1331E+00  1.9529E+00  5.6279E-01 -4.3533E+01 -8.4745E+00 -6.0794E+00 -1.7598E+01
             2.6982E+00
 GRADIENT:  -6.5542E-02  3.7269E-02  3.9933E+00 -4.1215E+00 -2.7780E-01  1.1905E-01  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             2.4154E+00

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2445
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.3064E-03  5.9259E-06  7.9806E-05 -1.8258E-04  7.4936E-06
 SE:             2.8628E-02  4.9551E-06  2.4155E-04  2.9844E-04  3.9770E-05
 N:                     100         100         100         100         100

 P VAL.:         9.3579E-01  2.3173E-01  7.4111E-01  5.4069E-01  8.5055E-01

 ETASHRINKSD(%)  4.0925E+00  9.9983E+01  9.9191E+01  9.9000E+01  9.9867E+01
 ETASHRINKVR(%)  8.0175E+00  1.0000E+02  9.9993E+01  9.9990E+01  1.0000E+02
 EBVSHRINKSD(%)  4.2212E+00  9.9977E+01  9.9134E+01  9.8912E+01  9.9826E+01
 EBVSHRINKVR(%)  8.2642E+00  1.0000E+02  9.9993E+01  9.9988E+01  1.0000E+02
 RELATIVEINF(%)  4.6917E+00  6.4401E-07  3.4612E-05  5.3926E-05  4.8718E-05
 EPSSHRINKSD(%)  5.6913E+00
 EPSSHRINKVR(%)  1.1059E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -583.87330744596920     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       151.27751911776897     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    32.70
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:     6.63
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -583.873       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.62E-01  1.00E-02  2.41E-02  2.91E-01  6.38E+00  1.59E+00  1.00E-02  1.00E-02  1.00E-02  1.00E-02  1.34E+01
 


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
+        7.72E+01
 
 TH 2
+       -1.65E+02  3.52E+02
 
 TH 3
+       -6.78E+03  1.45E+04  5.97E+05
 
 TH 4
+        6.82E+02 -1.46E+03 -6.00E+04  6.03E+03
 
 TH 5
+        3.44E+00 -7.36E+00 -3.03E+02  3.04E+01  1.54E-01
 
 TH 6
+       -4.42E+00  9.45E+00  3.89E+02 -3.91E+01 -1.97E-01  2.54E-01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.24E+00  4.78E+00  1.97E+02 -1.98E+01 -9.99E-02  1.28E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  6.49E-02
 
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
+        1.87E+03
 
 TH 2
+       -4.89E+02  1.19E+04
 
 TH 3
+       -9.77E+03  2.04E+04  8.50E+05
 
 TH 4
+       -3.99E+02 -2.20E+03 -8.53E+04  9.76E+03
 
 TH 5
+        5.83E+00 -2.87E+01 -4.31E+02  4.52E+01  6.89E-01
 
 TH 6
+        1.89E+00  3.53E+01  5.51E+02 -7.97E+01  1.54E-01  6.50E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.99E+01  1.24E+01  2.81E+02 -1.72E+01 -1.81E-01  9.53E-01  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.01E+00
 
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
+        1.89E+03
 
 TH 2
+       -3.51E+02  8.48E+02
 
 TH 3
+       -1.32E+04  2.54E+04  1.21E+06
 
 TH 4
+       -2.11E+02 -2.53E+03 -1.17E+05  1.26E+04
 
 TH 5
+        8.00E+00 -1.77E+01 -6.92E+02  6.70E+01  4.64E-01
 
 TH 6
+        2.00E+02 -2.77E+01 -2.01E+03  1.93E+01  1.71E+00  8.08E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -3.53E+01  1.62E+01  5.33E+02 -2.87E+01 -3.64E-01 -2.22E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  1.32E+01
 
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
 #CPUT: Total CPU Time in Seconds,       39.384
Stop Time:
Wed Sep 29 20:24:39 CDT 2021
