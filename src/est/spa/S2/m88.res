Wed Sep 29 17:45:19 CDT 2021
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
$DATA ../../../../data/spa/S2/dat88.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1648.32837231784        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1692E+02 -4.7694E+01 -4.5463E+01  1.3765E+01  8.1328E+01  4.8393E+01 -7.2020E+00 -2.7380E+00  6.1741E+00  2.4012E+01
            -5.3075E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1657.55669123286        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0105E+00  1.0910E+00  1.0797E+00  9.9760E-01  1.0215E+00  9.6821E-01  1.0500E+00  1.0408E+00  9.9226E-01  8.1316E-01
             1.2768E+00
 PARAMETER:  1.1048E-01  1.8708E-01  1.7671E-01  9.7595E-02  1.2123E-01  6.7697E-02  1.4874E-01  1.4000E-01  9.2230E-02 -1.0682E-01
             3.4437E-01
 GRADIENT:   1.0066E+01 -1.5918E+01 -1.2238E+01  7.3828E+00  3.5987E+01 -4.6576E+00 -4.7207E+00 -8.2324E+00  2.2240E+00  8.1539E+00
             4.4822E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1661.27655974381        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0076E+00  1.1493E+00  9.1413E-01  9.5766E-01  9.1939E-01  9.9494E-01  1.1747E+00  1.5814E+00  9.6874E-01  4.7034E-01
             1.2019E+00
 PARAMETER:  1.0755E-01  2.3912E-01  1.0221E-02  5.6733E-02  1.5958E-02  9.4928E-02  2.6102E-01  5.5829E-01  6.8240E-02 -6.5430E-01
             2.8386E-01
 GRADIENT:   2.6434E+00  9.2730E+00 -3.1464E+00  1.3704E+01 -2.5867E+01  5.1557E+00  5.8845E+00  1.2805E+01  7.2014E+00  2.0658E+00
             2.3491E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1663.65566499292        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0056E+00  1.1969E+00  9.5932E-01  9.2173E-01  9.7631E-01  9.8092E-01  1.0942E+00  1.3584E+00  9.8528E-01  5.9379E-01
             1.1426E+00
 PARAMETER:  1.0559E-01  2.7972E-01  5.8466E-02  1.8496E-02  7.6027E-02  8.0733E-02  1.9004E-01  4.0634E-01  8.5176E-02 -4.2124E-01
             2.3329E-01
 GRADIENT:  -6.2515E-01  2.7333E-01  5.3624E-01  2.9219E-01  1.9796E+00 -4.5535E-01  2.6015E-01 -2.8854E-01 -6.3390E-01 -1.2118E+00
            -5.5399E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1663.70972863761        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0066E+00  1.3043E+00  8.8015E-01  8.5421E-01  9.9093E-01  9.8282E-01  1.0173E+00  1.3607E+00  1.0478E+00  6.1652E-01
             1.1428E+00
 PARAMETER:  1.0658E-01  3.6564E-01 -2.7663E-02 -5.7580E-02  9.0891E-02  8.2674E-02  1.1720E-01  4.0799E-01  1.4671E-01 -3.8366E-01
             2.3344E-01
 GRADIENT:   2.6173E-01  4.2913E+00  9.1230E-01  3.8984E+00 -3.2725E+00  3.8297E-02 -3.9653E-01  2.3012E-01  1.7600E-01  1.3484E-01
             9.6238E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1663.76636891168        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      880
 NPARAMETR:  1.0072E+00  1.4755E+00  6.8221E-01  7.3699E-01  9.8555E-01  9.8290E-01  9.4827E-01  1.2684E+00  1.1333E+00  5.7442E-01
             1.1456E+00
 PARAMETER:  1.0717E-01  4.8897E-01 -2.8241E-01 -2.0518E-01  8.5440E-02  8.2754E-02  4.6888E-02  3.3774E-01  2.2512E-01 -4.5439E-01
             2.3589E-01
 GRADIENT:  -7.3255E-01  8.1210E+00  1.1926E+00  5.1841E+00 -2.7163E+00 -3.7683E-01  1.5172E-01  8.7962E-02 -1.0486E-01 -5.6166E-01
             1.4956E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1663.83995159433        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1056
 NPARAMETR:  1.0079E+00  1.6271E+00  4.8865E-01  6.2766E-01  9.6344E-01  9.8434E-01  8.9253E-01  1.0703E+00  1.2227E+00  5.3262E-01
             1.1389E+00
 PARAMETER:  1.0789E-01  5.8682E-01 -6.1611E-01 -3.6576E-01  6.2752E-02  8.4215E-02 -1.3690E-02  1.6793E-01  3.0103E-01 -5.2995E-01
             2.3007E-01
 GRADIENT:   9.7952E-02  1.3912E+01  3.2563E+00  3.4419E+00 -9.7423E+00 -1.0247E-01  1.1815E-01  1.3350E-02 -3.9476E-01 -2.0390E-01
            -2.5336E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1664.33777930509        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1235
 NPARAMETR:  1.0067E+00  1.8944E+00  1.8919E-01  4.2117E-01  8.9693E-01  9.8714E-01  7.9131E-01  5.2792E-01  1.4273E+00  4.3883E-01
             1.1026E+00
 PARAMETER:  1.0665E-01  7.3893E-01 -1.5650E+00 -7.6472E-01 -8.7821E-03  8.7054E-02 -1.3407E-01 -5.3881E-01  4.5575E-01 -7.2364E-01
             1.9769E-01
 GRADIENT:   9.3488E+00  5.7550E+01  1.9443E+01 -1.7071E+01 -6.7739E+01  1.6338E+00 -1.2132E+00 -5.6464E-01 -2.3633E+00  2.7288E+00
            -2.1740E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1666.21056737792        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1413
 NPARAMETR:  1.0006E+00  2.1569E+00  9.4265E-02  2.6097E-01  1.0283E+00  9.9076E-01  6.8970E-01  2.4922E-01  1.9304E+00  6.1133E-01
             1.0517E+00
 PARAMETER:  1.0060E-01  8.6867E-01 -2.2616E+00 -1.2434E+00  1.2788E-01  9.0712E-02 -2.7150E-01 -1.2894E+00  7.5775E-01 -3.9212E-01
             1.5042E-01
 GRADIENT:   1.5119E+00  3.6901E+01  2.4676E+00  4.8590E+00 -2.1317E+01  1.7550E+00 -1.3932E+01 -5.6761E-01 -6.4452E+00  2.0093E+00
            -7.8665E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1666.63876248372        NO. OF FUNC. EVALS.: 140
 CUMULATIVE NO. OF FUNC. EVALS.:     1553
 NPARAMETR:  1.0008E+00  2.1613E+00  9.1930E-02  2.5148E-01  1.0397E+00  9.6074E-01  7.3212E-01  2.3747E-01  1.9692E+00  6.2072E-01
             1.0517E+00
 PARAMETER:  1.0081E-01  8.7071E-01 -2.2867E+00 -1.2804E+00  1.3893E-01  5.9948E-02 -2.1181E-01 -1.3377E+00  7.7761E-01 -3.7688E-01
             1.5044E-01
 GRADIENT:   3.9438E+02  1.2449E+03  9.5527E+00  7.0128E+01 -3.9316E+00  2.7586E+01  2.0416E+01 -4.3636E-01  6.6106E+00  3.8915E+00
            -5.6716E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1667.08683470221        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1713
 NPARAMETR:  9.9704E-01  2.1434E+00  9.1565E-02  2.5027E-01  1.0407E+00  9.8439E-01  7.1832E-01  2.5457E-01  1.9759E+00  6.1697E-01
             1.0543E+00
 PARAMETER:  9.7036E-02  8.6238E-01 -2.2907E+00 -1.2852E+00  1.3985E-01  8.4270E-02 -2.3084E-01 -1.2682E+00  7.8105E-01 -3.8294E-01
             1.5289E-01
 GRADIENT:   3.8121E+02  1.1955E+03  8.5632E+00  6.7263E+01  5.4373E+00  3.7678E+01  1.5886E+01 -4.9028E-01  6.6370E+00  3.6574E+00
            -3.9286E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1667.14702642498        NO. OF FUNC. EVALS.: 170
 CUMULATIVE NO. OF FUNC. EVALS.:     1883
 NPARAMETR:  9.9895E-01  2.1344E+00  9.2610E-02  2.4867E-01  1.0447E+00  9.8647E-01  7.1774E-01  2.5866E-01  1.9708E+00  5.7081E-01
             1.0811E+00
 PARAMETER:  9.8945E-02  8.5820E-01 -2.2794E+00 -1.2916E+00  1.4371E-01  8.6382E-02 -2.3165E-01 -1.2523E+00  7.7845E-01 -4.6069E-01
             1.7801E-01
 GRADIENT:  -2.4896E+00 -3.7241E+01 -5.4005E-01 -3.1117E+00  1.0342E+01  2.4369E-01 -1.1373E+00 -4.7659E-01 -6.7895E+00 -1.2227E+00
             6.5750E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1667.27820647408        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2040
 NPARAMETR:  9.9985E-01  2.1352E+00  9.2506E-02  2.4931E-01  1.0404E+00  9.8609E-01  7.2032E-01  2.8085E-01  1.9833E+00  5.8005E-01
             1.0760E+00
 PARAMETER:  9.9854E-02  8.5854E-01 -2.2805E+00 -1.2891E+00  1.3957E-01  8.5996E-02 -2.2806E-01 -1.1699E+00  7.8479E-01 -4.4463E-01
             1.7322E-01
 GRADIENT:   3.7206E+02  1.1284E+03  7.8451E+00  6.2805E+01  1.0332E+01  3.6717E+01  1.5176E+01 -5.3475E-01  6.4144E+00  1.1967E+00
             1.9915E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1667.28938396290        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:     2132
 NPARAMETR:  9.9983E-01  2.1351E+00  9.2503E-02  2.4936E-01  1.0404E+00  9.8607E-01  7.2030E-01  2.8452E-01  1.9846E+00  5.8000E-01
             1.0759E+00
 PARAMETER:  9.9828E-02  8.5852E-01 -2.2805E+00 -1.2889E+00  1.3956E-01  8.5977E-02 -2.2809E-01 -1.1569E+00  7.8542E-01 -4.4472E-01
             1.7317E-01
 GRADIENT:  -3.4629E+04  8.0584E+03 -2.9983E+03  5.3498E+03 -4.9624E+04 -9.8914E-02  2.5826E-01  5.9756E+03  8.7001E+03 -1.5577E+04
            -4.0008E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2132
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2748E-04 -2.0591E-02 -1.3255E-03  3.8952E-02 -4.2806E-02
 SE:             2.9816E-02  2.7886E-02  2.5517E-03  2.3136E-02  1.6883E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9391E-01  4.6027E-01  6.0343E-01  9.2249E-02  1.1232E-02

 ETASHRINKSD(%)  1.1192E-01  6.5789E+00  9.1452E+01  2.2493E+01  4.3438E+01
 ETASHRINKVR(%)  2.2371E-01  1.2725E+01  9.9269E+01  3.9926E+01  6.8008E+01
 EBVSHRINKSD(%)  5.0836E-01  7.6786E+00  9.3220E+01  2.2730E+01  4.3685E+01
 EBVSHRINKVR(%)  1.0141E+00  1.4768E+01  9.9540E+01  4.0293E+01  6.8286E+01
 RELATIVEINF(%)  9.8392E+01  2.7496E+01  2.3647E-01  1.2692E+01  8.3926E+00
 EPSSHRINKSD(%)  4.4598E+01
 EPSSHRINKVR(%)  6.9306E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1667.2893839629007     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -932.13855739916255     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    27.92
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.34
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1667.289       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.14E+00  9.25E-02  2.49E-01  1.04E+00  9.86E-01  7.20E-01  2.85E-01  1.98E+00  5.80E-01  1.08E+00
 


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
+        1.73E+07
 
 TH 2
+        1.82E+02  5.19E+04
 
 TH 3
+       -1.76E+03 -4.42E+05  3.80E+06
 
 TH 4
+        1.11E+03 -2.02E+03  1.70E+04  1.67E+06
 
 TH 5
+       -2.59E+03 -1.15E+02 -2.27E+03  1.72E+03  8.21E+06
 
 TH 6
+       -5.13E+02  2.64E+01 -2.40E+02  1.50E+02 -3.52E+02  2.01E+02
 
 TH 7
+       -1.71E+02  1.69E+01 -1.30E+02  5.63E+01 -1.16E+02 -1.07E+07  6.42E+06
 
 TH 8
+       -5.25E+06 -1.01E+03  8.52E+03 -5.73E+03  7.78E+02  1.58E+02  5.28E+01  1.59E+06
 
 TH 9
+        2.31E+02 -1.58E+03  1.35E+04 -2.51E+03  1.46E+02  3.10E+01  1.87E+01 -1.16E+03  6.95E+04
 
 TH10
+        1.95E+03 -1.20E+02  1.02E+03 -6.33E+02  1.23E+03 -1.97E+02 -4.92E+01 -5.91E+02 -1.16E+02  2.60E+06
 
 TH11
+        5.91E+03 -3.39E+02  2.67E+03 -1.82E+03 -1.43E+03 -2.68E+02 -8.31E+01 -1.79E+03 -3.65E+02  1.08E+03  4.99E+06
 
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
 #CPUT: Total CPU Time in Seconds,       34.320
Stop Time:
Wed Sep 29 17:45:55 CDT 2021
