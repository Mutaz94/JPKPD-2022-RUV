Thu Sep 30 09:02:26 CDT 2021
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
$DATA ../../../../data/spa2/D/dat40.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      700
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

 TOT. NO. OF OBS RECS:      600
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
 RAW OUTPUT FILE (FILE): m40.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   17620.2570256976        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.8628E+02  3.0824E+02 -1.7960E+01  3.4456E+02  2.0860E+02 -1.8362E+03 -8.5002E+02 -3.4075E+01 -1.1613E+03 -5.3705E+02
            -3.5388E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -677.657020546387        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.4752E+00  1.2056E+00  1.0001E+00  1.2743E+00  1.0646E+00  2.0432E+00  1.8931E+00  9.9897E-01  1.5767E+00  1.1670E+00
             1.3884E+01
 PARAMETER:  4.8882E-01  2.8697E-01  1.0006E-01  3.4240E-01  1.6260E-01  8.1449E-01  7.3821E-01  9.8972E-02  5.5535E-01  2.5448E-01
             2.7307E+00
 GRADIENT:   5.1182E+01 -6.5130E+01 -2.6076E+01 -1.7011E+01  3.6845E+01  6.1857E+01 -2.8297E+01  3.3402E+00  1.2473E+01  1.6364E+01
             3.3657E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -779.483430630559        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  1.3572E+00  1.6539E+00  8.9650E+00  1.3434E+00  4.2488E+00  2.4525E+00  5.8780E+00  9.4667E-01  1.3265E+00  3.2791E+00
             1.2078E+01
 PARAMETER:  4.0545E-01  6.0316E-01  2.2933E+00  3.9522E-01  1.5466E+00  9.9711E-01  1.8712E+00  4.5195E-02  3.8255E-01  1.2876E+00
             2.5914E+00
 GRADIENT:   3.1786E+01  1.1575E+01  4.6568E-01 -3.9590E+01 -3.0828E+01  8.0245E+01  6.9835E+01  3.2907E-02  1.8744E+01  5.2461E+01
             2.9033E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -808.889257836815        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.2343E+00  8.4363E-01  1.1598E+01  1.7099E+00  2.5822E+00  2.1262E+00  5.7692E+00  2.8669E+00  1.0834E+00  2.0190E+00
             1.1625E+01
 PARAMETER:  3.1051E-01 -7.0036E-02  2.5509E+00  6.3641E-01  1.0487E+00  8.5435E-01  1.8525E+00  1.1532E+00  1.8013E-01  8.0260E-01
             2.5531E+00
 GRADIENT:  -1.0585E+01  2.8572E-01 -1.6201E+00  3.9292E+01 -1.9511E+00  4.3759E+01  2.1516E+01  5.1694E-01 -4.5210E+00  2.5266E+01
             2.8017E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -870.663281241869        NO. OF FUNC. EVALS.: 125
 CUMULATIVE NO. OF FUNC. EVALS.:      356
 NPARAMETR:  1.2165E+00  8.3783E-01  3.5402E+02  1.7072E+00  1.9480E+00  2.0662E+00  5.9147E+00  2.2434E-02  1.4230E+00  4.7614E-02
             9.1618E+00
 PARAMETER:  2.9601E-01 -7.6940E-02  5.9694E+00  6.3483E-01  7.6682E-01  8.2573E-01  1.8774E+00 -3.6972E+00  4.5278E-01 -2.9446E+00
             2.3150E+00
 GRADIENT:   3.1483E+01  1.1874E+01  7.0805E-02  8.0067E+01 -2.6756E+01  3.7483E+01  3.6495E+01  5.0411E-08 -8.2169E+00  2.0978E-02
             8.3712E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -874.739918089354        NO. OF FUNC. EVALS.:  92
 CUMULATIVE NO. OF FUNC. EVALS.:      448
 NPARAMETR:  1.2152E+00  8.3810E-01  5.9518E+02  1.7087E+00  2.1410E+00  2.0551E+00  5.9387E+00  1.0745E-02  1.5829E+00  2.7399E-02
             8.5636E+00
 PARAMETER:  2.9488E-01 -7.6613E-02  6.4889E+00  6.3575E-01  8.6125E-01  8.2033E-01  1.8815E+00 -4.4333E+00  5.5924E-01 -3.4973E+00
             2.2475E+00
 GRADIENT:   1.7665E+01  1.1182E+01  1.7641E-04  4.5736E+01 -4.9751E+00  9.9427E+00 -1.6212E+01 -5.1211E-09 -3.2797E+00  6.7161E-03
            -1.0764E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -878.706638067004        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      584
 NPARAMETR:  1.1564E+00  6.6743E-01  1.2376E+03  1.5000E+00  2.1561E+00  1.9583E+00  6.1901E+00  1.1524E-02  1.5463E+00  1.8913E-02
             8.5458E+00
 PARAMETER:  2.4530E-01 -3.0433E-01  7.2209E+00  5.0549E-01  8.6828E-01  7.7209E-01  1.9230E+00 -4.3633E+00  5.3589E-01 -3.8679E+00
             2.2454E+00
 GRADIENT:  -2.4004E+00 -3.7997E+00 -1.9566E-03 -1.6621E+01  7.5116E+00 -5.3213E+00 -4.8798E+00 -1.4943E-08  1.3860E+01  3.4045E-03
             1.5896E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -881.061951986683        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      762
 NPARAMETR:  1.1587E+00  7.1292E-01  5.6253E+01  1.4380E+00  2.1009E+00  1.9817E+00  6.4066E+00  1.1514E-02  1.1462E+00  1.0000E-02
             8.5016E+00
 PARAMETER:  2.4730E-01 -2.3839E-01  4.1299E+00  4.6325E-01  8.4238E-01  7.8393E-01  1.9573E+00 -4.3642E+00  2.3644E-01 -4.5716E+00
             2.2403E+00
 GRADIENT:  -6.1482E-01 -8.5587E-01  2.4119E-02  9.5969E-01  2.4814E+00 -5.7529E-01  1.5085E-01  2.8959E-07 -1.4733E-01  0.0000E+00
            -2.7358E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -881.158419692571        NO. OF FUNC. EVALS.: 147
 CUMULATIVE NO. OF FUNC. EVALS.:      909
 NPARAMETR:  1.1587E+00  7.4368E-01  3.6044E+01  1.4307E+00  2.0463E+00  1.9858E+00  6.5978E+00  1.1405E-02  1.1484E+00  1.0000E-02
             8.4977E+00
 PARAMETER:  2.4727E-01 -1.9614E-01  3.6847E+00  4.5813E-01  8.1601E-01  7.8601E-01  1.9867E+00 -4.3737E+00  2.3834E-01 -4.5755E+00
             2.2398E+00
 GRADIENT:  -1.7070E-01  1.8448E+00  9.6287E-02 -3.1018E-01 -8.0126E-01 -3.2356E-01  8.2333E+00  7.2861E-07  8.4182E-01  0.0000E+00
            -5.5545E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -881.239291233035        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1089
 NPARAMETR:  1.1600E+00  7.1587E-01  2.1069E+01  1.4265E+00  2.0153E+00  1.9917E+00  6.5337E+00  1.0728E-02  1.1200E+00  1.0000E-02
             8.5036E+00
 PARAMETER:  2.4841E-01 -2.3426E-01  3.1478E+00  4.5525E-01  8.0075E-01  7.8899E-01  1.9770E+00 -4.4349E+00  2.1333E-01 -4.5755E+00
             2.2405E+00
 GRADIENT:   4.7943E-01  9.2310E-02  4.4308E-02  1.0261E+00  9.8905E-02  9.6255E-01  4.4848E+00  1.9715E-06 -4.5206E-01  0.0000E+00
            -1.7024E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -881.259470166287        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1274             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1592E+00  7.0360E-01  1.8536E+01  1.4280E+00  2.0045E+00  1.9899E+00  6.6118E+00  1.1553E-02  1.1275E+00  1.0000E-02
             8.5021E+00
 PARAMETER:  2.4772E-01 -2.5154E-01  3.0197E+00  4.5631E-01  7.9540E-01  7.8811E-01  1.9889E+00 -4.3608E+00  2.2000E-01 -4.5755E+00
             2.2403E+00
 GRADIENT:   2.0348E+01  1.7403E+00  1.1382E-03  1.4513E+01  2.0696E+00  2.2544E+01  9.0759E+01  3.2305E-06  9.3500E-01  0.0000E+00
             2.7252E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -881.262604620462        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1456
 NPARAMETR:  1.1588E+00  6.9904E-01  1.8333E+01  1.4299E+00  1.9996E+00  1.9902E+00  6.6291E+00  1.0000E-02  1.1287E+00  1.0000E-02
             8.5011E+00
 PARAMETER:  2.4740E-01 -2.5805E-01  3.0087E+00  4.5763E-01  7.9296E-01  7.8822E-01  1.9915E+00 -4.5625E+00  2.2106E-01 -4.5755E+00
             2.2402E+00
 GRADIENT:   2.2392E-01 -3.8499E-02 -4.3939E-03 -1.2770E+00  4.3377E-01  6.3848E-01  6.6573E+00  0.0000E+00  1.8535E-01  0.0000E+00
             2.4966E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -881.266302015282        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     1645             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1585E+00  6.9111E-01  1.8024E+01  1.4355E+00  1.9948E+00  1.9897E+00  6.6555E+00  1.1340E-02  1.1323E+00  1.0000E-02
             8.5004E+00
 PARAMETER:  2.4712E-01 -2.6945E-01  2.9917E+00  4.6149E-01  7.9053E-01  7.8801E-01  1.9954E+00 -4.3794E+00  2.2429E-01 -4.5755E+00
             2.2401E+00
 GRADIENT:   2.0098E+01  1.8170E+00  1.5725E-02  1.5410E+01  1.4370E+00  2.2583E+01  9.1915E+01  3.3547E-06  7.5345E-01  0.0000E+00
             2.6749E+01

0ITERATION NO.:   62    OBJECTIVE VALUE:  -881.266302015282        NO. OF FUNC. EVALS.:  63
 CUMULATIVE NO. OF FUNC. EVALS.:     1708
 NPARAMETR:  1.1589E+00  6.9298E-01  1.8143E+01  1.4396E+00  1.9917E+00  1.9905E+00  6.6649E+00  1.1087E-02  1.1314E+00  1.0000E-02
             8.4991E+00
 PARAMETER:  2.4712E-01 -2.6945E-01  2.9917E+00  4.6149E-01  7.9053E-01  7.8801E-01  1.9954E+00 -4.3794E+00  2.2429E-01 -4.5755E+00
             2.2401E+00
 GRADIENT:  -7.7671E-02 -4.0386E-02 -2.6734E-03 -5.9751E-01  1.2989E-01 -4.7523E-02 -1.1769E-01  4.7808E-06  2.2483E-02  0.0000E+00
             8.0619E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1708
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0569E-02  4.4326E-02  2.4446E-06 -7.8778E-02  1.1770E-05
 SE:             2.8454E-02  2.3038E-02  3.2975E-06  1.4148E-02  6.8900E-05
 N:                     100         100         100         100         100

 P VAL.:         7.1032E-01  5.4344E-02  4.5848E-01  2.5816E-08  8.6436E-01

 ETASHRINKSD(%)  4.6763E+00  2.2821E+01  9.9989E+01  5.2602E+01  9.9769E+01
 ETASHRINKVR(%)  9.1340E+00  4.0434E+01  1.0000E+02  7.7535E+01  9.9999E+01
 EBVSHRINKSD(%)  5.2295E+00  1.8669E+01  9.9982E+01  5.2627E+01  9.9685E+01
 EBVSHRINKVR(%)  1.0186E+01  3.3853E+01  1.0000E+02  7.7558E+01  9.9999E+01
 RELATIVEINF(%)  8.8598E+01  3.3785E+01  4.4746E-07  1.0797E+01  1.3336E-04
 EPSSHRINKSD(%)  1.0075E+01
 EPSSHRINKVR(%)  1.9135E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -881.26630201528224     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       221.45993783032486     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.16
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0INVERSE COVARIANCE MATRIX SET TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S
 Elapsed covariance  time in seconds:    11.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -881.266       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.16E+00  6.91E-01  1.80E+01  1.44E+00  1.99E+00  1.99E+00  6.66E+00  1.13E-02  1.13E+00  1.00E-02  8.50E+00
 


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
 ********************                            STANDARD ERROR OF ESTIMATE                          ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         8.83E-02  5.00E-01  5.71E+01  4.85E-01  4.14E-01  1.79E-01  2.17E+00  9.85E-01  5.33E-01  0.00E+00  1.75E+00
 


 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS  ****


         EPS1     
 
 EPS1
+       .........
 
1


 OMEGA - CORR MATRIX FOR RANDOM EFFECTS - ETAS  *******


         ETA1      ETA2      ETA3      ETA4      ETA5     
 
 ETA1
+       .........
 
 ETA2
+       ......... .........
 
 ETA3
+       ......... ......... .........
 
 ETA4
+       ......... ......... ......... .........
 
 ETA5
+       ......... ......... ......... ......... .........
 


 SIGMA - CORR MATRIX FOR RANDOM EFFECTS - EPSILONS  ***


         EPS1     
 
 EPS1
+       .........
 
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                          COVARIANCE MATRIX OF ESTIMATE                         ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        7.80E-03
 
 TH 2
+        5.04E-04  2.50E-01
 
 TH 3
+        7.69E-01 -2.49E+01  3.26E+03
 
 TH 4
+        1.37E-02 -2.22E-01  2.54E+01  2.35E-01
 
 TH 5
+        6.85E-03 -1.67E-01  2.21E+01  1.74E-01  1.71E-01
 
 TH 6
+        1.03E-03  4.03E-02 -3.07E+00 -3.64E-02 -2.47E-02  3.20E-02
 
 TH 7
+        3.41E-02 -9.84E-01  1.13E+02  9.50E-01  8.16E-01 -1.28E-01  4.71E+00
 
 TH 8
+       -1.09E-02 -4.35E-01  3.40E+01  3.33E-01  2.38E-01 -8.80E-02  1.66E+00  9.71E-01
 
 TH 9
+        5.25E-03 -2.10E-01  2.63E+01  2.17E-01  1.60E-01 -1.64E-02  8.06E-01  2.19E-01  2.84E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        8.26E-02 -2.05E-01  1.70E+01  3.45E-01  2.07E-01 -1.24E-01  1.31E+00  5.58E-01 -5.54E-03  0.00E+00  3.07E+00
 
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
 ********************                          CORRELATION MATRIX OF ESTIMATE                        ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        8.83E-02
 
 TH 2
+        1.14E-02  5.00E-01
 
 TH 3
+        1.52E-01 -8.74E-01  5.71E+01
 
 TH 4
+        3.19E-01 -9.15E-01  9.18E-01  4.85E-01
 
 TH 5
+        1.87E-01 -8.07E-01  9.34E-01  8.69E-01  4.14E-01
 
 TH 6
+        6.50E-02  4.52E-01 -3.01E-01 -4.19E-01 -3.34E-01  1.79E-01
 
 TH 7
+        1.78E-01 -9.08E-01  9.11E-01  9.03E-01  9.08E-01 -3.29E-01  2.17E+00
 
 TH 8
+       -1.26E-01 -8.83E-01  6.04E-01  6.96E-01  5.82E-01 -4.99E-01  7.78E-01  9.85E-01
 
 TH 9
+        1.12E-01 -7.88E-01  8.65E-01  8.41E-01  7.23E-01 -1.72E-01  6.97E-01  4.18E-01  5.33E-01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+        5.34E-01 -2.34E-01  1.70E-01  4.07E-01  2.85E-01 -3.97E-01  3.44E-01  3.24E-01 -5.93E-03  0.00E+00  1.75E+00
 
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
 ********************                      INVERSE COVARIANCE MATRIX OF ESTIMATE                     ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 

            TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11      OM11  
             OM12      OM13      OM14      OM15      OM22      OM23      OM24      OM25      OM33      OM34      OM35      OM44  
            OM45      OM55      SG11  
 
 TH 1
+        4.62E+02
 
 TH 2
+       -1.05E+02  5.68E+01
 
 TH 3
+        6.59E-02 -6.76E-02  1.84E-04
 
 TH 4
+       -2.27E+02  7.71E+01  1.32E-02  2.26E+02
 
 TH 5
+        8.24E+00  2.86E+00 -2.83E-02 -2.26E+01  6.28E+00
 
 TH 6
+       -8.60E+01  1.08E+01  4.13E-02  5.65E+01 -1.06E+01  6.04E+01
 
 TH 7
+        1.61E+00  3.41E+00 -1.39E-02 -6.84E+00  2.66E+00 -3.14E+00  1.29E+00
 
 TH 8
+       -3.32E-01  2.43E-01 -4.60E-04  1.41E-01  5.46E-02 -2.71E-02  3.24E-02  1.34E-03
 
 TH 9
+        6.17E+01 -1.36E+01 -2.97E-02 -7.13E+01  1.09E+01 -2.11E+01  4.18E+00  2.40E-02  2.63E+01
 
 TH10
+        1.70E-41 -2.63E-41  6.82E-44  1.44E-42 -1.00E-41  1.81E-41 -5.04E-42 -1.73E-43 -9.94E-42  1.55E-75
 
 TH11
+        1.10E+00 -2.87E+00  8.93E-04 -7.60E+00  6.59E-01  8.78E-01  1.83E-01 -6.54E-03  2.29E+00  4.66E-35  8.74E-01
 
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
+        1.88E+02
 
 TH 2
+        5.35E+01  3.65E+01
 
 TH 3
+        5.35E-02  2.15E-02  4.46E-04
 
 TH 4
+        8.32E+01  3.42E+01  5.43E-02  1.33E+02
 
 TH 5
+       -2.05E+01 -2.54E+00 -6.72E-02 -1.96E+01  1.45E+01
 
 TH 6
+        3.14E+01  1.04E+01  4.74E-03 -1.10E+01 -2.40E+00  3.98E+01
 
 TH 7
+       -1.20E+00  5.47E+00 -3.87E-03 -7.23E+00  3.61E+00  2.29E+00  3.85E+00
 
 TH 8
+       -1.24E-03 -2.91E-04  4.30E-07 -7.07E-04 -6.44E-05  1.53E-04  2.67E-05  4.31E-07
 
 TH 9
+       -9.10E+00 -2.27E+00  1.70E-04 -3.09E+01  2.99E+00  5.59E+00  1.45E+00  3.83E-05  1.91E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.78E+01 -1.97E+00 -2.56E-02 -4.02E+01  9.33E+00  1.13E+01  6.79E+00 -1.51E-03  2.81E+00  0.00E+00  2.02E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.04
 #CPUT: Total CPU Time in Seconds,       50.771
Stop Time:
Thu Sep 30 09:03:18 CDT 2021
