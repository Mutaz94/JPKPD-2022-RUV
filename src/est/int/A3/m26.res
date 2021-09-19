Sat Sep 18 01:29:43 CDT 2021
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
$DATA ../../../../data/int/A3/dat26.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:     1000
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E21.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      900
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1526.56504811542        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.1852E+02  1.3624E+02  1.7993E+02 -8.9170E+01  1.0176E+02  4.5448E+01 -1.4524E+02 -1.7314E+02 -6.2314E+01 -1.6728E+02
            -4.0563E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2835.12049807367        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:       87
 NPARAMETR:  9.7465E-01  9.3731E-01  9.0256E-01  1.0719E+00  9.1227E-01  8.5725E-01  1.0249E+00  1.0247E+00  8.6487E-01  1.0988E+00
             2.0090E+00
 PARAMETER:  7.4320E-02  3.5257E-02 -2.5220E-03  1.6942E-01  8.1829E-03 -5.4024E-02  1.2461E-01  1.2437E-01 -4.5170E-02  1.9423E-01
             7.9762E-01
 GRADIENT:   1.2761E+01  2.2900E+01  2.7480E+01 -1.5385E+01 -4.7656E+00 -3.3813E+00 -1.3265E+01 -1.0307E+01 -9.7540E+00 -6.8835E+00
            -3.6524E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2840.84576800918        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.8224E-01  7.0725E-01  6.6285E-01  1.2244E+00  6.6579E-01  8.8887E-01  1.0957E+00  8.0539E-01  8.1447E-01  9.4075E-01
             1.9932E+00
 PARAMETER:  8.2082E-02 -2.4637E-01 -3.1120E-01  3.0242E-01 -3.0678E-01 -1.7805E-02  1.9141E-01 -1.1642E-01 -1.0522E-01  3.8918E-02
             7.8973E-01
 GRADIENT:   3.2181E+01  5.1384E+01  1.6971E+01  1.1159E+02  9.1783E+00  8.8252E+00 -1.8021E+01 -7.6503E+00 -3.1597E+01 -5.4205E+00
            -3.5265E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2877.64426895989        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.7960E-01  4.9495E-01  4.2480E-01  1.2948E+00  4.3428E-01  8.9039E-01  1.3563E+00  5.4563E-01  9.3281E-01  7.5667E-01
             2.3274E+00
 PARAMETER:  7.9386E-02 -6.0330E-01 -7.5614E-01  3.5835E-01 -7.3407E-01 -1.6093E-02  4.0475E-01 -5.0582E-01  3.0448E-02 -1.7883E-01
             9.4474E-01
 GRADIENT:   9.6025E+00  5.9781E+01  3.6209E+01  1.2526E+02 -3.6932E+01  1.0154E+01  5.7150E-01  4.7994E-01  4.2410E+00  2.7302E+00
             7.3008E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2885.32100596124        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.8436E-01  3.1545E-01  2.6192E-01  1.3623E+00  2.7460E-01  9.0398E-01  1.5717E+00  8.1153E-01  1.0511E+00  5.5225E-01
             2.2283E+00
 PARAMETER:  8.4240E-02 -1.0537E+00 -1.2397E+00  4.0915E-01 -1.1924E+00 -9.4697E-04  5.5215E-01 -1.0883E-01  1.4987E-01 -4.9376E-01
             9.0123E-01
 GRADIENT:   2.0052E+01  4.1678E+01  1.5056E+02  2.6420E+02 -1.8961E+02  1.2860E+01  1.2198E+01 -8.4834E+00 -3.8104E+01 -6.5516E+00
             1.1342E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2918.44516614476        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      374
 NPARAMETR:  9.7016E-01  1.7783E-01  1.2560E-01  1.1563E+00  1.7280E-01  8.5176E-01  1.5366E+00  3.7602E-01  1.5107E+00  6.5711E-01
             1.9458E+00
 PARAMETER:  6.9710E-02 -1.6269E+00 -1.9747E+00  2.4522E-01 -1.6556E+00 -6.0448E-02  5.2958E-01 -8.7811E-01  5.1256E-01 -3.1990E-01
             7.6570E-01
 GRADIENT:   1.9511E+00 -6.6234E+01  1.1078E+02  1.4405E+02 -8.4036E+01 -5.6353E+00 -7.2906E-01 -1.1847E+01 -2.6840E+01 -2.5204E+01
            -7.6708E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2930.48798204858        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      446
 NPARAMETR:  9.6609E-01  1.6655E-01  9.5209E-02  9.5180E-01  1.5165E-01  8.5309E-01  1.4861E+00  1.6035E-01  1.7356E+00  7.9333E-01
             1.9929E+00
 PARAMETER:  6.5497E-02 -1.6924E+00 -2.2517E+00  5.0597E-02 -1.7862E+00 -5.8889E-02  4.9614E-01 -1.7304E+00  6.5136E-01 -1.3151E-01
             7.8961E-01
 GRADIENT:  -1.0780E+01 -1.2092E+01 -8.4386E+00  5.0257E+01  2.0108E+01 -5.4466E+00 -3.3798E+00 -1.3497E+00  2.0260E+00  7.7066E+00
            -2.9559E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2930.82968124670        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.6780E-01  1.6639E-01  9.2991E-02  9.1119E-01  1.4937E-01  8.5885E-01  1.4953E+00  1.4842E-01  1.7359E+00  7.8843E-01
             2.0130E+00
 PARAMETER:  6.7272E-02 -1.6934E+00 -2.2753E+00  6.9920E-03 -1.8013E+00 -5.2159E-02  5.0235E-01 -1.8077E+00  6.5152E-01 -1.3771E-01
             7.9962E-01
 GRADIENT:  -7.0215E+00 -1.1882E+00 -7.9319E+00  2.8982E+01  1.0636E+01 -2.4893E+00 -1.1441E+00 -1.1679E+00  3.0224E+00  4.1579E+00
            -1.1313E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2930.82988876592        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      593
 NPARAMETR:  9.6791E-01  1.6632E-01  9.2930E-02  9.0947E-01  1.4929E-01  8.5901E-01  1.4957E+00  1.4836E-01  1.7350E+00  7.8811E-01
             2.0135E+00
 PARAMETER:  6.7382E-02 -1.6938E+00 -2.2759E+00  5.1117E-03 -1.8018E+00 -5.1979E-02  5.0257E-01 -1.8081E+00  6.5101E-01 -1.3812E-01
             7.9986E-01
 GRADIENT:  -6.7689E+00 -1.1440E+00 -7.5287E+00  2.8031E+01  1.0122E+01 -2.3993E+00 -1.0918E+00 -1.1674E+00  2.9233E+00  3.9757E+00
            -1.0912E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2930.83004245147        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      668
 NPARAMETR:  9.6800E-01  1.6626E-01  9.2876E-02  9.0796E-01  1.4922E-01  8.5914E-01  1.4959E+00  1.4852E-01  1.7342E+00  7.8781E-01
             2.0139E+00
 PARAMETER:  6.7481E-02 -1.6942E+00 -2.2765E+00  3.4464E-03 -1.8023E+00 -5.1822E-02  5.0276E-01 -1.8070E+00  6.5052E-01 -1.3850E-01
             8.0007E-01
 GRADIENT:  -6.5415E+00 -1.1215E+00 -7.1531E+00  2.7188E+01  9.6562E+00 -2.3205E+00 -1.0471E+00 -1.1701E+00  2.8267E+00  3.8101E+00
            -1.0572E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2930.83022701826        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      743
 NPARAMETR:  9.6810E-01  1.6619E-01  9.2820E-02  9.0637E-01  1.4915E-01  8.5928E-01  1.4962E+00  1.4894E-01  1.7333E+00  7.8748E-01
             2.0143E+00
 PARAMETER:  6.7585E-02 -1.6946E+00 -2.2771E+00  1.6971E-03 -1.8028E+00 -5.1656E-02  5.0296E-01 -1.8042E+00  6.5001E-01 -1.3892E-01
             8.0028E-01
 GRADIENT:  -6.3039E+00 -1.1000E+00 -6.7442E+00  2.6301E+01  9.1513E+00 -2.2377E+00 -9.9931E-01 -1.1763E+00  2.7229E+00  3.6315E+00
            -1.0213E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2930.83046882405        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      817
 NPARAMETR:  9.6822E-01  1.6611E-01  9.2756E-02  9.0454E-01  1.4906E-01  8.5945E-01  1.4966E+00  1.4974E-01  1.7323E+00  7.8706E-01
             2.0148E+00
 PARAMETER:  6.7708E-02 -1.6951E+00 -2.2778E+00 -3.3022E-04 -1.8034E+00 -5.1460E-02  5.0319E-01 -1.7988E+00  6.4942E-01 -1.3944E-01
             8.0052E-01
 GRADIENT:  -6.0211E+00 -1.0767E+00 -6.2516E+00  2.5270E+01  8.5458E+00 -2.1404E+00 -9.4230E-01 -1.1879E+00  2.6011E+00  3.4188E+00
            -9.7960E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2930.83076463554        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      891
 NPARAMETR:  9.6836E-01  1.6602E-01  9.2684E-02  9.0248E-01  1.4896E-01  8.5965E-01  1.4970E+00  1.5105E-01  1.7311E+00  7.8657E-01
             2.0153E+00
 PARAMETER:  6.7848E-02 -1.6956E+00 -2.2786E+00 -2.6052E-03 -1.8041E+00 -5.1234E-02  5.0345E-01 -1.7901E+00  6.4876E-01 -1.4007E-01
             8.0078E-01
 GRADIENT:  -5.6976E+00 -1.0524E+00 -5.6753E+00  2.4110E+01  7.8409E+00 -2.0296E+00 -8.7634E-01 -1.2066E+00  2.4627E+00  3.1728E+00
            -9.3247E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2930.83113419162        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      966
 NPARAMETR:  9.6851E-01  1.6592E-01  9.2603E-02  9.0018E-01  1.4885E-01  8.5987E-01  1.4974E+00  1.5303E-01  1.7298E+00  7.8598E-01
             2.0159E+00
 PARAMETER:  6.8008E-02 -1.6963E+00 -2.2794E+00 -5.1587E-03 -1.8048E+00 -5.0975E-02  5.0375E-01 -1.7771E+00  6.4802E-01 -1.4082E-01
             8.0107E-01
 GRADIENT:  -5.3309E+00 -1.0273E+00 -4.9984E+00  2.2805E+01  7.0163E+00 -1.9037E+00 -7.9983E-01 -1.2349E+00  2.3045E+00  2.8873E+00
            -8.7913E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2930.83166303189        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1041
 NPARAMETR:  9.6871E-01  1.6579E-01  9.2503E-02  8.9732E-01  1.4871E-01  8.6015E-01  1.4980E+00  1.5625E-01  1.7283E+00  7.8519E-01
             2.0166E+00
 PARAMETER:  6.8210E-02 -1.6970E+00 -2.2805E+00 -8.3373E-03 -1.8058E+00 -5.0644E-02  5.0413E-01 -1.7563E+00  6.4711E-01 -1.4183E-01
             8.0142E-01
 GRADIENT:  -4.8645E+00 -9.9933E-01 -4.1107E+00  2.1175E+01  5.9407E+00 -1.7445E+00 -7.0079E-01 -1.2809E+00  2.1037E+00  2.5179E+00
            -8.1207E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2930.83233877799        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1116
 NPARAMETR:  9.6894E-01  1.6564E-01  9.2386E-02  8.9397E-01  1.4854E-01  8.6050E-01  1.4987E+00  1.6111E-01  1.7264E+00  7.8417E-01
             2.0174E+00
 PARAMETER:  6.8452E-02 -1.6980E+00 -2.2818E+00 -1.2086E-02 -1.8069E+00 -5.0242E-02  5.0458E-01 -1.7257E+00  6.4605E-01 -1.4313E-01
             8.0182E-01
 GRADIENT:  -4.3036E+00 -9.7114E-01 -2.9989E+00  1.9244E+01  4.6016E+00 -1.5535E+00 -5.7858E-01 -1.3511E+00  1.8607E+00  2.0626E+00
            -7.3193E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2930.83325937230        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1191
 NPARAMETR:  9.6925E-01  1.6544E-01  9.2238E-02  8.8975E-01  1.4832E-01  8.6095E-01  1.4995E+00  1.6888E-01  1.7242E+00  7.8275E-01
             2.0184E+00
 PARAMETER:  6.8765E-02 -1.6992E+00 -2.2834E+00 -1.6820E-02 -1.8084E+00 -4.9717E-02  5.0516E-01 -1.6786E+00  6.4474E-01 -1.4494E-01
             8.0229E-01
 GRADIENT:  -3.5771E+00 -9.4199E-01 -1.4932E+00  1.6795E+01  2.7996E+00 -1.3069E+00 -4.1529E-01 -1.4647E+00  1.5431E+00  1.4573E+00
            -6.2882E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2930.83455581038        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:     1266
 NPARAMETR:  9.6966E-01  1.6516E-01  9.2041E-02  8.8418E-01  1.4802E-01  8.6158E-01  1.5007E+00  1.8191E-01  1.7213E+00  7.8064E-01
             2.0195E+00
 PARAMETER:  6.9192E-02 -1.7008E+00 -2.2855E+00 -2.3100E-02 -1.8104E+00 -4.8988E-02  5.0594E-01 -1.6043E+00  6.4307E-01 -1.4764E-01
             8.0286E-01
 GRADIENT:  -2.5840E+00 -9.1321E-01  6.8110E-01  1.3528E+01  2.1786E-01 -9.6993E-01 -1.8228E-01 -1.6596E+00  1.1007E+00  6.0463E-01
            -4.8804E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2931.13826547608        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:     1340
 NPARAMETR:  9.7571E-01  1.6103E-01  8.9157E-02  8.0958E-01  1.4359E-01  8.7109E-01  1.5181E+00  6.2220E-01  1.6855E+00  7.4599E-01
             2.0333E+00
 PARAMETER:  7.5409E-02 -1.7262E+00 -2.3174E+00 -1.1124E-01 -1.8408E+00 -3.8004E-02  5.1749E-01 -3.7450E-01  6.2205E-01 -1.9304E-01
             8.0964E-01
 GRADIENT:   1.2463E+01  1.9692E+00  4.4725E+01 -3.6088E+01 -3.2869E+01  4.1665E+00  4.7232E+00 -5.5600E+00 -1.0605E+01 -4.0205E+00
             3.5955E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2939.68202286090        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1523
 NPARAMETR:  9.8060E-01  1.8676E-01  1.1390E-01  8.8901E-01  1.6777E-01  8.7950E-01  1.4923E+00  7.0049E-01  1.5640E+00  6.7307E-01
             2.0100E+00
 PARAMETER:  8.0408E-02 -1.5779E+00 -2.0724E+00 -1.7645E-02 -1.6851E+00 -2.8407E-02  5.0033E-01 -2.5598E-01  5.4727E-01 -2.9590E-01
             7.9814E-01
 GRADIENT:   1.9371E+01  2.2687E+00  3.6707E+01 -5.9802E+01 -4.2952E+01  6.2410E+00  2.9585E+00 -7.4874E-01 -1.6541E+00 -6.4089E+00
             1.9323E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -2941.87232780130        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1700
 NPARAMETR:  9.7825E-01  2.0066E-01  1.2692E-01  9.8006E-01  1.8047E-01  8.7545E-01  1.4953E+00  6.4712E-01  1.3994E+00  6.5100E-01
             2.0137E+00
 PARAMETER:  7.8006E-02 -1.5062E+00 -1.9642E+00  7.9855E-02 -1.6122E+00 -3.3019E-02  5.0232E-01 -3.3523E-01  4.3604E-01 -3.2925E-01
             7.9999E-01
 GRADIENT:   9.9741E+00  2.9090E+00  1.9852E+01 -1.8619E+01 -2.7735E+01  3.5284E+00  2.6007E+00 -3.5568E+00 -1.1799E+01 -8.8636E+00
             1.1947E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -2943.06096349779        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1875
 NPARAMETR:  9.7361E-01  1.9553E-01  1.2072E-01  9.8468E-01  1.7685E-01  8.6642E-01  1.4775E+00  7.1338E-01  1.5093E+00  6.7680E-01
             1.9913E+00
 PARAMETER:  7.3256E-02 -1.5320E+00 -2.0143E+00  8.4559E-02 -1.6325E+00 -4.3382E-02  4.9037E-01 -2.3774E-01  5.1164E-01 -2.9038E-01
             7.8878E-01
 GRADIENT:  -2.1250E-01  7.9748E-01  8.5846E-02 -3.5297E-01  2.4597E-01  9.8513E-03  1.5327E-01  7.2680E-02  4.2734E-01  1.7312E-01
             1.3572E+00

0ITERATION NO.:  110    OBJECTIVE VALUE:  -2943.06269990465        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2052            RESET HESSIAN, TYPE II
 NPARAMETR:  9.7367E-01  1.9506E-01  1.2042E-01  9.8410E-01  1.7657E-01  8.6642E-01  1.4772E+00  7.1343E-01  1.5089E+00  6.7675E-01
             1.9899E+00
 PARAMETER:  7.3312E-02 -1.5344E+00 -2.0168E+00  8.3971E-02 -1.6340E+00 -4.3389E-02  4.9012E-01 -2.3768E-01  5.1136E-01 -2.9046E-01
             7.8809E-01
 GRADIENT:   7.1936E+00  1.0790E+01  1.2117E+01  1.1966E+00  7.1823E+01  4.4792E-01  7.2309E-01  6.7626E-02  1.2595E+00  5.0764E-01
             1.4320E+00

0ITERATION NO.:  112    OBJECTIVE VALUE:  -2943.06269990465        NO. OF FUNC. EVALS.:  58
 CUMULATIVE NO. OF FUNC. EVALS.:     2110
 NPARAMETR:  9.7367E-01  1.9506E-01  1.2042E-01  9.8410E-01  1.7657E-01  8.6642E-01  1.4772E+00  7.1343E-01  1.5089E+00  6.7675E-01
             1.9899E+00
 PARAMETER:  7.3312E-02 -1.5344E+00 -2.0168E+00  8.3971E-02 -1.6340E+00 -4.3389E-02  4.9012E-01 -2.3768E-01  5.1136E-01 -2.9046E-01
             7.8809E-01
 GRADIENT:   8.6699E-03  2.1942E-02  1.4998E-02  2.6554E-03 -6.2411E-02  8.2984E-03 -1.6120E-04  2.3419E-03 -1.4063E-04  4.4126E-03
             1.3444E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2110
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5425E-03  1.2486E-02  1.5241E-02 -2.2043E-03  1.6382E-02
 SE:             2.9451E-02  2.6163E-02  1.6047E-02  2.7399E-02  2.4287E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5823E-01  6.3319E-01  3.4224E-01  9.3588E-01  4.9999E-01

 ETASHRINKSD(%)  1.3358E+00  1.2349E+01  4.6239E+01  8.2099E+00  1.8636E+01
 ETASHRINKVR(%)  2.6538E+00  2.3173E+01  7.1098E+01  1.5746E+01  3.3799E+01
 EBVSHRINKSD(%)  1.5559E+00  1.1509E+01  4.6299E+01  6.1875E+00  1.9197E+01
 EBVSHRINKVR(%)  3.0875E+00  2.1694E+01  7.1162E+01  1.1992E+01  3.4708E+01
 RELATIVEINF(%)  9.6862E+01  3.0231E+01  4.2279E+00  5.1185E+01  9.4042E+00
 EPSSHRINKSD(%)  2.1449E+01
 EPSSHRINKVR(%)  3.8297E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2943.0626999046503     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1288.9733401362396     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    46.16
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    14.53
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2943.063       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.74E-01  1.95E-01  1.20E-01  9.84E-01  1.77E-01  8.66E-01  1.48E+00  7.13E-01  1.51E+00  6.77E-01  1.99E+00
 


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
+        1.54E+03
 
 TH 2
+       -1.93E+01  1.05E+04
 
 TH 3
+        6.48E+00 -3.74E+03  3.95E+04
 
 TH 4
+        1.25E+01 -1.40E+02 -6.24E+02  4.36E+02
 
 TH 5
+       -2.55E+01 -6.14E+03 -3.41E+04 -6.21E+02  5.15E+04
 
 TH 6
+        2.52E+01 -1.37E+01  2.96E+04  5.30E+01 -2.63E+01  2.80E+02
 
 TH 7
+        7.27E+00  4.77E+01  4.18E+01 -4.69E+00 -3.83E+01 -4.61E+00  5.47E+01
 
 TH 8
+        2.94E+00  1.45E+01  5.76E+01  7.24E-01  1.78E+01 -4.87E+00  3.78E+00  2.78E+01
 
 TH 9
+        1.42E+01  2.50E+01  1.92E+02 -7.51E+00  1.45E+02  5.97E+00  2.06E+00  1.17E+00  5.10E+01
 
 TH10
+        1.03E+01  2.41E+01  1.53E+02  1.22E+01  1.68E+02 -5.14E+00  1.15E+01  3.78E+01  4.97E+00  2.12E+02
 
 TH11
+       -2.25E+01 -2.89E+01 -9.92E+01  3.83E+00  7.18E+01  5.98E-01  5.57E+00  1.63E+01  9.81E+00  2.54E+00  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       60.806
Stop Time:
Sat Sep 18 01:30:45 CDT 2021
