Sat Sep 25 12:56:30 CDT 2021
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
$DATA ../../../../data/spa/TD1/dat50.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m50.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1689.67915742881        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.6608E+01 -4.2227E+01 -2.4959E+01 -5.3540E+01 -5.6669E+00 -1.2489E-02  7.4330E+00  1.8105E+01  5.9939E+00  2.2973E+01
             6.4952E+00

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1696.01380098936        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9385E-01  9.7739E-01  1.1813E+00  1.0807E+00  1.0751E+00  9.9976E-01  8.9527E-01  7.8619E-01  9.9592E-01  9.0424E-01
             1.0442E+00
 PARAMETER:  9.3829E-02  7.7128E-02  2.6662E-01  1.7763E-01  1.7242E-01  9.9763E-02 -1.0632E-02 -1.4056E-01  9.5915E-02 -6.5873E-04
             1.4323E-01
 GRADIENT:   3.1174E+01  1.9414E+01 -5.3134E+00  4.1426E+01  2.3600E+01 -3.9819E-02  2.4453E+00  4.0096E+00  1.6741E+00 -1.2163E+01
             1.3043E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1697.71715550951        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      164
 NPARAMETR:  1.0021E+00  9.3162E-01  1.2125E+00  1.1159E+00  1.0465E+00  1.0024E+00  7.6456E-01  5.7441E-01  9.8994E-01  1.0410E+00
             9.9336E-01
 PARAMETER:  1.0209E-01  2.9170E-02  2.9271E-01  2.0962E-01  1.4548E-01  1.0238E-01 -1.6846E-01 -4.5441E-01  8.9886E-02  1.4021E-01
             9.3333E-02
 GRADIENT:   5.6142E+01  4.2982E+01  1.5501E+01  5.7873E+01 -2.9672E+01  1.5851E+00 -5.7193E-01  1.2454E+00 -3.5284E+00  2.3461E+00
            -5.2533E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.77858246563        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      235
 NPARAMETR:  9.8520E-01  9.2083E-01  1.0491E+00  1.0916E+00  9.8054E-01  1.0008E+00  9.2396E-01  4.0764E-01  9.5820E-01  9.5391E-01
             9.9110E-01
 PARAMETER:  8.5085E-02  1.7520E-02  1.4796E-01  1.8765E-01  8.0346E-02  1.0082E-01  2.0913E-02 -7.9737E-01  5.7301E-02  5.2813E-02
             9.1060E-02
 GRADIENT:   1.5822E+01  9.2639E+00  3.0902E-01  1.4153E+01 -6.4563E+00  4.1895E-01 -2.0970E-01  1.2320E+00 -1.0507E+00  1.9877E+00
            -2.2296E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.79456403344        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  9.8138E-01  9.2025E-01  1.0246E+00  1.0856E+00  9.7124E-01  1.0008E+00  9.7252E-01  3.3919E-01  9.4985E-01  9.3951E-01
             9.9332E-01
 PARAMETER:  8.1201E-02  1.6891E-02  1.2426E-01  1.8216E-01  7.0816E-02  1.0082E-01  7.2140E-02 -9.8121E-01  4.8549E-02  3.7608E-02
             9.3293E-02
 GRADIENT:   6.6943E+00  2.9601E+00 -8.8743E-01  5.1525E+00 -2.2122E+00  1.2505E-01  1.7180E-01  8.8180E-01 -1.4401E-01  1.6661E+00
            -8.2842E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.79893125723        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      375
 NPARAMETR:  9.7973E-01  9.2184E-01  1.0026E+00  1.0813E+00  9.6197E-01  1.0010E+00  1.0015E+00  2.8182E-01  9.4306E-01  9.2530E-01
             9.9436E-01
 PARAMETER:  7.9519E-02  1.8619E-02  1.0264E-01  1.7818E-01  6.1225E-02  1.0103E-01  1.0149E-01 -1.1665E+00  4.1371E-02  2.2364E-02
             9.4339E-02
 GRADIENT:   2.4849E+00  4.0833E-01 -1.1731E+00  1.3105E+00 -4.4814E-01  6.8970E-03  2.0774E-01  6.2304E-01  5.0415E-02  1.2222E+00
            -2.1173E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.80893116002        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      445
 NPARAMETR:  9.7850E-01  9.2586E-01  9.7963E-01  1.0760E+00  9.5264E-01  1.0013E+00  1.0243E+00  2.0912E-01  9.3764E-01  9.1093E-01
             9.9526E-01
 PARAMETER:  7.8270E-02  2.2970E-02  7.9419E-02  1.7327E-01  5.1477E-02  1.0132E-01  1.2400E-01 -1.4649E+00  3.5610E-02  6.7131E-03
             9.5249E-02
 GRADIENT:  -8.2212E-01 -1.2917E+00 -1.2215E+00 -1.4830E+00  7.7588E-01 -6.9981E-02  2.7110E-01  3.4959E-01  1.7218E-01  7.5379E-01
             2.8379E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.85452520229        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.7801E-01  9.3376E-01  9.6050E-01  1.0696E+00  9.4619E-01  1.0016E+00  1.0254E+00  8.5058E-02  9.3641E-01  8.9984E-01
             9.9605E-01
 PARAMETER:  7.7765E-02  3.1462E-02  5.9701E-02  1.6724E-01  4.4689E-02  1.0162E-01  1.2508E-01 -2.3644E+00  3.4298E-02 -5.5405E-03
             9.6044E-02
 GRADIENT:  -2.4549E+00 -1.7497E+00 -6.5767E-01 -2.6485E+00  1.1266E+00 -8.9427E-02 -2.1618E-01  5.5812E-02 -1.9323E-01 -1.0783E-01
             3.9111E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.97651216563        NO. OF FUNC. EVALS.:  96
 CUMULATIVE NO. OF FUNC. EVALS.:      614
 NPARAMETR:  9.7938E-01  9.2895E-01  9.6004E-01  1.0742E+00  9.4232E-01  1.0018E+00  1.0559E+00  1.6964E-02  9.2931E-01  8.9589E-01
             9.9435E-01
 PARAMETER:  7.9161E-02  2.6304E-02  5.9216E-02  1.7155E-01  4.0591E-02  1.0175E-01  1.5442E-01 -3.9767E+00  2.6689E-02 -9.9322E-03
             9.4329E-02
 GRADIENT:  -4.0969E+01 -1.4348E+00  1.7848E+00 -1.3367E+01 -2.8509E+00 -4.4467E+00  4.2221E-01  2.0252E-03 -4.2132E-01 -3.8312E-01
            -5.2252E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1699.99529849900        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      792
 NPARAMETR:  9.9425E-01  7.3015E-01  1.0531E+00  1.2134E+00  9.0488E-01  1.0053E+00  1.1366E+00  1.1386E-02  8.7451E-01  9.3338E-01
             9.9568E-01
 PARAMETER:  9.4234E-02 -2.1451E-01  1.5172E-01  2.9340E-01  4.9891E-05  1.0525E-01  2.2805E-01 -4.3753E+00 -3.4095E-02  3.1060E-02
             9.5670E-02
 GRADIENT:  -3.1398E+00  7.9373E+00  4.1681E+00  9.6570E+00 -9.5656E+00 -1.1917E+00 -4.3547E-01  4.6676E-04 -6.3441E-01 -2.6358E-01
            -9.7935E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1700.37648833871        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      973
 NPARAMETR:  9.9713E-01  5.9801E-01  1.2140E+00  1.3031E+00  9.4282E-01  1.0073E+00  1.0484E+00  2.9234E-01  8.6734E-01  1.0246E+00
             9.9659E-01
 PARAMETER:  9.7129E-02 -4.1414E-01  2.9394E-01  3.6476E-01  4.1115E-02  1.0732E-01  1.4722E-01 -1.1298E+00 -4.2323E-02  1.2431E-01
             9.6588E-02
 GRADIENT:   8.8312E+00  6.0919E+00  3.0337E-01  1.3742E+01 -2.3789E+00  7.3030E-01  1.1317E+00  2.0433E-01  2.0345E+00  2.6740E+00
            -6.9508E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1701.08796944475        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1148
 NPARAMETR:  9.8918E-01  4.0308E-01  1.2364E+00  1.4146E+00  8.8949E-01  1.0016E+00  8.9304E-01  9.6065E-02  8.1379E-01  1.0253E+00
             1.0028E+00
 PARAMETER:  8.9116E-02 -8.0862E-01  3.1221E-01  4.4686E-01 -1.7112E-02  1.0156E-01 -1.3127E-02 -2.2427E+00 -1.0605E-01  1.2494E-01
             1.0282E-01
 GRADIENT:  -2.6306E+00  2.6273E+00  2.9668E+00  5.6716E+00 -5.2278E+00 -1.7565E-01  3.9196E-02 -3.1316E-03 -9.0277E-01 -9.9997E-02
             6.5984E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1701.31885073207        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1323
 NPARAMETR:  9.8621E-01  2.3204E-01  1.2600E+00  1.5167E+00  8.5542E-01  9.9820E-01  6.0681E-01  2.6395E-02  7.6874E-01  1.0308E+00
             1.0001E+00
 PARAMETER:  8.6115E-02 -1.3609E+00  3.3112E-01  5.1655E-01 -5.6165E-02  9.8196E-02 -3.9953E-01 -3.5346E+00 -1.6300E-01  1.3036E-01
             1.0012E-01
 GRADIENT:  -2.7397E+00  6.6089E-01 -1.2514E+00  6.2796E+00  1.7807E+00 -2.9695E-01 -7.6529E-03 -5.2747E-04 -7.8485E-01 -6.5162E-01
            -2.9597E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1701.34397991172        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1498
 NPARAMETR:  9.8636E-01  1.7706E-01  1.2613E+00  1.5467E+00  8.4030E-01  9.9780E-01  4.7744E-01  1.1143E-02  7.5448E-01  1.0310E+00
             1.0000E+00
 PARAMETER:  8.6265E-02 -1.6312E+00  3.3217E-01  5.3610E-01 -7.4002E-02  9.7795E-02 -6.3932E-01 -4.3969E+00 -1.8172E-01  1.3055E-01
             1.0003E-01
 GRADIENT:  -9.4723E-02  6.3065E-02 -1.6287E-01  7.0775E-01  3.8247E-02 -2.9562E-02  2.0614E-03 -7.2589E-05 -8.0300E-02  1.4439E-02
            -7.1997E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1701.34406355454        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1678
 NPARAMETR:  9.8628E-01  1.7097E-01  1.2634E+00  1.5501E+00  8.3962E-01  9.9775E-01  4.6156E-01  1.0141E-02  7.5290E-01  1.0314E+00
             1.0003E+00
 PARAMETER:  8.6186E-02 -1.6663E+00  3.3383E-01  5.3830E-01 -7.4804E-02  9.7750E-02 -6.7315E-01 -4.4911E+00 -1.8382E-01  1.3093E-01
             1.0028E-01
 GRADIENT:   1.7205E-03 -1.4935E-02 -6.2872E-02 -1.4067E-01  1.1427E-01  1.0447E-03  2.6299E-03 -6.4427E-05  1.8641E-02 -5.4571E-04
             1.0735E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1701.34459181441        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1858
 NPARAMETR:  9.8624E-01  1.6936E-01  1.2631E+00  1.5509E+00  8.3895E-01  9.9771E-01  2.4478E-01  1.3227E-02  7.5313E-01  1.0314E+00
             1.0002E+00
 PARAMETER:  8.6143E-02 -1.6757E+00  3.3354E-01  5.3886E-01 -7.5603E-02  9.7711E-02 -1.3074E+00 -4.2255E+00 -1.8351E-01  1.3087E-01
             1.0019E-01
 GRADIENT:  -2.9328E-02  1.3546E-03 -1.9405E-03 -4.8428E-03 -8.8266E-02 -2.0444E-03  2.7616E-04 -1.0867E-04  8.5453E-02 -7.9562E-03
            -1.6889E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1701.34585511283        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2037
 NPARAMETR:  9.8607E-01  1.5937E-01  1.2657E+00  1.5571E+00  8.3744E-01  9.9752E-01  1.0000E-02  8.8677E-02  7.5034E-01  1.0319E+00
             1.0003E+00
 PARAMETER:  8.5969E-02 -1.7365E+00  3.3560E-01  5.4282E-01 -7.7402E-02  9.7521E-02 -5.2379E+00 -2.3228E+00 -1.8722E-01  1.3135E-01
             1.0030E-01
 GRADIENT:   1.9320E-02 -1.1631E-02 -2.0917E-01  1.4658E-01 -5.2492E-02 -2.4934E-03  0.0000E+00 -3.2179E-03  1.1894E-01  3.2759E-01
             1.9813E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1701.35159924050        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     2215
 NPARAMETR:  9.8597E-01  1.5544E-01  1.2701E+00  1.5597E+00  8.3746E-01  9.9751E-01  1.0000E-02  2.3604E-01  7.4848E-01  1.0223E+00
             9.9915E-01
 PARAMETER:  8.5867E-02 -1.7615E+00  3.3910E-01  5.4450E-01 -7.7384E-02  9.7511E-02 -7.0571E+00 -1.3438E+00 -1.8971E-01  1.2210E-01
             9.9150E-02
 GRADIENT:  -9.2803E-04 -8.1498E-03 -5.0528E-01 -1.4520E-01 -4.8873E-02  1.2885E-03  0.0000E+00  1.0253E-02 -1.9137E-01  1.4714E-01
             8.6559E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1701.35352258753        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2391
 NPARAMETR:  9.8597E-01  1.5669E-01  1.2813E+00  1.5600E+00  8.4255E-01  9.9754E-01  1.0000E-02  2.5111E-01  7.4907E-01  1.0246E+00
             9.9931E-01
 PARAMETER:  8.5870E-02 -1.7535E+00  3.4786E-01  5.4470E-01 -7.1326E-02  9.7537E-02 -6.9201E+00 -1.2819E+00 -1.8892E-01  1.2434E-01
             9.9311E-02
 GRADIENT:  -8.0417E-03 -6.4933E-04  2.9495E-02  8.5399E-02 -1.7261E-02 -6.4412E-04  0.0000E+00 -3.1950E-04 -7.0260E-03 -2.6460E-02
            -1.6715E-02

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1701.35402878799        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2573
 NPARAMETR:  9.8609E-01  1.6340E-01  1.2801E+00  1.5561E+00  8.4377E-01  9.9765E-01  1.0000E-02  2.4904E-01  7.5117E-01  1.0248E+00
             9.9932E-01
 PARAMETER:  8.5991E-02 -1.7115E+00  3.4693E-01  5.4216E-01 -6.9870E-02  9.7651E-02 -5.7547E+00 -1.2901E+00 -1.8612E-01  1.2450E-01
             9.9323E-02
 GRADIENT:  -3.3713E-02  3.1650E-02  5.4339E-02  3.4336E-01 -1.3596E-01 -5.9910E-03  0.0000E+00  5.7760E-04  1.0981E-02 -2.1301E-02
            -1.2078E-02

0ITERATION NO.:   96    OBJECTIVE VALUE:  -1701.35402878799        NO. OF FUNC. EVALS.:  29
 CUMULATIVE NO. OF FUNC. EVALS.:     2602
 NPARAMETR:  9.8611E-01  1.6312E-01  1.2800E+00  1.5559E+00  8.4385E-01  9.9767E-01  1.0000E-02  2.4873E-01  7.5115E-01  1.0249E+00
             9.9934E-01
 PARAMETER:  8.5991E-02 -1.7115E+00  3.4693E-01  5.4216E-01 -6.9870E-02  9.7651E-02 -5.7547E+00 -1.2901E+00 -1.8612E-01  1.2450E-01
             9.9323E-02
 GRADIENT:  -4.6092E-02  2.0404E-02  5.4195E-02  4.8225E-01 -9.0765E-02 -8.6729E-03  0.0000E+00  4.7003E-04  9.8788E-03 -1.3059E-02
            -1.0139E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2602
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.0745E-06 -5.9726E-05 -5.2462E-03 -5.8517E-03 -2.3492E-02
 SE:             2.9830E-02  3.0831E-05  4.3569E-03  2.9207E-02  2.4996E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9986E-01  5.2717E-02  2.2855E-01  8.4120E-01  3.4731E-01

 ETASHRINKSD(%)  6.5868E-02  9.9897E+01  8.5404E+01  2.1535E+00  1.6261E+01
 ETASHRINKVR(%)  1.3169E-01  1.0000E+02  9.7870E+01  4.2606E+00  2.9877E+01
 EBVSHRINKSD(%)  4.0925E-01  9.9904E+01  8.5575E+01  2.2890E+00  1.3377E+01
 EBVSHRINKVR(%)  8.1683E-01  1.0000E+02  9.7919E+01  4.5256E+00  2.4965E+01
 RELATIVEINF(%)  9.7046E+01  4.4532E-06  2.4952E-01  6.0322E+00  5.8722E+00
 EPSSHRINKSD(%)  4.1766E+01
 EPSSHRINKVR(%)  6.6089E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1701.3540287879948     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -966.20320222425664     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    28.06
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.36
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1701.354       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.86E-01  1.63E-01  1.28E+00  1.56E+00  8.44E-01  9.98E-01  1.00E-02  2.49E-01  7.51E-01  1.02E+00  9.99E-01
 


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
+        1.15E+03
 
 TH 2
+       -2.76E+01  3.99E+02
 
 TH 3
+        2.06E+00  9.58E+01  2.33E+02
 
 TH 4
+       -1.10E+01  4.69E+02 -4.61E+01  7.77E+02
 
 TH 5
+        2.74E+00 -3.19E+02 -4.62E+02 -3.93E+01  1.14E+03
 
 TH 6
+        3.00E+00 -5.77E+00 -7.12E-01 -3.01E+00  2.16E+00  2.10E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        2.27E+00 -2.25E+00 -7.65E+00 -8.31E-01  3.85E+00  2.18E+00  0.00E+00  3.32E+00
 
 TH 9
+        1.75E+00 -9.97E+01  7.90E+00 -2.53E+00 -3.10E+00  2.60E+00  0.00E+00  5.31E-01  3.24E+02
 
 TH10
+        5.67E-02  9.09E+00 -1.06E+01 -1.70E+00 -7.44E+01 -1.58E+00  0.00E+00  1.02E+01  1.73E+00  1.03E+02
 
 TH11
+       -5.93E+00 -1.50E+01 -3.03E+01 -1.10E+01  1.66E+01  3.79E+00  0.00E+00  5.21E+00  2.29E+01  1.83E+01  2.32E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.492
Stop Time:
Sat Sep 25 12:57:05 CDT 2021
