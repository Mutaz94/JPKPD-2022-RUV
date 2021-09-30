Wed Sep 29 14:55:45 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1688.02945519678        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.2112E+02  3.9078E+01 -4.0549E+01  1.2811E+02  2.3022E+00  6.7852E+01  7.6852E+00  1.5745E+01  3.2330E+01  1.6359E+01
            -5.4708E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.52174909862        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      192
 NPARAMETR:  9.9236E-01  1.1383E+00  1.3187E+00  9.2331E-01  1.1981E+00  9.5031E-01  9.7162E-01  8.7112E-01  8.3007E-01  8.9964E-01
             1.2214E+00
 PARAMETER:  9.2329E-02  2.2952E-01  3.7667E-01  2.0205E-02  2.8076E-01  4.9033E-02  7.1212E-02 -3.7970E-02 -8.6250E-02 -5.7573E-03
             2.9996E-01
 GRADIENT:  -6.9402E+01  2.0579E+01 -9.8766E+00  3.6150E+01  3.4007E+01 -7.2609E+00 -4.1180E+00  1.0821E+00 -9.8759E+00 -2.0726E+01
             1.8592E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1702.14454997464        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0042E+00  9.8560E-01  1.6735E+00  1.0439E+00  1.2259E+00  9.5955E-01  7.9335E-01  5.7575E-01  9.0535E-01  1.1101E+00
             1.1892E+00
 PARAMETER:  1.0421E-01  8.5493E-02  6.1490E-01  1.4299E-01  3.0371E-01  5.8705E-02 -1.3149E-01 -4.5208E-01  5.7002E-04  2.0442E-01
             2.7326E-01
 GRADIENT:  -3.6005E+01  4.3599E+01 -2.1291E+00  6.6256E+01  1.1821E+00 -1.4941E+00 -1.5526E+00 -1.1108E-01 -9.4568E+00 -1.9712E+00
             1.1805E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1704.39647623434        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.0184E+00  1.0555E+00  1.6993E+00  9.5780E-01  1.2603E+00  9.6203E-01  5.4785E-01  5.6905E-01  1.0490E+00  1.1543E+00
             1.1586E+00
 PARAMETER:  1.1828E-01  1.5399E-01  6.3023E-01  5.6883E-02  3.3134E-01  6.1288E-02 -5.0176E-01 -4.6380E-01  1.4780E-01  2.4347E-01
             2.4720E-01
 GRADIENT:  -9.6527E-01 -3.8542E-01 -3.2089E-01  1.6938E+00  1.8240E+00 -1.6048E-01  9.4429E-02 -1.7954E-01  3.9532E-01 -1.0274E+00
             1.9366E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1704.53795673333        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.0197E+00  1.1729E+00  1.6843E+00  8.7933E-01  1.2949E+00  9.6388E-01  5.1049E-01  9.2569E-01  1.1319E+00  1.1646E+00
             1.1487E+00
 PARAMETER:  1.1951E-01  2.5948E-01  6.2134E-01 -2.8598E-02  3.5842E-01  6.3216E-02 -5.7239E-01  2.2783E-02  2.2391E-01  2.5235E-01
             2.3863E-01
 GRADIENT:   4.2905E-01  7.3410E-01  9.2296E-02  1.2859E+00 -9.8943E-01  1.8585E-01 -1.7608E-01  1.7619E-01 -2.3860E-01  1.7725E-01
             2.1410E-02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1704.57825723481        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.0204E+00  1.3216E+00  1.4833E+00  7.8455E-01  1.3108E+00  9.6489E-01  5.4635E-01  6.7843E-01  1.2222E+00  1.1640E+00
             1.1514E+00
 PARAMETER:  1.2019E-01  3.7885E-01  4.9429E-01 -1.4265E-01  3.7065E-01  6.4260E-02 -5.0450E-01 -2.8797E-01  3.0063E-01  2.5190E-01
             2.4095E-01
 GRADIENT:  -3.0515E-01  7.2911E+00 -4.4136E-01  7.7227E+00 -2.0968E-01  1.1246E-01 -1.5078E-01 -1.3225E-02 -5.7755E-01  2.0745E-01
            -1.9533E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1704.59199703313        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1077
 NPARAMETR:  1.0210E+00  1.4157E+00  1.3857E+00  7.2180E-01  1.3276E+00  9.6582E-01  5.5111E-01  5.2502E-01  1.2931E+00  1.1667E+00
             1.1549E+00
 PARAMETER:  1.2081E-01  4.4761E-01  4.2623E-01 -2.2600E-01  3.8341E-01  6.5227E-02 -4.9582E-01 -5.4432E-01  3.5703E-01  2.5414E-01
             2.4404E-01
 GRADIENT:  -9.3417E-02  8.0719E+00 -4.3219E-01  7.3759E+00  1.4852E-01  2.6483E-01 -2.7734E-01 -1.4201E-02 -1.1362E+00  2.6702E-01
             2.6772E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1704.62121071702        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1259
 NPARAMETR:  1.0221E+00  1.4768E+00  1.3223E+00  6.7152E-01  1.3400E+00  9.6613E-01  5.5813E-01  4.2770E-01  1.3604E+00  1.1670E+00
             1.1564E+00
 PARAMETER:  1.2182E-01  4.8990E-01  3.7937E-01 -2.9820E-01  3.9266E-01  6.5546E-02 -4.8316E-01 -7.4932E-01  4.0776E-01  2.5444E-01
             2.4528E-01
 GRADIENT:   1.9612E+00 -9.5672E+00 -1.2407E-01 -4.0484E+00  9.3807E-01  2.7758E-01  8.6074E-01 -1.4048E-03  3.6866E-01  2.8747E-01
             8.7247E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1704.64028972710        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1437             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0229E+00  1.4830E+00  1.3255E+00  6.6895E-01  1.3415E+00  9.6618E-01  5.4257E-01  5.1444E-01  1.3746E+00  1.1661E+00
             1.1545E+00
 PARAMETER:  1.2263E-01  4.9406E-01  3.8175E-01 -3.0204E-01  3.9377E-01  6.5592E-02 -5.1143E-01 -5.6467E-01  4.1814E-01  2.5369E-01
             2.4369E-01
 GRADIENT:   4.4891E+02  3.4668E+02  5.4384E-01  7.9317E+01  1.4060E+01  4.3410E+01  9.6467E+00  2.9490E-02  1.7319E+01  1.3467E+00
             2.0958E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1704.64205223225        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1617             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0222E+00  1.4839E+00  1.3275E+00  6.6919E-01  1.3418E+00  9.6591E-01  5.4143E-01  5.0917E-01  1.3747E+00  1.1670E+00
             1.1546E+00
 PARAMETER:  1.2199E-01  4.9465E-01  3.8331E-01 -3.0168E-01  3.9404E-01  6.5315E-02 -5.1354E-01 -5.7498E-01  4.1823E-01  2.5446E-01
             2.4374E-01
 GRADIENT:   4.4497E+02  3.4918E+02  6.4733E-01  8.0105E+01  1.3669E+01  4.3381E+01  9.5692E+00  2.6078E-02  1.7177E+01  1.4073E+00
             2.0023E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1704.64214503615        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1798
 NPARAMETR:  1.0222E+00  1.4839E+00  1.3287E+00  6.6922E-01  1.3420E+00  9.6591E-01  5.4104E-01  5.1140E-01  1.3749E+00  1.1673E+00
             1.1546E+00
 PARAMETER:  1.2199E-01  4.9467E-01  3.8421E-01 -3.0164E-01  3.9419E-01  6.5311E-02 -5.1426E-01 -5.7060E-01  4.1838E-01  2.5467E-01
             2.4372E-01
 GRADIENT:   2.2240E+00 -3.2692E+00 -6.2230E-02 -2.1244E-01 -1.6512E-01  1.5724E-01  7.1502E-02  1.1356E-03  4.5566E-02  5.1296E-04
             3.3655E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1704.64222864184        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1981             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0222E+00  1.4839E+00  1.3302E+00  6.6924E-01  1.3424E+00  9.6590E-01  5.4057E-01  5.0925E-01  1.3752E+00  1.1675E+00
             1.1545E+00
 PARAMETER:  1.2199E-01  4.9467E-01  3.8536E-01 -3.0161E-01  3.9442E-01  6.5307E-02 -5.1513E-01 -5.7482E-01  4.1862E-01  2.5488E-01
             2.4369E-01
 GRADIENT:   4.4503E+02  3.4954E+02  7.2232E-01  8.0165E+01  1.3590E+01  4.3389E+01  9.5618E+00  2.3766E-02  1.7200E+01  1.3998E+00
             1.9169E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1704.64225872971        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2161
 NPARAMETR:  1.0222E+00  1.4839E+00  1.3312E+00  6.6926E-01  1.3425E+00  9.6590E-01  5.4028E-01  5.1015E-01  1.3754E+00  1.1677E+00
             1.1545E+00
 PARAMETER:  1.2199E-01  4.9466E-01  3.8607E-01 -3.0158E-01  3.9456E-01  6.5303E-02 -5.1567E-01 -5.7306E-01  4.1876E-01  2.5500E-01
             2.4367E-01
 GRADIENT:   2.2289E+00 -3.1198E+00  3.5509E-03 -1.6120E-01 -2.0746E-01  1.5783E-01  6.3983E-02 -1.2937E-03  4.8422E-02 -2.3550E-02
            -5.8209E-02

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1704.64227331422        NO. OF FUNC. EVALS.:  98
 CUMULATIVE NO. OF FUNC. EVALS.:     2259
 NPARAMETR:  1.0222E+00  1.4838E+00  1.3311E+00  6.6925E-01  1.3426E+00  9.6590E-01  5.4023E-01  5.1502E-01  1.3755E+00  1.1678E+00
             1.1546E+00
 PARAMETER:  1.2199E-01  4.9463E-01  3.8604E-01 -3.0160E-01  3.9462E-01  6.5302E-02 -5.1576E-01 -5.6355E-01  4.1881E-01  2.5508E-01
             2.4372E-01
 GRADIENT:  -6.2181E-06 -8.3178E-03 -3.0913E-02 -1.6759E-02 -4.7102E-02  1.5068E-04  7.8554E-03  3.9675E-04 -1.0448E-02 -5.1211E-03
             1.4569E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2259
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -5.0339E-04 -3.2736E-02 -4.3273E-03  1.0672E-02 -3.5121E-02
 SE:             2.9798E-02  1.5026E-02  3.3235E-03  2.4988E-02  2.3435E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8652E-01  2.9357E-02  1.9290E-01  6.6932E-01  1.3396E-01

 ETASHRINKSD(%)  1.7288E-01  4.9662E+01  8.8866E+01  1.6286E+01  2.1489E+01
 ETASHRINKVR(%)  3.4547E-01  7.4660E+01  9.8760E+01  2.9920E+01  3.8361E+01
 EBVSHRINKSD(%)  5.6875E-01  4.9337E+01  8.9338E+01  1.6432E+01  1.9638E+01
 EBVSHRINKVR(%)  1.1343E+00  7.4332E+01  9.8863E+01  3.0165E+01  3.5420E+01
 RELATIVEINF(%)  9.8688E+01  1.5907E+00  4.0393E-01  4.7473E+00  2.2205E+01
 EPSSHRINKSD(%)  3.9986E+01
 EPSSHRINKVR(%)  6.3983E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.6422733142233     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -969.49144675048512     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    30.80
 Elapsed covariance  time in seconds:     6.13
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.642       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.48E+00  1.33E+00  6.69E-01  1.34E+00  9.66E-01  5.40E-01  5.15E-01  1.38E+00  1.17E+00  1.15E+00
 


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
 
         3.07E-02  7.92E-01  8.55E-01  5.45E-01  1.84E-01  6.66E-02  1.64E-01  3.59E+00  7.84E-01  1.47E-01  8.82E-02
 


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
+        9.44E-04
 
 TH 2
+        6.28E-03  6.28E-01
 
 TH 3
+       -6.71E-03 -6.23E-01  7.31E-01
 
 TH 4
+       -4.14E-03 -4.31E-01  4.28E-01  2.97E-01
 
 TH 5
+        1.41E-03  1.22E-01 -9.51E-02 -8.37E-02  3.38E-02
 
 TH 6
+        1.94E-05  8.92E-04 -1.64E-03 -8.76E-04  1.19E-03  4.44E-03
 
 TH 7
+       -4.74E-04  6.82E-03 -2.72E-02 -4.11E-03 -1.20E-03  1.57E-03  2.70E-02
 
 TH 8
+       -2.44E-02 -2.23E+00  2.22E+00  1.52E+00 -4.55E-01 -1.90E-02 -6.15E-03  1.29E+01
 
 TH 9
+        6.22E-03  6.10E-01 -5.89E-01 -4.19E-01  1.22E-01  1.27E-03 -6.20E-03 -2.22E+00  6.14E-01
 
 TH10
+        2.30E-04  1.97E-02  2.52E-03 -1.37E-02  1.14E-02  1.52E-03 -5.52E-03 -1.45E-01  2.41E-02  2.17E-02
 
 TH11
+        3.48E-04  3.35E-02 -2.84E-02 -2.29E-02  8.67E-03 -1.81E-04 -1.35E-03 -2.03E-01  3.29E-02  1.98E-03  7.79E-03
 
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
+        3.07E-02
 
 TH 2
+        2.58E-01  7.92E-01
 
 TH 3
+       -2.55E-01 -9.20E-01  8.55E-01
 
 TH 4
+       -2.47E-01 -9.98E-01  9.20E-01  5.45E-01
 
 TH 5
+        2.50E-01  8.39E-01 -6.05E-01 -8.36E-01  1.84E-01
 
 TH 6
+        9.48E-03  1.69E-02 -2.89E-02 -2.41E-02  9.69E-02  6.66E-02
 
 TH 7
+       -9.38E-02  5.23E-02 -1.94E-01 -4.58E-02 -3.97E-02  1.43E-01  1.64E-01
 
 TH 8
+       -2.21E-01 -7.83E-01  7.22E-01  7.78E-01 -6.89E-01 -7.92E-02 -1.04E-02  3.59E+00
 
 TH 9
+        2.58E-01  9.82E-01 -8.79E-01 -9.83E-01  8.47E-01  2.43E-02 -4.81E-02 -7.87E-01  7.84E-01
 
 TH10
+        5.08E-02  1.69E-01  2.00E-02 -1.71E-01  4.20E-01  1.55E-01 -2.28E-01 -2.75E-01  2.09E-01  1.47E-01
 
 TH11
+        1.28E-01  4.80E-01 -3.76E-01 -4.77E-01  5.34E-01 -3.09E-02 -9.30E-02 -6.42E-01  4.77E-01  1.53E-01  8.82E-02
 
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
+        1.28E+03
 
 TH 2
+       -1.16E+02  5.66E+02
 
 TH 3
+        4.92E+01  2.69E+01  3.17E+01
 
 TH 4
+       -2.63E+02  7.00E+02 -3.29E+01  1.15E+03
 
 TH 5
+       -1.48E+02 -1.23E+02 -7.88E+01  1.70E+01  3.62E+02
 
 TH 6
+        3.47E-01  8.62E+01  1.61E+01  6.85E+01 -6.79E+01  2.63E+02
 
 TH 7
+        5.81E+01 -2.22E+01  1.48E+01 -2.05E+01 -4.51E+01 -1.57E+01  6.52E+01
 
 TH 8
+       -1.40E-01  2.05E+00 -4.64E-01  3.87E+00 -7.36E-01  9.83E-01  4.02E-01  3.23E-01
 
 TH 9
+       -2.04E+00 -2.74E+01 -5.14E+00  6.57E+01 -9.03E+00 -8.32E+00  3.10E+01  1.16E+00  7.35E+01
 
 TH10
+        1.20E+01  2.05E+01 -3.06E+00  2.41E+01 -5.08E+01 -1.08E+01  1.23E+01  1.48E+00  3.94E+00  7.60E+01
 
 TH11
+        2.75E+01  2.48E+01  1.62E+00  4.45E+01 -8.15E+01  3.13E+01  2.40E+01  4.97E+00  2.58E+01  3.21E+01  2.65E+02
 
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
 
 Elapsed finaloutput time in seconds:     0.03
 #CPUT: Total CPU Time in Seconds,       36.978
Stop Time:
Wed Sep 29 14:56:23 CDT 2021
