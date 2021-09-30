Wed Sep 29 09:36:21 CDT 2021
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
$DATA ../../../../data/int/D/dat72.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26181.6585565760        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   2.7361E+02  4.5055E+02 -2.1754E+00  1.0224E+02  2.7963E+02 -2.6995E+03 -1.1388E+03 -7.7344E+01 -1.9883E+03 -8.5615E+02
            -5.2826E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1006.31923159104        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  2.1722E+00  1.4986E+00  7.6688E-01  2.8491E+00  8.4681E-01  5.4265E+00  4.7686E+00  1.0443E+00  5.0156E+00  2.6684E+00
             1.0468E+01
 PARAMETER:  8.7573E-01  5.0453E-01 -1.6543E-01  1.1470E+00 -6.6275E-02  1.7913E+00  1.6620E+00  1.4339E-01  1.7126E+00  1.0815E+00
             2.4483E+00
 GRADIENT:   8.2285E+01  7.3524E+00 -4.4183E+01  8.2809E+01 -3.7104E+01  2.2274E+02  9.0997E+01  5.3830E+00  1.1548E+02  5.6797E+01
             3.9639E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1007.56498521454        NO. OF FUNC. EVALS.: 113
 CUMULATIVE NO. OF FUNC. EVALS.:      196
 NPARAMETR:  2.1680E+00  1.5400E+00  7.8903E-01  2.9038E+00  8.7715E-01  5.3422E+00  4.7650E+00  1.0424E+00  5.0524E+00  2.7041E+00
             1.0487E+01
 PARAMETER:  8.7379E-01  5.3180E-01 -1.3694E-01  1.1660E+00 -3.1078E-02  1.7756E+00  1.6613E+00  1.4152E-01  1.7199E+00  1.0948E+00
             2.4501E+00
 GRADIENT:   2.9386E+01 -2.6291E+00 -4.4252E+01  6.7408E+01 -3.7151E+01  1.1381E+02  3.8207E+01  5.0773E+00  8.5724E+01  5.1239E+01
             3.6325E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1147.97888072654        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  3.1861E+00  2.2792E-01  6.0187E+01  2.9580E+00  2.5779E+00  8.0231E+00  1.8126E+01  5.1283E+00  3.4034E+00  1.7406E+00
             9.4238E+00
 PARAMETER:  1.2588E+00 -1.3788E+00  4.1975E+00  1.1845E+00  1.0470E+00  2.1823E+00  2.9973E+00  1.7348E+00  1.3248E+00  6.5424E-01
             2.3432E+00
 GRADIENT:   3.1273E+01  4.5045E+00 -3.6909E+00  3.3315E+01 -4.5892E-01  1.4015E+02  3.3081E+01  2.5591E+00  6.2226E+01  4.3976E+01
             3.1822E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1339.16902258004        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      551
 NPARAMETR:  1.3524E+00  2.4568E-01  6.7161E+01  1.9389E+00  2.4929E+00  3.1570E+00  9.5190E+00  8.0225E-01  2.0938E+00  9.6044E-01
             8.1075E+00
 PARAMETER:  4.0186E-01 -1.3037E+00  4.3071E+00  7.6214E-01  1.0134E+00  1.2496E+00  2.3533E+00 -1.2034E-01  8.3899E-01  5.9637E-02
             2.1928E+00
 GRADIENT:   4.2055E+00 -3.9290E+00 -7.8766E-01  1.4748E+00 -7.6138E+00  1.5705E+01  1.6116E+01  4.0827E-03  1.0764E+00  1.1159E+01
             8.9940E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1347.39161155725        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      727
 NPARAMETR:  1.2553E+00  5.8495E-01  3.5105E+01  1.7314E+00  2.5760E+00  3.0701E+00  6.8671E+00  2.0124E-01  2.2831E+00  4.7742E-01
             7.7163E+00
 PARAMETER:  3.2735E-01 -4.3623E-01  3.6583E+00  6.4895E-01  1.0462E+00  1.2217E+00  2.0267E+00 -1.5033E+00  9.2553E-01 -6.3935E-01
             2.1433E+00
 GRADIENT:  -9.3457E+00 -4.4886E+00 -2.2582E+00  9.9501E+00  1.3220E+01  3.1724E+00 -4.1276E+00 -6.2891E-04  8.4912E+00  1.9299E+00
            -4.9543E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1352.76127007220        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      902
 NPARAMETR:  1.3080E+00  1.3920E+00  1.7795E+01  1.0779E+00  2.5206E+00  3.0383E+00  5.4275E+00  4.1403E-01  1.4981E+00  2.9977E-01
             7.7633E+00
 PARAMETER:  3.6853E-01  4.3073E-01  2.9789E+00  1.7501E-01  1.0245E+00  1.2113E+00  1.7915E+00 -7.8183E-01  5.0423E-01 -1.1047E+00
             2.1494E+00
 GRADIENT:   8.4035E-02  7.9694E-01  1.3242E+00 -1.8743E+00 -2.3477E+00 -5.3917E-01  1.4600E-01 -4.8512E-03  1.8902E+00  5.8924E-01
            -2.8643E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1354.99851202059        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1080
 NPARAMETR:  1.3065E+00  1.8647E+00  3.2939E+00  7.6598E-01  2.2160E+00  3.0558E+00  4.6901E+00  1.0000E-02  7.9067E-01  1.1496E-01
             7.8395E+00
 PARAMETER:  3.6733E-01  7.2312E-01  1.2921E+00 -1.6660E-01  8.9569E-01  1.2170E+00  1.6455E+00 -5.3075E+00 -1.3487E-01 -2.0631E+00
             2.1592E+00
 GRADIENT:  -6.3555E-01 -1.1817E+00  2.7806E-01 -2.3981E+00  8.0703E+00  2.7186E+00 -1.3068E-01  0.0000E+00 -4.8149E-01  1.4125E-01
             1.1239E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1355.27575143662        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1255
 NPARAMETR:  1.3077E+00  1.8403E+00  2.8312E+00  7.8121E-01  2.0993E+00  3.0325E+00  4.7297E+00  1.0000E-02  8.2068E-01  1.0764E-01
             7.7977E+00
 PARAMETER:  3.6827E-01  7.0995E-01  1.1407E+00 -1.4691E-01  8.4161E-01  1.2094E+00  1.6539E+00 -5.8119E+00 -9.7627E-02 -2.1290E+00
             2.1538E+00
 GRADIENT:  -3.9881E-02  7.5079E-02 -2.0646E-02  6.1128E-01  4.4723E-03  5.1471E-02 -1.7351E-01  0.0000E+00 -3.3767E-01  1.2129E-01
             5.2634E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1355.28894916291        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1432
 NPARAMETR:  1.3075E+00  1.8322E+00  2.8860E+00  7.8480E-01  2.1043E+00  3.0318E+00  4.7362E+00  1.0000E-02  8.6271E-01  1.0194E-01
             7.7900E+00
 PARAMETER:  3.6809E-01  7.0552E-01  1.1599E+00 -1.4233E-01  8.4398E-01  1.2092E+00  1.6552E+00 -5.9043E+00 -4.7678E-02 -2.1833E+00
             2.1528E+00
 GRADIENT:  -3.5201E-02 -7.1631E-02 -8.8103E-02 -8.1549E-02  2.5224E-01 -7.8555E-02 -3.1858E-02  0.0000E+00  7.9109E-02  1.0927E-01
            -4.1305E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1355.34473859850        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1610
 NPARAMETR:  1.3076E+00  1.8309E+00  2.8890E+00  7.8599E-01  2.1046E+00  3.0325E+00  4.7400E+00  1.0000E-02  8.5648E-01  1.5196E-02
             7.7947E+00
 PARAMETER:  3.6818E-01  7.0480E-01  1.1609E+00 -1.4081E-01  8.4412E-01  1.2094E+00  1.6560E+00 -8.8776E+00 -5.4927E-02 -4.0867E+00
             2.1534E+00
 GRADIENT:  -4.3596E-02 -6.5238E-03 -5.6922E-02  2.5309E-01  2.1114E-01  3.2713E-02 -1.0552E-01  0.0000E+00 -3.5463E-02  2.4009E-03
             1.9723E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1355.42309035478        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1769
 NPARAMETR:  1.3062E+00  1.8235E+00  2.8887E+00  7.8516E-01  2.1041E+00  3.1237E+00  4.7385E+00  1.0000E-02  8.5837E-01  1.0000E-02
             7.7908E+00
 PARAMETER:  3.6709E-01  7.0076E-01  1.1608E+00 -1.4186E-01  8.4391E-01  1.2390E+00  1.6557E+00 -9.5013E+00 -5.2721E-02 -4.5657E+00
             2.1529E+00
 GRADIENT:   3.6047E+01  2.2436E+01  3.5976E-01  2.1339E+00  4.1673E+00  1.2312E+02  1.1293E+02  0.0000E+00  9.5123E-02  0.0000E+00
             3.9581E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1355.43958037310        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     1957             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.8248E+00  2.8967E+00  7.8627E-01  2.1036E+00  3.1002E+00  4.7653E+00  1.0000E-02  8.5988E-01  1.0000E-02
             7.7895E+00
 PARAMETER:  3.6981E-01  7.0149E-01  1.1636E+00 -1.4045E-01  8.4363E-01  1.2315E+00  1.6614E+00 -9.5013E+00 -5.0958E-02 -4.5654E+00
             2.1528E+00
 GRADIENT:   3.6913E+01  2.2725E+01  3.5034E-01  1.6347E+00  4.1140E+00  1.1916E+02  1.1472E+02  0.0000E+00  1.3321E-01  0.0000E+00
             3.8982E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1355.44090785785        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     2149
 NPARAMETR:  1.3097E+00  1.8213E+00  2.9054E+00  7.8792E-01  2.1033E+00  3.1002E+00  4.7706E+00  1.0000E-02  8.6204E-01  1.0000E-02
             7.7893E+00
 PARAMETER:  3.6982E-01  6.9956E-01  1.1666E+00 -1.3836E-01  8.4349E-01  1.2315E+00  1.6625E+00 -9.5013E+00 -4.8450E-02 -4.5654E+00
             2.1528E+00
 GRADIENT:   5.1082E-01 -1.2605E-01 -3.4592E-02 -5.6645E-01 -2.1562E-01  8.5415E+00  1.2188E+00  0.0000E+00  3.8292E-02  0.0000E+00
             1.5913E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1355.44189098665        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2342
 NPARAMETR:  1.3097E+00  1.8181E+00  2.9146E+00  7.8931E-01  2.1032E+00  3.1002E+00  4.7752E+00  1.0000E-02  8.6443E-01  1.0000E-02
             7.7893E+00
 PARAMETER:  3.6981E-01  6.9782E-01  1.1697E+00 -1.3660E-01  8.4348E-01  1.2315E+00  1.6634E+00 -9.5013E+00 -4.5688E-02 -4.5654E+00
             2.1527E+00
 GRADIENT:   5.0999E-01 -1.3786E-01 -3.5796E-02 -5.2140E-01 -2.4630E-01  8.5408E+00  1.2408E+00  0.0000E+00  3.7990E-02  0.0000E+00
             1.5100E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1355.44300009807        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     2539             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.8151E+00  2.9264E+00  7.9228E-01  2.1041E+00  3.1002E+00  4.7791E+00  1.0000E-02  8.6279E-01  1.0000E-02
             7.7880E+00
 PARAMETER:  3.6981E-01  6.9615E-01  1.1738E+00 -1.3284E-01  8.4389E-01  1.2315E+00  1.6643E+00 -9.5013E+00 -4.7589E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6909E+01  2.2527E+01  2.4556E-01  2.4704E+00  4.3977E+00  1.1913E+02  1.1489E+02  0.0000E+00  1.3578E-02  0.0000E+00
             3.8083E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1355.44384700644        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     2707
 NPARAMETR:  1.3097E+00  1.8124E+00  2.9370E+00  7.9298E-01  2.1044E+00  3.1001E+00  4.7830E+00  1.0000E-02  8.6571E-01  1.0000E-02
             7.7880E+00
 PARAMETER:  3.6979E-01  6.9465E-01  1.1774E+00 -1.3195E-01  8.4404E-01  1.2314E+00  1.6651E+00 -9.5013E+00 -4.4200E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6911E+01  2.2408E+01  2.9137E-01  2.2779E+00  4.2565E+00  1.1912E+02  1.1516E+02  0.0000E+00  3.6278E-02  0.0000E+00
             3.8221E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1355.44436352319        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2896
 NPARAMETR:  1.3097E+00  1.8106E+00  2.9446E+00  7.9373E-01  2.1048E+00  3.1002E+00  4.7862E+00  1.0000E-02  8.6855E-01  1.0000E-02
             7.7883E+00
 PARAMETER:  3.6981E-01  6.9366E-01  1.1800E+00 -1.3102E-01  8.4420E-01  1.2315E+00  1.6657E+00 -9.5013E+00 -4.0929E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   5.0852E-01 -7.5089E-02 -8.5220E-02  6.3294E-03 -5.8892E-02  8.5357E+00  1.0738E+00  0.0000E+00 -2.7551E-02  0.0000E+00
            -3.5031E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1355.44495459294        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     3091             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.8077E+00  2.9566E+00  7.9461E-01  2.1054E+00  3.1002E+00  4.7907E+00  1.0000E-02  8.7244E-01  1.0000E-02
             7.7886E+00
 PARAMETER:  3.6981E-01  6.9208E-01  1.1840E+00 -1.2990E-01  8.4449E-01  1.2315E+00  1.6667E+00 -9.5013E+00 -3.6464E-02 -4.5654E+00
             2.1527E+00
 GRADIENT:   3.6910E+01  2.2218E+01  3.1964E-01  2.0172E+00  4.1760E+00  1.1913E+02  1.1565E+02  0.0000E+00  9.6508E-02  0.0000E+00
             3.8615E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1355.44565952179        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3280
 NPARAMETR:  1.3097E+00  1.8042E+00  2.9743E+00  7.9614E-01  2.1068E+00  3.1003E+00  4.7961E+00  1.0000E-02  8.7584E-01  1.0000E-02
             7.7882E+00
 PARAMETER:  3.6979E-01  6.9014E-01  1.1900E+00 -1.2799E-01  8.4519E-01  1.2315E+00  1.6678E+00 -9.5013E+00 -3.2572E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   5.0681E-01 -1.4480E-01 -5.6334E-02 -1.5994E-01 -6.4065E-02  8.5332E+00  1.2306E+00  0.0000E+00  1.1172E-02  0.0000E+00
            -1.7767E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1355.44593509971        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3472             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.8030E+00  2.9842E+00  7.9714E-01  2.1073E+00  3.1003E+00  4.7979E+00  1.0000E-02  8.7523E-01  1.0000E-02
             7.7879E+00
 PARAMETER:  3.6982E-01  6.8948E-01  1.1933E+00 -1.2672E-01  8.4543E-01  1.2315E+00  1.6682E+00 -9.5013E+00 -3.3272E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6920E+01  2.2097E+01  3.6223E-01  2.1329E+00  4.1494E+00  1.1913E+02  1.1584E+02  0.0000E+00  5.4109E-02  0.0000E+00
             3.8340E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1355.44598358052        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3661
 NPARAMETR:  1.3097E+00  1.8022E+00  2.9902E+00  7.9775E-01  2.1076E+00  3.1003E+00  4.7990E+00  1.0000E-02  8.7493E-01  1.0000E-02
             7.7876E+00
 PARAMETER:  3.6983E-01  6.8902E-01  1.1953E+00 -1.2596E-01  8.4556E-01  1.2315E+00  1.6684E+00 -9.5013E+00 -3.3614E-02 -4.5654E+00
             2.1525E+00
 GRADIENT:   5.1869E-01 -7.1590E-02  4.6480E-03  7.7900E-02 -2.0851E-01  8.5316E+00  1.0500E+00  0.0000E+00 -6.4741E-02  0.0000E+00
            -5.6596E-01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1355.44623997536        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3855             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.8004E+00  2.9945E+00  7.9804E-01  2.1084E+00  3.1003E+00  4.8023E+00  1.0000E-02  8.8027E-01  1.0000E-02
             7.7883E+00
 PARAMETER:  3.6980E-01  6.8798E-01  1.1968E+00 -1.2559E-01  8.4591E-01  1.2315E+00  1.6691E+00 -9.5013E+00 -2.7525E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6903E+01  2.1967E+01  3.2636E-01  1.9909E+00  4.2920E+00  1.1912E+02  1.1617E+02  0.0000E+00  1.1476E-01  0.0000E+00
             3.8661E+01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1355.44635613406        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     4043
 NPARAMETR:  1.3097E+00  1.7999E+00  2.9983E+00  7.9840E-01  2.1085E+00  3.1003E+00  4.8029E+00  1.0000E-02  8.8014E-01  1.0000E-02
             7.7881E+00
 PARAMETER:  3.6981E-01  6.8772E-01  1.1981E+00 -1.2514E-01  8.4599E-01  1.2315E+00  1.6692E+00 -9.5013E+00 -2.7671E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   5.0769E-01 -1.3784E-01 -3.8150E-02 -8.7219E-02 -5.1639E-02  8.5327E+00  1.2260E+00  0.0000E+00  6.0566E-03  0.0000E+00
            -1.8579E-01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1355.44659563878        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4237             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.7973E+00  3.0152E+00  7.9992E-01  2.1094E+00  3.1003E+00  4.8069E+00  1.0000E-02  8.8086E-01  1.0000E-02
             7.7877E+00
 PARAMETER:  3.6983E-01  6.8629E-01  1.2037E+00 -1.2324E-01  8.4641E-01  1.2315E+00  1.6700E+00 -9.5013E+00 -2.6859E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6925E+01  2.1919E+01  3.9964E-01  2.1166E+00  4.1171E+00  1.1912E+02  1.1620E+02  0.0000E+00  5.3894E-02  0.0000E+00
             3.8336E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1355.44667476160        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4422
 NPARAMETR:  1.3097E+00  1.7963E+00  3.0201E+00  8.0037E-01  2.1098E+00  3.1003E+00  4.8082E+00  1.0000E-02  8.8217E-01  1.0000E-02
             7.7876E+00
 PARAMETER:  3.6982E-01  6.8575E-01  1.2053E+00 -1.2268E-01  8.4657E-01  1.2315E+00  1.6703E+00 -9.5013E+00 -2.5370E-02 -4.5654E+00
             2.1525E+00
 GRADIENT:   5.1501E-01 -1.0754E-01  1.8074E-02  2.3695E-02 -1.7303E-01  8.5307E+00  1.1407E+00  0.0000E+00 -3.3089E-02  0.0000E+00
            -4.2231E-01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -1355.44675205032        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     4616             RESET HESSIAN, TYPE I
 NPARAMETR:  1.3097E+00  1.7951E+00  3.0235E+00  8.0071E-01  2.1103E+00  3.1003E+00  4.8104E+00  1.0000E-02  8.8458E-01  1.0000E-02
             7.7880E+00
 PARAMETER:  3.6981E-01  6.8509E-01  1.2064E+00 -1.2226E-01  8.4683E-01  1.2315E+00  1.6708E+00 -9.5013E+00 -2.2639E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:   3.6909E+01  2.1818E+01  3.6294E-01  2.0420E+00  4.2633E+00  1.1911E+02  1.1644E+02  0.0000E+00  9.6230E-02  0.0000E+00
             3.8543E+01

0ITERATION NO.:  132    OBJECTIVE VALUE:  -1355.44675205032        NO. OF FUNC. EVALS.:  65
 CUMULATIVE NO. OF FUNC. EVALS.:     4681
 NPARAMETR:  1.3097E+00  1.7946E+00  3.0294E+00  8.0121E-01  2.1104E+00  3.1003E+00  4.8109E+00  1.0000E-02  8.8370E-01  1.0000E-02
             7.7876E+00
 PARAMETER:  3.6981E-01  6.8509E-01  1.2064E+00 -1.2226E-01  8.4683E-01  1.2315E+00  1.6708E+00 -9.5013E+00 -2.2639E-02 -4.5654E+00
             2.1526E+00
 GRADIENT:  -2.3319E-04  6.2983E-03 -1.8760E-02 -3.9677E-02 -5.7214E-03  1.7295E-04 -6.9719E-03  0.0000E+00  3.0038E-03  0.0000E+00
             2.6716E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4681
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         5.9801E-03  1.0717E-02 -4.2387E-05 -4.2204E-02  2.3409E-05
 SE:             2.9235E-02  2.7459E-02  3.3398E-05  1.0295E-02  1.5298E-04
 N:                     100         100         100         100         100

 P VAL.:         8.3792E-01  6.9632E-01  2.0439E-01  4.1459E-05  8.7838E-01

 ETASHRINKSD(%)  2.0606E+00  8.0095E+00  9.9888E+01  6.5510E+01  9.9487E+01
 ETASHRINKVR(%)  4.0788E+00  1.5377E+01  1.0000E+02  8.8104E+01  9.9997E+01
 EBVSHRINKSD(%)  1.5588E+00  4.4500E+00  9.9878E+01  7.1919E+01  9.9459E+01
 EBVSHRINKVR(%)  3.0933E+00  8.7019E+00  1.0000E+02  9.2114E+01  9.9997E+01
 RELATIVEINF(%)  9.6814E+01  4.9360E+01  3.3855E-05  3.2399E+00  7.8095E-04
 EPSSHRINKSD(%)  7.6421E+00
 EPSSHRINKVR(%)  1.4700E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1355.4467520503249     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       298.64260771808586     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   166.88
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    15.62
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1355.447       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.31E+00  1.80E+00  3.02E+00  8.01E-01  2.11E+00  3.10E+00  4.81E+00  1.00E-02  8.85E-01  1.00E-02  7.79E+00
 


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
+        2.81E+01
 
 TH 2
+        3.53E+00  2.36E+00
 
 TH 3
+        4.15E+00 -6.82E-01  1.37E+00
 
 TH 4
+       -1.69E+01  1.96E+01 -1.62E+01  2.58E+02
 
 TH 5
+       -2.55E+01  2.08E+00 -7.09E+00  7.54E+01  3.78E+01
 
 TH 6
+        3.31E+00  2.11E-01  6.17E-01 -4.32E+00 -3.57E+00  4.15E-01
 
 TH 7
+       -1.75E-01 -1.53E+00  9.22E-01 -1.70E+01 -4.00E+00  1.41E-01  1.19E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.10E+00 -1.56E+00  7.31E-01 -1.55E+01 -2.92E+00  2.78E-02  1.13E+00  0.00E+00  1.11E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -5.02E-01 -8.21E-01  3.87E-01 -8.25E+00 -1.58E+00  5.40E-02  5.97E-01  0.00E+00  6.35E-01  0.00E+00  7.02E-01
 
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
+        6.44E+01
 
 TH 2
+       -2.35E-01  1.37E+01
 
 TH 3
+        2.44E-01  7.20E-01  2.36E+00
 
 TH 4
+       -2.25E+00  2.13E+01 -1.09E+01  2.23E+02
 
 TH 5
+       -1.43E+00 -4.94E+00 -1.06E+01  4.36E+01  6.66E+01
 
 TH 6
+        5.96E-01 -2.94E-02  5.68E-03  4.42E-01 -6.16E-01  1.90E+01
 
 TH 7
+        8.64E-02  1.22E+00 -5.52E-01 -1.61E+01  2.09E+00 -1.50E-01  6.39E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.95E-01 -1.15E+00 -1.15E+00 -1.54E+01  3.16E+00  4.61E-02  1.67E+00  0.00E+00  8.73E+00
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  2.15E+01
 
 TH11
+       -3.05E+00 -1.75E+00  1.09E-02 -1.15E+01 -1.64E-01  8.38E-01  8.39E-01  0.00E+00  3.03E+00  0.00E+00  1.71E+01
 
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
+        6.54E+01
 
 TH 2
+        2.52E+01  1.37E+01
 
 TH 3
+        1.25E+00  1.03E+00  1.59E+00
 
 TH 4
+        4.61E+01  2.26E+01 -2.47E+00  2.35E+02
 
 TH 5
+       -1.22E+01 -4.46E+00 -8.11E+00  9.44E+00  6.09E+01
 
 TH 6
+        1.37E+01  4.50E+00 -4.43E-01 -2.55E+01 -9.12E-01  1.89E+01
 
 TH 7
+        1.51E+00  2.70E+00 -8.99E-01 -2.07E+01  9.20E+00  7.39E+00  1.04E+01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -6.09E+00 -2.52E+00 -1.13E+00 -3.83E+01  6.95E+00  6.39E+00  3.87E+00  0.00E+00  1.42E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.46E+01 -2.38E+01 -2.23E+00 -1.77E+02  2.11E+01  1.83E+01  1.68E+01  0.00E+00  4.65E+01  0.00E+00  6.99E+02
 
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
 #CPUT: Total CPU Time in Seconds,      182.530
Stop Time:
Wed Sep 29 09:39:26 CDT 2021
