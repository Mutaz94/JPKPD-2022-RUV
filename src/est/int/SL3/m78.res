Sat Sep 25 02:40:34 CDT 2021
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
$DATA ../../../../data/int/SL3/dat78.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      985
 NO. OF DATA ITEMS IN DATA SET:   6
 ID DATA ITEM IS DATA ITEM NO.:   1
 DEP VARIABLE IS DATA ITEM NO.:   3
 MDV DATA ITEM IS DATA ITEM NO.:  5
0INDICES PASSED TO SUBROUTINE PRED:
   6   2   4   0   0   0   0   0   0   0   0
0LABELS FOR DATA ITEMS:
 ID TIME DV AMT MDV EVID
0FORMAT FOR DATA:
 (2E4.0,E20.0,E4.0,2E2.0)

 TOT. NO. OF OBS RECS:      885
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
 RAW OUTPUT FILE (FILE): m78.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1018.31873138934        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.4392E+01 -3.7120E+01  9.7063E+01 -1.1208E+01  1.1897E+02  5.2138E+00 -1.5889E+02 -2.0181E+02 -9.6099E+01 -7.3724E+00
            -5.1545E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2777.20868883365        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0023E+00  1.4109E+00  1.1036E+00  8.3294E-01  1.2195E+00  9.2132E-01  1.2967E+00  9.2241E-01  1.0067E+00  9.3858E-01
             2.3077E+00
 PARAMETER:  1.0232E-01  4.4423E-01  1.9855E-01 -8.2795E-02  2.9848E-01  1.8050E-02  3.5983E-01  1.9236E-02  1.0668E-01  3.6608E-02
             9.3624E-01
 GRADIENT:   3.8653E+01  2.9658E+01 -3.0797E-01 -1.5108E+01 -2.1566E+01 -2.2722E+01  2.7089E+01  3.4305E-01 -9.1557E+00 -3.0545E+01
            -2.0350E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2785.47657319735        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0307E+00  1.5962E+00  2.2598E+00  7.9910E-01  1.7028E+00  8.6075E-01  1.2112E+00  4.5176E-01  9.3368E-01  1.3442E+00
             2.4089E+00
 PARAMETER:  1.3023E-01  5.6762E-01  9.1526E-01 -1.2427E-01  6.3228E-01 -4.9954E-02  2.9159E-01 -6.9460E-01  3.1374E-02  3.9577E-01
             9.7918E-01
 GRADIENT:   1.1736E+02  6.4240E+01 -1.4938E+01  7.7997E+01  4.8713E+01 -5.4414E+01  2.3993E+01 -3.2000E-01 -1.5889E+00 -2.3765E+01
            -1.0030E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2803.76397727528        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.8313E-01  1.8073E+00  3.4678E+00  6.1218E-01  2.0628E+00  9.8470E-01  9.2845E-01  3.4980E-01  1.0838E+00  1.7388E+00
             2.4364E+00
 PARAMETER:  8.2983E-02  6.9181E-01  1.3435E+00 -3.9074E-01  8.2408E-01  8.4581E-02  2.5764E-02 -9.5039E-01  1.8047E-01  6.5321E-01
             9.9051E-01
 GRADIENT:  -1.3853E+01  1.0753E+01 -1.1519E+01  2.3677E+01  4.8383E+01  3.1252E+00 -2.6280E+00 -6.8974E-02 -7.8373E-01  8.6932E+00
            -3.8484E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2809.18301670787        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  9.8608E-01  1.9194E+00  7.3872E+00  5.1610E-01  2.0232E+00  9.6902E-01  7.6610E-01  7.9955E-01  1.7056E+00  1.7003E+00
             2.4622E+00
 PARAMETER:  8.5985E-02  7.5200E-01  2.0997E+00 -5.6145E-01  8.0468E-01  6.8528E-02 -1.6644E-01 -1.2370E-01  6.3392E-01  6.3082E-01
             1.0011E+00
 GRADIENT:  -7.8067E+00 -3.5443E+00 -3.5929E-01 -3.9804E+00 -6.3320E+00 -2.5958E+00  5.0332E+00 -3.7557E-02  3.9385E+00 -2.8628E+00
            -7.5330E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2809.91350912226        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      421
 NPARAMETR:  9.9133E-01  1.9454E+00  1.0394E+01  5.0598E-01  2.0859E+00  9.7464E-01  7.1042E-01  1.2758E+00  1.7898E+00  1.7379E+00
             2.4685E+00
 PARAMETER:  9.1294E-02  7.6548E-01  2.4412E+00 -5.8126E-01  8.3522E-01  7.4313E-02 -2.4190E-01  3.4357E-01  6.8210E-01  6.5266E-01
             1.0036E+00
 GRADIENT:  -3.5383E+00 -1.2575E+01 -6.5262E-01 -2.0075E+00  6.1052E-01 -1.1154E+00  4.5838E-01 -4.5579E-02  1.6689E+00 -7.7289E-01
            -2.8133E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2810.38851374567        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      598
 NPARAMETR:  9.9266E-01  2.0208E+00  1.9346E+01  4.6546E-01  2.1182E+00  9.7438E-01  7.4541E-01  3.8452E+00  1.7271E+00  1.7668E+00
             2.4828E+00
 PARAMETER:  9.2628E-02  8.0351E-01  3.0625E+00 -6.6473E-01  8.5059E-01  7.4046E-02 -1.9381E-01  1.4468E+00  6.4642E-01  6.6917E-01
             1.0094E+00
 GRADIENT:  -1.7423E+00 -3.0436E+00 -3.5047E-01 -8.6390E-01 -3.7161E-01 -1.2679E+00 -2.2361E+00 -7.0254E-03 -7.2430E-01  1.4383E+00
             1.2139E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2810.71245942902        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      775
 NPARAMETR:  9.9313E-01  2.0662E+00  3.2513E+01  4.4099E-01  2.1356E+00  9.7915E-01  7.9348E-01  7.9797E+00  1.6321E+00  1.7666E+00
             2.4724E+00
 PARAMETER:  9.3103E-02  8.2570E-01  3.5816E+00 -7.1873E-01  8.5874E-01  7.8925E-02 -1.3133E-01  2.1769E+00  5.8986E-01  6.6904E-01
             1.0052E+00
 GRADIENT:  -4.3362E-01  2.0272E+00 -1.7790E-01  1.2523E+00  9.6646E-01  3.1461E-01 -2.9404E-01 -2.9536E-02  5.5056E-01  4.0566E-01
             3.3348E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2810.80516338730        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      952
 NPARAMETR:  9.9346E-01  2.0878E+00  4.8979E+01  4.2587E-01  2.1390E+00  9.7914E-01  8.1392E-01  1.2767E+01  1.5582E+00  1.7666E+00
             2.4680E+00
 PARAMETER:  9.3443E-02  8.3613E-01  3.9914E+00 -7.5361E-01  8.6036E-01  7.8918E-02 -1.0590E-01  2.6469E+00  5.4355E-01  6.6905E-01
             1.0034E+00
 GRADIENT:   4.7701E-01  4.7429E-01 -3.8625E-01  1.3481E+00  4.3419E-01  2.6615E-01  4.5997E-01  4.7197E-01  1.1525E+00 -2.5116E-01
            -5.7870E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2810.85546078970        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1132
 NPARAMETR:  9.9345E-01  2.1127E+00  7.7965E+01  4.0770E-01  2.1403E+00  9.7846E-01  8.2195E-01  1.9249E+01  1.4866E+00  1.7673E+00
             2.4668E+00
 PARAMETER:  9.3433E-02  8.4797E-01  4.4563E+00 -7.9721E-01  8.6097E-01  7.8228E-02 -9.6078E-02  3.0575E+00  4.9651E-01  6.6944E-01
             1.0029E+00
 GRADIENT:   4.5271E-01 -1.3325E+00 -6.9703E-01  1.7361E+00 -3.2949E-01  8.4745E-03 -6.3839E-01  1.1685E+00  1.4023E+00 -7.7170E-01
            -1.8009E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2810.85927261786        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1309
 NPARAMETR:  9.9343E-01  2.1157E+00  8.1896E+01  4.0568E-01  2.1408E+00  9.7851E-01  8.2169E-01  2.0012E+01  1.4816E+00  1.7674E+00
             2.4671E+00
 PARAMETER:  9.3409E-02  8.4936E-01  4.5055E+00 -8.0220E-01  8.6116E-01  7.8271E-02 -9.6390E-02  3.0963E+00  4.9312E-01  6.6949E-01
             1.0031E+00
 GRADIENT:   1.6657E-01 -2.7431E+01 -1.3473E+01  3.2515E+01 -1.0734E+01  6.9063E-02 -2.2053E+01  2.4192E+01  2.2983E+01 -1.6673E+01
            -8.4966E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2810.86322070482        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1491            RESET HESSIAN, TYPE II
 NPARAMETR:  9.9330E-01  2.1163E+00  8.3806E+01  4.0446E-01  2.1401E+00  9.7843E-01  8.2116E-01  1.9982E+01  1.4677E+00  1.7673E+00
             2.4683E+00
 PARAMETER:  9.3275E-02  8.4965E-01  4.5285E+00 -8.0521E-01  8.6084E-01  7.8189E-02 -9.7037E-02  3.0948E+00  4.8367E-01  6.6944E-01
             1.0035E+00
 GRADIENT:   7.0569E+00  2.7516E+01 -2.4108E-02  2.1789E+00  3.6383E+00  6.1001E-01 -1.4743E-01 -1.9216E-02  3.1539E-01  9.0155E-01
             1.8208E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2810.99011115378        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1674
 NPARAMETR:  9.9320E-01  2.1306E+00  1.2213E+02  3.9290E-01  2.1404E+00  9.7847E-01  8.6278E-01  2.2228E+01  1.0869E+00  1.7666E+00
             2.4683E+00
 PARAMETER:  9.3172E-02  8.5642E-01  4.9051E+00 -8.3421E-01  8.6098E-01  7.8231E-02 -4.7596E-02  3.2013E+00  1.8334E-01  6.6903E-01
             1.0035E+00
 GRADIENT:  -1.4898E-01  1.0111E+00  7.0344E-04  3.0815E-01 -2.6086E-01  6.0540E-03  5.7706E-01 -4.5144E-02  5.2451E-01 -2.3723E-01
             2.7692E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2811.31894091414        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1849
 NPARAMETR:  9.9330E-01  2.1704E+00  5.7240E+02  3.5856E-01  2.1425E+00  9.7838E-01  8.8853E-01  3.6078E+01  2.0163E-01  1.7703E+00
             2.4685E+00
 PARAMETER:  9.3279E-02  8.7490E-01  6.4498E+00 -9.2567E-01  8.6197E-01  7.8147E-02 -1.8185E-02  3.6857E+00 -1.5013E+00  6.7117E-01
             1.0036E+00
 GRADIENT:   1.2534E-01  5.8247E-01 -1.1661E-03  8.2712E-02 -2.6677E-01 -1.2303E-02 -6.1309E-01 -5.0750E-03  3.1020E-02  1.0806E-02
             2.4130E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2811.32925387380        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2024
 NPARAMETR:  9.9321E-01  2.1849E+00  1.3569E+03  3.4839E-01  2.1437E+00  9.7843E-01  8.8448E-01  4.7780E+01  7.3965E-02  1.7703E+00
             2.4681E+00
 PARAMETER:  9.3184E-02  8.8157E-01  7.3130E+00 -9.5444E-01  8.6251E-01  7.8193E-02 -2.2759E-02  3.9666E+00 -2.5042E+00  6.7117E-01
             1.0035E+00
 GRADIENT:  -8.7916E-02  4.3680E-01 -9.7812E-04 -4.1716E-02 -1.1590E-02 -6.1065E-03 -3.8869E-01 -1.5458E-03  4.3378E-03 -2.4453E-02
            -2.7020E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -2811.35415445833        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2203
 NPARAMETR:  9.9398E-01  2.1144E+00  7.6155E+05  3.9934E-01  2.1444E+00  9.7893E-01  9.1722E-01  3.8962E+02  1.0000E-02  1.7720E+00
             2.4684E+00
 PARAMETER:  9.3963E-02  8.4879E-01  1.3643E+01 -8.1794E-01  8.6286E-01  7.8707E-02  1.3587E-02  6.0652E+00 -9.9455E+00  6.7211E-01
             1.0036E+00
 GRADIENT:  -8.8502E+00  9.9783E+00 -3.8800E-06  2.7079E+00 -7.1499E-01  1.5265E-01 -2.8679E-01  2.5116E-05  0.0000E+00  2.3323E-01
             1.8345E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -2811.36053471180        NO. OF FUNC. EVALS.: 212
 CUMULATIVE NO. OF FUNC. EVALS.:     2415             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9367E-01  2.1094E+00  1.1174E+06  3.9955E-01  2.1441E+00  9.7567E-01  9.1820E-01  4.3382E+02  1.0000E-02  1.7709E+00
             2.4681E+00
 PARAMETER:  9.3651E-02  8.4642E-01  1.4027E+01 -8.1742E-01  8.6271E-01  7.5372E-02  1.4665E-02  6.1726E+00 -1.0386E+01  6.7146E-01
             1.0035E+00
 GRADIENT:  -1.8635E+02  8.1756E+01  6.4499E-04  1.6471E+01  3.9568E-01 -1.5109E+01  2.6611E+00  5.6893E-04  0.0000E+00  8.2598E-01
             4.8490E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -2811.36121634302        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     2584
 NPARAMETR:  9.9346E-01  2.1078E+00  1.0601E+06  4.0009E-01  2.1436E+00  9.7580E-01  9.1809E-01  4.3408E+02  1.0000E-02  1.7705E+00
             2.4680E+00
 PARAMETER:  9.3441E-02  8.4564E-01  1.3974E+01 -8.1606E-01  8.6250E-01  7.5500E-02  1.4538E-02  6.1732E+00 -1.0386E+01  6.7124E-01
             1.0034E+00
 GRADIENT:  -1.1737E+02  6.8031E+01  1.5939E-04  1.2449E+01  1.8129E+00  1.0546E+01  4.3760E+00  2.5933E-04  0.0000E+00  1.9795E-01
             4.9921E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -2811.36169084706        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2773
 NPARAMETR:  9.9346E-01  2.1076E+00  1.0584E+06  4.0050E-01  2.1437E+00  9.7614E-01  9.1834E-01  4.3094E+02  1.0000E-02  1.7708E+00
             2.4682E+00
 PARAMETER:  9.3434E-02  8.4556E-01  1.3972E+01 -8.1503E-01  8.6251E-01  7.5847E-02  1.4812E-02  6.1660E+00 -1.0386E+01  6.7142E-01
             1.0035E+00
 GRADIENT:  -1.6072E+01  7.1135E+00  6.7281E-07  2.3935E+00 -7.7190E-01  2.3057E-01 -1.7353E+00  3.0585E-06  0.0000E+00  4.9254E-02
             2.5134E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -2811.36261113324        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:     2872
 NPARAMETR:  9.9355E-01  2.1079E+00  1.0563E+06  4.0066E-01  2.1442E+00  9.7810E-01  9.1845E-01  4.3138E+02  1.0000E-02  1.7716E+00
             2.4682E+00
 PARAMETER:  9.3533E-02  8.4569E-01  1.3970E+01 -8.1464E-01  8.6275E-01  7.7861E-02  1.4929E-02  6.1670E+00 -1.0386E+01  6.7190E-01
             1.0035E+00
 GRADIENT:  -2.9229E+01  2.5688E+01  9.7735E-05  2.2355E+00  2.8648E+00  8.9407E+00 -1.6371E+00  2.4422E-04  0.0000E+00  1.0920E+00
             4.7593E+00

0ITERATION NO.:   97    OBJECTIVE VALUE:  -2811.36261113324        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     2936
 NPARAMETR:  9.9354E-01  2.1074E+00  1.0511E+06  4.0068E-01  2.1440E+00  9.7811E-01  9.1843E-01  4.3405E+02  1.0000E-02  1.7714E+00
             2.4682E+00
 PARAMETER:  9.3533E-02  8.4569E-01  1.3970E+01 -8.1464E-01  8.6275E-01  7.7861E-02  1.4929E-02  6.1670E+00 -1.0386E+01  6.7190E-01
             1.0035E+00
 GRADIENT:   2.8695E-01  6.7841E-01  5.4685E-03 -3.2019E-02  5.8883E-02 -2.0990E-01  3.4233E-01 -6.3363E-03  0.0000E+00  1.2831E-01
            -7.1275E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2936
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         9.6937E-04 -6.3104E-03  1.3031E-06 -2.5247E-04 -1.3779E-02
 SE:             2.9347E-02  2.8698E-02  2.0675E-06  9.9651E-05  2.6951E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7365E-01  8.2596E-01  5.2851E-01  1.1293E-02  6.0918E-01

 ETASHRINKSD(%)  1.6823E+00  3.8590E+00  9.9993E+01  9.9666E+01  9.7095E+00
 ETASHRINKVR(%)  3.3363E+00  7.5691E+00  1.0000E+02  9.9999E+01  1.8476E+01
 EBVSHRINKSD(%)  1.6340E+00  3.9744E+00  9.9995E+01  9.9694E+01  7.9253E+00
 EBVSHRINKVR(%)  3.2413E+00  7.7909E+00  1.0000E+02  9.9999E+01  1.5223E+01
 RELATIVEINF(%)  9.6692E+01  6.1585E+00  1.4883E-07  6.2200E-05  5.2631E+01
 EPSSHRINKSD(%)  1.5921E+01
 EPSSHRINKVR(%)  2.9307E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          885
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1626.5212037722706     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2811.3626111332360     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1184.8414073609654     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    78.56
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.83
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2811.363       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.94E-01  2.11E+00  1.06E+06  4.01E-01  2.14E+00  9.78E-01  9.18E-01  4.31E+02  1.00E-02  1.77E+00  2.47E+00
 


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
+        3.65E+03
 
 TH 2
+       -7.44E+01  2.87E+02
 
 TH 3
+        2.09E-05 -1.27E-06  8.08E-14
 
 TH 4
+       -3.03E+01  4.23E+02 -2.23E-06  1.00E+03
 
 TH 5
+        2.82E+01 -1.08E+01  2.02E-07 -5.80E+00  7.26E+01
 
 TH 6
+       -1.14E+03 -1.12E+02  1.94E-06  1.07E+03  1.11E+02  5.87E+03
 
 TH 7
+        7.54E+02 -4.56E+01 -2.13E-07  2.37E+02  1.11E+02 -2.08E+03  4.12E+02
 
 TH 8
+       -2.99E-02 -3.68E-03  1.05E-09  4.47E-02  3.88E-04 -8.11E-02  2.64E-02  2.51E-06
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        1.20E+02  3.52E+00 -1.66E-06  4.69E+01 -1.23E+00  4.79E+01  7.44E+01  1.86E-02  0.00E+00  5.53E+01
 
 TH11
+        5.83E+01 -5.25E+00  9.30E-07 -2.58E+01  2.06E+00  4.23E+01 -1.57E+02 -3.83E-03  0.00E+00  1.28E+01  1.97E+02
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       92.469
Stop Time:
Sat Sep 25 02:42:09 CDT 2021
