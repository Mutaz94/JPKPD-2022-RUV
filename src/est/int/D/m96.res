Wed Sep 29 10:10:02 CDT 2021
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
$DATA ../../../../data/int/D/dat96.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   54838.4244971869        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.0646E+03  9.6773E+02  5.6398E+01  9.5141E+02 -4.4427E+01 -3.7063E+03 -2.2514E+03 -8.3945E+01 -2.6289E+03 -1.1920E+03
            -1.0605E+05

0ITERATION NO.:    5    OBJECTIVE VALUE:  -438.022077655620        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.5457E-01  2.7614E+00  8.5119E-01  2.1028E+00  9.0114E-01  7.8059E+00  7.6083E+00  9.8479E-01  2.3339E+00  1.3240E+00
             1.2067E+01
 PARAMETER:  5.3504E-02  1.1157E+00 -6.1118E-02  8.4329E-01 -4.0948E-03  2.1549E+00  2.1292E+00  8.4675E-02  9.4755E-01  3.8064E-01
             2.5905E+00
 GRADIENT:  -8.1364E+00  4.6947E+01 -6.7039E+01  1.7924E+02  1.3580E+00  2.6021E+02  1.0034E+02  4.7228E+00  2.3272E-01  2.1848E+01
            -5.3203E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -562.365127312566        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      161
 NPARAMETR:  1.1168E+00  3.1611E+00  1.9112E+01  4.3185E+00  2.8411E+00  5.2255E+00  9.3005E+00  9.3645E-01  5.4905E+00  9.0148E-01
             1.2312E+01
 PARAMETER:  2.1043E-01  1.2509E+00  3.0503E+00  1.5629E+00  1.1442E+00  1.7536E+00  2.3301E+00  3.4340E-02  1.8030E+00 -3.7197E-03
             2.6106E+00
 GRADIENT:   7.7737E+00  5.0145E+01 -2.0110E+01  1.2004E+02  5.0526E+01  1.8156E+02  1.5277E+02 -7.1158E-03  7.9395E+01  1.0430E+01
             1.8361E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -691.802699913000        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.3924E-01  2.1245E+00  9.6569E+00  1.2312E+00  2.1627E+00  2.1951E+00  4.2857E+00  1.4718E+00  2.1385E+00  7.8422E-01
             1.3225E+01
 PARAMETER:  3.7321E-02  8.5353E-01  2.3677E+00  3.0796E-01  8.7134E-01  8.8624E-01  1.5553E+00  4.8649E-01  8.6012E-01 -1.4306E-01
             2.6821E+00
 GRADIENT:  -4.6444E+01  2.2827E+01 -1.7133E+00  5.0018E+01 -1.3708E+01  2.3544E+01 -8.1576E+01  2.4867E-01  6.6401E+00  9.3743E+00
             1.2607E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -728.334466294086        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      308
 NPARAMETR:  1.0477E+00  1.4741E+00  1.8097E+01  1.0143E+00  2.4832E+00  2.3433E+00  5.1017E+00  2.2490E+00  1.3162E+00  2.6697E-01
             1.2672E+01
 PARAMETER:  1.4661E-01  4.8801E-01  2.9957E+00  1.1422E-01  1.0095E+00  9.5154E-01  1.7296E+00  9.1047E-01  3.7478E-01 -1.2206E+00
             2.6394E+00
 GRADIENT:  -5.1247E-01 -7.7625E+00 -1.1452E+00 -5.1295E+00  2.4263E+01  3.4179E+01  2.0474E+01  1.4496E-01  7.6134E+00  1.0844E+00
             9.3084E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -734.345967654115        NO. OF FUNC. EVALS.: 122
 CUMULATIVE NO. OF FUNC. EVALS.:      430
 NPARAMETR:  1.0550E+00  1.5053E+00  9.3296E+00  9.7724E-01  2.2151E+00  2.1482E+00  5.1793E+00  4.2231E-01  1.2857E+00  8.0497E-02
             1.2480E+01
 PARAMETER:  1.5351E-01  5.0901E-01  2.3332E+00  7.6975E-02  8.9531E-01  8.6461E-01  1.7447E+00 -7.6201E-01  3.5134E-01 -2.4195E+00
             2.6241E+00
 GRADIENT:   1.2212E+01 -4.3081E+00 -6.3522E-01 -6.2722E+00  3.4140E+00  1.8683E+01  3.4619E+01  1.7991E-02  8.2650E+00  1.0308E-01
             6.0195E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -744.710655537653        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      613
 NPARAMETR:  1.0540E+00  1.5119E+00  1.1644E+01  9.7801E-01  2.2360E+00  2.1417E+00  6.4293E+00  4.9524E-02  1.2813E+00  2.0222E-02
             1.2344E+01
 PARAMETER:  1.5259E-01  5.1339E-01  2.5548E+00  7.7761E-02  9.0468E-01  8.6161E-01  1.9609E+00 -2.9053E+00  3.4789E-01 -3.8010E+00
             2.6132E+00
 GRADIENT:   9.2406E+00 -8.7989E-01 -4.3579E-01 -3.2614E+01  1.2910E+00  1.6725E-01  1.3669E+00  1.5111E-04  1.1743E+01  6.3499E-03
             2.4266E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -744.860211163213        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      793
 NPARAMETR:  1.0501E+00  1.5177E+00  1.2021E+01  9.8287E-01  2.3012E+00  2.1361E+00  6.3603E+00  1.2822E-02  1.2537E+00  1.0000E-02
             1.2269E+01
 PARAMETER:  1.4892E-01  5.1719E-01  2.5866E+00  8.2717E-02  9.3342E-01  8.5896E-01  1.9501E+00 -4.2566E+00  3.2613E-01 -5.1214E+00
             2.6071E+00
 GRADIENT:   8.6385E+00 -8.8466E-01 -1.0598E+00 -2.8995E+01  9.1994E+00 -8.1253E-01 -2.5329E+00  1.0312E-05  1.0855E+01  0.0000E+00
            -8.9818E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -746.908792049943        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      971
 NPARAMETR:  1.0261E+00  1.5584E+00  4.2893E+03  1.0143E+00  2.5013E+00  2.1025E+00  6.5125E+00  1.6192E-01  1.0936E+00  1.0000E-02
             1.2394E+01
 PARAMETER:  1.2572E-01  5.4368E-01  8.4639E+00  1.1415E-01  1.0168E+00  8.4313E-01  1.9737E+00 -1.7207E+00  1.8944E-01 -1.3968E+01
             2.6172E+00
 GRADIENT:  -1.5150E+00  2.0300E+00  1.2915E-03 -2.4744E+01  6.3508E+00 -3.3739E+00 -1.7337E+00 -3.1029E-07  5.8271E+00  0.0000E+00
             1.1978E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -748.829252785265        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1157             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0141E+00  1.5568E+00  4.5991E+06  1.1210E+00  2.4815E+00  2.0884E+00  6.7478E+00  3.4994E+02  9.9780E-01  1.0000E-02
             1.2446E+01
 PARAMETER:  1.1404E-01  5.4265E-01  1.5441E+01  2.1421E-01  1.0089E+00  8.3638E-01  2.0092E+00  5.9578E+00  9.7797E-02 -2.0157E+01
             2.6214E+00
 GRADIENT:  -9.5975E+00  1.1353E+01  3.5579E-05  4.6433E+00  6.4060E-01  9.4228E+00  1.2736E+02 -7.0281E-05 -2.4465E+00  0.0000E+00
             3.7821E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -751.306762851395        NO. OF FUNC. EVALS.: 167
 CUMULATIVE NO. OF FUNC. EVALS.:     1324
 NPARAMETR:  1.0273E+00  1.3271E+00  4.3285E+06  1.1710E+00  2.4673E+00  2.1259E+00  7.5084E+00  3.4130E+02  1.0012E+00  1.0000E-02
             1.2480E+01
 PARAMETER:  1.2690E-01  3.8297E-01  1.5381E+01  2.5789E-01  1.0031E+00  8.5419E-01  2.1160E+00  5.9328E+00  1.0116E-01 -2.0157E+01
             2.6241E+00
 GRADIENT:  -9.9094E+00  1.0115E+01 -4.2664E-04  7.9155E+00 -5.7189E+00  1.9318E-02  1.1609E+01 -1.7096E-03 -2.1773E+00  0.0000E+00
            -1.1493E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -752.440902145395        NO. OF FUNC. EVALS.: 199
 CUMULATIVE NO. OF FUNC. EVALS.:     1523
 NPARAMETR:  1.0379E+00  9.5437E-01  5.4770E+06  1.2183E+00  2.5048E+00  2.1164E+00  7.7740E+00  1.1136E+01  1.0215E+00  1.0000E-02
             1.2575E+01
 PARAMETER:  1.3723E-01  5.3294E-02  1.5616E+01  2.9749E-01  1.0182E+00  8.4972E-01  2.1508E+00  2.5102E+00  1.2124E-01 -2.0157E+01
             2.6317E+00
 GRADIENT:  -2.9416E+01  2.2718E+00  1.0571E-06 -9.3823E+00  3.9328E+00 -5.1106E+00  1.5519E+01  2.2083E-05 -4.1362E+00  0.0000E+00
             1.6775E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -753.553567213422        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1701
 NPARAMETR:  1.0577E+00  7.7321E-01  5.6756E+06  1.4171E+00  2.4516E+00  2.1685E+00  7.7708E+00  3.7697E+00  1.2292E+00  1.0000E-02
             1.2567E+01
 PARAMETER:  1.5610E-01 -1.5721E-01  1.5652E+01  4.4859E-01  9.9673E-01  8.7404E-01  2.1504E+00  1.4270E+00  3.0639E-01 -2.0157E+01
             2.6311E+00
 GRADIENT:  -9.8672E+01  7.3637E+00  7.6012E-06  4.7821E+01 -7.9291E+00  1.9036E+01 -9.1233E+00 -5.1103E-05 -1.8656E+01  0.0000E+00
             1.4857E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -754.398340005382        NO. OF FUNC. EVALS.: 149
 CUMULATIVE NO. OF FUNC. EVALS.:     1850
 NPARAMETR:  1.0797E+00  7.9775E-01  4.2913E+08  1.4765E+00  2.4832E+00  2.1517E+00  7.8000E+00  2.2704E+00  1.4599E+00  1.0000E-02
             1.2581E+01
 PARAMETER:  1.7665E-01 -1.2595E-01  1.9977E+01  4.8969E-01  1.0095E+00  8.6627E-01  2.1541E+00  9.1994E-01  4.7836E-01 -2.0157E+01
             2.6322E+00
 GRADIENT:   1.4355E+02 -2.1948E+01  2.0353E-04 -3.7164E+01  1.0059E+01  3.4145E+01  1.7837E+02  2.5830E-02 -4.1921E+00  0.0000E+00
             2.9828E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -754.674729029713        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2034
 NPARAMETR:  1.0614E+00  7.9628E-01  1.7297E+04  1.4202E+00  2.4764E+00  2.1131E+00  7.8495E+00  2.2713E+00  1.4744E+00  1.0000E-02
             1.2526E+01
 PARAMETER:  1.5954E-01 -1.2781E-01  9.8583E+00  4.5083E-01  1.0068E+00  8.4815E-01  2.1604E+00  9.2033E-01  4.8822E-01 -2.0157E+01
             2.6278E+00
 GRADIENT:   8.8441E+00 -6.0679E-01 -1.0976E-04 -8.6086E+00  4.9614E-01 -1.3814E-01  3.1197E+00 -1.0223E-06  2.6654E+00  0.0000E+00
            -5.4639E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -755.682327085476        NO. OF FUNC. EVALS.: 188
 CUMULATIVE NO. OF FUNC. EVALS.:     2222             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0230E+00  8.0364E-01  7.2240E+05  1.4715E+00  2.4785E+00  2.1175E+00  9.1042E+00  2.2514E+00  1.4040E+00  1.0000E-02
             1.2549E+01
 PARAMETER:  1.2273E-01 -1.1861E-01  1.3590E+01  4.8628E-01  1.0076E+00  8.5025E-01  2.3087E+00  9.1155E-01  4.3929E-01 -2.0157E+01
             2.6296E+00
 GRADIENT:  -1.1472E+01 -7.4343E-01 -8.5203E-05  5.9458E-01  6.8034E+00  1.4572E+01  2.4836E+02  1.4234E-04  8.6396E+00  0.0000E+00
             5.9124E+01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -756.351877182120        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     2345
 NPARAMETR:  1.0512E+00  7.3026E-01  7.1804E+05  1.5098E+00  2.4876E+00  2.1070E+00  9.0273E+00  2.2359E+00  1.4765E+00  1.0000E-02
             1.2540E+01
 PARAMETER:  1.4989E-01 -2.1435E-01  1.3584E+01  5.1197E-01  1.0113E+00  8.4527E-01  2.3003E+00  9.0463E-01  4.8971E-01 -2.0157E+01
             2.6289E+00
 GRADIENT:  -2.0000E+01  6.0401E+00 -1.9800E-06  1.5362E+01  3.4539E+00  1.4293E+01  2.4236E+02 -2.1940E-05  5.6435E-01  0.0000E+00
             5.2117E+01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -756.969614012890        NO. OF FUNC. EVALS.: 203
 CUMULATIVE NO. OF FUNC. EVALS.:     2548
 NPARAMETR:  1.0391E+00  6.0058E-01  6.6077E+05  1.5453E+00  2.4816E+00  2.1084E+00  9.0603E+00  2.2624E+00  1.4879E+00  1.0000E-02
             1.2544E+01
 PARAMETER:  1.3832E-01 -4.0986E-01  1.3501E+01  5.3523E-01  1.0089E+00  8.4594E-01  2.3039E+00  9.1642E-01  4.9736E-01 -2.0157E+01
             2.6293E+00
 GRADIENT:  -1.9415E+01  1.9876E+00 -9.8695E-06  5.8158E-01  1.4715E-01  2.6385E+00  2.4764E+01 -9.1783E-07 -2.9413E+00  0.0000E+00
             4.5781E+00

0ITERATION NO.:   90    OBJECTIVE VALUE:  -757.436771190049        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     2742             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0218E+00  5.1250E-01  8.1372E+05  1.6397E+00  2.4817E+00  2.1116E+00  9.3968E+00  2.2925E+00  1.5894E+00  1.0000E-02
             1.2584E+01
 PARAMETER:  1.2153E-01 -5.6846E-01  1.3709E+01  5.9451E-01  1.0089E+00  8.4747E-01  2.3404E+00  9.2965E-01  5.6339E-01 -2.0157E+01
             2.6324E+00
 GRADIENT:  -3.1106E+01  5.5878E+00 -2.5215E-06  2.9052E+01  1.5302E+00  1.9102E+01  2.6993E+02  6.2695E-05 -1.2541E+00  0.0000E+00
             5.9237E+01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -757.616776664765        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2924             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0388E+00  4.8365E-01  5.8102E+05  1.6356E+00  2.4839E+00  2.1036E+00  9.7605E+00  2.2902E+00  1.6185E+00  1.0000E-02
             1.2564E+01
 PARAMETER:  1.3809E-01 -6.2640E-01  1.3373E+01  5.9200E-01  1.0098E+00  8.4366E-01  2.3783E+00  9.2865E-01  5.8148E-01 -2.0157E+01
             2.6308E+00
 GRADIENT:  -3.3710E+01  3.9241E+00 -4.0840E-05  1.1153E+01  4.4565E+00  2.4866E+01  2.9035E+02  6.6638E-05  2.9728E+00  0.0000E+00
             6.5848E+01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -757.727119608600        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     3113
 NPARAMETR:  1.0382E+00  4.3209E-01  5.4358E+05  1.6586E+00  2.4820E+00  2.1039E+00  9.7856E+00  2.3062E+00  1.6260E+00  1.0000E-02
             1.2563E+01
 PARAMETER:  1.3752E-01 -7.3911E-01  1.3306E+01  6.0600E-01  1.0091E+00  8.4377E-01  2.3809E+00  9.3559E-01  5.8611E-01 -2.0157E+01
             2.6308E+00
 GRADIENT:  -3.5152E+01  1.5532E+00 -2.1955E-05  7.7957E+00  1.5128E-01  1.1111E+00  3.2637E+01 -6.6327E-06 -1.5399E+00  0.0000E+00
             1.0388E+01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -757.805055324058        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     3307             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0378E+00  4.2397E-01  5.1759E+05  1.6867E+00  2.4794E+00  2.1033E+00  9.9617E+00  2.3003E+00  1.6349E+00  1.0000E-02
             1.2558E+01
 PARAMETER:  1.3706E-01 -7.5809E-01  1.3257E+01  6.2278E-01  1.0080E+00  8.4353E-01  2.3988E+00  9.3303E-01  5.9161E-01 -2.0157E+01
             2.6303E+00
 GRADIENT:  -1.9485E+01  3.1189E+00  7.5381E-06  2.0790E+01  4.1046E+00  1.2964E+01  3.0383E+02  2.2678E-06  4.0893E+00  0.0000E+00
             5.5593E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -757.828577462177        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:     3494
 NPARAMETR:  1.0380E+00  4.1053E-01  5.9748E+05  1.6958E+00  2.4795E+00  2.1026E+00  1.0058E+01  2.3349E+00  1.6440E+00  1.0000E-02
             1.2560E+01
 PARAMETER:  1.3734E-01 -7.9030E-01  1.3400E+01  6.2813E-01  1.0080E+00  8.4318E-01  2.4084E+00  9.4795E-01  5.9713E-01 -2.0157E+01
             2.6305E+00
 GRADIENT:  -2.6741E+01  1.3104E+00 -9.2563E-06  9.3087E+00 -7.5518E-01  2.5631E+00  3.6416E+01 -1.1767E-05 -2.2909E+00  0.0000E+00
             4.9752E+00

0ITERATION NO.:  115    OBJECTIVE VALUE:  -757.837805799557        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3677
 NPARAMETR:  1.0394E+00  4.0277E-01  6.0421E+05  1.6991E+00  2.4819E+00  2.0994E+00  1.0156E+01  2.3407E+00  1.6571E+00  1.0000E-02
             1.2573E+01
 PARAMETER:  1.3868E-01 -8.0939E-01  1.3412E+01  6.3011E-01  1.0090E+00  8.4165E-01  2.4180E+00  9.5044E-01  6.0504E-01 -2.0157E+01
             2.6315E+00
 GRADIENT:  -1.2381E+01 -1.2162E+00 -1.8370E-05 -5.1772E+00  1.9000E+00  1.7743E+00  3.8997E+01  9.3892E-06  1.9169E+00  0.0000E+00
             8.1208E+00

0ITERATION NO.:  120    OBJECTIVE VALUE:  -757.847917652961        NO. OF FUNC. EVALS.: 197
 CUMULATIVE NO. OF FUNC. EVALS.:     3874             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0380E+00  3.9226E-01  7.4494E+05  1.7088E+00  2.4801E+00  2.1018E+00  1.0191E+01  1.9858E+00  1.6585E+00  1.0000E-02
             1.2566E+01
 PARAMETER:  1.3728E-01 -8.3583E-01  1.3621E+01  6.3582E-01  1.0083E+00  8.4278E-01  2.4215E+00  7.8600E-01  6.0591E-01 -2.0157E+01
             2.6310E+00
 GRADIENT:  -4.9794E+02  1.5282E+01 -4.0678E-05  1.2768E+02  7.2235E+00  1.0088E+02  3.0914E+02  6.8563E-05 -2.4046E+01  0.0000E+00
             1.9636E+02

0ITERATION NO.:  125    OBJECTIVE VALUE:  -757.849575342419        NO. OF FUNC. EVALS.: 191
 CUMULATIVE NO. OF FUNC. EVALS.:     4065
 NPARAMETR:  1.0382E+00  3.8609E-01  9.3591E+05  1.7121E+00  2.4797E+00  2.1016E+00  1.0217E+01  2.0042E+00  1.6601E+00  1.0000E-02
             1.2565E+01
 PARAMETER:  1.3748E-01 -8.5167E-01  1.3849E+01  6.3769E-01  1.0081E+00  8.4271E-01  2.4240E+00  7.9526E-01  6.0690E-01 -2.0157E+01
             2.6309E+00
 GRADIENT:  -1.0680E+01  1.7681E-01 -1.1275E-05  2.0036E+00  7.1608E-02  5.9048E-01  3.9445E+01  1.0247E-05 -1.1495E+00  0.0000E+00
             2.0172E+00

0ITERATION NO.:  130    OBJECTIVE VALUE:  -757.850471727860        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     4247
 NPARAMETR:  1.0380E+00  3.8457E-01  9.5178E+05  1.7126E+00  2.4811E+00  2.1015E+00  1.0257E+01  2.0419E+00  1.6658E+00  1.0000E-02
             1.2570E+01
 PARAMETER:  1.3725E-01 -8.5528E-01  1.3848E+01  6.3879E-01  1.0082E+00  8.4265E-01  2.4263E+00  8.0581E-01  6.0789E-01 -2.0157E+01
             2.6310E+00
 GRADIENT:  -2.2508E-02  4.8532E-03 -2.3512E-04  3.4215E-01 -1.8685E-01 -5.1906E-04 -3.2865E-01 -3.4000E-02 -2.8071E-01  0.0000E+00
            -5.0759E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4247
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5508E-02  6.4436E-02  1.3408E-08 -9.8626E-02 -1.0527E-05
 SE:             2.7836E-02  2.1219E-02  3.3495E-08  1.6478E-02  8.8843E-05
 N:                     100         100         100         100         100

 P VAL.:         5.7744E-01  2.3920E-03  6.8895E-01  2.1683E-09  9.0568E-01

 ETASHRINKSD(%)  6.7470E+00  2.8913E+01  1.0000E+02  4.4796E+01  9.9702E+01
 ETASHRINKVR(%)  1.3039E+01  4.9467E+01  1.0000E+02  6.9525E+01  9.9999E+01
 EBVSHRINKSD(%)  8.5802E+00  2.6887E+01  1.0000E+02  3.7610E+01  9.9604E+01
 EBVSHRINKVR(%)  1.6424E+01  4.6545E+01  1.0000E+02  6.1074E+01  9.9998E+01
 RELATIVEINF(%)  8.3329E+01  3.3983E+01  0.0000E+00  2.4199E+01  3.5014E-04
 EPSSHRINKSD(%)  4.4831E+00
 EPSSHRINKVR(%)  8.7653E+00

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -757.85047172786028     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       896.23888804055048     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   161.77
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    18.58
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -757.850       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.04E+00  3.85E-01  9.35E+05  1.71E+00  2.48E+00  2.10E+00  1.02E+01  2.03E+00  1.66E+00  1.00E-02  1.26E+01
 


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
+        1.79E+02
 
 TH 2
+       -9.78E+00  5.03E+01
 
 TH 3
+        2.76E-07 -7.01E-08  8.50E-16
 
 TH 4
+        1.98E+00  1.62E+01 -6.13E-09  7.66E+01
 
 TH 5
+       -3.62E+00 -3.02E+00  6.24E-10 -8.78E+00  2.65E+01
 
 TH 6
+       -4.28E+00  2.04E+00 -3.94E-09  1.12E+00  4.18E-01  2.85E+01
 
 TH 7
+        3.38E-01  3.97E+00  1.04E-09 -2.81E+00  1.76E-01  6.30E-02  9.40E-01
 
 TH 8
+        5.43E+00 -1.91E-01 -3.32E-09 -2.85E-01 -6.64E-04 -1.70E-01 -7.41E-03  5.97E-02
 
 TH 9
+        1.72E+00 -4.97E-01 -1.85E-09 -2.22E+01  3.25E+00 -1.14E+00  5.39E-01 -5.44E-01  2.09E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.23E+00 -1.59E+00  7.91E-11 -5.64E+00  2.46E-01  1.69E+00  7.14E-02  8.31E-04  2.28E+00  0.00E+00  5.08E+00
 
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
 #CPUT: Total CPU Time in Seconds,      180.434
Stop Time:
Wed Sep 29 10:13:04 CDT 2021
