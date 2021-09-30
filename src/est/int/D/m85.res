Wed Sep 29 09:52:46 CDT 2021
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
$DATA ../../../../data/int/D/dat85.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m85.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   42573.9432930132        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.2715E+02  7.2317E+02 -5.4138E+00  6.0203E+02  1.0242E+02 -3.0380E+03 -1.7620E+03 -4.5534E+01 -2.5997E+03 -7.2467E+02
            -8.3505E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -578.399467848238        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  8.6945E-01  2.3103E+00  8.8503E-01  2.1670E+00  9.5313E-01  4.0003E+00  4.5246E+00  9.7069E-01  2.0749E+00  1.3303E+00
             1.2891E+01
 PARAMETER: -3.9898E-02  9.3737E-01 -2.2137E-02  8.7334E-01  5.1999E-02  1.4864E+00  1.6095E+00  7.0248E-02  8.2993E-01  3.8544E-01
             2.6566E+00
 GRADIENT:  -4.8458E+01  6.6163E+01 -4.4664E+01  1.9177E+02 -1.1903E+01  1.5810E+02  2.9744E+01  4.7955E+00 -1.7144E+01  2.6707E+01
             2.0588E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -681.275443437228        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  7.7792E-01  3.2564E+00  1.3238E+01  1.8436E+00  2.0243E+00  2.9479E+00  8.6895E+00  5.5888E-01  1.3772E+00  1.2623E+00
             1.3140E+01
 PARAMETER: -1.5113E-01  1.2806E+00  2.6831E+00  7.1172E-01  8.0523E-01  1.1811E+00  2.2621E+00 -4.8182E-01  4.2008E-01  3.3296E-01
             2.6757E+00
 GRADIENT:  -9.0312E+01  5.3218E+01  3.0385E-01  7.6101E+01 -3.3913E+01  9.5561E+01  4.7464E+01  9.8790E-02  6.0356E+00  2.4406E+01
             3.0379E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -786.203885796685        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.1768E+00  1.0997E+00  8.9345E+00  1.2207E+00  2.0772E+00  1.9191E+00  4.1583E+00  3.8946E-01  1.1630E+00  7.7077E-01
             1.2124E+01
 PARAMETER:  2.6277E-01  1.9505E-01  2.2899E+00  2.9940E-01  8.3101E-01  7.5184E-01  1.5251E+00 -8.4300E-01  2.5102E-01 -1.6037E-01
             2.5952E+00
 GRADIENT:   2.1966E+01 -2.2598E+01 -6.9498E-01 -9.5803E+00 -5.3174E+00 -4.7302E+00 -5.7244E+00  5.4037E-02  9.4115E+00  1.1006E+01
             2.2329E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -803.012659006414        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      297
 NPARAMETR:  1.0369E+00  1.2433E+00  1.1124E+01  1.0544E+00  2.2897E+00  1.9948E+00  4.0557E+00  2.2257E-01  1.2216E+00  3.8186E-01
             1.0689E+01
 PARAMETER:  1.3623E-01  3.1775E-01  2.5091E+00  1.5293E-01  9.2844E-01  7.9052E-01  1.5001E+00 -1.4025E+00  3.0020E-01 -8.6269E-01
             2.4692E+00
 GRADIENT:  -1.8195E+01 -1.6979E+01 -1.4163E+00 -1.1967E+01  1.2565E+01  1.2010E-02  1.1480E+01  7.1723E-03  8.5568E+00  2.8692E+00
             6.7497E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -803.680537598143        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  1.0527E+00  1.5288E+00  8.5745E+00  9.2290E-01  2.2420E+00  2.0025E+00  3.7593E+00  1.9657E-01  1.0635E+00  2.5920E-01
             1.0481E+01
 PARAMETER:  1.5131E-01  5.2451E-01  2.2488E+00  1.9770E-02  9.0736E-01  7.9438E-01  1.4242E+00 -1.5267E+00  1.6153E-01 -1.2502E+00
             2.4496E+00
 GRADIENT:  -6.7831E+00 -5.0918E+00 -3.5738E-01 -5.2459E+00  3.8636E+00  4.1158E-01  3.8079E+00  6.4729E-03  4.7456E+00  1.3401E+00
             2.4942E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -803.738618307336        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      437
 NPARAMETR:  1.0603E+00  1.6140E+00  7.8500E+00  8.7874E-01  2.2287E+00  2.0046E+00  3.6821E+00  1.6209E-01  9.4544E-01  2.0215E-01
             1.0420E+01
 PARAMETER:  1.5851E-01  5.7875E-01  2.1605E+00 -2.9267E-02  9.0142E-01  7.9542E-01  1.4035E+00 -1.7196E+00  4.3891E-02 -1.4987E+00
             2.4438E+00
 GRADIENT:  -1.9118E+00 -1.4870E+00  1.0611E-01 -3.5979E+00  1.1451E+00  5.1418E-01  1.5188E+00  4.7458E-03  2.9446E+00  8.1570E-01
             1.0968E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -803.847550635744        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  1.0659E+00  1.6706E+00  7.0123E+00  8.4494E-01  2.2094E+00  2.0050E+00  3.6351E+00  1.0621E-01  7.7542E-01  1.3003E-01
             1.0382E+01
 PARAMETER:  1.6384E-01  6.1318E-01  2.0477E+00 -6.8493E-02  8.9273E-01  7.9565E-01  1.3906E+00 -2.1424E+00 -1.5435E-01 -1.9400E+00
             2.4401E+00
 GRADIENT:   1.5954E+00  9.7117E-01  5.1617E-01 -1.4694E+00 -9.1552E-01  4.6671E-01 -3.9610E-01  2.4776E-03  9.6620E-01  3.3717E-01
             3.7225E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -803.910216672615        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      579
 NPARAMETR:  1.0659E+00  1.6842E+00  6.1403E+00  8.2897E-01  2.1832E+00  2.0036E+00  3.6238E+00  6.6126E-02  6.4207E-01  8.1158E-02
             1.0385E+01
 PARAMETER:  1.6385E-01  6.2128E-01  1.9149E+00 -8.7566E-02  8.8078E-01  7.9496E-01  1.3875E+00 -2.6162E+00 -3.4306E-01 -2.4114E+00
             2.4404E+00
 GRADIENT:   1.5510E+00  1.0257E+00  5.8949E-01 -5.3938E-02 -1.1166E+00  3.1637E-01 -8.0956E-01  1.2927E-03 -9.7706E-02  1.3210E-01
            -1.5014E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -806.759066302797        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      742
 NPARAMETR:  1.0665E+00  1.6988E+00  6.4047E+00  8.4344E-01  2.2126E+00  2.0223E+00  3.8138E+00  1.3641E-02  5.8369E-01  1.2041E-02
             1.0672E+01
 PARAMETER:  1.6442E-01  6.2991E-01  1.9570E+00 -7.0272E-02  8.9418E-01  8.0424E-01  1.4386E+00 -4.1947E+00 -4.3838E-01 -4.3195E+00
             2.4676E+00
 GRADIENT:  -2.8060E+00  2.8857E+00 -2.8683E-03 -6.1404E+00  3.4273E+00  4.8473E+00  1.7536E+01  5.6428E-05  5.9930E-02  2.9820E-03
             5.0026E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -808.236149402406        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      920
 NPARAMETR:  1.0658E+00  1.7125E+00  6.4067E+00  8.4405E-01  2.2214E+00  2.0316E+00  4.1604E+00  1.0000E-02  5.2414E-01  1.0000E-02
             1.0614E+01
 PARAMETER:  1.6371E-01  6.3798E-01  1.9573E+00 -6.9542E-02  8.9813E-01  8.0883E-01  1.5256E+00 -5.7980E+00 -5.4600E-01 -5.8577E+00
             2.4621E+00
 GRADIENT:  -9.6845E+00 -2.7996E+00 -5.3244E-01 -1.5473E+01  2.7122E+00 -1.5376E+01 -4.9175E-01  0.0000E+00 -6.4503E-02  0.0000E+00
             7.3788E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -809.233849916701        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1106             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0674E+00  1.7225E+00  6.1184E+00  9.0589E-01  2.1833E+00  2.1136E+00  4.2064E+00  1.0000E-02  5.4146E-01  1.0000E-02
             1.0637E+01
 PARAMETER:  1.6525E-01  6.4378E-01  1.9113E+00  1.1623E-03  8.8082E-01  8.4840E-01  1.5366E+00 -5.0184E+00 -5.1349E-01 -4.7937E+00
             2.4643E+00
 GRADIENT:  -1.0437E+00  1.5187E+01 -6.2034E-01  2.5944E+00  2.2041E+00  2.0366E+01  4.0556E+01  0.0000E+00 -1.2038E+00  0.0000E+00
             3.6606E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -809.497836693842        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1244
 NPARAMETR:  1.0676E+00  1.7201E+00  7.3613E+00  9.4110E-01  2.1852E+00  2.1152E+00  4.2162E+00  1.0000E-02  6.9179E-01  1.0000E-02
             1.0712E+01
 PARAMETER:  1.6539E-01  6.4236E-01  2.0962E+00  3.9295E-02  8.8171E-01  8.4913E-01  1.5389E+00 -5.0184E+00 -2.6847E-01 -4.7937E+00
             2.4714E+00
 GRADIENT:  -1.0542E+01  8.0712E+00  2.9830E-01  2.0279E+00 -7.9609E+00 -9.4420E-01 -1.0916E+01  0.0000E+00 -4.0232E-01  0.0000E+00
             3.2650E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -810.479200042660        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1422
 NPARAMETR:  1.0706E+00  1.6643E+00  1.0692E+01  1.0232E+00  2.2660E+00  2.1470E+00  4.6829E+00  1.0000E-02  7.8324E-01  1.0000E-02
             1.0693E+01
 PARAMETER:  1.6824E-01  6.0940E-01  2.4695E+00  1.2294E-01  9.1802E-01  8.6407E-01  1.6439E+00 -5.0184E+00 -1.4432E-01 -4.7937E+00
             2.4696E+00
 GRADIENT:  -8.5089E+00  1.3287E+01  5.8967E-02  1.2974E+00 -4.9950E+00  3.3800E+00  4.7751E+00  0.0000E+00 -1.4791E-01  0.0000E+00
            -1.1012E+00

0ITERATION NO.:   70    OBJECTIVE VALUE:  -811.553899141605        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1601
 NPARAMETR:  1.0740E+00  1.4713E+00  7.6541E+01  1.1076E+00  2.5400E+00  2.1098E+00  4.8095E+00  1.0000E-02  8.9153E-01  1.0000E-02
             1.0718E+01
 PARAMETER:  1.7136E-01  4.8618E-01  4.4378E+00  2.0217E-01  1.0322E+00  8.4661E-01  1.6706E+00 -5.0184E+00 -1.4819E-02 -4.7937E+00
             2.4719E+00
 GRADIENT:  -7.8992E+00  8.1150E+00 -7.0857E-03  3.2317E-01  8.3096E+00 -1.7748E+00 -1.7194E-01  0.0000E+00 -6.2394E-03  0.0000E+00
             1.6381E+00

0ITERATION NO.:   75    OBJECTIVE VALUE:  -812.532043777116        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1738
 NPARAMETR:  1.0961E+00  1.2737E+00  2.8950E+01  1.1341E+00  2.4483E+00  2.1193E+00  4.8491E+00  1.0000E-02  8.8321E-01  1.0000E-02
             1.0664E+01
 PARAMETER:  1.9173E-01  3.4190E-01  3.4656E+00  2.2587E-01  9.9538E-01  8.5107E-01  1.6788E+00 -5.0184E+00 -2.4188E-02 -4.7937E+00
             2.4669E+00
 GRADIENT:   1.6417E+00  2.6471E-01 -9.2423E-02 -2.6831E+00  6.3760E+00 -2.8666E-01 -6.3713E+00  0.0000E+00 -7.8720E-01  0.0000E+00
            -6.6452E+00

0ITERATION NO.:   80    OBJECTIVE VALUE:  -813.594100010397        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1914
 NPARAMETR:  1.0979E+00  1.0444E+00  2.1716E+01  1.2730E+00  2.3109E+00  2.1327E+00  5.4353E+00  1.0000E-02  1.0465E+00  1.0000E-02
             1.0741E+01
 PARAMETER:  1.9339E-01  1.4346E-01  3.1780E+00  3.4135E-01  9.3764E-01  8.5738E-01  1.7929E+00 -5.0184E+00  1.4541E-01 -4.7937E+00
             2.4741E+00
 GRADIENT:   1.3628E+00  1.3879E+00  1.8124E-01 -8.0846E-01 -6.6191E+00  2.7702E+00  1.3548E+00  0.0000E+00  1.0984E+00  0.0000E+00
             5.3417E+00

0ITERATION NO.:   85    OBJECTIVE VALUE:  -813.813763247383        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     2091
 NPARAMETR:  1.0936E+00  9.8461E-01  3.4206E+01  1.3021E+00  2.4087E+00  2.1148E+00  5.5218E+00  1.0000E-02  1.0474E+00  1.0000E-02
             1.0714E+01
 PARAMETER:  1.8946E-01  8.4486E-02  3.6324E+00  3.6400E-01  9.7910E-01  8.4895E-01  1.8087E+00 -5.0184E+00  1.4628E-01 -4.7937E+00
             2.4715E+00
 GRADIENT:   6.9872E-02  9.5273E-03 -3.6367E-02  1.8156E-01  4.7083E-01  1.7111E-02 -1.1078E-01  0.0000E+00  8.1976E-02  0.0000E+00
             1.8119E-02

0ITERATION NO.:   90    OBJECTIVE VALUE:  -813.878661473778        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2267
 NPARAMETR:  1.0943E+00  9.6162E-01  3.5666E+01  1.3119E+00  2.4100E+00  2.1204E+00  5.6196E+00  1.0000E-02  1.0371E+00  1.0000E-02
             1.0725E+01
 PARAMETER:  1.9016E-01  6.0866E-02  3.6742E+00  3.7146E-01  9.7963E-01  8.5161E-01  1.8263E+00 -5.0184E+00  1.3645E-01 -4.7937E+00
             2.4725E+00
 GRADIENT:   2.6437E-01  3.6119E-02 -1.6282E-02 -2.6500E-01  1.5263E-01  9.6991E-01  1.8001E+00  0.0000E+00 -3.4475E-01  0.0000E+00
             1.5020E+00

0ITERATION NO.:   95    OBJECTIVE VALUE:  -813.918907914737        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2449
 NPARAMETR:  1.0943E+00  9.4336E-01  3.7176E+01  1.3217E+00  2.4101E+00  2.1237E+00  5.7190E+00  1.0000E-02  1.0389E+00  1.0000E-02
             1.0720E+01
 PARAMETER:  1.9014E-01  4.1695E-02  3.7157E+00  3.7892E-01  9.7968E-01  8.5318E-01  1.8438E+00 -5.0184E+00  1.3814E-01 -4.7937E+00
             2.4721E+00
 GRADIENT:   3.9835E-01  3.4928E-01  1.8897E-03 -7.6694E-01 -3.0029E-01  1.4297E+00  4.0317E+00  0.0000E+00 -4.3538E-01  0.0000E+00
             7.7583E-01

0ITERATION NO.:  100    OBJECTIVE VALUE:  -813.938078744093        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     2629
 NPARAMETR:  1.0939E+00  9.3246E-01  3.8625E+01  1.3320E+00  2.4132E+00  2.1205E+00  5.7768E+00  1.0000E-02  1.0645E+00  1.0000E-02
             1.0714E+01
 PARAMETER:  1.8976E-01  3.0075E-02  3.7539E+00  3.8665E-01  9.8095E-01  8.5163E-01  1.8539E+00 -5.0184E+00  1.6251E-01 -4.7937E+00
             2.4715E+00
 GRADIENT:   3.4854E-01  6.3865E-01 -1.5700E-02 -1.5699E+00 -9.0725E-02  8.7560E-01  5.4118E+00  0.0000E+00  2.8274E-01  0.0000E+00
             6.5225E-01

0ITERATION NO.:  105    OBJECTIVE VALUE:  -813.954980784766        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2811             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0938E+00  9.1446E-01  4.0415E+01  1.3410E+00  2.4170E+00  2.1205E+00  5.8029E+00  1.0000E-02  1.0698E+00  1.0000E-02
             1.0714E+01
 PARAMETER:  1.8964E-01  1.0575E-02  3.7992E+00  3.9344E-01  9.8255E-01  8.5164E-01  1.8584E+00 -5.0184E+00  1.6745E-01 -4.7937E+00
             2.4715E+00
 GRADIENT:   9.8881E+00  1.2436E+00  8.1563E-03  8.7564E+00  2.8576E+00  2.2507E+01  8.9415E+01  0.0000E+00  4.6726E-01  0.0000E+00
             4.4848E+01

0ITERATION NO.:  110    OBJECTIVE VALUE:  -813.962375646664        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     2971
 NPARAMETR:  1.0930E+00  9.0332E-01  4.1518E+01  1.3459E+00  2.4191E+00  2.1197E+00  5.8296E+00  1.0000E-02  1.0750E+00  1.0000E-02
             1.0706E+01
 PARAMETER:  1.8897E-01 -1.6753E-03  3.8261E+00  3.9704E-01  9.8339E-01  8.5128E-01  1.8629E+00 -5.0184E+00  1.7234E-01 -4.7937E+00
             2.4708E+00
 GRADIENT:   1.3875E-01  1.7840E-01 -1.4955E-02 -1.2307E+00  7.3375E-04  7.7532E-01  5.0843E+00  0.0000E+00  1.4995E-01  0.0000E+00
            -6.0469E-01

0ITERATION NO.:  115    OBJECTIVE VALUE:  -813.967918165777        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     3140
 NPARAMETR:  1.0931E+00  8.9530E-01  4.2833E+01  1.3515E+00  2.4211E+00  2.1200E+00  5.8432E+00  1.0000E-02  1.0787E+00  1.0000E-02
             1.0709E+01
 PARAMETER:  1.8904E-01 -1.0595E-02  3.8573E+00  4.0120E-01  9.8422E-01  8.5143E-01  1.8653E+00 -5.0184E+00  1.7573E-01 -4.7937E+00
             2.4710E+00
 GRADIENT:   9.6872E+00  1.0275E+00  1.2056E-02  9.4572E+00  2.7866E+00  2.2464E+01  9.0609E+01  0.0000E+00  4.2286E-01  0.0000E+00
             4.3997E+01

0ITERATION NO.:  120    OBJECTIVE VALUE:  -813.972692628933        NO. OF FUNC. EVALS.: 192
 CUMULATIVE NO. OF FUNC. EVALS.:     3332             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0932E+00  8.8701E-01  4.4742E+01  1.3581E+00  2.4243E+00  2.1194E+00  5.8856E+00  1.0000E-02  1.0838E+00  1.0000E-02
             1.0706E+01
 PARAMETER:  1.8915E-01 -1.9905E-02  3.9009E+00  4.0610E-01  9.8555E-01  8.5113E-01  1.8725E+00 -5.0184E+00  1.8045E-01 -4.7937E+00
             2.4708E+00
 GRADIENT:   9.7863E+00  1.2039E+00  1.5084E-02  9.8392E+00  2.7151E+00  2.2348E+01  9.2111E+01  0.0000E+00  4.3067E-01  0.0000E+00
             4.3445E+01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -813.975016688866        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     3522             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0931E+00  8.8225E-01  4.5713E+01  1.3619E+00  2.4261E+00  2.1209E+00  5.8931E+00  1.0000E-02  1.0848E+00  1.0000E-02
             1.0710E+01
 PARAMETER:  1.8901E-01 -2.5282E-02  3.9224E+00  4.0887E-01  9.8628E-01  8.5183E-01  1.8738E+00 -5.0184E+00  1.8141E-01 -4.7937E+00
             2.4712E+00
 GRADIENT:   9.6472E+00  1.1711E+00  1.7802E-02  1.0455E+01  2.6468E+00  2.2635E+01  9.2101E+01  0.0000E+00  3.0463E-01  0.0000E+00
             4.3779E+01

0ITERATION NO.:  130    OBJECTIVE VALUE:  -813.976452549100        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     3707             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0931E+00  8.7762E-01  4.6544E+01  1.3642E+00  2.4274E+00  2.1209E+00  5.9012E+00  1.0000E-02  1.0868E+00  1.0000E-02
             1.0710E+01
 PARAMETER:  1.8901E-01 -3.0539E-02  3.9404E+00  4.1057E-01  9.8682E-01  8.5183E-01  1.8752E+00 -5.0184E+00  1.8327E-01 -4.7937E+00
             2.4712E+00
 GRADIENT:   9.6469E+00  1.0861E+00  1.8734E-02  1.0551E+01  2.6525E+00  2.2645E+01  9.2341E+01  0.0000E+00  3.0255E-01  0.0000E+00
             4.3786E+01

0ITERATION NO.:  135    OBJECTIVE VALUE:  -813.976985820858        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     3885
 NPARAMETR:  1.0933E+00  8.7297E-01  4.6784E+01  1.3651E+00  2.4290E+00  2.1209E+00  5.9027E+00  1.0000E-02  1.0893E+00  1.0000E-02
             1.0714E+01
 PARAMETER:  1.8923E-01 -3.5857E-02  3.9455E+00  4.1125E-01  9.8747E-01  8.5184E-01  1.8754E+00 -5.0184E+00  1.8550E-01 -4.7937E+00
             2.4715E+00
 GRADIENT:   1.2335E-01 -4.6428E-02 -8.3806E-03 -4.6945E-01  2.1979E-02  1.0807E+00  5.0677E+00  0.0000E+00  2.0073E-02  0.0000E+00
             2.0840E-01

0ITERATION NO.:  140    OBJECTIVE VALUE:  -813.977726065358        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4070
 NPARAMETR:  1.0932E+00  8.7032E-01  4.7455E+01  1.3667E+00  2.4302E+00  2.1206E+00  5.9122E+00  1.0000E-02  1.0913E+00  1.0000E-02
             1.0712E+01
 PARAMETER:  1.8914E-01 -3.8892E-02  3.9598E+00  4.1238E-01  9.8797E-01  8.5172E-01  1.8770E+00 -5.0184E+00  1.8735E-01 -4.7937E+00
             2.4714E+00
 GRADIENT:   1.1381E-01 -4.9039E-02 -8.6275E-03 -5.4074E-01  6.2744E-02  1.0319E+00  5.2036E+00  0.0000E+00  5.1798E-02  0.0000E+00
             2.2505E-02

0ITERATION NO.:  145    OBJECTIVE VALUE:  -813.978276954856        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     4255
 NPARAMETR:  1.0932E+00  8.6872E-01  4.8072E+01  1.3683E+00  2.4312E+00  2.1206E+00  5.9185E+00  1.0000E-02  1.0925E+00  1.0000E-02
             1.0712E+01
 PARAMETER:  1.8913E-01 -4.0737E-02  3.9727E+00  4.1357E-01  9.8839E-01  8.5171E-01  1.8781E+00 -5.0184E+00  1.8850E-01 -4.7937E+00
             2.4713E+00
 GRADIENT:   1.0949E-01 -1.2389E-02 -7.1825E-03 -3.8017E-01  3.9448E-02  1.0331E+00  5.2416E+00  0.0000E+00  2.9912E-02  0.0000E+00
            -1.0512E-01

0ITERATION NO.:  146    OBJECTIVE VALUE:  -813.978276954856        NO. OF FUNC. EVALS.:  30
 CUMULATIVE NO. OF FUNC. EVALS.:     4285
 NPARAMETR:  1.0928E+00  8.6948E-01  4.9807E+01  1.3716E+00  2.4300E+00  2.1211E+00  5.9256E+00  1.0000E-02  1.0905E+00  1.0000E-02
             1.0707E+01
 PARAMETER:  1.8913E-01 -4.0737E-02  3.9727E+00  4.1357E-01  9.8839E-01  8.5171E-01  1.8781E+00 -5.0184E+00  1.8850E-01 -4.7937E+00
             2.4713E+00
 GRADIENT:   4.5040E-02 -1.2751E-02 -7.0376E-03 -3.4189E-01  5.7513E-02 -2.2109E-02 -7.6256E-02  0.0000E+00  2.9127E-02  0.0000E+00
             2.7292E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     4285
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5357E-02  3.7535E-02 -5.3973E-07 -7.4356E-02  3.4531E-06
 SE:             2.8348E-02  2.3690E-02  2.7501E-06  1.3938E-02  9.7542E-05
 N:                     100         100         100         100         100

 P VAL.:         5.8800E-01  1.1310E-01  8.4441E-01  9.5982E-08  9.7176E-01

 ETASHRINKSD(%)  5.0295E+00  2.0635E+01  9.9991E+01  5.3304E+01  9.9673E+01
 ETASHRINKVR(%)  9.8061E+00  3.7011E+01  1.0000E+02  7.8195E+01  9.9999E+01
 EBVSHRINKSD(%)  5.4892E+00  1.5705E+01  9.9987E+01  5.5861E+01  9.9567E+01
 EBVSHRINKVR(%)  1.0677E+01  2.8943E+01  1.0000E+02  8.0517E+01  9.9998E+01
 RELATIVEINF(%)  8.9127E+01  3.3866E+01  3.4073E-07  9.3089E+00  3.4421E-04
 EPSSHRINKSD(%)  5.1972E+00
 EPSSHRINKVR(%)  1.0124E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -813.97827695485580     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       840.11108281355496     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:   131.36
0R MATRIX ALGORITHMICALLY SINGULAR
0COVARIANCE MATRIX UNOBTAINABLE
0R MATRIX IS OUTPUT
0S MATRIX ALGORITHMICALLY SINGULAR
0S MATRIX IS OUTPUT
0T MATRIX - EQUAL TO RS*RMAT, WHERE S* IS A PSEUDO INVERSE OF S - IS OUTPUT
 Elapsed covariance  time in seconds:    15.61
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -813.978       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.09E+00  8.69E-01  4.81E+01  1.37E+00  2.43E+00  2.12E+00  5.92E+00  1.00E-02  1.09E+00  1.00E-02  1.07E+01
 


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
+        1.54E+02
 
 TH 2
+        1.08E+01  1.31E+01
 
 TH 3
+       -3.00E-02 -8.02E-03  3.16E-05
 
 TH 4
+       -4.24E+01  2.26E+01  3.45E-02  1.30E+02
 
 TH 5
+        1.86E+01  2.78E+00 -1.90E-02 -2.68E+01  1.18E+01
 
 TH 6
+       -1.40E+01 -1.49E+01  2.85E-02  6.91E+00 -1.56E+01  3.30E+01
 
 TH 7
+        7.72E+00  3.01E-01 -5.99E-03 -1.04E+01  3.86E+00 -4.25E+00  1.32E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        2.25E+01 -2.44E+00 -1.88E-02 -4.21E+01  1.27E+01 -1.10E+01  4.50E+00  0.00E+00  1.63E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -1.85E+00 -2.26E+00 -4.01E-04 -7.11E+00  6.72E-01  1.27E+00  3.27E-01  0.00E+00  1.80E+00  0.00E+00  8.02E-01
 
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
+        1.75E+02
 
 TH 2
+       -9.52E-01  3.17E+01
 
 TH 3
+        2.74E-03  5.90E-03  1.44E-04
 
 TH 4
+       -4.66E+00  3.14E+01  6.56E-03  1.26E+02
 
 TH 5
+       -9.62E-01 -5.07E+00 -5.72E-02 -1.22E+01  3.35E+01
 
 TH 6
+       -1.99E-01 -3.50E-01  8.45E-04  1.66E+00 -5.64E-01  3.40E+01
 
 TH 7
+        4.71E-01  3.90E+00 -1.17E-03 -8.14E+00  5.81E-01 -1.79E-01  3.04E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.37E-01 -2.99E+00 -8.13E-03 -3.21E+01  3.89E+00 -7.50E-01  2.10E+00  0.00E+00  2.12E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -6.70E+00 -2.73E+00 -4.66E-04 -9.97E+00  7.77E-01  1.76E+00  4.26E-01  0.00E+00  3.59E+00  0.00E+00  8.16E+00
 
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
+        1.86E+02
 
 TH 2
+        5.37E+01  3.09E+01
 
 TH 3
+        1.29E-02  8.47E-03  5.58E-05
 
 TH 4
+        8.26E+01  3.32E+01  2.79E-02  1.35E+02
 
 TH 5
+       -1.01E+01 -3.72E+00 -3.08E-02 -2.28E+01  1.97E+01
 
 TH 6
+        3.99E+01  1.03E+01 -3.23E-03 -1.51E+01  3.37E+00  4.26E+01
 
 TH 7
+        3.02E-01  3.75E+00 -2.38E-03 -8.90E+00  3.11E+00  3.57E+00  3.10E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -1.37E+01 -3.56E+00 -6.46E-03 -4.07E+01  6.61E+00  9.30E+00  2.50E+00  0.00E+00  2.73E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -8.69E+01 -2.39E+01 -1.71E-02 -7.78E+01  1.22E+01  1.18E+01  3.67E+00  0.00E+00  2.26E+01  0.00E+00  2.27E+02
 
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
 #CPUT: Total CPU Time in Seconds,      147.053
Stop Time:
Wed Sep 29 09:55:15 CDT 2021
