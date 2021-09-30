Thu Sep 30 08:06:55 CDT 2021
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
$DATA ../../../../data/spa2/TD2/dat45.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m45.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2047.04285293665        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.1535E+02 -2.4461E+01  4.6634E+01  9.0309E+01  1.6739E+02 -7.0546E+00 -2.4076E+01 -4.6298E+02 -1.2197E+02 -1.7768E+00
            -1.3039E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2214.29549121785        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  9.4562E-01  1.0202E+00  1.0565E+00  8.9357E-01  9.4501E-01  1.0181E+00  9.8461E-01  1.5004E+00  8.6423E-01  8.9952E-01
             9.8314E-01
 PARAMETER:  4.4088E-02  1.1997E-01  1.5499E-01 -1.2526E-02  4.3440E-02  1.1792E-01  8.4490E-02  5.0571E-01 -4.5921E-02 -5.8947E-03
             8.2998E-02
 GRADIENT:   4.2804E+02 -5.1254E+01  5.8761E+01 -1.3437E+02 -4.5105E+01  4.1274E+01 -6.3736E+00 -2.6218E+02 -8.5793E+01  1.4504E+01
            -1.0612E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2225.12664292079        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      236
 NPARAMETR:  9.4559E-01  1.0438E+00  1.0802E+00  9.0393E-01  9.5897E-01  1.0538E+00  1.0129E+00  1.5270E+00  8.3831E-01  8.6743E-01
             9.9324E-01
 PARAMETER:  4.4054E-02  1.4291E-01  1.7711E-01 -1.0019E-03  5.8102E-02  1.5243E-01  1.1278E-01  5.2331E-01 -7.6370E-02 -4.2215E-02
             9.3218E-02
 GRADIENT:  -3.7496E+01 -9.9942E+01  5.0196E+01 -1.7873E+02 -7.1628E+01 -4.5336E+01 -6.6679E+00 -2.9068E+02 -9.0573E+01  5.7456E+00
            -9.6509E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2230.24981587459        NO. OF FUNC. EVALS.: 141
 CUMULATIVE NO. OF FUNC. EVALS.:      377
 NPARAMETR:  9.4559E-01  1.0631E+00  1.0802E+00  9.0444E-01  9.7152E-01  1.1684E+00  1.0186E+00  1.5271E+00  8.3832E-01  8.5559E-01
             9.9380E-01
 PARAMETER:  4.4054E-02  1.6116E-01  1.7711E-01 -4.3449E-04  7.1108E-02  2.5563E-01  1.1839E-01  5.2337E-01 -7.6354E-02 -5.5966E-02
             9.3784E-02
 GRADIENT:  -3.0285E+01 -8.3661E+01  4.6588E+01 -1.5512E+02 -6.0844E+01  1.3567E+00 -3.4061E+00 -2.9121E+02 -8.8783E+01  2.0484E+00
            -9.5942E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2230.43342623858        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:      522             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4559E-01  1.0627E+00  1.0802E+00  9.0464E-01  9.7471E-01  1.1697E+00  1.0179E+00  1.5268E+00  8.3850E-01  8.5541E-01
             9.9402E-01
 PARAMETER:  4.4053E-02  1.6082E-01  1.7711E-01 -2.2117E-04  7.4389E-02  2.5671E-01  1.1776E-01  5.2320E-01 -7.6142E-02 -5.6178E-02
             9.3997E-02
 GRADIENT:   4.2763E+02  3.1274E+01  5.2025E+01 -7.6190E+01 -3.0162E+01  2.2527E+02  5.2742E+00 -2.5896E+02 -8.2291E+01  3.4493E+00
            -9.4436E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2230.61492669926        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      703
 NPARAMETR:  9.4559E-01  1.0628E+00  1.0801E+00  9.0473E-01  9.7472E-01  1.1339E+00  1.0179E+00  1.5287E+00  8.3856E-01  8.5539E-01
             9.9409E-01
 PARAMETER:  4.4058E-02  1.6086E-01  1.7709E-01 -1.1378E-04  7.4391E-02  2.2570E-01  1.1775E-01  5.2441E-01 -7.6074E-02 -5.6194E-02
             9.4069E-02
 GRADIENT:  -3.2287E+01 -8.5759E+01  4.4538E+01 -1.5307E+02 -5.4902E+01 -1.1357E+01 -3.6661E+00 -2.9062E+02 -8.8279E+01  1.8772E+00
            -9.5642E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2231.01748830195        NO. OF FUNC. EVALS.: 166
 CUMULATIVE NO. OF FUNC. EVALS.:      869            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4560E-01  1.0623E+00  1.0801E+00  9.0498E-01  9.7645E-01  1.1695E+00  1.0190E+00  1.5289E+00  8.3877E-01  8.5432E-01
             9.9435E-01
 PARAMETER:  4.4059E-02  1.6048E-01  1.7709E-01  1.5457E-04  7.6166E-02  2.5658E-01  1.1886E-01  5.2455E-01 -7.5815E-02 -5.7449E-02
             9.4329E-02
 GRADIENT:   4.2733E+02  3.0073E+01  5.0727E+01 -7.4855E+01 -2.6781E+01  2.2499E+02  5.4399E+00 -2.5817E+02 -8.1853E+01  3.2622E+00
            -9.3932E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2231.02563872170        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:     1015
 NPARAMETR:  9.4560E-01  1.0621E+00  1.0801E+00  9.0497E-01  9.7689E-01  1.1651E+00  1.0183E+00  1.5289E+00  8.3877E-01  8.5368E-01
             9.9447E-01
 PARAMETER:  4.4059E-02  1.6029E-01  1.7709E-01  1.5001E-04  7.6623E-02  2.5283E-01  1.1815E-01  5.2454E-01 -7.5819E-02 -5.8195E-02
             9.4450E-02
 GRADIENT:  -3.0479E+01 -8.7453E+01  4.3251E+01 -1.5211E+02 -5.0581E+01  1.8681E-01 -3.8045E+00 -2.9056E+02 -8.8045E+01  1.5510E+00
            -9.5126E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2231.10926546797        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1205
 NPARAMETR:  9.4560E-01  1.0618E+00  1.0801E+00  9.0518E-01  9.7765E-01  1.1653E+00  1.0188E+00  1.5290E+00  8.3877E-01  8.5313E-01
             9.9469E-01
 PARAMETER:  4.4059E-02  1.5995E-01  1.7709E-01  3.7936E-04  7.7397E-02  2.5301E-01  1.1865E-01  5.2459E-01 -7.5819E-02 -5.8838E-02
             9.4678E-02
 GRADIENT:  -3.0471E+01 -8.7998E+01  4.2725E+01 -1.5161E+02 -4.8944E+01  2.6249E-01 -3.8065E+00 -2.9053E+02 -8.7922E+01  1.4555E+00
            -9.4870E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2231.35299987561        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1388
 NPARAMETR:  9.4560E-01  1.0614E+00  1.0801E+00  9.0525E-01  9.7740E-01  1.1692E+00  1.0185E+00  1.5303E+00  8.3878E-01  8.5291E-01
             9.9499E-01
 PARAMETER:  4.4063E-02  1.5958E-01  1.7707E-01  4.5450E-04  7.7144E-02  2.5630E-01  1.1834E-01  5.2543E-01 -7.5810E-02 -5.9101E-02
             9.4980E-02
 GRADIENT:  -3.0253E+01 -8.8284E+01  4.2649E+01 -1.5195E+02 -4.9081E+01  1.6265E+00 -3.9027E+00 -2.9005E+02 -8.7860E+01  1.4929E+00
            -9.4435E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2231.55537014367        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1572
 NPARAMETR:  9.4560E-01  1.0610E+00  1.0801E+00  9.0549E-01  9.7715E-01  1.1694E+00  1.0182E+00  1.5312E+00  8.3878E-01  8.5268E-01
             9.9529E-01
 PARAMETER:  4.4065E-02  1.5919E-01  1.7706E-01  7.1693E-04  7.6885E-02  2.5649E-01  1.1802E-01  5.2606E-01 -7.5803E-02 -5.9368E-02
             9.5276E-02
 GRADIENT:  -3.0240E+01 -8.8402E+01  4.2553E+01 -1.5194E+02 -4.9123E+01  1.7041E+00 -3.9923E+00 -2.8968E+02 -8.7771E+01  1.5326E+00
            -9.4042E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2231.88610199645        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1756
 NPARAMETR:  9.4560E-01  1.0606E+00  1.0801E+00  9.0580E-01  9.7690E-01  1.1691E+00  1.0179E+00  1.5328E+00  8.3879E-01  8.5245E-01
             9.9560E-01
 PARAMETER:  4.4069E-02  1.5882E-01  1.7704E-01  1.0589E-03  7.6631E-02  2.5623E-01  1.1771E-01  5.2708E-01 -7.5792E-02 -5.9635E-02
             9.5589E-02
 GRADIENT:  -3.0254E+01 -8.8423E+01  4.2377E+01 -1.5173E+02 -4.9135E+01  1.5968E+00 -4.0725E+00 -2.8910E+02 -8.7620E+01  1.5813E+00
            -9.3597E+01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2232.20461665562        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1940
 NPARAMETR:  9.4561E-01  1.0602E+00  1.0801E+00  9.0610E-01  9.7665E-01  1.1691E+00  1.0175E+00  1.5343E+00  8.3880E-01  8.5223E-01
             9.9591E-01
 PARAMETER:  4.4073E-02  1.5845E-01  1.7703E-01  1.3976E-03  7.6377E-02  2.5623E-01  1.1740E-01  5.2805E-01 -7.5781E-02 -5.9902E-02
             9.5900E-02
 GRADIENT:  -3.0251E+01 -8.8448E+01  4.2209E+01 -1.5154E+02 -4.9149E+01  1.5984E+00 -4.1529E+00 -2.8854E+02 -8.7474E+01  1.6291E+00
            -9.3158E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -2232.52317568964        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2124
 NPARAMETR:  9.4561E-01  1.0598E+00  1.0801E+00  9.0641E-01  9.7641E-01  1.1691E+00  1.0172E+00  1.5358E+00  8.3881E-01  8.5200E-01
             9.9622E-01
 PARAMETER:  4.4076E-02  1.5809E-01  1.7701E-01  1.7363E-03  7.6124E-02  2.5623E-01  1.1709E-01  5.2903E-01 -7.5771E-02 -6.0168E-02
             9.6210E-02
 GRADIENT:  -3.0249E+01 -8.8471E+01  4.2041E+01 -1.5135E+02 -4.9164E+01  1.5998E+00 -4.2330E+00 -2.8799E+02 -8.7328E+01  1.6769E+00
            -9.2719E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -2232.74558714332        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     2308
 NPARAMETR:  9.4561E-01  1.0590E+00  1.0800E+00  9.0694E-01  9.7591E-01  1.1996E+00  1.0166E+00  1.5373E+00  8.3882E-01  8.5156E-01
             9.9678E-01
 PARAMETER:  4.4080E-02  1.5733E-01  1.7699E-01  2.3259E-03  7.5620E-02  2.8196E-01  1.1648E-01  5.3004E-01 -7.5760E-02 -6.0685E-02
             9.6771E-02
 GRADIENT:  -2.8629E+01 -8.8694E+01  4.1911E+01 -1.5126E+02 -4.9244E+01  1.1999E+01 -4.4159E+00 -2.8739E+02 -8.7173E+01  1.7537E+00
            -9.1918E+01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -2232.74558714332        NO. OF FUNC. EVALS.:  25
 CUMULATIVE NO. OF FUNC. EVALS.:     2333
 NPARAMETR:  9.4561E-01  1.0590E+00  1.0800E+00  9.0694E-01  9.7591E-01  1.1996E+00  1.0166E+00  1.5373E+00  8.3882E-01  8.5156E-01
             9.9678E-01
 PARAMETER:  4.4080E-02  1.5733E-01  1.7699E-01  2.3259E-03  7.5620E-02  2.8196E-01  1.1648E-01  5.3004E-01 -7.5760E-02 -6.0685E-02
             9.6771E-02
 GRADIENT:  -2.9807E+01 -2.3701E+04  9.0255E+01  7.4121E+04 -3.7193E+04 -2.6336E+04 -3.1893E+04  1.3608E+04  7.4191E+04 -7.4285E+04
             7.4180E+04

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2333
 NO. OF SIG. DIGITS IN FINAL EST.:  2.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6843E-02  2.6493E-02 -5.9678E-02  6.5600E-02 -1.3722E-02
 SE:             2.9153E-02  2.2591E-02  4.0400E-02  3.0712E-02  2.1627E-02
 N:                     100         100         100         100         100

 P VAL.:         5.6345E-01  2.4091E-01  1.3963E-01  3.2682E-02  5.2577E-01

 ETASHRINKSD(%)  2.3325E+00  2.4318E+01  1.0000E-10  1.0000E-10  2.7547E+01
 ETASHRINKVR(%)  4.6105E+00  4.2722E+01  1.0000E-10  1.0000E-10  4.7506E+01
 EBVSHRINKSD(%)  2.3183E-01  2.4696E+01  3.8633E+01  2.0710E+01  2.5165E+01
 EBVSHRINKVR(%)  4.6311E-01  4.3293E+01  6.2341E+01  3.7131E+01  4.3997E+01
 RELATIVEINF(%)  9.9530E+01  1.7995E+01  2.1514E+01  2.2229E+01  2.1422E+01
 EPSSHRINKSD(%)  2.5602E+01
 EPSSHRINKVR(%)  4.4649E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          600
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1102.7262398456071     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2232.7455871433240     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1130.0193472977169     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    51.73
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    10.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2232.746       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  1.06E+00  1.08E+00  9.07E-01  9.76E-01  1.20E+00  1.02E+00  1.54E+00  8.39E-01  8.52E-01  9.97E-01
 


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
+        4.15E+07
 
 TH 2
+       -3.44E+00  6.69E+06
 
 TH 3
+        9.14E+01  5.84E+06  1.02E+07
 
 TH 4
+       -4.59E+00 -3.69E+03 -1.07E+07  2.26E+07
 
 TH 5
+       -1.58E+00  5.92E+02  9.96E+06 -1.37E+03  1.95E+07
 
 TH 6
+        2.64E+01 -1.22E+02  6.18E-01  2.21E+02 -2.08E+02  1.62E+06
 
 TH 7
+        4.65E-02  1.25E+03  8.21E+06 -2.12E+03  1.23E+03 -1.71E+02  1.32E+07
 
 TH 8
+       -8.31E-02 -4.67E+02 -6.97E+01 -8.91E+02 -1.88E+02  2.46E+01 -2.38E+02  5.60E+05
 
 TH 9
+        8.37E-01 -4.49E+03 -1.16E+07  2.44E+07 -2.27E+07  2.40E+02 -1.87E+07  5.97E+01  5.28E+07
 
 TH10
+       -5.94E-02  2.42E+02  1.14E+07 -4.98E+02  1.70E+03 -2.37E+02  2.29E+03 -4.98E+01 -2.60E+07  2.56E+07
 
 TH11
+       -4.08E+00 -3.81E+03 -9.75E+06 -7.55E+03 -1.49E+03  2.04E+02 -1.94E+03 -8.01E+02  2.22E+07 -4.60E+02  1.87E+07
 
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
 #CPUT: Total CPU Time in Seconds,       61.963
Stop Time:
Thu Sep 30 08:07:59 CDT 2021
