Sat Sep 18 08:43:49 CDT 2021
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
$DATA ../../../../data/spa/B/dat78.csv ignore=@
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1672.11352374190        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   9.8105E+01 -1.0873E+02 -6.3616E+01 -8.2230E+01  3.4944E+01  3.5812E+00 -1.7981E+01  1.8088E+01  6.9965E+00  5.7693E+00
             1.9965E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1691.75573166084        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.6260E-01  1.3565E+00  2.0134E+00  8.4180E-01  1.4271E+00  9.7225E-01  1.2371E+00  7.0573E-01  8.5434E-01  1.1669E+00
             9.6646E-01
 PARAMETER:  6.1884E-02  4.0490E-01  7.9983E-01 -7.2218E-02  4.5567E-01  7.1854E-02  3.1279E-01 -2.4852E-01 -5.7421E-02  2.5438E-01
             6.5883E-02
 GRADIENT:   1.8121E+01  7.8221E+00  9.0352E+00 -3.3280E+01  7.8217E+00 -6.5864E+00  2.2644E+01 -8.4059E-01 -2.2208E+00 -3.3177E+01
            -1.5565E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1695.93271523400        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      159
 NPARAMETR:  9.8486E-01  1.2195E+00  2.2759E+00  9.3890E-01  1.4595E+00  1.0245E+00  1.0891E+00  4.2152E-01  1.0364E+00  1.2520E+00
             9.7109E-01
 PARAMETER:  8.4742E-02  2.9845E-01  9.2240E-01  3.6950E-02  4.7809E-01  1.2422E-01  1.8538E-01 -7.6389E-01  1.3578E-01  3.2472E-01
             7.0660E-02
 GRADIENT:   6.9031E+01 -1.0281E+01  1.7739E+00 -2.3434E+01  3.1978E+01  1.4415E+01  8.8492E+00 -2.3685E-01  3.6156E+00 -1.9562E+01
            -6.7850E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1698.21436043273        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.6633E-01  1.1041E+00  1.7783E+00  1.0285E+00  1.3019E+00  9.9654E-01  1.1580E+00  1.6423E-01  9.6686E-01  1.2377E+00
             9.6760E-01
 PARAMETER:  6.5754E-02  1.9902E-01  6.7566E-01  1.2809E-01  3.6385E-01  9.6530E-02  2.4667E-01 -1.7065E+00  6.6298E-02  3.1327E-01
             6.7060E-02
 GRADIENT:   2.9942E+01  5.3176E+00 -7.3306E-01  5.3978E+00  7.3988E+00  3.8816E+00  4.2937E+00 -3.4167E-03  2.4265E+00 -4.4677E+00
            -6.4844E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1698.21472938769        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      303
 NPARAMETR:  9.6499E-01  1.1020E+00  1.7792E+00  1.0294E+00  1.3005E+00  9.9556E-01  1.1540E+00  1.6220E-01  9.6678E-01  1.2395E+00
             9.6762E-01
 PARAMETER:  6.4362E-02  1.9713E-01  6.7615E-01  1.2896E-01  3.6278E-01  9.5551E-02  2.4326E-01 -1.7189E+00  6.6214E-02  3.1473E-01
             6.7083E-02
 GRADIENT:   2.6842E+01  4.8338E+00 -6.5310E-01  4.8439E+00  6.6272E+00  3.4874E+00  3.9323E+00 -2.5820E-03  2.2235E+00 -4.0244E+00
            -5.7855E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1698.43307425289        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  9.7354E-01  1.1038E+00  1.9325E+00  1.0335E+00  1.3254E+00  9.9523E-01  1.0642E+00  1.6770E-01  9.8903E-01  1.2885E+00
             9.7200E-01
 PARAMETER:  7.3184E-02  1.9877E-01  7.5882E-01  1.3295E-01  3.8173E-01  9.5220E-02  1.6221E-01 -1.6856E+00  8.8965E-02  3.5349E-01
             7.1605E-02
 GRADIENT:   2.9342E+00 -1.3375E-01  9.6969E-01 -1.6268E+00 -2.4800E+00 -3.7235E-01 -2.0430E-02 -1.4376E-02 -8.4241E-01  3.7205E-01
            -3.8876E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1698.45505671285        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      638
 NPARAMETR:  9.7269E-01  1.2141E+00  1.8705E+00  9.6440E-01  1.3572E+00  9.9655E-01  9.8246E-01  1.9894E-01  1.0511E+00  1.2973E+00
             9.7340E-01
 PARAMETER:  7.2307E-02  2.9397E-01  7.2621E-01  6.3751E-02  4.0543E-01  9.6540E-02  8.2300E-02 -1.5147E+00  1.4981E-01  3.6028E-01
             7.3040E-02
 GRADIENT:  -4.2977E-01  2.8065E+00 -2.0417E-01  4.0925E+00  1.1670E+00 -1.3526E-01 -5.5519E-01 -1.9284E-02 -5.1683E-01 -1.0416E-01
             1.5678E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1698.49731963183        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      815
 NPARAMETR:  9.7315E-01  1.3423E+00  1.7231E+00  8.7906E-01  1.3739E+00  9.9811E-01  9.4741E-01  2.3840E-01  1.1057E+00  1.2958E+00
             9.7074E-01
 PARAMETER:  7.2781E-02  3.9441E-01  6.4415E-01 -2.8902E-02  4.1762E-01  9.8106E-02  4.5973E-02 -1.3338E+00  2.0049E-01  3.5916E-01
             7.0306E-02
 GRADIENT:  -6.9667E-01  4.3743E+00  1.2100E-01  4.3749E+00 -7.8703E-01  2.4089E-01 -3.0930E-01 -1.0363E-02 -4.0889E-01 -4.1968E-01
            -4.1631E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1698.57283348403        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      992
 NPARAMETR:  9.7407E-01  1.5227E+00  1.5771E+00  7.5729E-01  1.4176E+00  9.9871E-01  8.8655E-01  3.0920E-01  1.2161E+00  1.3150E+00
             9.7034E-01
 PARAMETER:  7.3732E-02  5.2051E-01  5.5562E-01 -1.7801E-01  4.4895E-01  9.8710E-02 -2.0413E-02 -1.0738E+00  2.9568E-01  3.7384E-01
             6.9891E-02
 GRADIENT:  -3.0505E-01  3.3043E+00  3.1575E-02  2.6334E+00 -1.0943E+00  1.7830E-01 -1.5567E-01  2.2291E-02 -2.3067E-01  9.6026E-02
            -1.8687E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1698.64980248472        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1167
 NPARAMETR:  9.7513E-01  1.7321E+00  1.3957E+00  6.1655E-01  1.4695E+00  9.9927E-01  8.3240E-01  3.7893E-01  1.3865E+00  1.3375E+00
             9.7165E-01
 PARAMETER:  7.4817E-02  6.4936E-01  4.3336E-01 -3.8361E-01  4.8495E-01  9.9270E-02 -8.3438E-02 -8.7040E-01  4.2677E-01  3.9082E-01
             7.1237E-02
 GRADIENT:   3.5382E-01  3.4317E+00 -9.1897E-02  2.1350E+00 -1.9094E+00  1.0124E-01  1.8155E-01  9.8545E-02  3.0380E-01  2.1960E-01
             5.2644E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1698.73749211889        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1343
 NPARAMETR:  9.7549E-01  1.9378E+00  1.1687E+00  4.7948E-01  1.5228E+00  9.9947E-01  7.9635E-01  2.8374E-01  1.5882E+00  1.3621E+00
             9.7080E-01
 PARAMETER:  7.5187E-02  7.6154E-01  2.5590E-01 -6.3505E-01  5.2054E-01  9.9472E-02 -1.2771E-01 -1.1597E+00  5.6262E-01  4.0902E-01
             7.0366E-02
 GRADIENT:  -3.4372E-01  6.0170E+00 -5.4957E-01  2.6726E+00  1.6882E-01 -8.5567E-02 -2.5564E-01  7.5530E-02 -3.0471E-01  4.3642E-02
             1.9631E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1698.80136406955        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1518
 NPARAMETR:  9.7581E-01  2.0739E+00  1.0668E+00  3.8315E-01  1.5669E+00  1.0002E+00  7.6973E-01  1.4438E-01  1.8219E+00  1.3879E+00
             9.7251E-01
 PARAMETER:  7.5510E-02  8.2941E-01  1.6465E-01 -8.5934E-01  5.4909E-01  1.0025E-01 -1.6172E-01 -1.8353E+00  6.9989E-01  4.2782E-01
             7.2122E-02
 GRADIENT:  -4.8042E-01 -1.6140E+00 -2.9658E-01 -2.4004E-01  4.2282E-01  8.5122E-02 -5.0707E-02  2.3821E-02  8.5164E-02  1.4890E-01
             1.6560E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1698.80285673294        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1694
 NPARAMETR:  9.7618E-01  2.1068E+00  1.0740E+00  3.6147E-01  1.5801E+00  1.0001E+00  7.6220E-01  1.2085E-01  1.8927E+00  1.3966E+00
             9.7397E-01
 PARAMETER:  7.5887E-02  8.4517E-01  1.7135E-01 -9.1758E-01  5.5747E-01  1.0012E-01 -1.7155E-01 -2.0132E+00  7.3799E-01  4.3401E-01
             7.3623E-02
 GRADIENT:   9.8615E-02 -6.3724E-01 -1.0583E-01 -1.8218E-01  1.6667E-01 -5.3174E-03 -3.9143E-02  1.7401E-02 -1.6277E-02  7.7384E-02
             1.3787E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1698.80297529665        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1870
 NPARAMETR:  9.7614E-01  2.1120E+00  1.0811E+00  3.5811E-01  1.5823E+00  1.0002E+00  7.6073E-01  1.1548E-01  1.9075E+00  1.3982E+00
             9.7424E-01
 PARAMETER:  7.5851E-02  8.4765E-01  1.7802E-01 -9.2692E-01  5.5891E-01  1.0020E-01 -1.7347E-01 -2.0586E+00  7.4578E-01  4.3519E-01
             7.3899E-02
 GRADIENT:  -2.3206E-02 -2.9702E-01 -1.6246E-02 -1.2334E-01 -4.7595E-02  1.8104E-02  2.0447E-02  1.5757E-02  1.7865E-02  6.4589E-02
             9.8115E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1698.81370426618        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     2052
 NPARAMETR:  9.7598E-01  2.0662E+00  1.1059E+00  3.8970E-01  1.5674E+00  1.0000E+00  7.6746E-01  8.5317E-02  1.8206E+00  1.3884E+00
             9.7320E-01
 PARAMETER:  7.5688E-02  8.2569E-01  2.0067E-01 -8.4238E-01  5.4939E-01  1.0003E-01 -1.6466E-01 -2.3614E+00  6.9916E-01  4.2816E-01
             7.2831E-02
 GRADIENT:  -6.5990E-02  1.0304E+00 -2.5417E-02  3.8072E-01  3.2958E-02 -1.9540E-03 -7.1997E-02  8.2739E-03  1.0592E-01  1.0568E-02
             7.4378E-02

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1698.81757629434        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2228
 NPARAMETR:  9.7602E-01  2.0711E+00  1.0997E+00  3.8579E-01  1.5687E+00  1.0000E+00  7.6727E-01  1.8389E-02  1.8277E+00  1.3892E+00
             9.7306E-01
 PARAMETER:  7.5726E-02  8.2808E-01  1.9499E-01 -8.5246E-01  5.5026E-01  1.0005E-01 -1.6492E-01 -3.8960E+00  7.0308E-01  4.2872E-01
             7.2688E-02
 GRADIENT:  -5.2929E-03  5.8609E-02 -5.8793E-03  2.4397E-02  1.6711E-02  3.9072E-04  2.3159E-03  3.8703E-04  9.6085E-03 -1.4084E-03
             4.2994E-03

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1698.81840105592        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2417
 NPARAMETR:  9.7598E-01  2.0621E+00  1.1063E+00  3.9176E-01  1.5663E+00  1.0000E+00  7.6843E-01  1.0000E-02  1.8128E+00  1.3880E+00
             9.7294E-01
 PARAMETER:  7.5684E-02  8.2374E-01  2.0102E-01 -8.3712E-01  5.4873E-01  1.0004E-01 -1.6340E-01 -4.6438E+00  6.9490E-01  4.2787E-01
             7.2568E-02
 GRADIENT:  -3.7527E-02 -1.3097E-01 -2.4632E-02  3.7396E-02  1.3691E-01  9.2479E-03 -1.8603E-02  0.0000E+00  4.3815E-02  4.6224E-02
             3.4734E-02

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1698.81843143183        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     2600
 NPARAMETR:  9.7600E-01  2.0621E+00  1.1062E+00  3.9167E-01  1.5660E+00  1.0000E+00  7.6859E-01  1.0000E-02  1.8121E+00  1.3876E+00
             9.7289E-01
 PARAMETER:  7.5708E-02  8.2375E-01  2.0096E-01 -8.3733E-01  5.4852E-01  1.0002E-01 -1.6320E-01 -4.6438E+00  6.9448E-01  4.2756E-01
             7.2520E-02
 GRADIENT:   1.9299E-02 -1.9412E-01  8.3528E-03 -4.1959E-02 -8.1326E-03  2.8497E-03  8.0599E-03  0.0000E+00 -1.5539E-03  2.1366E-03
            -4.0846E-03

0ITERATION NO.:   87    OBJECTIVE VALUE:  -1698.81843339251        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     2657
 NPARAMETR:  9.7600E-01  2.0622E+00  1.1060E+00  3.9168E-01  1.5660E+00  1.0000E+00  7.6858E-01  1.0000E-02  1.8121E+00  1.3876E+00
             9.7290E-01
 PARAMETER:  7.5707E-02  8.2376E-01  2.0076E-01 -8.3730E-01  5.4852E-01  1.0002E-01 -1.6321E-01 -4.6438E+00  6.9448E-01  4.2756E-01
             7.2521E-02
 GRADIENT:   1.5090E-02 -1.5562E-01  2.7969E-03 -2.3991E-02  5.7618E-03  1.9207E-03  2.7351E-03  0.0000E+00  3.2524E-03  3.5563E-03
            -8.1733E-04

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     2657
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -6.7508E-04 -2.3804E-02 -6.4704E-05  1.7487E-02 -4.1103E-02
 SE:             2.9763E-02  2.3279E-02  3.2933E-05  1.8623E-02  2.3720E-02
 N:                     100         100         100         100         100

 P VAL.:         9.8190E-01  3.0653E-01  4.9448E-02  3.4775E-01  8.3121E-02

 ETASHRINKSD(%)  2.8971E-01  2.2012E+01  9.9890E+01  3.7609E+01  2.0536E+01
 ETASHRINKVR(%)  5.7857E-01  3.9178E+01  1.0000E+02  6.1074E+01  3.6855E+01
 EBVSHRINKSD(%)  4.2397E-01  1.8779E+01  9.9894E+01  4.2813E+01  1.6931E+01
 EBVSHRINKVR(%)  8.4614E-01  3.4031E+01  1.0000E+02  6.7296E+01  3.0996E+01
 RELATIVEINF(%)  9.9065E+01  3.0702E+00  4.5696E-05  1.5114E+00  3.3978E+01
 EPSSHRINKSD(%)  4.1086E+01
 EPSSHRINKVR(%)  6.5292E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1698.8184333925087     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -963.66760682877054     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    34.69
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.70
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1698.818       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.76E-01  2.06E+00  1.11E+00  3.92E-01  1.57E+00  1.00E+00  7.69E-01  1.00E-02  1.81E+00  1.39E+00  9.73E-01
 


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
+        1.16E+03
 
 TH 2
+       -7.26E+00  3.01E+02
 
 TH 3
+       -4.76E-01  7.90E+00  1.11E+01
 
 TH 4
+       -5.43E+00  4.07E+02 -2.40E+01  7.45E+02
 
 TH 5
+       -2.82E+00 -3.31E+01 -1.89E+01  4.39E+01  1.64E+02
 
 TH 6
+        3.98E+00 -1.84E+00  6.04E-01 -2.54E+00 -9.72E-01  2.02E+02
 
 TH 7
+        1.20E+00 -3.45E+00  9.55E+00 -1.54E+01 -8.74E+00  5.29E+00  1.60E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        4.31E-02 -6.00E+00 -4.66E+00  3.32E+01  1.69E+00 -4.27E-01  1.65E+01  0.00E+00  1.50E+01
 
 TH10
+       -6.54E-01 -5.28E+00 -3.57E+00  2.42E+00 -3.03E+01 -7.94E-02 -3.59E+00  0.00E+00  1.68E+00  5.51E+01
 
 TH11
+       -7.71E+00 -1.43E+01 -8.87E+00 -7.42E-01  6.59E+00 -1.50E+00  1.12E+01  0.00E+00  1.75E+00  1.65E+01  2.57E+02
 
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
 #CPUT: Total CPU Time in Seconds,       41.459
Stop Time:
Sat Sep 18 08:44:33 CDT 2021
