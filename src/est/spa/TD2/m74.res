Sat Sep 25 13:45:14 CDT 2021
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
$DATA ../../../../data/spa/TD2/dat74.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m74.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1691.59872131163        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4771E+01 -5.4862E+01 -1.7458E+01 -7.0266E+01 -1.5839E+01  2.3490E+01 -1.0757E+01  1.4421E+01 -5.8828E+00  1.6169E+01
             2.2866E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.66474619049        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      112
 NPARAMETR:  9.9181E-01  1.0727E+00  1.1270E+00  1.0135E+00  1.0991E+00  9.4263E-01  1.1041E+00  8.8728E-01  1.0270E+00  9.6224E-01
             9.8165E-01
 PARAMETER:  9.1772E-02  1.7019E-01  2.1954E-01  1.1345E-01  1.9448E-01  4.0914E-02  1.9901E-01 -1.9593E-02  1.2664E-01  6.1514E-02
             8.1477E-02
 GRADIENT:  -4.5604E+00  1.6179E+00  3.5644E-01  8.3818E-01  4.2855E+00 -2.2758E+00 -8.6288E-01  3.7059E+00  1.8494E+00 -7.3924E+00
             8.6706E+00

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1699.66234477021        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      290
 NPARAMETR:  9.9115E-01  9.5728E-01  1.1462E+00  1.0865E+00  1.0618E+00  9.5612E-01  1.1831E+00  6.4743E-01  9.9259E-01  1.0386E+00
             9.5547E-01
 PARAMETER:  9.1107E-02  5.6336E-02  2.3646E-01  1.8293E-01  1.5992E-01  5.5126E-02  2.6810E-01 -3.3474E-01  9.2567E-02  1.3785E-01
             5.4451E-02
 GRADIENT:  -3.2957E+00  1.5590E+00  4.8156E+00  1.3907E+00 -1.7236E+00  3.6965E+00 -1.3305E+00 -7.2253E-01  4.1139E+00  2.1864E+00
            -1.7927E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1700.31122363374        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      467
 NPARAMETR:  9.9291E-01  8.0486E-01  1.0373E+00  1.1716E+00  9.3647E-01  9.4456E-01  1.4147E+00  4.7499E-01  8.9715E-01  9.2709E-01
             9.5722E-01
 PARAMETER:  9.2884E-02 -1.1708E-01  1.3666E-01  2.5837E-01  3.4364E-02  4.2963E-02  4.4689E-01 -6.4446E-01 -8.5287E-03  2.4293E-02
             5.6279E-02
 GRADIENT:   1.6920E+00  2.3320E+00  3.7376E+00 -2.2052E-01 -4.5012E+00 -8.0828E-01 -3.7721E-01 -5.0052E-01 -4.0777E-01 -2.7836E-01
            -3.0917E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1700.55766602006        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      644
 NPARAMETR:  9.8856E-01  5.7668E-01  1.1695E+00  1.3205E+00  9.1156E-01  9.4335E-01  1.7644E+00  6.5059E-01  8.4134E-01  9.4193E-01
             9.5763E-01
 PARAMETER:  8.8496E-02 -4.5047E-01  2.5658E-01  3.7804E-01  7.4028E-03  4.1677E-02  6.6780E-01 -3.2988E-01 -7.2763E-02  4.0175E-02
             5.6704E-02
 GRADIENT:  -9.4876E-01  4.4609E+00  5.4702E+00  6.3505E+00 -1.0442E+01  8.8613E-02 -1.5818E-01  2.6865E-01  1.5617E-01  1.0392E+00
             5.8702E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1700.75091299811        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      822
 NPARAMETR:  9.8692E-01  4.5352E-01  1.1879E+00  1.3918E+00  8.8638E-01  9.4217E-01  2.0613E+00  6.4046E-01  8.1531E-01  9.3442E-01
             9.5688E-01
 PARAMETER:  8.6831E-02 -6.9071E-01  2.7221E-01  4.3060E-01 -2.0604E-02  4.0436E-02  8.2332E-01 -3.4557E-01 -1.0419E-01  3.2172E-02
             5.5921E-02
 GRADIENT:  -4.6999E-01  1.5198E+00  7.1766E-01  4.4495E+00 -1.3691E-01  3.6606E-01  1.5232E-01 -4.3700E-01 -1.7769E-01 -3.3600E-01
            -1.8720E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1701.02708170285        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      999
 NPARAMETR:  9.8322E-01  2.4801E-01  1.2704E+00  1.5133E+00  8.6326E-01  9.3857E-01  2.9052E+00  7.5293E-01  7.8004E-01  9.4014E-01
             9.5617E-01
 PARAMETER:  8.3073E-02 -1.2943E+00  3.3930E-01  5.1428E-01 -4.7039E-02  3.6599E-02  1.1665E+00 -1.8378E-01 -1.4840E-01  3.8278E-02
             5.5178E-02
 GRADIENT:  -2.6514E-01  1.0282E+00 -1.2543E+00 -4.3891E+00  1.3344E+00  3.2286E-01  1.8747E+00  1.5649E-01 -2.6235E+00  2.7956E-01
             1.2330E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1701.46319882616        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1176
 NPARAMETR:  9.8029E-01  9.3797E-02  1.4339E+00  1.6227E+00  8.8164E-01  9.3524E-01  4.6191E+00  9.5284E-01  7.5638E-01  9.6057E-01
             9.5624E-01
 PARAMETER:  8.0093E-02 -2.2666E+00  4.6037E-01  5.8409E-01 -2.5968E-02  3.3045E-02  1.6302E+00  5.1687E-02 -1.7922E-01  5.9770E-02
             5.5249E-02
 GRADIENT:  -5.0447E-01  1.1945E+00  1.7667E+00  1.4271E+01 -6.2755E+00 -1.1828E-01  2.5423E-01  7.1455E-01 -6.9836E-01  2.7348E-01
             9.1453E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1701.94866579360        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1355
 NPARAMETR:  9.7930E-01  3.6576E-02  1.5741E+00  1.6696E+00  9.1788E-01  9.3427E-01  6.9109E+00  1.0175E+00  7.3927E-01  1.0020E+00
             9.5790E-01
 PARAMETER:  7.9085E-02 -3.2084E+00  5.5368E-01  6.1256E-01  1.4311E-02  3.2005E-02  2.0331E+00  1.1731E-01 -2.0209E-01  1.0200E-01
             5.6992E-02
 GRADIENT:  -1.0421E+00 -1.2839E+00  8.3731E+00  3.8665E+01 -9.6394E+00 -5.4890E-02 -6.0298E+00 -1.6030E+00  6.7286E-01  2.8851E+00
            -6.7623E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.12030054831        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1538
 NPARAMETR:  9.7928E-01  1.0000E-02  1.5499E+00  1.6620E+00  9.0697E-01  9.3539E-01  1.2500E+01  1.0541E+00  7.3087E-01  9.5428E-01
             9.5702E-01
 PARAMETER:  7.9059E-02 -4.6913E+00  5.3819E-01  6.0805E-01  2.3583E-03  3.3212E-02  2.6257E+00  1.5271E-01 -2.1352E-01  5.3205E-02
             5.6065E-02
 GRADIENT:   4.5832E-01  0.0000E+00  3.1516E+00 -5.8506E+00  1.8732E+00  6.3812E-01 -1.8316E+00  3.3948E-01 -1.2231E+00 -4.5003E-01
            -4.6074E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1703.30985025873        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1717
 NPARAMETR:  9.7928E-01  1.0000E-02  1.4359E+00  1.6619E+00  8.6615E-01  9.3317E-01  1.2247E+01  9.2336E-01  7.3732E-01  9.4077E-01
             9.5805E-01
 PARAMETER:  7.9058E-02 -4.6867E+00  4.6179E-01  6.0793E-01 -4.3701E-02  3.0834E-02  2.6053E+00  2.0262E-02 -2.0474E-01  3.8947E-02
             5.7147E-02
 GRADIENT:   7.6312E-01  0.0000E+00  2.9305E+00  1.0001E+01 -1.0109E+00 -3.7634E-01  5.1244E-01 -1.0056E+00 -5.9045E-01 -7.1888E-01
            -4.2955E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1703.37230334566        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1894
 NPARAMETR:  9.7884E-01  1.0000E-02  1.3631E+00  1.6510E+00  8.3944E-01  9.3410E-01  1.1983E+01  8.7186E-01  7.4000E-01  9.2685E-01
             9.5748E-01
 PARAMETER:  7.8612E-02 -4.6933E+00  4.0973E-01  6.0135E-01 -7.5022E-02  3.1833E-02  2.5835E+00 -3.7124E-02 -2.0111E-01  2.4032E-02
             5.6551E-02
 GRADIENT:  -3.4440E-02  0.0000E+00  4.4921E-01 -1.0552E+00 -5.1463E-01  1.2512E-02 -1.1384E-01 -5.6853E-02 -1.9140E-01  2.3304E-01
             1.5816E-01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1703.37294607060        NO. OF FUNC. EVALS.:  64
 CUMULATIVE NO. OF FUNC. EVALS.:     1958
 NPARAMETR:  9.7884E-01  1.0000E-02  1.3640E+00  1.6525E+00  8.3966E-01  9.3407E-01  1.1957E+01  8.7520E-01  7.4056E-01  9.2556E-01
             9.5714E-01
 PARAMETER:  7.8640E-02 -4.6931E+00  4.1000E-01  6.0169E-01 -7.4723E-02  3.1782E-02  2.5839E+00 -3.3402E-02 -2.0051E-01  2.2577E-02
             5.6161E-02
 GRADIENT:   4.0227E-02  0.0000E+00 -2.0449E+03 -1.3880E+03  1.6225E-01 -7.9150E-03  3.2113E+02 -6.0942E-02 -6.1696E-01 -2.3191E-01
            -2.3200E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1958
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.1580E-04  9.4092E-03 -2.1732E-02 -8.0693E-03 -2.9025E-02
 SE:             2.9834E-02  6.3140E-03  1.4920E-02  2.8941E-02  2.1807E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9690E-01  1.3617E-01  1.4522E-01  7.8039E-01  1.8319E-01

 ETASHRINKSD(%)  5.2308E-02  7.8847E+01  5.0018E+01  3.0429E+00  2.6945E+01
 ETASHRINKVR(%)  1.0459E-01  9.5526E+01  7.5018E+01  5.9932E+00  4.6629E+01
 EBVSHRINKSD(%)  4.3233E-01  8.4440E+01  5.1743E+01  3.0659E+00  2.3875E+01
 EBVSHRINKVR(%)  8.6280E-01  9.7579E+01  7.6713E+01  6.0378E+00  4.2049E+01
 RELATIVEINF(%)  9.9086E+01  1.9839E+00  3.7813E+00  7.1005E+01  9.4132E+00
 EPSSHRINKSD(%)  4.3889E+01
 EPSSHRINKVR(%)  6.8515E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1703.3729460706047     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.22211950686653     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.50
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.89
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1703.373       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.79E-01  1.00E-02  1.36E+00  1.65E+00  8.40E-01  9.34E-01  1.20E+01  8.75E-01  7.40E-01  9.25E-01  9.57E-01
 


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
+        1.41E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+        3.39E+01  0.00E+00  6.71E+05
 
 TH 4
+        7.35E+00  0.00E+00 -3.39E+02  2.12E+05
 
 TH 5
+        1.35E+02  0.00E+00  1.35E+03  7.62E+02  2.97E+07
 
 TH 6
+        8.36E+00  0.00E+00 -2.67E+01 -1.56E+01 -6.72E+01  3.31E+02
 
 TH 7
+       -4.77E-01  0.00E+00  8.60E+00  3.18E+01 -2.34E+01  4.13E-01  2.18E+02
 
 TH 8
+       -3.05E+01  0.00E+00 -4.23E+02 -1.80E+02 -1.48E+03  7.42E+00  5.01E+00  4.71E+02
 
 TH 9
+       -3.81E+02  0.00E+00  2.53E+06 -1.71E+02 -1.68E+07  2.94E+02  2.99E+00  1.61E+07  9.51E+06
 
 TH10
+       -8.09E+01  0.00E+00 -1.32E+03 -5.97E+02 -2.70E+07  7.57E+01  1.68E+01  1.22E+03  1.53E+07  2.45E+07
 
 TH11
+       -6.69E+01  0.00E+00 -2.90E+02 -1.30E+02 -1.02E+03  5.15E+01  3.33E+00  2.46E+02  1.48E+07  7.94E+02  4.53E+02
 
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
 #CPUT: Total CPU Time in Seconds,       33.402
Stop Time:
Sat Sep 25 13:45:52 CDT 2021
