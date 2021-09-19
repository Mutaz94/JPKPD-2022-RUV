Sat Sep 18 10:45:30 CDT 2021
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
$DATA ../../../../data/spa/A3/dat85.csv ignore=@
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
 (E4.0,E3.0,E21.0,E4.0,2E2.0)

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


0ITERATION NO.:    0    OBJECTIVE VALUE:   446.879066341616        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4203E+02  3.4282E+01  9.2753E+01 -8.6023E+01  1.3667E+02  3.2394E+01 -8.8111E+01 -3.8740E+01 -1.8461E+02 -1.6076E+02
            -3.7095E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1180.10840598557        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0262E+00  1.0196E+00  9.5231E-01  1.2019E+00  1.0097E+00  7.5881E-01  1.0160E+00  9.8575E-01  1.1098E+00  9.8344E-01
             5.3461E+00
 PARAMETER:  1.2587E-01  1.1937E-01  5.1137E-02  2.8392E-01  1.0965E-01 -1.7600E-01  1.1583E-01  8.5645E-02  2.0420E-01  8.3305E-02
             1.7764E+00
 GRADIENT:   2.0896E+01  2.7560E+00 -2.0477E+01  2.7523E+01 -2.0657E+00 -1.6850E+01  8.7024E+00  6.7425E+00  1.7559E+01  1.7581E+01
             1.6037E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1200.55991271954        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  9.9662E-01  8.4515E-01  4.4328E-01  1.1847E+00  5.4638E-01  8.8113E-01  4.1302E-01  1.6225E-01  1.1719E+00  3.3737E-01
             4.7874E+00
 PARAMETER:  9.6619E-02 -6.8235E-02 -7.1355E-01  2.6951E-01 -5.0444E-01 -2.6552E-02 -7.8425E-01 -1.7186E+00  2.5859E-01 -9.8658E-01
             1.6660E+00
 GRADIENT:  -3.6280E+01  3.2566E+01  2.4884E+00  3.9302E+01 -3.5839E+01  1.2163E+01  7.3641E-02  3.9999E-01  1.1595E+01  3.6158E+00
             1.0242E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1218.23096698193        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      225
 NPARAMETR:  9.9449E-01  6.4991E-01  9.4379E-01  1.3382E+00  8.3255E-01  7.9839E-01  5.0078E-01  2.2739E-01  1.0154E+00  2.2888E-01
             4.2170E+00
 PARAMETER:  9.4472E-02 -3.3092E-01  4.2152E-02  3.9131E-01 -8.3256E-02 -1.2516E-01 -5.9160E-01 -1.3811E+00  1.1529E-01 -1.3745E+00
             1.5391E+00
 GRADIENT:   1.7688E+01  5.9572E-01 -7.8853E+00  1.0832E+01  7.8933E+00 -1.0818E+01  4.1852E-01  6.5707E-01  2.1023E+00  9.6210E-01
             5.4557E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1219.50971498561        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      295
 NPARAMETR:  9.8394E-01  4.9844E-01  9.8640E-01  1.4105E+00  7.9807E-01  8.2304E-01  2.2973E-01  5.1748E-02  9.6141E-01  7.4496E-02
             4.1832E+00
 PARAMETER:  8.3813E-02 -5.9627E-01  8.6310E-02  4.4397E-01 -1.2556E-01 -9.4752E-02 -1.3708E+00 -2.8614E+00  6.0645E-02 -2.4970E+00
             1.5311E+00
 GRADIENT:   1.7874E-01 -1.4666E+00 -1.1746E-01 -2.6336E+00  1.6981E+00 -2.3665E-01  5.9551E-02  3.2933E-02  2.0787E+00  9.0293E-02
             1.6592E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1219.59242755462        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      367
 NPARAMETR:  9.8351E-01  4.3262E-01  9.4787E-01  1.4459E+00  7.5682E-01  8.2469E-01  1.3805E-01  2.1562E-02  9.2887E-01  3.4725E-02
             4.1771E+00
 PARAMETER:  8.3376E-02 -7.3790E-01  4.6458E-02  4.6873E-01 -1.7863E-01 -9.2747E-02 -1.8801E+00 -3.7368E+00  2.6210E-02 -3.2603E+00
             1.5296E+00
 GRADIENT:   8.4410E-01 -2.3787E-01 -1.3891E-01 -1.4848E-01  2.1009E-01  6.6377E-02  1.2596E-02  6.1773E-03 -1.2218E-01  1.9750E-02
            -9.2764E-02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1219.59358041911        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  9.8265E-01  4.2875E-01  9.4682E-01  1.4471E+00  7.5531E-01  8.2402E-01  1.1565E-01  1.6016E-02  9.3046E-01  2.9658E-02
             4.1772E+00
 PARAMETER:  8.2493E-02 -7.4688E-01  4.5359E-02  4.6955E-01 -1.8063E-01 -9.3556E-02 -2.0572E+00 -4.0342E+00  2.7920E-02 -3.4180E+00
             1.5296E+00
 GRADIENT:  -1.2074E+00 -5.0055E-01 -2.9509E-01 -8.9847E-01  6.0423E-01 -9.8134E-02  9.2612E-03  3.4246E-03  3.1272E-01  1.4495E-02
             4.0238E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1219.60484621365        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      518
 NPARAMETR:  9.8400E-01  5.2564E-01  9.7354E-01  1.3967E+00  7.9832E-01  8.2424E-01  3.0910E-02  1.0000E-02  9.6535E-01  1.6133E-02
             4.1825E+00
 PARAMETER:  8.3872E-02 -5.4314E-01  7.3180E-02  4.3408E-01 -1.2525E-01 -9.3294E-02 -3.3767E+00 -6.4572E+00  6.4735E-02 -4.0269E+00
             1.5309E+00
 GRADIENT:  -2.1052E+00  5.9104E-01  5.9587E-01  1.9505E+00 -1.0837E+00 -1.9705E-01  9.8457E-04  0.0000E+00  7.9317E-02  4.1317E-03
             3.1672E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1219.61782375118        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      676
 NPARAMETR:  9.8480E-01  4.7637E-01  9.7100E-01  1.4262E+00  7.8294E-01  8.2502E-01  2.9304E-02  1.0000E-02  9.4323E-01  1.2438E-02
             4.1878E+00
 PARAMETER:  8.4682E-02 -6.4157E-01  7.0574E-02  4.5499E-01 -1.4470E-01 -9.2349E-02 -3.4300E+00 -6.5058E+00  4.1552E-02 -4.2870E+00
             1.5322E+00
 GRADIENT:  -1.9181E-01 -4.5480E-02 -1.9696E-02 -2.2036E-01  6.8030E-02 -1.4287E-03  6.3489E-04  0.0000E+00  1.7150E-03  2.4585E-03
            -5.5656E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1219.61851236494        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      851
 NPARAMETR:  9.8483E-01  4.7143E-01  9.6932E-01  1.4291E+00  7.8048E-01  8.2505E-01  1.9251E-02  1.0000E-02  9.4135E-01  1.0000E-02
             4.1878E+00
 PARAMETER:  8.4709E-02 -6.5199E-01  6.8836E-02  4.5704E-01 -1.4784E-01 -9.2309E-02 -3.8502E+00 -7.2365E+00  3.9556E-02 -4.6511E+00
             1.5322E+00
 GRADIENT:   1.0202E-02  1.7212E-02  1.7118E-03  4.8468E-02 -1.5980E-02  2.0064E-03  2.6662E-04  0.0000E+00 -2.0070E-02  0.0000E+00
            -5.0059E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1219.61861269220        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1028
 NPARAMETR:  9.8481E-01  4.7109E-01  9.6935E-01  1.4293E+00  7.8040E-01  8.2503E-01  1.0000E-02  1.0000E-02  9.4131E-01  1.0000E-02
             4.1879E+00
 PARAMETER:  8.4698E-02 -6.5272E-01  6.8870E-02  4.5715E-01 -1.4795E-01 -9.2332E-02 -4.5602E+00 -8.4760E+00  3.9513E-02 -5.2271E+00
             1.5322E+00
 GRADIENT:   2.1856E-03  2.1397E-04 -8.9898E-04  7.6539E-04  4.0004E-04  2.2466E-04  0.0000E+00  0.0000E+00 -1.7900E-03  0.0000E+00
            -4.3694E-03

0ITERATION NO.:   51    OBJECTIVE VALUE:  -1219.61861269220        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1050
 NPARAMETR:  9.8481E-01  4.7109E-01  9.6935E-01  1.4293E+00  7.8040E-01  8.2503E-01  1.0000E-02  1.0000E-02  9.4131E-01  1.0000E-02
             4.1879E+00
 PARAMETER:  8.4698E-02 -6.5272E-01  6.8870E-02  4.5715E-01 -1.4795E-01 -9.2332E-02 -4.5602E+00 -8.4760E+00  3.9513E-02 -5.2271E+00
             1.5322E+00
 GRADIENT:   2.1856E-03  2.1397E-04 -8.9898E-04  7.6539E-04  4.0004E-04  2.2466E-04  0.0000E+00  0.0000E+00 -1.7900E-03  0.0000E+00
            -4.3694E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1050
 NO. OF SIG. DIGITS IN FINAL EST.:  3.2
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5500E-03 -1.5839E-04  1.0084E-04 -1.4295E-02 -1.4332E-05
 SE:             2.7723E-02  6.8946E-05  1.0067E-04  2.3865E-02  1.5740E-04
 N:                     100         100         100         100         100

 P VAL.:         9.5541E-01  2.1601E-02  3.1654E-01  5.4916E-01  9.2745E-01

 ETASHRINKSD(%)  7.1246E+00  9.9769E+01  9.9663E+01  2.0051E+01  9.9473E+01
 ETASHRINKVR(%)  1.3742E+01  9.9999E+01  9.9999E+01  3.6081E+01  9.9997E+01
 EBVSHRINKSD(%)  7.0130E+00  9.9778E+01  9.9577E+01  1.9681E+01  9.9403E+01
 EBVSHRINKVR(%)  1.3534E+01  1.0000E+02  9.9998E+01  3.5488E+01  9.9996E+01
 RELATIVEINF(%)  8.0190E+01  1.7282E-05  1.1513E-04  4.0910E+00  1.1930E-04
 EPSSHRINKSD(%)  1.8601E+01
 EPSSHRINKVR(%)  3.3742E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1219.6186126921959     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -484.46778612845776     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    11.45
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1219.619       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.85E-01  4.71E-01  9.69E-01  1.43E+00  7.80E-01  8.25E-01  1.00E-02  1.00E-02  9.41E-01  1.00E-02  4.19E+00
 


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
+        1.46E+03
 
 TH 2
+       -1.06E+02  2.60E+02
 
 TH 3
+        9.18E+00  1.10E+02  1.85E+02
 
 TH 4
+       -1.18E+02  2.66E+02  2.58E+01  3.87E+02
 
 TH 5
+        3.37E+01 -2.77E+02 -3.28E+02 -1.27E+02  6.32E+02
 
 TH 6
+       -1.08E+01 -1.68E+01  1.21E+01 -2.42E+01 -8.67E+00  2.11E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.02E+01 -4.44E+01 -2.63E+00 -5.75E+00  2.18E+01  4.34E+00  0.00E+00  0.00E+00  8.25E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.35E+01 -1.29E+01 -3.18E+00 -9.98E+00  7.69E+00  5.30E+00  0.00E+00  0.00E+00  1.05E+01  0.00E+00  2.59E+01
 
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
 #CPUT: Total CPU Time in Seconds,       18.176
Stop Time:
Sat Sep 18 10:45:49 CDT 2021
