Sat Sep 25 02:50:51 CDT 2021
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
$DATA ../../../../data/int/SL3/dat92.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      979
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

 TOT. NO. OF OBS RECS:      879
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
 RAW OUTPUT FILE (FILE): m92.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -617.524636050453        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   6.9679E+01 -4.3278E+01  1.4998E+02  2.7674E+01  1.4379E+02  2.4818E+01 -1.3705E+02 -2.3738E+02 -1.4071E+02 -2.2786E+01
            -5.8302E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2690.35476864519        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0303E+00  1.2619E+00  9.3652E-01  8.8811E-01  1.0702E+00  8.8872E-01  1.2343E+00  1.0191E+00  1.1475E+00  9.6565E-01
             2.6127E+00
 PARAMETER:  1.2990E-01  3.3262E-01  3.4414E-02 -1.8664E-02  1.6784E-01 -1.7968E-02  3.1050E-01  1.1889E-01  2.3761E-01  6.5041E-02
             1.0604E+00
 GRADIENT:   7.7549E+01 -4.3715E+00 -2.2793E+00 -1.1414E+01 -1.7604E+01 -1.7961E+01  7.4484E+00  1.0271E+00 -2.6911E+00  8.5853E-02
            -1.8921E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2692.73498485776        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0212E+00  1.6617E+00  1.1123E+00  7.4199E-01  1.3442E+00  9.7403E-01  1.0472E+00  1.0604E+00  1.3032E+00  1.2632E+00
             2.6104E+00
 PARAMETER:  1.2101E-01  6.0783E-01  2.0643E-01 -1.9841E-01  3.9578E-01  7.3685E-02  1.4614E-01  1.5866E-01  3.6480E-01  3.3368E-01
             1.0595E+00
 GRADIENT:   4.4204E+01  9.7884E+01  1.2591E+01  5.3975E+01 -3.0417E+01  1.6428E+01  4.7958E+00 -4.0945E+00 -2.8782E+00  3.9936E+00
            -2.0645E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2700.24485243665        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.9410E-01  1.8109E+00  1.1029E+00  5.7599E-01  1.5785E+00  9.3400E-01  8.7763E-01  1.5915E+00  1.5723E+00  1.3841E+00
             2.6275E+00
 PARAMETER:  9.4080E-02  6.9382E-01  1.9793E-01 -4.5166E-01  5.5647E-01  3.1721E-02 -3.0530E-02  5.6468E-01  5.5254E-01  4.2505E-01
             1.0660E+00
 GRADIENT:  -2.0031E+01  7.6616E+00 -4.9214E-01  1.6862E+01  1.9513E+01  2.1057E+00 -4.6991E+00 -1.2596E+00  2.7858E+00  1.8182E+00
             1.3025E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2702.95904815534        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      349
 NPARAMETR:  1.0054E+00  2.0636E+00  8.2835E-01  4.2510E-01  1.6827E+00  9.2855E-01  8.6673E-01  1.6570E+00  1.7031E+00  1.4881E+00
             2.6032E+00
 PARAMETER:  1.0540E-01  8.2443E-01 -8.8321E-02 -7.5544E-01  6.2037E-01  2.5870E-02 -4.3028E-02  6.0503E-01  6.3245E-01  4.9753E-01
             1.0567E+00
 GRADIENT:   1.3497E+00  2.8735E+01  1.0660E+00  1.2871E+01  8.8788E-01 -5.6153E-01 -2.0883E+00 -1.4009E-01 -3.6328E+00  2.6458E+00
            -9.3769E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2706.25678265847        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      524
 NPARAMETR:  1.0040E+00  2.3585E+00  3.3379E-01  2.0203E-01  1.8582E+00  9.2546E-01  8.1166E-01  7.7857E-01  2.5815E+00  1.6188E+00
             2.6060E+00
 PARAMETER:  1.0395E-01  9.5803E-01 -9.9724E-01 -1.4994E+00  7.1963E-01  2.2533E-02 -1.0868E-01 -1.5029E-01  1.0484E+00  5.8167E-01
             1.0578E+00
 GRADIENT:  -1.9093E+00 -7.0656E+00 -2.6211E+00  3.3958E+00  4.3313E+00 -1.8011E+00 -7.8513E-01  7.3661E-01  1.4174E+00  2.7371E+00
             3.0208E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2706.72966561277        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      700
 NPARAMETR:  1.0044E+00  2.4186E+00  3.3978E-01  1.6285E-01  1.9062E+00  9.3034E-01  8.0695E-01  4.0490E-01  2.7646E+00  1.6401E+00
             2.6064E+00
 PARAMETER:  1.0439E-01  9.8319E-01 -9.7946E-01 -1.7149E+00  7.4513E-01  2.7797E-02 -1.1450E-01 -8.0412E-01  1.1169E+00  5.9477E-01
             1.0580E+00
 GRADIENT:  -2.1088E-01 -1.0628E+00 -2.8983E-02 -4.3971E-01 -1.0883E+00  3.4225E-01 -3.5916E-01  2.0164E-01  1.4656E-01 -2.0887E-01
             9.6143E-02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2707.79978500584        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      879
 NPARAMETR:  1.0041E+00  2.5137E+00  3.2046E-01  1.0528E-01  1.9950E+00  9.3217E-01  8.0300E-01  9.2468E-02  3.4278E+00  1.7081E+00
             2.5959E+00
 PARAMETER:  1.0404E-01  1.0218E+00 -1.0380E+00 -2.1512E+00  7.9064E-01  2.9755E-02 -1.1941E-01 -2.2809E+00  1.3319E+00  6.3541E-01
             1.0539E+00
 GRADIENT:  -4.8510E-01  1.0057E+01 -1.1135E+00  1.4639E+00  1.1239E+00  9.7926E-01  2.7009E+00  1.0321E-02 -1.2027E+00  5.9565E-01
            -7.0992E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2707.99571803161        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1054
 NPARAMETR:  1.0045E+00  2.5414E+00  2.8100E-01  7.8180E-02  2.0171E+00  9.2797E-01  7.9339E-01  3.1701E-02  4.0050E+00  1.7199E+00
             2.6051E+00
 PARAMETER:  1.0453E-01  1.0327E+00 -1.1694E+00 -2.4487E+00  8.0166E-01  2.5245E-02 -1.3144E-01 -3.3514E+00  1.4876E+00  6.4228E-01
             1.0575E+00
 GRADIENT:   6.9492E-01 -6.3234E+00 -1.4444E+00  1.7997E+00 -3.6861E-01 -4.4214E-01 -1.4864E+00  1.1404E-03  2.6090E+00 -5.9603E-01
             3.7881E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2708.03393407923        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1229
 NPARAMETR:  1.0042E+00  2.5439E+00  2.9718E-01  7.7576E-02  2.0172E+00  9.2935E-01  7.9551E-01  2.8042E-02  4.0587E+00  1.7181E+00
             2.5990E+00
 PARAMETER:  1.0419E-01  1.0337E+00 -1.1134E+00 -2.4565E+00  8.0171E-01  2.6736E-02 -1.2877E-01 -3.4741E+00  1.5009E+00  6.4124E-01
             1.0551E+00
 GRADIENT:   6.5055E-03 -1.1523E+00 -1.3349E-01  6.4443E-01 -4.6280E-01  3.3811E-02 -7.2166E-01  8.8764E-04  1.2202E+00 -3.2339E-01
            -3.4522E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2708.03442842934        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1406
 NPARAMETR:  1.0043E+00  2.5440E+00  2.9633E-01  7.7538E-02  2.0179E+00  9.2928E-01  7.9572E-01  1.8355E-02  4.0576E+00  1.7184E+00
             2.5994E+00
 PARAMETER:  1.0424E-01  1.0337E+00 -1.1163E+00 -2.4570E+00  8.0203E-01  2.6652E-02 -1.2851E-01 -3.8979E+00  1.5006E+00  6.4140E-01
             1.0553E+00
 GRADIENT:   8.2558E-02  8.6590E-02 -1.1219E-02 -2.9747E-02  1.6386E-01 -1.7112E-02  1.4865E-01  3.7944E-04 -8.2247E-02  2.8787E-02
             1.8615E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2708.03459851616        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1581
 NPARAMETR:  1.0042E+00  2.5440E+00  2.9641E-01  7.7512E-02  2.0174E+00  9.2931E-01  7.9538E-01  1.0000E-02  4.0573E+00  1.7182E+00
             2.5993E+00
 PARAMETER:  1.0422E-01  1.0337E+00 -1.1160E+00 -2.4573E+00  8.0183E-01  2.6690E-02 -1.2893E-01 -4.5088E+00  1.5005E+00  6.4127E-01
             1.0552E+00
 GRADIENT:   3.1196E-02  1.2487E-01 -1.0718E-03 -2.2906E-02  2.1514E-02 -1.2441E-03 -2.7173E-02  1.0973E-05 -5.8894E-02 -6.0708E-03
             3.2075E-02

0ITERATION NO.:   56    OBJECTIVE VALUE:  -2708.03459851616        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     1603
 NPARAMETR:  1.0042E+00  2.5440E+00  2.9641E-01  7.7512E-02  2.0174E+00  9.2931E-01  7.9538E-01  1.0000E-02  4.0573E+00  1.7182E+00
             2.5993E+00
 PARAMETER:  1.0422E-01  1.0337E+00 -1.1160E+00 -2.4573E+00  8.0183E-01  2.6690E-02 -1.2893E-01 -4.5088E+00  1.5005E+00  6.4127E-01
             1.0552E+00
 GRADIENT:   3.1196E-02  1.2487E-01 -1.0718E-03 -2.2906E-02  2.1514E-02 -1.2441E-03 -2.7173E-02  1.0973E-05 -5.8894E-02 -6.0708E-03
             3.2075E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1603
 NO. OF SIG. DIGITS IN FINAL EST.:  3.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3278E-03 -2.0757E-02 -9.6833E-06  3.2378E-02 -2.3202E-02
 SE:             2.9354E-02  2.7176E-02  1.2137E-05  1.4243E-02  2.5868E-02
 N:                     100         100         100         100         100

 P VAL.:         9.6392E-01  4.4499E-01  4.2496E-01  2.3016E-02  3.6975E-01

 ETASHRINKSD(%)  1.6616E+00  8.9569E+00  9.9959E+01  5.2283E+01  1.3339E+01
 ETASHRINKVR(%)  3.2955E+00  1.7112E+01  1.0000E+02  7.7231E+01  2.4898E+01
 EBVSHRINKSD(%)  1.8357E+00  6.7782E+00  9.9925E+01  6.6889E+01  9.1608E+00
 EBVSHRINKVR(%)  3.6376E+00  1.3097E+01  1.0000E+02  8.9037E+01  1.7482E+01
 RELATIVEINF(%)  9.6312E+01  3.9737E+01  4.8460E-05  4.6667E+00  6.7936E+01
 EPSSHRINKSD(%)  1.6569E+01
 EPSSHRINKVR(%)  3.0393E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          879
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1615.4939413738146     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2708.0345985161575     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1092.5406571423430     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    40.92
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.75
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2708.035       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  2.54E+00  2.96E-01  7.75E-02  2.02E+00  9.29E-01  7.95E-01  1.00E-02  4.06E+00  1.72E+00  2.60E+00
 


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
+        1.24E+03
 
 TH 2
+       -2.98E+01  5.39E+02
 
 TH 3
+       -1.73E+01  2.64E+02  3.47E+02
 
 TH 4
+        3.26E+02 -5.73E+03 -4.61E+03  2.91E+06
 
 TH 5
+       -9.32E+00  7.86E+01  7.90E+01 -1.54E+03  1.07E+02
 
 TH 6
+        9.87E-01 -1.24E+01 -8.02E+00  1.61E+02 -4.48E+00  2.19E+02
 
 TH 7
+       -1.72E+01  3.92E+02  3.46E+02 -6.99E+03  1.33E+02 -1.75E+01  8.67E+02
 
 TH 8
+        4.17E+00  3.57E-01  8.63E+00 -9.99E+00  8.77E-01  2.06E+01  1.13E+00  5.04E+01
 
 TH 9
+        1.35E+01 -6.68E+03 -1.75E+02  9.17E+04 -6.30E+01  6.34E+00 -2.75E+02 -2.64E-01  2.89E+03
 
 TH10
+       -7.19E+00  8.13E+01  8.28E+01 -1.44E+03  2.21E+01 -2.96E+00  1.18E+02  2.08E-01 -5.73E+01  7.09E+01
 
 TH11
+       -1.96E+01  2.07E+01  3.71E+01 -5.77E+02  1.08E+01  1.75E+00  4.99E+01 -3.96E-01 -2.10E+01  1.54E+01  1.73E+02
 
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
 #CPUT: Total CPU Time in Seconds,       54.774
Stop Time:
Sat Sep 25 02:51:47 CDT 2021
