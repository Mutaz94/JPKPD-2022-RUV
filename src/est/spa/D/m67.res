Wed Sep 29 20:14:59 CDT 2021
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
$DATA ../../../../data/spa/D/dat67.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m67.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   15227.1989436768        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.0094E+02  2.8144E+02 -7.6928E+01 -7.1816E+01  4.0136E+02 -2.5079E+03 -7.2280E+02 -4.4213E+01 -1.7706E+03 -8.5888E+02
            -2.7517E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -566.801082602098        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.3029E+00  9.3914E-01  8.3956E-01  1.7736E+00  1.5329E+00  2.4735E+00  1.2174E+00  9.6763E-01  1.3558E+00  1.0962E+00
             1.3827E+01
 PARAMETER:  3.6459E-01  3.7206E-02 -7.4877E-02  6.7300E-01  5.2717E-01  1.0056E+00  2.9671E-01  6.7099E-02  4.0442E-01  1.9182E-01
             2.7266E+00
 GRADIENT:   1.5559E+00  1.7996E+01 -7.0199E+00  1.4495E+01 -7.4840E+00  6.3566E+01 -8.8470E-01  6.9222E+00 -1.6706E+00  3.8117E-01
             7.8762E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -595.164343037097        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.3432E+00  9.0312E-01  1.0108E+00  2.1960E+00  6.3617E+00  2.4042E+00  2.8574E+00  1.8517E-01  2.7808E+00  7.9709E+00
             1.1241E+01
 PARAMETER:  3.9505E-01 -1.8977E-03  1.1072E-01  8.8662E-01  1.9503E+00  9.7721E-01  1.1499E+00 -1.5865E+00  1.1228E+00  2.1758E+00
             2.5195E+00
 GRADIENT:   1.8347E+01  2.5084E+01  2.7739E+00  6.2591E+01 -7.6246E+00 -2.3202E+01  1.0482E+01 -1.4525E-02  3.2056E+01  3.0425E+00
             4.6666E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -606.616095609131        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  1.2429E+00  7.5829E-01  8.0306E-01  1.8582E+00  3.9631E+00  2.4305E+00  1.4260E+00  9.0112E-02  2.3704E+00  8.0150E+00
             1.1460E+01
 PARAMETER:  3.1742E-01 -1.7669E-01 -1.1933E-01  7.1961E-01  1.4770E+00  9.8812E-01  4.5490E-01 -2.3067E+00  9.6305E-01  2.1813E+00
             2.5389E+00
 GRADIENT:  -6.9643E+00  2.0229E+01  2.3333E-01  2.4059E+01 -1.4677E+01 -2.9483E+00  5.4282E+00  9.4172E-03  1.9362E+01  7.0065E-01
             7.4625E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -726.762064356151        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      305
 NPARAMETR:  7.6319E-01  5.0848E-02  5.7652E-02  6.2753E-01  6.5158E+02  1.9578E+00  1.0000E-02  1.0000E-02  9.6978E-01  9.4348E+00
             9.6834E+00
 PARAMETER: -1.7025E-01 -2.8789E+00 -2.7533E+00 -3.6596E-01  6.5794E+00  7.7183E-01 -4.5728E+00 -4.8271E+00  6.9314E-02  2.3444E+00
             2.3704E+00
 GRADIENT:   9.0533E+01  1.5239E+01 -5.7131E+01  1.0091E+02 -8.0662E-03 -1.5845E+01  0.0000E+00  0.0000E+00 -1.4400E+00  2.2386E-03
            -2.5453E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -736.282710025855        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      464
 NPARAMETR:  6.5450E-01  3.8776E-02  4.9368E-02  5.2015E-01  9.8796E+02  2.0618E+00  1.0285E-02  1.0000E-02  7.8544E-01  9.5648E+00
             9.1922E+00
 PARAMETER: -3.2388E-01 -3.1499E+00 -2.9085E+00 -5.5363E-01  6.9956E+00  8.2357E-01 -4.4770E+00 -5.0290E+00 -1.4151E-01  2.3581E+00
             2.3184E+00
 GRADIENT:   3.4802E+01  6.2094E+00 -2.9785E+01  3.5766E+01 -1.8489E-03 -1.6177E+00  4.9264E-04  0.0000E+00 -2.3455E+01  3.9675E-04
            -6.7177E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -746.521607544922        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      639
 NPARAMETR:  4.7734E-01  1.5531E-02  2.5164E-02  3.0598E-01  5.0575E+03  2.0315E+00  1.0000E-02  1.0000E-02  8.8478E-01  1.1870E+01
             9.5902E+00
 PARAMETER: -6.3953E-01 -4.0649E+00 -3.5824E+00 -1.0842E+00  8.6286E+00  8.0879E-01 -8.5444E+00 -6.3784E+00 -2.2411E-02  2.5740E+00
             2.3607E+00
 GRADIENT:   7.7743E-01  1.7182E+00 -3.9821E+00  1.0163E+00 -1.0873E-04  1.1268E+01  0.0000E+00  0.0000E+00 -4.5767E-01  8.0469E-07
            -2.0405E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -746.849401269118        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:      823
 NPARAMETR:  4.7943E-01  1.1729E-02  2.5400E-02  3.0836E-01  3.2649E+04  1.9732E+00  1.0000E-02  1.0000E-02  8.7947E-01  1.1566E+01
             9.6222E+00
 PARAMETER: -6.3516E-01 -4.3457E+00 -3.5730E+00 -1.0765E+00  1.0494E+01  7.7967E-01 -8.6111E+00 -6.0741E+00 -2.8436E-02  2.5481E+00
             2.3641E+00
 GRADIENT:   2.1832E+00  1.0290E-01 -4.4366E+00  3.8562E+00  2.3613E-05  4.8787E-01  0.0000E+00  0.0000E+00 -8.3532E-01  6.6925E-11
            -1.6899E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -746.890007669955        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1005
 NPARAMETR:  4.7766E-01  1.0655E-02  2.5474E-02  3.0841E-01  1.3215E+01  1.9704E+00  1.0000E-02  1.0000E-02  8.8293E-01  7.2349E+00
             9.6351E+00
 PARAMETER: -6.3886E-01 -4.4418E+00 -3.5701E+00 -1.0763E+00  2.6813E+00  7.7824E-01 -8.6111E+00 -6.0741E+00 -2.4511E-02  2.0789E+00
             2.3654E+00
 GRADIENT:  -2.0959E+00  1.9246E-01  1.3597E+00 -1.5739E+00 -5.8596E-02 -3.0613E-01  0.0000E+00  0.0000E+00  2.2700E-01  2.0088E-02
             8.0781E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -746.912480711608        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1184
 NPARAMETR:  4.7883E-01  1.0000E-02  2.5421E-02  3.0821E-01  1.2298E+01  1.9720E+00  1.0000E-02  1.0000E-02  8.8205E-01  2.9133E+00
             9.6314E+00
 PARAMETER: -6.3640E-01 -4.5440E+00 -3.5722E+00 -1.0770E+00  2.6094E+00  7.7905E-01 -8.6111E+00 -6.0741E+00 -2.5502E-02  1.1693E+00
             2.3650E+00
 GRADIENT:  -2.1536E-01  5.3486E-03 -3.7028E-01 -2.3315E-01  6.2675E-04 -8.3872E-02  0.0000E+00  0.0000E+00  4.7265E-03  1.6419E-03
             2.3212E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -746.915792075482        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1370
 NPARAMETR:  4.7933E-01  1.0000E-02  2.5393E-02  3.0810E-01  1.1988E+01  1.9739E+00  1.0000E-02  1.0000E-02  8.8229E-01  1.3294E+00
             9.6304E+00
 PARAMETER: -6.3537E-01 -4.5718E+00 -3.5733E+00 -1.0773E+00  2.5839E+00  7.8000E-01 -8.6111E+00 -6.0741E+00 -2.5236E-02  3.8473E-01
             2.3649E+00
 GRADIENT:   5.1335E-01  0.0000E+00 -1.1793E+00  3.8434E-01  2.3471E-03  2.5669E-01  0.0000E+00  0.0000E+00 -5.9865E-03  4.3066E-04
            -2.0180E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -746.923558918027        NO. OF FUNC. EVALS.: 129
 CUMULATIVE NO. OF FUNC. EVALS.:     1499
 NPARAMETR:  4.7634E-01  1.0000E-02  2.5153E-02  3.0650E-01  1.1808E+01  1.9694E+00  1.0000E-02  1.0000E-02  8.8292E-01  1.0344E+00
             9.6149E+00
 PARAMETER: -6.4163E-01 -4.5718E+00 -3.5828E+00 -1.0826E+00  2.5688E+00  7.7772E-01 -8.6111E+00 -6.0741E+00 -2.4525E-02  1.3383E-01
             2.3633E+00
 GRADIENT:  -8.2398E-01  0.0000E+00 -3.6346E+00  4.0125E+00 -7.0027E-03 -4.7861E-01  0.0000E+00  0.0000E+00 -1.3010E-01  3.3994E-04
            -1.1696E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -746.930840208390        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1693             RESET HESSIAN, TYPE I
 NPARAMETR:  4.7795E-01  1.0000E-02  2.5178E-02  3.0605E-01  1.1970E+01  1.9732E+00  1.0000E-02  1.0000E-02  8.8397E-01  1.0163E+00
             9.6311E+00
 PARAMETER: -6.3825E-01 -4.5718E+00 -3.5818E+00 -1.0840E+00  2.5824E+00  7.7967E-01 -8.6111E+00 -6.0741E+00 -2.3335E-02  1.1621E-01
             2.3650E+00
 GRADIENT:   4.8333E+01  0.0000E+00  6.3441E+01  2.4327E+01  1.0266E-02  2.0243E+01  0.0000E+00  0.0000E+00  4.7310E-01  2.9087E-04
             1.8922E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -746.933558644151        NO. OF FUNC. EVALS.: 194
 CUMULATIVE NO. OF FUNC. EVALS.:     1887
 NPARAMETR:  4.7761E-01  1.0000E-02  2.5143E-02  3.0561E-01  1.1923E+01  1.9731E+00  1.0000E-02  1.0000E-02  8.8354E-01  1.0129E+00
             9.6317E+00
 PARAMETER: -6.3896E-01 -4.5718E+00 -3.5832E+00 -1.0855E+00  2.5785E+00  7.7959E-01 -8.6111E+00 -6.0741E+00 -2.3822E-02  1.1284E-01
             2.3651E+00
 GRADIENT:   5.9585E-01  0.0000E+00 -8.5589E-01 -2.9910E-01 -3.8744E-03  3.5392E-01  0.0000E+00  0.0000E+00  2.3107E-01  2.9281E-04
             9.5409E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -746.938689411439        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     2082
 NPARAMETR:  4.7751E-01  1.0000E-02  2.5016E-02  3.0488E-01  1.2056E+01  1.9741E+00  1.0000E-02  1.0000E-02  8.8315E-01  1.0024E+00
             9.6284E+00
 PARAMETER: -6.3917E-01 -4.5718E+00 -3.5882E+00 -1.0878E+00  2.5896E+00  7.8012E-01 -8.6111E+00 -6.0741E+00 -2.4256E-02  1.0238E-01
             2.3647E+00
 GRADIENT:   1.5457E+00  0.0000E+00 -3.2289E+00  2.1288E+00  8.0875E-03  6.0848E-01  0.0000E+00  0.0000E+00  3.0533E-02  2.5686E-04
            -4.8237E-01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -746.941350397399        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:     2240
 NPARAMETR:  4.7512E-01  1.0000E-02  2.4982E-02  3.0425E-01  1.1903E+01  1.9681E+00  1.0000E-02  1.0000E-02  8.8258E-01  8.9325E-01
             9.6222E+00
 PARAMETER: -6.4419E-01 -4.5718E+00 -3.5896E+00 -1.0899E+00  2.5768E+00  7.7705E-01 -8.6111E+00 -6.0741E+00 -2.4903E-02 -1.2888E-02
             2.3641E+00
 GRADIENT:  -8.1909E-01  0.0000E+00 -1.3203E+00  9.1600E-01 -1.0900E-02 -4.6426E-01  0.0000E+00  0.0000E+00  1.4596E-02  2.4160E-04
            -3.9450E-01

0ITERATION NO.:   80    OBJECTIVE VALUE:  -746.947651207255        NO. OF FUNC. EVALS.: 193
 CUMULATIVE NO. OF FUNC. EVALS.:     2433
 NPARAMETR:  4.7584E-01  1.0000E-02  2.4931E-02  3.0355E-01  1.1910E+01  1.9714E+00  1.0000E-02  1.0000E-02  8.8236E-01  8.9157E-01
             9.6320E+00
 PARAMETER: -6.4268E-01 -4.5718E+00 -3.5916E+00 -1.0922E+00  2.5774E+00  7.7875E-01 -8.6111E+00 -6.0741E+00 -2.5158E-02 -1.4772E-02
             2.3651E+00
 GRADIENT:   3.3159E-01  0.0000E+00 -8.3165E-01 -3.7850E-01 -8.5983E-03  2.6348E-01  0.0000E+00  0.0000E+00  7.0985E-02  2.4559E-04
             1.8118E-01

0ITERATION NO.:   81    OBJECTIVE VALUE:  -746.947651207255        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     2455
 NPARAMETR:  4.7584E-01  1.0000E-02  2.4931E-02  3.0355E-01  1.1910E+01  1.9714E+00  1.0000E-02  1.0000E-02  8.8236E-01  8.9157E-01
             9.6320E+00
 PARAMETER: -6.4268E-01 -4.5718E+00 -3.5916E+00 -1.0922E+00  2.5774E+00  7.7875E-01 -8.6111E+00 -6.0741E+00 -2.5158E-02 -1.4772E-02
             2.3651E+00
 GRADIENT:   3.3159E-01  0.0000E+00 -8.3165E-01 -3.7850E-01 -8.5983E-03  2.6348E-01  0.0000E+00  0.0000E+00  7.0985E-02  2.4559E-04
             1.8118E-01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     2455
 NO. OF SIG. DIGITS IN FINAL EST.:  2.1
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.6404E-03  1.2107E-06  1.2195E-04 -2.2646E-02 -5.7433E-06
 SE:             2.9139E-02  8.1509E-07  2.2464E-04  2.2807E-02  7.4517E-05
 N:                     100         100         100         100         100

 P VAL.:         9.5511E-01  1.3744E-01  5.8722E-01  3.2072E-01  9.3857E-01

 ETASHRINKSD(%)  2.3795E+00  9.9997E+01  9.9247E+01  2.3595E+01  9.9750E+01
 ETASHRINKVR(%)  4.7024E+00  1.0000E+02  9.9994E+01  4.1623E+01  9.9999E+01
 EBVSHRINKSD(%)  2.1204E+00  9.9996E+01  9.9302E+01  2.4376E+01  9.9734E+01
 EBVSHRINKVR(%)  4.1959E+00  1.0000E+02  9.9995E+01  4.2810E+01  9.9999E+01
 RELATIVEINF(%)  5.8181E+00  3.6517E-08  4.6865E-05  4.8857E-01  5.8511E-05
 EPSSHRINKSD(%)  1.3261E+01
 EPSSHRINKVR(%)  2.4763E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -746.94765120725481     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -11.796824643516629     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    41.19
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.19
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -746.948       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.76E-01  1.00E-02  2.49E-02  3.04E-01  1.19E+01  1.97E+00  1.00E-02  1.00E-02  8.82E-01  8.92E-01  9.63E+00
 


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
+        1.20E+03
 
 TH 2
+        0.00E+00  1.69E+03
 
 TH 3
+       -7.74E+03  0.00E+00  8.16E+05
 
 TH 4
+       -1.04E+02  0.00E+00 -7.84E+04  8.29E+03
 
 TH 5
+        1.96E-01  0.00E+00 -4.45E+00  3.86E-01  1.91E-03
 
 TH 6
+        2.03E+00  0.00E+00 -3.57E+00 -2.63E+01  4.46E-03  4.59E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+       -4.73E+00  0.00E+00  1.26E+03 -1.30E+02 -1.18E-02 -1.60E+00  0.00E+00  0.00E+00  8.79E+01
 
 TH10
+        2.63E-03  0.00E+00  1.42E-03 -2.53E-03 -9.20E-05  1.30E-03  0.00E+00  0.00E+00  4.71E-02  2.30E-02
 
 TH11
+       -1.34E+01  0.00E+00  2.43E+02 -1.76E+01 -3.22E-03  9.05E-01  0.00E+00  0.00E+00  3.41E+00  2.22E-04  3.80E+00
 
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
 #CPUT: Total CPU Time in Seconds,       49.443
Stop Time:
Wed Sep 29 20:15:50 CDT 2021
