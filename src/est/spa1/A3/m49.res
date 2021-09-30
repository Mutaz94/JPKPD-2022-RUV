Thu Sep 30 00:18:40 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat49.csv ignore=@
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
 NO. OF DATA RECS IN DATA SET:      600
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

 TOT. NO. OF OBS RECS:      500
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
 RAW OUTPUT FILE (FILE): m49.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   490.954027631496        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.4921E+02  5.9256E+01  2.0673E+02 -8.1936E+00  2.2733E+02  2.4430E+01 -1.3917E+02 -3.0001E+02 -1.5953E+02 -1.8147E+02
            -4.2660E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1452.14449748682        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9822E-01  1.0599E+00  1.0567E+00  1.0740E+00  9.9991E-01  9.3600E-01  1.2805E+00  9.7812E-01  1.3235E+00  8.7351E-01
             3.2157E+00
 PARAMETER:  9.8220E-02  1.5815E-01  1.5519E-01  1.7139E-01  9.9911E-02  3.3863E-02  3.4722E-01  7.7879E-02  3.8029E-01 -3.5240E-02
             1.2680E+00
 GRADIENT:   1.6143E+02  1.3292E+01 -1.1279E+00  3.0345E+01 -3.2077E+00 -5.6714E+00  1.3509E+00  4.9596E+00  2.3714E+01  5.7459E+00
            -4.3709E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1460.22467203649        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  9.9252E-01  1.0331E+00  1.1830E+00  1.1049E+00  1.0985E+00  9.6180E-01  1.2335E+00  5.4560E-01  1.1083E+00  3.4862E-01
             3.5319E+00
 PARAMETER:  9.2493E-02  1.3256E-01  2.6808E-01  1.9977E-01  1.9392E-01  6.1055E-02  3.0983E-01 -5.0587E-01  2.0287E-01 -9.5378E-01
             1.3618E+00
 GRADIENT:   1.2778E+02  8.8728E+00 -4.6847E+00  2.8396E+01  7.0967E+00  7.3734E+00 -5.0882E+00  9.5062E-01 -2.0391E+00  8.5380E-01
             1.8799E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1463.26694966008        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  9.4618E-01  9.5278E-01  7.4039E-01  1.0915E+00  8.0678E-01  9.2496E-01  1.4945E+00  3.4374E-01  1.0626E+00  3.0828E-01
             3.3708E+00
 PARAMETER:  4.4679E-02  5.1633E-02 -2.0058E-01  1.8758E-01 -1.1470E-01  2.1991E-02  5.0178E-01 -9.6786E-01  1.6070E-01 -1.0767E+00
             1.3151E+00
 GRADIENT:   2.2459E+01  1.4867E+00 -3.9096E+00  1.6593E+01  7.4860E+00 -3.8091E+00 -6.0535E-01  6.0863E-01 -5.4036E-01  5.1309E-01
             3.9820E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1464.27895959821        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      407
 NPARAMETR:  9.4662E-01  7.3038E-01  7.4162E-01  1.2214E+00  7.0659E-01  9.4602E-01  1.7984E+00  8.8178E-02  1.0215E+00  3.6357E-01
             3.3723E+00
 PARAMETER:  4.5139E-02 -2.1419E-01 -1.9892E-01  3.0001E-01 -2.4731E-01  4.4511E-02  6.8690E-01 -2.3284E+00  1.2131E-01 -9.1177E-01
             1.3156E+00
 GRADIENT:  -4.8743E+00  6.0223E+00  3.9980E+00  6.2489E+00 -6.0728E+00  2.4968E-01 -7.7797E-01  2.9761E-02 -1.4944E+00 -7.1552E-01
            -2.9061E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1486.93898716008        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  9.0980E-01  4.3986E-01  2.3630E-01  1.2184E+00  2.8048E-01  9.4431E-01  2.1831E+00  1.0000E-02  1.2765E+00  6.3742E-01
             2.8901E+00
 PARAMETER:  5.4662E-03 -7.2129E-01 -1.3427E+00  2.9751E-01 -1.1713E+00  4.2700E-02  8.8076E-01 -7.2205E+00  3.4411E-01 -3.5033E-01
             1.1613E+00
 GRADIENT:  -9.5954E+01  3.8243E+01  1.5939E-01  6.8601E+01 -3.7812E+01 -1.4992E+01  2.0576E+01  0.0000E+00 -2.1666E+01 -2.4902E+00
             3.5547E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1498.35799366050        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:      771             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4345E-01  3.7345E-01  2.0424E-01  1.1496E+00  2.5648E-01  9.8122E-01  1.6149E+00  1.0000E-02  1.4779E+00  7.0282E-01
             2.7280E+00
 PARAMETER:  4.1788E-02 -8.8497E-01 -1.4885E+00  2.3943E-01 -1.2607E+00  8.1038E-02  5.7930E-01 -8.5253E+00  4.9065E-01 -2.5265E-01
             1.1036E+00
 GRADIENT:   3.9862E+01  2.1562E-01  1.9082E+01  4.6481E+01  1.3869E+02  6.0945E+00  1.6306E+00  0.0000E+00  2.3697E+01 -3.1830E+00
             1.2620E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1498.70685897819        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      947
 NPARAMETR:  9.4432E-01  3.9327E-01  2.0302E-01  1.1496E+00  2.5642E-01  9.7559E-01  1.4932E+00  1.0000E-02  1.4779E+00  7.3995E-01
             2.7280E+00
 PARAMETER:  4.2711E-02 -8.3325E-01 -1.4945E+00  2.3941E-01 -1.2609E+00  7.5284E-02  5.0092E-01 -8.5253E+00  4.9062E-01 -2.0117E-01
             1.1036E+00
 GRADIENT:  -1.0471E-01 -1.6465E-01 -1.1552E+00  2.3374E+01  1.9213E+01 -6.1517E-01  6.7335E-02  0.0000E+00  3.2210E+00  9.9737E-02
             3.0348E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1499.44848891587        NO. OF FUNC. EVALS.: 134
 CUMULATIVE NO. OF FUNC. EVALS.:     1081
 NPARAMETR:  9.3303E-01  3.7400E-01  1.9238E-01  1.0971E+00  2.4263E-01  9.6846E-01  1.4700E+00  1.0000E-02  1.4193E+00  7.4665E-01
             2.6927E+00
 PARAMETER:  3.0684E-02 -8.8349E-01 -1.5483E+00  1.9268E-01 -1.3162E+00  6.7950E-02  4.8525E-01 -8.5253E+00  4.5016E-01 -1.9216E-01
             1.0905E+00
 GRADIENT:   1.3259E+01  1.4438E+01  4.7410E+01  1.2981E+01  1.0394E+02  1.2815E+00  1.3537E+00  0.0000E+00  1.8926E+00 -1.7781E+00
             6.4546E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1500.26472780783        NO. OF FUNC. EVALS.: 160
 CUMULATIVE NO. OF FUNC. EVALS.:     1241            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4435E-01  3.7075E-01  1.8874E-01  1.1043E+00  2.4030E-01  9.7740E-01  1.4643E+00  1.0000E-02  1.4905E+00  7.6407E-01
             2.7011E+00
 PARAMETER:  4.2744E-02 -8.9222E-01 -1.5674E+00  1.9919E-01 -1.3259E+00  7.7136E-02  4.8139E-01 -8.5253E+00  4.9912E-01 -1.6910E-01
             1.0937E+00
 GRADIENT:   4.0936E+01  1.2528E+01  4.4717E+01  2.1878E+01  1.1563E+02  4.5872E+00  3.1961E+00  0.0000E+00  1.7209E+01  1.0248E+00
             1.1750E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1500.67728878098        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1376             RESET HESSIAN, TYPE I
 NPARAMETR:  9.4438E-01  3.6310E-01  1.8204E-01  1.1005E+00  2.3677E-01  9.7740E-01  1.3498E+00  1.0000E-02  1.5158E+00  7.7227E-01
             2.7061E+00
 PARAMETER:  4.2772E-02 -9.1307E-01 -1.6035E+00  1.9574E-01 -1.3407E+00  7.7142E-02  3.9996E-01 -8.5253E+00  5.1597E-01 -1.5843E-01
             1.0955E+00
 GRADIENT:   4.0809E+01  3.2666E+00  3.3169E+01  2.3737E+01  1.4471E+02  4.3548E+00  1.0811E+00  0.0000E+00  1.8811E+01 -1.5704E+00
             1.5863E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1501.14298450895        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1552
 NPARAMETR:  9.4466E-01  3.9445E-01  1.7531E-01  1.0776E+00  2.3676E-01  9.7821E-01  9.4540E-01  1.0000E-02  1.5157E+00  8.4895E-01
             2.7065E+00
 PARAMETER:  4.3065E-02 -8.3026E-01 -1.6412E+00  1.7471E-01 -1.3407E+00  7.7968E-02  4.3853E-02 -8.5253E+00  5.1588E-01 -6.3758E-02
             1.0956E+00
 GRADIENT:  -4.2961E-01  5.1445E+00  1.0715E+00 -2.4203E+00  2.2042E+01  1.4959E-01  1.0779E+00  0.0000E+00 -7.0597E+00  1.5036E+00
             6.3822E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1501.55977372793        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1728
 NPARAMETR:  9.4500E-01  3.9157E-01  1.7278E-01  1.0707E+00  2.3638E-01  9.7777E-01  2.3598E-01  1.0000E-02  1.5155E+00  9.0205E-01
             2.7071E+00
 PARAMETER:  4.3432E-02 -8.3759E-01 -1.6557E+00  1.6829E-01 -1.3423E+00  7.7519E-02 -1.3440E+00 -8.5253E+00  5.1577E-01 -3.0828E-03
             1.0959E+00
 GRADIENT:  -6.9810E-01  4.1379E-01  2.9651E+00 -9.9580E-01  3.3746E+01  5.7187E-01  4.1143E-02  0.0000E+00 -1.1140E+01  4.1477E-01
            -3.2540E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1504.33479721818        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1906
 NPARAMETR:  9.4216E-01  3.3582E-01  1.4121E-01  1.0343E+00  2.0135E-01  9.8446E-01  1.0000E-02  1.0000E-02  1.5069E+00  8.8426E-01
             2.7750E+00
 PARAMETER:  4.0423E-02 -9.9118E-01 -1.8575E+00  1.3372E-01 -1.5027E+00  8.4333E-02 -9.6973E+01 -8.5253E+00  5.1004E-01 -2.3006E-02
             1.1207E+00
 GRADIENT:  -9.4518E+00  1.4553E+00 -3.2895E+00  9.3878E+00  1.8165E+00  2.7430E+00  0.0000E+00  0.0000E+00 -2.8911E+01 -5.6706E+00
             2.9004E+01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1505.37266618874        NO. OF FUNC. EVALS.: 157
 CUMULATIVE NO. OF FUNC. EVALS.:     2063
 NPARAMETR:  9.4666E-01  3.3769E-01  1.4325E-01  1.0258E+00  2.0305E-01  9.7587E-01  1.0000E-02  1.0000E-02  1.5486E+00  9.0421E-01
             2.7428E+00
 PARAMETER:  4.5185E-02 -9.8562E-01 -1.8432E+00  1.2551E-01 -1.4943E+00  7.5571E-02 -9.1387E+01 -8.5253E+00  5.3737E-01 -6.9056E-04
             1.1090E+00
 GRADIENT:   3.9695E+01  9.7750E+00  4.4732E+01  1.0160E+01  1.5171E+02  3.8524E+00  0.0000E+00  0.0000E+00 -7.4143E-01  6.3768E-01
             3.0641E+01

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1505.47586408834        NO. OF FUNC. EVALS.: 161
 CUMULATIVE NO. OF FUNC. EVALS.:     2224
 NPARAMETR:  9.4658E-01  3.3419E-01  1.4560E-01  1.0239E+00  2.0616E-01  9.7656E-01  1.0000E-02  1.0000E-02  1.5489E+00  9.0506E-01
             2.7739E+00
 PARAMETER:  4.4097E-02 -9.8619E-01 -1.8453E+00  1.2485E-01 -1.4939E+00  7.5285E-02 -9.1387E+01 -8.5253E+00  5.4298E-01 -7.5673E-04
             1.1092E+00
 GRADIENT:  -4.2169E+04  4.2727E+03 -2.2687E+03  3.3767E+04 -2.8136E+03 -2.1084E+04  0.0000E+00  0.0000E+00  7.7429E+03 -4.2170E+04
            -3.8635E+03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2224
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.2209E-03 -6.9776E-05  2.5451E-04 -1.2094E-02 -5.4433E-04
 SE:             2.8963E-02  1.2159E-04  1.7154E-04  2.7758E-02  2.6860E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3888E-01  5.6608E-01  1.3790E-01  6.6305E-01  9.8383E-01

 ETASHRINKSD(%)  2.9708E+00  9.9593E+01  9.9425E+01  7.0066E+00  1.0015E+01
 ETASHRINKVR(%)  5.8534E+00  9.9998E+01  9.9997E+01  1.3522E+01  1.9028E+01
 EBVSHRINKSD(%)  2.7291E+00  9.9574E+01  9.9422E+01  8.2527E+00  1.0727E+01
 EBVSHRINKVR(%)  5.3837E+00  9.9998E+01  9.9997E+01  1.5824E+01  2.0303E+01
 RELATIVEINF(%)  9.4359E+01  3.8362E-04  3.6613E-04  4.6129E+01  5.1161E+00
 EPSSHRINKSD(%)  2.9609E+01
 EPSSHRINKVR(%)  5.0452E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1505.4758640883433     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -586.53733088367062     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    39.23
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.66
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1505.476       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.46E-01  3.37E-01  1.43E-01  1.03E+00  2.03E-01  9.76E-01  1.00E-02  1.00E-02  1.56E+00  9.04E-01  2.74E+00
 


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
+        1.18E+07
 
 TH 2
+        6.97E+02  9.53E+05
 
 TH 3
+       -1.06E+03  3.19E+03  1.51E+06
 
 TH 4
+        1.85E+03 -2.44E+03  1.40E+03  6.43E+06
 
 TH 5
+       -6.44E+02 -3.66E+03 -1.34E+06  1.47E+03  1.17E+06
 
 TH 6
+       -3.39E+03  9.43E+02 -1.17E+03  2.50E+03 -1.07E+03  1.11E+07
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        3.06E+02 -4.15E+02  5.83E+03 -6.22E+02  4.12E+05  3.93E+02  0.00E+00  0.00E+00  1.47E+05
 
 TH10
+       -2.63E+03  7.29E+02 -8.61E+02  1.96E+03 -8.64E+02 -2.55E+03  0.00E+00  0.00E+00  3.11E+02  1.29E+07
 
 TH11
+       -1.04E+02  1.13E+02 -1.69E+03  1.72E+02 -1.12E+05 -1.09E+02  0.00E+00  0.00E+00  9.53E+01 -7.91E+01  1.19E+04
 
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
 #CPUT: Total CPU Time in Seconds,       46.055
Stop Time:
Thu Sep 30 00:19:32 CDT 2021
