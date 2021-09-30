Wed Sep 29 11:31:18 CDT 2021
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
$DATA ../../../../data/spa/B/dat77.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m77.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1692.08863680293        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.9680E+02  3.9289E+01  1.3444E+01  7.2925E+01 -4.7773E+01  4.0926E+01  6.1898E+00 -1.1457E+00  4.6296E+01  1.5507E+01
            -1.4145E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1698.80897355146        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      194
 NPARAMETR:  1.0109E+00  1.1035E+00  1.0498E+00  9.5413E-01  1.1091E+00  1.0626E+00  9.9904E-01  1.0193E+00  7.4296E-01  9.1099E-01
             1.0842E+00
 PARAMETER:  1.1083E-01  1.9852E-01  1.4860E-01  5.3049E-02  2.0351E-01  1.6071E-01  9.9044E-02  1.1910E-01 -1.9711E-01  6.7788E-03
             1.8080E-01
 GRADIENT:  -2.5359E+01  2.3240E+01  8.6923E+00  2.5011E+01  6.8214E+00  1.2486E+01 -1.0910E+01 -1.0100E+01 -1.0068E+01 -4.1338E+00
             1.0905E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1700.03576387163        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      371
 NPARAMETR:  1.0115E+00  1.0219E+00  1.0353E+00  9.9265E-01  1.0512E+00  1.0391E+00  1.2359E+00  1.2995E+00  6.4838E-01  7.9132E-01
             1.0630E+00
 PARAMETER:  1.1145E-01  1.2167E-01  1.3473E-01  9.2625E-02  1.4989E-01  1.3832E-01  3.1183E-01  3.6197E-01 -3.3328E-01 -1.3406E-01
             1.6113E-01
 GRADIENT:  -2.3920E+01  1.2929E+01 -3.8226E+00  7.2728E+00  4.8637E+00  4.0111E+00  5.1642E+00  3.5128E+00 -6.5006E+00 -4.4198E+00
             3.5841E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1702.41713259352        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      548
 NPARAMETR:  1.0191E+00  7.9063E-01  2.0647E+00  1.1853E+00  1.3078E+00  1.0241E+00  1.3201E+00  2.0205E+00  6.9788E-01  1.1566E+00
             1.0422E+00
 PARAMETER:  1.1894E-01 -1.3492E-01  8.2499E-01  2.6997E-01  3.6832E-01  1.2377E-01  3.7771E-01  8.0333E-01 -2.5971E-01  2.4551E-01
             1.4138E-01
 GRADIENT:  -6.4054E-01  1.8581E+01 -4.7285E+00  3.8014E+01  4.1703E+00 -6.2083E-01  2.1909E+00  1.5275E+00  1.7230E+00  2.6910E+00
            -2.0207E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1703.06452253640        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      724
 NPARAMETR:  1.0185E+00  6.4197E-01  2.6406E+00  1.2773E+00  1.3387E+00  1.0243E+00  1.3952E+00  2.3123E+00  6.7398E-01  1.1613E+00
             1.0525E+00
 PARAMETER:  1.1833E-01 -3.4322E-01  1.0710E+00  3.4473E-01  3.9170E-01  1.2406E-01  4.3304E-01  9.3823E-01 -2.9455E-01  2.4953E-01
             1.5120E-01
 GRADIENT:  -7.5073E-02  5.7735E+00  1.5040E+00  9.6778E+00 -1.9200E+00 -4.6354E-02 -1.2503E+00 -1.4369E-01 -1.0977E+00 -5.0161E-01
             1.6282E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1703.12446026910        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      899
 NPARAMETR:  1.0183E+00  5.4321E-01  3.0550E+00  1.3476E+00  1.3557E+00  1.0241E+00  1.5339E+00  2.4904E+00  6.5888E-01  1.1833E+00
             1.0526E+00
 PARAMETER:  1.1818E-01 -5.1025E-01  1.2168E+00  3.9829E-01  4.0432E-01  1.2379E-01  5.2782E-01  1.0124E+00 -3.1722E-01  2.6834E-01
             1.5122E-01
 GRADIENT:   5.2073E-01  6.1564E+00  3.3293E+00  1.1285E+01 -5.6034E+00  8.4794E-03 -7.2764E-01 -1.5239E+00 -7.6694E-01  1.0255E-01
            -1.9314E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1703.19900317083        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1074
 NPARAMETR:  1.0180E+00  3.9300E-01  3.8166E+00  1.4569E+00  1.3853E+00  1.0237E+00  1.8836E+00  2.8191E+00  6.3314E-01  1.2160E+00
             1.0521E+00
 PARAMETER:  1.1786E-01 -8.3395E-01  1.4394E+00  4.7630E-01  4.2589E-01  1.2342E-01  7.3317E-01  1.1364E+00 -3.5707E-01  2.9558E-01
             1.5077E-01
 GRADIENT:   8.9130E-01  7.1457E+00  3.7526E+00  1.9372E+01 -8.3391E+00  5.1503E-02  3.0751E-01 -2.5663E+00 -6.4792E-01  9.4005E-01
            -8.8420E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1703.27022347985        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1249
 NPARAMETR:  1.0178E+00  3.0202E-01  4.4098E+00  1.5243E+00  1.4057E+00  1.0236E+00  2.2446E+00  3.0630E+00  6.1820E-01  1.2376E+00
             1.0513E+00
 PARAMETER:  1.1769E-01 -1.0973E+00  1.5838E+00  5.2156E-01  4.4057E-01  1.2329E-01  9.0854E-01  1.2194E+00 -3.8095E-01  3.1313E-01
             1.5006E-01
 GRADIENT:   9.7619E-01  7.4287E+00  3.3171E+00  2.6737E+01 -9.7680E+00  6.9390E-02  1.3523E+00 -2.7583E+00 -5.9114E-01  1.4843E+00
            -1.5939E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1703.40840700637        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  1.0175E+00  2.2310E-01  4.9801E+00  1.5791E+00  1.4237E+00  1.0236E+00  2.7228E+00  3.2915E+00  6.0883E-01  1.2511E+00
             1.0520E+00
 PARAMETER:  1.1739E-01 -1.4001E+00  1.7055E+00  5.5687E-01  4.5328E-01  1.2332E-01  1.1017E+00  1.2913E+00 -3.9621E-01  3.2405E-01
             1.5070E-01
 GRADIENT:   7.5114E-01  5.7718E+00  2.5589E+00  2.1579E+01 -8.8527E+00  4.9054E-02  2.8892E+00 -2.1691E+00  1.9799E-02  1.6717E+00
            -1.7284E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1703.57952650919        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1599
 NPARAMETR:  1.0173E+00  1.4258E-01  5.6965E+00  1.6394E+00  1.4454E+00  1.0240E+00  3.2089E+00  3.5589E+00  6.1400E-01  1.2601E+00
             1.0554E+00
 PARAMETER:  1.1712E-01 -1.8479E+00  1.8398E+00  5.9436E-01  4.6839E-01  1.2369E-01  1.2659E+00  1.3694E+00 -3.8776E-01  3.3116E-01
             1.5393E-01
 GRADIENT:   3.4972E-02  3.8867E+00  1.0337E+00  2.3553E+01 -6.1266E+00 -9.7400E-02  2.2357E+00 -7.8983E-01  1.7752E+00  1.7094E+00
            -9.2096E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1703.97965127497        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1782             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0182E+00  5.9650E-02  6.2556E+00  1.6876E+00  1.4542E+00  1.0241E+00  4.6230E+00  3.7846E+00  5.9323E-01  1.2429E+00
             1.0608E+00
 PARAMETER:  1.1807E-01 -2.7193E+00  1.9335E+00  6.2328E-01  4.7443E-01  1.2383E-01  1.6310E+00  1.4309E+00 -4.2217E-01  3.1742E-01
             1.5905E-01
 GRADIENT:   4.8797E+02  8.8172E+00  5.7517E+00  1.1769E+03  1.0494E+01  5.8230E+01  5.2107E+00  1.7706E+01  2.8414E+01  1.4702E+00
             8.2514E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1704.03243876365        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1964
 NPARAMETR:  1.0179E+00  5.8897E-02  6.2323E+00  1.6797E+00  1.4638E+00  1.0238E+00  4.5983E+00  3.7636E+00  5.9052E-01  1.2478E+00
             1.0620E+00
 PARAMETER:  1.1770E-01 -2.7320E+00  1.9298E+00  6.1863E-01  4.8103E-01  1.2355E-01  1.6257E+00  1.4254E+00 -4.2675E-01  3.2136E-01
             1.6019E-01
 GRADIENT:   2.5780E+00 -1.6398E+00  4.2277E-01 -2.0565E+01  1.3116E+00  3.0919E-01 -5.0921E+00 -1.0329E-01  1.0347E+00 -5.4547E-01
             1.3976E+00

0ITERATION NO.:   56    OBJECTIVE VALUE:  -1704.03243876365        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1992
 NPARAMETR:  1.0178E+00  6.0529E-02  6.2194E+00  1.6796E+00  1.4634E+00  1.0238E+00  4.6637E+00  3.7656E+00  5.9020E-01  1.2490E+00
             1.0612E+00
 PARAMETER:  1.1770E-01 -2.7320E+00  1.9298E+00  6.1863E-01  4.8103E-01  1.2355E-01  1.6257E+00  1.4254E+00 -4.2675E-01  3.2136E-01
             1.6019E-01
 GRADIENT:   3.6380E-01 -1.6381E+00  2.9845E-01  1.9573E+00  1.4226E+00  1.0544E-01 -4.5204E+00 -7.7790E-01  5.1207E-01 -4.3044E-01
             1.2013E+00
 NUMSIGDIG:         2.8         1.2         2.2         3.3         2.5         2.7         1.3         2.7         2.1         1.8
                    1.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1992
 NO. OF SIG. DIGITS IN FINAL EST.:  1.2

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5724E-03  6.7415E-03 -5.2425E-02 -1.9214E-02 -6.1911E-02
 SE:             2.9725E-02  8.5041E-03  1.7839E-02  2.7464E-02  2.0283E-02
 N:                     100         100         100         100         100

 P VAL.:         9.3104E-01  4.2793E-01  3.2944E-03  4.8416E-01  2.2709E-03

 ETASHRINKSD(%)  4.1603E-01  7.1510E+01  4.0238E+01  7.9926E+00  3.2049E+01
 ETASHRINKVR(%)  8.3034E-01  9.1883E+01  6.4285E+01  1.5346E+01  5.3826E+01
 EBVSHRINKSD(%)  4.8027E-01  8.0135E+01  5.3039E+01  6.3832E+00  2.6450E+01
 EBVSHRINKVR(%)  9.5823E-01  9.6054E+01  7.7947E+01  1.2359E+01  4.5904E+01
 RELATIVEINF(%)  9.8608E+01  9.5575E-01  1.2610E+01  2.1324E+01  3.0934E+01
 EPSSHRINKSD(%)  4.4055E+01
 EPSSHRINKVR(%)  6.8702E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1704.0324387636510     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -968.88161219991287     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    29.11
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.01
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1704.032       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  5.89E-02  6.23E+00  1.68E+00  1.46E+00  1.02E+00  4.60E+00  3.76E+00  5.91E-01  1.25E+00  1.06E+00
 


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
+        1.01E+03
 
 TH 2
+       -2.25E+01  4.26E+04
 
 TH 3
+       -7.55E-01 -2.03E+00  6.37E-01
 
 TH 4
+       -2.10E+01 -5.83E+03 -3.19E+00  9.05E+04
 
 TH 5
+       -7.20E-01 -3.44E+01 -2.29E+00 -2.33E+01  2.21E+02
 
 TH 6
+       -1.62E+00 -2.89E+00 -1.02E-01 -5.97E+00  8.80E-01  1.86E+02
 
 TH 7
+        4.13E-01  2.61E+03 -2.71E-02 -1.41E+02  1.63E+00  3.98E-01  5.70E+01
 
 TH 8
+       -5.09E-02  8.83E-01 -1.87E+00 -1.74E+04 -1.68E+01  2.48E-01 -6.53E-02  8.55E+00
 
 TH 9
+       -5.01E+00 -6.06E+02 -5.76E-01 -8.74E+00  7.03E+00 -7.49E+00 -5.78E+02  5.49E+00  4.68E+02
 
 TH10
+        1.10E+00 -2.18E+01  9.11E-01  4.49E+00 -3.32E+01  1.22E+00 -2.71E+00 -3.77E+00  8.67E+00  5.00E+01
 
 TH11
+       -7.72E+00 -2.59E+01  1.35E+00 -1.09E+01  2.35E+01  2.70E+00 -8.89E-01 -9.99E+00  1.14E+01  2.52E+01  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       37.169
Stop Time:
Wed Sep 29 11:31:57 CDT 2021
