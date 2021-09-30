Wed Sep 29 08:12:20 CDT 2021
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
$DATA ../../../../data/int/D/dat20.csv ignore=@
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
 (2E4.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m20.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -2409.54312308821        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4810E+02  9.2685E-01 -3.6636E+01 -2.3639E+02  4.2299E+02 -6.4973E+02 -1.8556E+02 -3.0516E+02 -8.0407E+02 -3.2481E+02
            -8.5171E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -3354.33598517672        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:      150
 NPARAMETR:  8.7342E-01  1.1020E+00  9.2056E-01  1.1713E+00  6.8174E-01  1.9865E+00  1.1674E+00  1.9177E+00  2.3991E+00  1.8248E+00
             1.5143E+00
 PARAMETER: -3.5338E-02  1.9709E-01  1.7222E-02  2.5812E-01 -2.8311E-01  7.8637E-01  2.5479E-01  7.5115E-01  9.7509E-01  7.0146E-01
             5.1494E-01
 GRADIENT:  -1.4743E+02  1.1245E+02  5.1155E+00  1.5884E+00 -1.2060E+02 -6.6784E+01  7.4857E-01 -9.5850E+01 -7.5832E+01 -3.6647E+01
             5.9052E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -3366.96002344874        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  8.6281E-01  1.1671E+00  1.0392E+00  1.1057E+00  7.1219E-01  1.9822E+00  1.0376E+00  2.2564E+00  2.3834E+00  2.0999E+00
             1.5091E+00
 PARAMETER: -4.7566E-02  2.5455E-01  1.3843E-01  2.0048E-01 -2.3941E-01  7.8420E-01  1.3695E-01  9.1375E-01  9.6852E-01  8.4190E-01
             5.1151E-01
 GRADIENT:  -1.5512E+02  1.2010E+02  1.8790E+01 -3.0556E-02 -1.1090E+02 -7.0174E+01 -8.9515E+00 -7.3704E+01 -7.8850E+01 -7.3634E+00
             5.8862E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -3399.01179644527        NO. OF FUNC. EVALS.: 187
 CUMULATIVE NO. OF FUNC. EVALS.:      515
 NPARAMETR:  8.5096E-01  1.1317E+00  9.9029E-01  1.0444E+00  7.0822E-01  2.0649E+00  8.0973E-01  2.2694E+00  2.6182E+00  2.1654E+00
             1.4513E+00
 PARAMETER: -6.1389E-02  2.2368E-01  9.0240E-02  1.4345E-01 -2.4500E-01  8.2509E-01 -1.1105E-01  9.1951E-01  1.0625E+00  8.7262E-01
             4.7249E-01
 GRADIENT:  -1.4981E+02  1.1342E+02  1.6130E+01 -1.1690E+01 -1.0057E+02 -4.4447E+01 -9.0700E+00 -6.3378E+01 -5.9720E+01 -7.7941E+00
             5.4459E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -3409.98588353229        NO. OF FUNC. EVALS.: 114
 CUMULATIVE NO. OF FUNC. EVALS.:      629
 NPARAMETR:  8.5085E-01  1.1320E+00  9.9015E-01  1.0442E+00  7.0799E-01  2.1813E+00  8.0985E-01  2.2724E+00  2.8584E+00  2.1659E+00
             1.0981E+00
 PARAMETER: -6.1525E-02  2.2399E-01  9.0103E-02  1.4326E-01 -2.4533E-01  8.7990E-01 -1.1090E-01  9.2084E-01  1.1503E+00  8.7282E-01
             1.9361E-01
 GRADIENT:   2.1083E+02  3.2120E+02  5.4004E+01  6.2880E+01  9.4862E+01  1.5908E+03 -1.6589E+01 -1.0756E+02  3.8389E+02  4.7141E+02
             1.4300E+02

0ITERATION NO.:   25    OBJECTIVE VALUE:  -3461.02559212322        NO. OF FUNC. EVALS.: 156
 CUMULATIVE NO. OF FUNC. EVALS.:      785
 NPARAMETR:  8.5084E-01  1.1320E+00  9.9016E-01  1.0355E+00  7.0798E-01  2.1752E+00  8.0986E-01  2.3853E+00  2.8594E+00  2.1656E+00
             1.0937E+00
 PARAMETER: -6.1530E-02  2.2400E-01  9.0111E-02  1.3493E-01 -2.4534E-01  8.7711E-01 -1.1090E-01  9.6931E-01  1.1506E+00  8.7271E-01
             1.8953E-01
 GRADIENT:  -1.3303E+02  1.0329E+02  2.8132E+01 -2.5121E+01 -1.0569E+02  2.0529E+01 -1.4404E+01 -1.2479E+02 -4.9223E+01 -1.2925E+01
             1.5520E+02

0ITERATION NO.:   30    OBJECTIVE VALUE:  -3500.02368745685        NO. OF FUNC. EVALS.: 146
 CUMULATIVE NO. OF FUNC. EVALS.:      931             RESET HESSIAN, TYPE I
 NPARAMETR:  9.3311E-01  1.0950E+00  9.9024E-01  1.0687E+00  7.6853E-01  2.4377E+00  8.8411E-01  2.7700E+00  3.0329E+00  2.1754E+00
             1.0936E+00
 PARAMETER:  3.0764E-02  1.9079E-01  9.0191E-02  1.6648E-01 -1.6327E-01  9.9107E-01 -2.3175E-02  1.1188E+00  1.2095E+00  8.7721E-01
             1.8945E-01
 GRADIENT:   2.8196E+02  1.8185E+02  2.0267E+01  9.8480E+01  5.8671E+01  1.7951E+03  4.0367E+00 -7.9180E+00  5.0548E+02  4.7224E+02
             1.7328E+02

0ITERATION NO.:   35    OBJECTIVE VALUE:  -3508.07848776101        NO. OF FUNC. EVALS.:  94
 CUMULATIVE NO. OF FUNC. EVALS.:     1025
 NPARAMETR:  9.5575E-01  1.0853E+00  9.9027E-01  1.0712E+00  7.9885E-01  2.4212E+00  9.0702E-01  2.9150E+00  3.0003E+00  2.1232E+00
             1.0936E+00
 PARAMETER:  5.4739E-02  1.8182E-01  9.0225E-02  1.6879E-01 -1.2458E-01  9.8425E-01  2.4116E-03  1.1699E+00  1.1987E+00  8.5292E-01
             1.8951E-01
 GRADIENT:  -6.3252E+01  8.1703E+00  8.2333E+00 -1.6276E+01 -5.4421E+01  8.7240E+01 -1.8504E-01 -6.5012E+01 -3.9571E-01 -3.6313E+00
             1.7403E+02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -3509.58787827109        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1209
 NPARAMETR:  9.6095E-01  1.0879E+00  9.9028E-01  1.0721E+00  8.0655E-01  2.4164E+00  9.1328E-01  2.9437E+00  2.9909E+00  2.1178E+00
             1.0937E+00
 PARAMETER:  6.0167E-02  1.8422E-01  9.0231E-02  1.6964E-01 -1.1499E-01  9.8228E-01  9.2844E-03  1.1797E+00  1.1956E+00  8.5036E-01
             1.8952E-01
 GRADIENT:  -6.1453E+01  7.0908E+00  7.1193E+00 -1.5733E+01 -5.1728E+01  8.6172E+01 -2.4816E-01 -6.3144E+01 -2.3349E-01 -3.7679E+00
             1.7457E+02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -3521.56374059346        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1374
 NPARAMETR:  1.0306E+00  1.0755E+00  9.9045E-01  1.0720E+00  8.0649E-01  2.4120E+00  8.9783E-01  3.4621E+00  2.9848E+00  2.1208E+00
             1.0940E+00
 PARAMETER:  1.3011E-01  1.7281E-01  9.0400E-02  1.6954E-01 -1.1506E-01  9.8047E-01 -7.7762E-03  1.3419E+00  1.1935E+00  8.5179E-01
             1.8984E-01
 GRADIENT:  -3.4182E+01 -2.7516E+00  8.6394E-01 -1.3559E+01 -5.5499E+01  8.1204E+01  6.4315E+00 -2.3614E+01  1.1399E+01 -9.0009E+00
             1.8091E+02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -3527.00984104238        NO. OF FUNC. EVALS.: 148
 CUMULATIVE NO. OF FUNC. EVALS.:     1522
 NPARAMETR:  1.1054E+00  1.0758E+00  9.8587E-01  1.0928E+00  8.5710E-01  2.2075E+00  8.5213E-01  3.4620E+00  2.7935E+00  2.0514E+00
             1.0942E+00
 PARAMETER:  2.0021E-01  1.7304E-01  8.5766E-02  1.8874E-01 -5.4198E-02  8.9186E-01 -6.0019E-02  1.3418E+00  1.1273E+00  8.1853E-01
             1.9001E-01
 GRADIENT:  -6.3629E+00 -4.0452E+00 -1.0051E+00 -8.9160E+00 -2.8936E+01  3.1760E+01 -1.2095E+00 -2.9559E+01 -8.5820E+00 -1.6777E+01
             1.7967E+02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -3533.10736697841        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1660
 NPARAMETR:  1.1011E+00  1.0788E+00  9.9931E-01  1.1146E+00  8.9487E-01  2.1800E+00  8.6518E-01  3.5610E+00  2.8500E+00  2.1146E+00
             1.0630E+00
 PARAMETER:  1.9634E-01  1.7581E-01  9.9307E-02  2.0848E-01 -1.1082E-02  8.7931E-01 -4.4813E-02  1.3700E+00  1.1473E+00  8.4887E-01
             1.6112E-01
 GRADIENT:   7.3114E+02  1.3591E+02  5.9794E+00  1.5373E+02  7.1864E+01  1.3806E+03  4.0729E+00  6.1610E+01  5.2769E+02  5.9877E+02
             1.3352E+02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -3537.61971885600        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:     1730
 NPARAMETR:  1.1010E+00  1.0862E+00  1.0088E+00  1.1150E+00  9.1030E-01  2.1799E+00  9.0669E-01  3.6899E+00  2.8298E+00  2.0886E+00
             1.0263E+00
 PARAMETER:  1.9626E-01  1.8268E-01  1.0873E-01  2.0881E-01  6.0163E-03  8.7927E-01  2.0467E-03  1.4056E+00  1.1402E+00  8.3649E-01
             1.2593E-01
 GRADIENT:   7.8562E+02  1.5400E+02  7.0182E+00  1.6606E+02  8.1736E+01  1.4821E+03  4.0834E+00  7.6686E+01  5.6623E+02  6.4042E+02
             7.0748E+01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -3538.01237490656        NO. OF FUNC. EVALS.: 130
 CUMULATIVE NO. OF FUNC. EVALS.:     1860
 NPARAMETR:  1.1010E+00  1.0865E+00  1.0090E+00  1.1149E+00  9.1253E-01  2.1792E+00  9.1267E-01  3.7075E+00  2.8259E+00  2.0859E+00
             1.0215E+00
 PARAMETER:  1.9622E-01  1.8295E-01  1.0899E-01  2.0874E-01  8.4650E-03  8.7896E-01  8.6168E-03  1.4104E+00  1.1388E+00  8.3519E-01
             1.2123E-01
 GRADIENT:  -1.2372E+01 -9.7801E+03 -3.3715E+03 -8.0490E+00  7.3424E+03 -1.5901E+01 -2.2924E+00  4.9719E+02  9.7502E-01  3.3976E+03
             5.8898E+01

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1860
 NO. OF SIG. DIGITS IN FINAL EST.:  2.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.3467E-02 -6.6177E-02 -3.1095E-02  1.7825E-02 -3.6954E-02
 SE:             3.0955E-02  1.9888E-02  2.6591E-02  2.8063E-02  2.5846E-02
 N:                     100         100         100         100         100

 P VAL.:         6.6352E-01  8.7659E-04  2.4224E-01  5.2530E-01  1.5278E-01

 ETASHRINKSD(%)  1.0000E-10  3.3372E+01  1.0918E+01  5.9868E+00  1.3413E+01
 ETASHRINKVR(%)  1.0000E-10  5.5607E+01  2.0643E+01  1.1615E+01  2.5027E+01
 EBVSHRINKSD(%)  5.9400E-02  4.3652E+01  1.3496E+01  3.0229E+00  8.1257E+00
 EBVSHRINKVR(%)  1.1876E-01  6.8249E+01  2.5170E+01  5.9544E+00  1.5591E+01
 RELATIVEINF(%)  9.9881E+01  2.3703E+01  7.2266E+01  8.7077E+01  6.4740E+01
 EPSSHRINKSD(%)  2.6567E+01
 EPSSHRINKVR(%)  4.6076E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -3538.0123749065551     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1883.9230151381444     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    60.21
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    17.03
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -3538.012       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.10E+00  1.09E+00  1.01E+00  1.11E+00  9.13E-01  2.18E+00  9.13E-01  3.71E+00  2.83E+00  2.09E+00  1.02E+00
 


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
+        3.94E+05
 
 TH 2
+       -2.32E-01  9.30E+05
 
 TH 3
+       -7.74E+05 -1.08E+03  3.04E+06
 
 TH 4
+        2.64E-01  1.14E+06  1.23E+03  1.99E+06
 
 TH 5
+       -1.28E-01 -1.24E+02  7.08E+06  8.64E+05  1.29E+07
 
 TH 6
+       -2.16E+05 -4.82E+04  1.74E+05 -2.68E-01 -3.01E+05  1.01E+04
 
 TH 7
+       -2.67E+06  1.01E+06 -1.17E+01 -7.60E+00 -1.07E+07 -1.05E+05  4.41E+06
 
 TH 8
+        6.25E+04  8.67E+02 -3.27E+04  1.55E+04  2.48E+05 -1.88E+03  3.95E+04  1.36E+03
 
 TH 9
+        1.02E+05 -7.90E+01  5.69E+01  1.57E+00  1.79E+05 -1.84E-01  6.25E+04 -1.12E+03  8.65E+03
 
 TH10
+        1.88E+05 -1.51E+05  1.27E+02 -1.75E+05 -2.14E+03 -5.50E+03 -2.75E+02  9.83E+03 -9.42E+03  6.09E+03
 
 TH11
+       -3.34E+06 -7.46E+05 -7.61E+00  6.37E+05 -4.67E+06  5.45E-02  2.74E+01 -2.91E+04 -1.78E+05 -8.60E+04  9.05E+02
 
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
 #CPUT: Total CPU Time in Seconds,       77.353
Stop Time:
Wed Sep 29 08:13:39 CDT 2021
