Wed Sep 29 23:33:59 CDT 2021
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
$DATA ../../../../data/spa1/A2/dat64.csv ignore=@
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
 (E4.0,E3.0,E20.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m64.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1387.89773533076        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.5168E+02 -1.7903E+01 -2.6342E+01  3.8482E+01  1.0264E+02  2.1762E+01 -4.3102E+00  2.8253E+00  1.6291E+01 -1.3259E+01
            -1.4233E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1798.26079257367        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  1.0738E+00  1.0944E+00  1.1378E+00  1.0185E+00  9.9997E-01  1.1376E+00  9.7032E-01  9.2693E-01  8.9102E-01  8.3452E-01
             2.1628E+00
 PARAMETER:  1.7125E-01  1.9019E-01  2.2906E-01  1.1831E-01  9.9970E-02  2.2892E-01  6.9866E-02  2.4118E-02 -1.5391E-02 -8.0904E-02
             8.7142E-01
 GRADIENT:   2.7102E+02  3.6159E+01  1.0493E+01  2.9355E+01 -2.7148E+01  4.1599E+01 -4.8721E-01  1.4967E+00 -1.8314E+00  5.9576E+00
             2.2653E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1802.25891677933        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0326E+00  1.0873E+00  7.5666E-01  9.9498E-01  8.5194E-01  1.0689E+00  1.0025E+00  5.1240E-01  9.9656E-01  5.9690E-01
             2.1264E+00
 PARAMETER:  1.3204E-01  1.8371E-01 -1.7885E-01  9.4967E-02 -6.0239E-02  1.6663E-01  1.0252E-01 -5.6865E-01  9.6553E-02 -4.1601E-01
             8.5444E-01
 GRADIENT:   1.7209E+02  9.2938E+00 -1.4688E+01  3.8952E+01  2.3236E+01  1.9552E+01  1.7005E+00  1.1781E+00  1.8459E+01  9.8749E-01
             7.0610E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1804.19756010975        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      293
 NPARAMETR:  1.0077E+00  1.0215E+00  8.8227E-01  1.0350E+00  8.8151E-01  1.0562E+00  1.0849E+00  4.7105E-01  8.6777E-01  6.7677E-01
             2.1213E+00
 PARAMETER:  1.0765E-01  1.2131E-01 -2.5259E-02  1.3437E-01 -2.6121E-02  1.5470E-01  1.8146E-01 -6.5279E-01 -4.1833E-02 -2.9043E-01
             8.5204E-01
 GRADIENT:   1.4573E+00  6.3633E-01 -4.7019E-01  9.0993E-02  5.4419E-01  1.0305E-01  1.4050E-01  2.8296E-01  2.1665E-01  3.2862E-01
            -4.2267E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1804.44838605915        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      469
 NPARAMETR:  1.0035E+00  8.2392E-01  9.5588E-01  1.1639E+00  8.3532E-01  1.0533E+00  1.2556E+00  1.1877E-01  8.0968E-01  7.2421E-01
             2.1366E+00
 PARAMETER:  1.0346E-01 -9.3682E-02  5.4874E-02  2.5182E-01 -7.9939E-02  1.5196E-01  3.2765E-01 -2.0306E+00 -1.1111E-01 -2.2267E-01
             8.5921E-01
 GRADIENT:  -3.2241E+00  7.8485E+00  1.6517E+00  1.1903E+01 -5.5317E+00  2.6863E-01 -8.2542E-02  2.2666E-02 -6.3313E-01 -1.7886E-01
             2.8023E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1804.76733647122        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      647
 NPARAMETR:  1.0015E+00  5.9798E-01  1.1516E+00  1.3091E+00  8.4718E-01  1.0457E+00  1.4977E+00  1.4600E-02  7.6837E-01  8.1871E-01
             2.1224E+00
 PARAMETER:  1.0152E-01 -4.1420E-01  2.4118E-01  3.6935E-01 -6.5844E-02  1.4471E-01  5.0395E-01 -4.1267E+00 -1.6349E-01 -1.0003E-01
             8.5256E-01
 GRADIENT:   4.8626E-01  3.8878E+00  2.4157E+00  7.5097E+00 -4.6555E+00 -8.4348E-01 -1.8582E-01  6.7096E-04 -6.3462E-01 -7.6760E-02
            -1.6564E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1804.85926351832        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      825
 NPARAMETR:  9.9815E-01  4.5866E-01  1.2257E+00  1.3948E+00  8.3659E-01  1.0469E+00  1.7773E+00  1.0000E-02  7.4416E-01  8.4952E-01
             2.1247E+00
 PARAMETER:  9.8153E-02 -6.7944E-01  3.0349E-01  4.3277E-01 -7.8421E-02  1.4581E-01  6.7508E-01 -6.1182E+00 -1.9550E-01 -6.3087E-02
             8.5365E-01
 GRADIENT:  -1.7269E+00  2.3057E+00  1.5342E+00  4.8713E+00 -3.3137E+00  4.7131E-01  6.1218E-01  0.0000E+00 -2.2696E-02  2.7171E-01
             1.0431E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1804.88821134831        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1001
 NPARAMETR:  9.9733E-01  3.4381E-01  1.3241E+00  1.4699E+00  8.4177E-01  1.0437E+00  2.0202E+00  1.0000E-02  7.3188E-01  8.8843E-01
             2.1193E+00
 PARAMETER:  9.7323E-02 -9.6765E-01  3.8072E-01  4.8519E-01 -7.2250E-02  1.4277E-01  8.0319E-01 -8.2892E+00 -2.1214E-01 -1.8294E-02
             8.5108E-01
 GRADIENT:   4.9468E-01  1.8105E+00  1.4596E+00  4.2950E+00 -3.1987E+00 -3.1120E-02  8.7148E-01  0.0000E+00  2.0802E-01  2.7314E-01
            -8.9920E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1804.90067296828        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1176
 NPARAMETR:  9.9643E-01  2.7540E-01  1.3951E+00  1.5183E+00  8.4828E-01  1.0422E+00  2.1768E+00  1.0000E-02  7.2444E-01  9.1022E-01
             2.1193E+00
 PARAMETER:  9.6419E-02 -1.1895E+00  4.3298E-01  5.1756E-01 -6.4540E-02  1.4135E-01  8.7787E-01 -9.9786E+00 -2.2236E-01  5.9281E-03
             8.5108E-01
 GRADIENT:   7.0148E-01  2.0838E+00  1.4938E+00  9.9146E+00 -3.7026E+00 -2.7745E-01  4.9343E-01  0.0000E+00  2.9096E-02  2.1162E-01
            -1.2734E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1804.92890702823        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1353
 NPARAMETR:  9.9509E-01  2.0671E-01  1.4515E+00  1.5619E+00  8.5095E-01  1.0413E+00  2.3697E+00  1.0000E-02  7.1624E-01  9.2457E-01
             2.1211E+00
 PARAMETER:  9.5078E-02 -1.4764E+00  4.7260E-01  5.4592E-01 -6.1405E-02  1.4050E-01  9.6278E-01 -1.2223E+01 -2.3374E-01  2.1575E-02
             8.5193E-01
 GRADIENT:   1.4736E-01  1.0155E+00  4.0921E-01  7.7358E+00 -1.4207E+00 -2.7383E-01 -1.9437E-01  0.0000E+00 -2.3382E-01  1.9893E-02
            -5.6207E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1804.93443924494        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1528
 NPARAMETR:  9.9436E-01  1.6879E-01  1.4767E+00  1.5851E+00  8.4989E-01  1.0411E+00  2.5349E+00  1.0000E-02  7.1166E-01  9.3102E-01
             2.1215E+00
 PARAMETER:  9.4344E-02 -1.6791E+00  4.8984E-01  5.6067E-01 -6.2646E-02  1.4030E-01  1.0302E+00 -1.3840E+01 -2.4015E-01  2.8528E-02
             8.5210E-01
 GRADIENT:  -9.2695E-02  6.1108E-01  1.9491E-01  5.3810E+00 -8.3350E-01 -1.5650E-01 -2.7885E-01  0.0000E+00 -2.6354E-01 -5.0384E-02
            -2.5985E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1804.95407180971        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:     1690
 NPARAMETR:  9.9432E-01  1.6399E-01  1.4747E+00  1.5851E+00  8.4937E-01  1.0414E+00  2.6950E+00  1.0000E-02  7.1010E-01  9.3094E-01
             2.1208E+00
 PARAMETER:  9.4302E-02 -1.7079E+00  4.8847E-01  5.6067E-01 -6.3265E-02  1.4054E-01  1.0914E+00 -1.3840E+01 -2.4234E-01  2.8441E-02
             8.5181E-01
 GRADIENT:   1.8040E-01  2.4515E-01 -4.1034E-01 -1.0816E+00  1.1529E+00  1.6157E-02  4.8581E-02  0.0000E+00 -2.2328E-02 -1.1615E-01
            -1.9897E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1804.95615344069        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1871             RESET HESSIAN, TYPE I
 NPARAMETR:  9.9419E-01  1.4747E-01  1.4732E+00  1.5915E+00  8.4455E-01  1.0413E+00  2.8266E+00  1.0000E-02  7.0967E-01  9.3183E-01
             2.1208E+00
 PARAMETER:  9.4173E-02 -1.8141E+00  4.8743E-01  5.6469E-01 -6.8957E-02  1.4047E-01  1.1391E+00 -1.3840E+01 -2.4295E-01  2.9400E-02
             8.5181E-01
 GRADIENT:   9.6996E+01  4.8449E+00  2.1211E+00  2.4579E+02  2.8546E+00  1.5167E+01  1.3118E+00  0.0000E+00  6.1443E+00  1.0531E-01
             8.8626E+00

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1804.95929028449        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2047
 NPARAMETR:  9.9370E-01  1.3746E-01  1.4746E+00  1.5993E+00  8.4288E-01  1.0409E+00  2.8849E+00  1.0000E-02  7.0801E-01  9.3262E-01
             2.1204E+00
 PARAMETER:  9.3679E-02 -1.8844E+00  4.8842E-01  5.6957E-01 -7.0932E-02  1.4006E-01  1.1595E+00 -1.3840E+01 -2.4529E-01  3.0246E-02
             8.5159E-01
 GRADIENT:  -4.2619E-02 -4.1070E-02 -6.6545E-01 -4.2173E+00  1.4925E+00 -6.0256E-03 -1.0197E-01  0.0000E+00  2.8282E-01 -9.8234E-02
             1.0217E-01

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1804.96225614581        NO. OF FUNC. EVALS.: 189
 CUMULATIVE NO. OF FUNC. EVALS.:     2236
 NPARAMETR:  9.9335E-01  1.1506E-01  1.4822E+00  1.6133E+00  8.3832E-01  1.0406E+00  3.1891E+00  1.0000E-02  7.0496E-01  9.3551E-01
             2.1198E+00
 PARAMETER:  9.3327E-02 -2.0623E+00  4.9351E-01  5.7831E-01 -7.6357E-02  1.3982E-01  1.2597E+00 -1.3840E+01 -2.4962E-01  3.3341E-02
             8.5134E-01
 GRADIENT:  -1.8914E-02  1.4932E-01  2.9725E-01 -3.2831E+00 -5.9344E-01  2.6904E-02  5.3823E-02  0.0000E+00  3.6347E-01  4.7198E-02
             1.2376E-01

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1804.96225614581        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     2264
 NPARAMETR:  9.9371E-01  1.1335E-01  1.4800E+00  1.6125E+00  8.3916E-01  1.0408E+00  3.1627E+00  1.0000E-02  7.0340E-01  9.3464E-01
             2.1192E+00
 PARAMETER:  9.3327E-02 -2.0623E+00  4.9351E-01  5.7831E-01 -7.6357E-02  1.3982E-01  1.2597E+00 -1.3840E+01 -2.4962E-01  3.3341E-02
             8.5134E-01
 GRADIENT:  -3.7131E-01  8.7602E-02  3.0123E-01  1.0096E+00 -2.9989E-01 -3.0025E-02  6.0150E-02  0.0000E+00  3.3560E-01  4.7198E-02
             1.5303E-01

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     2264
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0659E-03  1.8093E-03  2.2828E-05 -9.0258E-03 -2.4078E-02
 SE:             2.9479E-02  6.4341E-03  1.1844E-04  2.7486E-02  2.0447E-02
 N:                     100         100         100         100         100

 P VAL.:         9.7116E-01  7.7855E-01  8.4717E-01  7.4262E-01  2.3896E-01

 ETASHRINKSD(%)  1.2412E+00  7.8445E+01  9.9603E+01  7.9189E+00  3.1501E+01
 ETASHRINKVR(%)  2.4669E+00  9.5354E+01  9.9998E+01  1.5211E+01  5.3079E+01
 EBVSHRINKSD(%)  1.2957E+00  7.9892E+01  9.9579E+01  7.5956E+00  3.1781E+01
 EBVSHRINKVR(%)  2.5747E+00  9.5957E+01  9.9998E+01  1.4614E+01  5.3462E+01
 RELATIVEINF(%)  9.3188E+01  7.3383E-02  1.1712E-04  1.7373E+00  3.4677E+00
 EPSSHRINKSD(%)  2.6145E+01
 EPSSHRINKVR(%)  4.5455E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1804.9622561458052     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -886.02372294113252     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    35.64
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.78
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1804.962       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.93E-01  1.15E-01  1.48E+00  1.61E+00  8.38E-01  1.04E+00  3.19E+00  1.00E-02  7.05E-01  9.36E-01  2.12E+00
 


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
+       -3.38E+01  4.37E+02
 
 TH 3
+       -6.26E+00  4.14E+01  9.00E+01
 
 TH 4
+       -2.02E+01  4.32E+02 -4.16E+01  7.37E+02
 
 TH 5
+        1.50E+01 -2.37E+02 -2.54E+02 -5.96E+01  8.30E+02
 
 TH 6
+        3.04E+00 -5.89E+00  5.36E-01 -7.40E+00 -1.91E+00  1.73E+02
 
 TH 7
+        7.98E-02  9.43E+00  4.57E-02 -9.18E-01  2.38E-01  4.14E-02  8.19E-01
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        7.13E+00 -2.30E+01  1.06E+01 -5.31E+00 -5.17E+00 -7.83E-01  1.03E+00  0.00E+00  2.97E+02
 
 TH10
+       -6.83E-01  1.08E+01  2.46E+00  3.89E+00 -4.12E+01  1.28E+00 -3.38E-01  0.00E+00  5.26E-01  5.66E+01
 
 TH11
+       -1.24E+01 -9.86E+00 -3.90E+00 -1.52E+01 -6.07E+00  2.01E+00  5.18E-02  0.00E+00  1.13E+01  2.10E+01  1.06E+02
 
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
 #CPUT: Total CPU Time in Seconds,       43.315
Stop Time:
Wed Sep 29 23:34:47 CDT 2021
