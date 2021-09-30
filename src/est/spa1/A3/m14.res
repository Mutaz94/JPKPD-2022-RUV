Thu Sep 30 00:00:50 CDT 2021
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
$DATA ../../../../data/spa1/A3/dat14.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m14.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -107.166248565747        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1086E+02 -4.9472E+01  1.4080E+02 -6.6101E+01  2.0495E+02  3.9657E+01 -1.1155E+01 -2.5435E+02 -8.2705E+00 -7.5279E+01
            -3.5920E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1478.33924877611        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.0574E+00  1.1404E+00  9.0169E-01  1.1421E+00  9.1171E-01  8.4644E-01  8.7490E-01  1.1249E+00  7.4496E-01  9.3636E-01
             5.3298E+00
 PARAMETER:  1.5581E-01  2.3135E-01 -3.4862E-03  2.3290E-01  7.5629E-03 -6.6712E-02 -3.3647E-02  2.1766E-01 -1.9442E-01  3.4247E-02
             1.7733E+00
 GRADIENT:   4.5202E+01  2.7413E+01 -9.7393E+00  4.9838E+01 -2.1537E+01 -1.6542E+01  1.4861E+01  9.0301E+00  2.1905E+01  2.0785E+01
             3.6340E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1511.05188480125        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      156
 NPARAMETR:  1.0427E+00  9.1549E-01  2.9597E-01  1.1469E+00  4.2257E-01  9.0682E-01  4.2801E-01  1.2674E+00  7.4326E-01  4.9791E-01
             4.7227E+00
 PARAMETER:  1.4182E-01  1.1701E-02 -1.1175E+00  2.3709E-01 -7.6141E-01  2.1928E-03 -7.4862E-01  3.3698E-01 -1.9671E-01 -5.9734E-01
             1.6524E+00
 GRADIENT:   5.6263E+00  9.2395E+01  1.6917E+01  1.1042E+02 -8.5233E+01 -1.3015E+01  9.1145E-01  1.9545E+01  4.9282E+00  1.0918E+01
             2.9078E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1599.30371125224        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      230
 NPARAMETR:  9.6649E-01  5.8928E-01  2.1044E-01  1.1082E+00  2.9382E-01  9.4033E-01  9.0829E-02  3.9084E-01  9.6464E-01  6.9165E-01
             2.8475E+00
 PARAMETER:  6.5919E-02 -4.2886E-01 -1.4586E+00  2.0271E-01 -1.1248E+00  3.8472E-02 -2.2988E+00 -8.3947E-01  6.3999E-02 -2.6867E-01
             1.1464E+00
 GRADIENT:  -6.6241E+01  7.1994E+01  2.1647E+01  4.6750E+01  2.7755E+01 -1.3661E+01  5.9317E-02 -2.6812E+00 -2.2199E+01  4.6445E+00
             1.1508E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1605.97380188323        NO. OF FUNC. EVALS.: 136
 CUMULATIVE NO. OF FUNC. EVALS.:      366
 NPARAMETR:  9.9851E-01  5.1800E-01  2.0123E-01  1.1252E+00  2.7137E-01  9.8526E-01  7.2961E-02  3.7901E-01  1.0995E+00  6.8360E-01
             2.6587E+00
 PARAMETER:  9.8509E-02 -5.5778E-01 -1.5033E+00  2.1799E-01 -1.2043E+00  8.5145E-02 -2.5178E+00 -8.7020E-01  1.9483E-01 -2.8038E-01
             1.0779E+00
 GRADIENT:  -1.4425E+01  4.0432E+01  2.5386E+01  3.9743E+01 -2.9521E+01  4.6387E+00  4.8304E-02 -4.3460E+00  1.4355E+00  1.0888E+00
            -4.4191E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1622.89456277513        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      544
 NPARAMETR:  1.0081E+00  4.0978E-01  1.6281E-01  1.0628E+00  2.2657E-01  9.8857E-01  1.1907E-01  1.1924E+00  1.1734E+00  5.3206E-01
             2.5253E+00
 PARAMETER:  1.0809E-01 -7.9214E-01 -1.7151E+00  1.6095E-01 -1.3847E+00  8.8509E-02 -2.0281E+00  2.7600E-01  2.5993E-01 -5.3099E-01
             1.0264E+00
 GRADIENT:   1.2584E+01 -4.3617E-01  7.8590E+00 -1.1862E+01  9.4952E+00  6.9790E+00  2.2517E-01 -6.7782E+00  3.8851E+00  6.9995E+00
            -6.5596E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1624.62877008673        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      719
 NPARAMETR:  1.0034E+00  3.9029E-01  1.5718E-01  1.0695E+00  2.2054E-01  9.6990E-01  1.2281E-01  1.3906E+00  1.1896E+00  4.2038E-01
             2.5263E+00
 PARAMETER:  1.0335E-01 -8.4088E-01 -1.7503E+00  1.6715E-01 -1.4117E+00  6.9437E-02 -1.9971E+00  4.2971E-01  2.7362E-01 -7.6660E-01
             1.0267E+00
 GRADIENT:   3.6056E+00 -3.6942E+00 -2.8141E+00 -2.1279E+00  2.1423E+00  2.1306E-01  2.2631E-01  2.8759E-01 -9.1967E-01  5.4178E-01
             4.9291E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1624.91909510757        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      896            RESET HESSIAN, TYPE II
 NPARAMETR:  1.0018E+00  4.1290E-01  1.7229E-01  1.0938E+00  2.3530E-01  9.6972E-01  2.6373E-02  1.4645E+00  1.1844E+00  3.1588E-01
             2.5202E+00
 PARAMETER:  1.0181E-01 -7.8454E-01 -1.6586E+00  1.8963E-01 -1.3469E+00  6.9255E-02 -3.5354E+00  4.8149E-01  2.6921E-01 -1.0524E+00
             1.0243E+00
 GRADIENT:   4.4802E+01  1.0722E+01  3.0942E+01  1.6033E+01  1.2512E+02  3.4525E+00  1.2910E-02  2.0135E+00  3.1198E+00  8.8552E-01
             1.0294E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1624.92272587196        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      969
 NPARAMETR:  1.0022E+00  4.1251E-01  1.7154E-01  1.0933E+00  2.3479E-01  9.6995E-01  1.0000E-02  1.4570E+00  1.1823E+00  3.2482E-01
             2.5222E+00
 PARAMETER:  1.0218E-01 -7.8549E-01 -1.6629E+00  1.8923E-01 -1.3491E+00  6.9488E-02 -8.9899E+00  4.7640E-01  2.6750E-01 -1.0245E+00
             1.0251E+00
 GRADIENT:   4.5573E+01  1.0447E+01  3.0263E+01  1.6543E+01  1.2685E+02  3.5185E+00  0.0000E+00  1.8294E+00  2.8034E+00  1.0209E+00
             1.0842E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1624.92277329144        NO. OF FUNC. EVALS.: 101
 CUMULATIVE NO. OF FUNC. EVALS.:     1070
 NPARAMETR:  1.0023E+00  4.1207E-01  1.7104E-01  1.0926E+00  2.3438E-01  9.7002E-01  1.0000E-02  1.4532E+00  1.1822E+00  3.2961E-01
             2.5225E+00
 PARAMETER:  1.0227E-01 -7.8655E-01 -1.6659E+00  1.8860E-01 -1.3508E+00  6.9557E-02 -1.3412E+01  4.7379E-01  2.6742E-01 -1.0099E+00
             1.0253E+00
 GRADIENT:   1.3795E+00 -5.5712E-01 -1.5045E+00  1.2057E+00  2.9631E+00  1.1518E-01  0.0000E+00 -9.6784E-02 -4.9623E-01  1.4468E-01
             7.3853E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1624.92425295078        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1250
 NPARAMETR:  1.0019E+00  4.1279E-01  1.7185E-01  1.0926E+00  2.3493E-01  9.6980E-01  1.0000E-02  1.4568E+00  1.1839E+00  3.2406E-01
             2.5214E+00
 PARAMETER:  1.0186E-01 -7.8482E-01 -1.6611E+00  1.8854E-01 -1.3484E+00  6.9336E-02 -7.4810E+00  4.7622E-01  2.6877E-01 -1.0268E+00
             1.0248E+00
 GRADIENT:   5.4047E-01 -2.4287E-01 -4.4710E-01  7.0321E-02  1.8004E+00  5.4649E-02  0.0000E+00 -1.4917E-01 -1.0419E-01  1.0825E-01
             3.5203E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1624.92869158432        NO. OF FUNC. EVALS.: 174
 CUMULATIVE NO. OF FUNC. EVALS.:     1424
 NPARAMETR:  1.0018E+00  4.1113E-01  1.7110E-01  1.0915E+00  2.3403E-01  9.6975E-01  1.0000E-02  1.4592E+00  1.1867E+00  3.2778E-01
             2.5195E+00
 PARAMETER:  1.0155E-01 -7.8858E-01 -1.6653E+00  1.8761E-01 -1.3521E+00  6.9184E-02 -1.6632E+01  4.7745E-01  2.7093E-01 -1.0257E+00
             1.0238E+00
 GRADIENT:  -1.4073E-01  6.2585E-02  3.7615E-02  2.7651E-02  1.9578E-01 -9.8897E-03  0.0000E+00 -2.0992E-02 -1.6227E-02 -4.3757E-02
            -5.4265E-02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1424
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         2.0570E-03 -3.5833E-04  1.2222E-02 -4.7354E-03  1.7821E-02
 SE:             2.9241E-02  1.6446E-04  2.4017E-02  2.6830E-02  1.2296E-02
 N:                     100         100         100         100         100

 P VAL.:         9.4392E-01  2.9342E-02  6.1084E-01  8.5990E-01  1.4726E-01

 ETASHRINKSD(%)  2.0384E+00  9.9449E+01  1.9540E+01  1.0116E+01  5.8805E+01
 ETASHRINKVR(%)  4.0352E+00  9.9997E+01  3.5261E+01  1.9210E+01  8.3030E+01
 EBVSHRINKSD(%)  2.1083E+00  9.9396E+01  1.8820E+01  8.9308E+00  6.0094E+01
 EBVSHRINKVR(%)  4.1721E+00  9.9996E+01  3.4098E+01  1.7064E+01  8.4075E+01
 RELATIVEINF(%)  9.5733E+01  5.1082E-04  1.1297E+01  5.6206E+01  9.4708E-01
 EPSSHRINKSD(%)  2.8947E+01
 EPSSHRINKVR(%)  4.9514E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          500
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    918.93853320467269     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1624.9286915843184     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -705.99015837964566     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    20.98
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     7.86
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1624.929       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  4.11E-01  1.71E-01  1.09E+00  2.34E-01  9.70E-01  1.00E-02  1.46E+00  1.19E+00  3.24E-01  2.52E+00
 


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
+        1.13E+03
 
 TH 2
+       -3.37E+01  2.55E+03
 
 TH 3
+       -7.91E+01  2.53E+03  1.16E+04
 
 TH 4
+       -1.30E+01  2.25E+02 -1.55E+02  5.48E+02
 
 TH 5
+        1.50E+02 -7.60E+03 -1.50E+04 -8.87E+02  3.29E+04
 
 TH 6
+        3.40E+00 -2.02E+01  1.06E+01 -8.31E+00  4.35E+01  1.93E+02
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        1.78E+00 -2.12E+00 -8.62E+00 -6.03E-01  7.97E+00  9.83E-01  0.00E+00  3.85E+01
 
 TH 9
+        9.82E+00 -7.25E+01  8.77E+01 -9.44E+00  2.56E+02  2.01E+00  0.00E+00 -3.06E+00  8.95E+01
 
 TH10
+       -3.51E-01 -1.42E+02  9.85E-01  4.12E+00  6.31E+02  2.29E+00  0.00E+00  2.39E+01  2.06E+01  7.29E+01
 
 TH11
+       -1.80E+01 -3.53E+01 -3.06E+01 -3.81E+00  6.25E+01  1.54E+00  0.00E+00  9.98E+00  9.28E+00  8.16E+00  6.67E+01
 
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
1THERE ARE ERROR MESSAGES IN FILE PRDERR                                                                  
 #CPUT: Total CPU Time in Seconds,       28.910
Stop Time:
Thu Sep 30 00:01:20 CDT 2021
