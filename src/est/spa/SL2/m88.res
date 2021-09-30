Wed Sep 29 16:12:03 CDT 2021
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
$DATA ../../../../data/spa/SL2/dat88.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m88.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1633.39811543434        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.1518E+02 -5.1472E+01 -4.6652E+01  1.1670E+01  8.6194E+01  4.6736E+01 -8.5679E+00 -4.4029E+00  6.3484E+00  2.1999E+01
            -7.6032E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1644.47158663692        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:      172
 NPARAMETR:  1.0129E+00  1.0867E+00  1.0685E+00  9.9967E-01  1.0275E+00  9.8153E-01  1.0533E+00  1.0484E+00  9.8924E-01  8.2859E-01
             1.3445E+00
 PARAMETER:  1.1279E-01  1.8313E-01  1.6627E-01  9.9672E-02  1.2709E-01  8.1359E-02  1.5192E-01  1.4725E-01  8.9180E-02 -8.8033E-02
             3.9601E-01
 GRADIENT:   1.2543E+01 -2.7814E+01 -2.3394E+01  6.5104E+00  5.3302E+01 -1.4321E+00 -5.1311E+00 -7.4436E+00  3.5828E+00  9.1750E+00
             4.9621E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1649.06125545923        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      350
 NPARAMETR:  1.0021E+00  1.1730E+00  9.3785E-01  9.4119E-01  9.2055E-01  1.0011E+00  1.1756E+00  1.7677E+00  9.5096E-01  3.8051E-01
             1.2525E+00
 PARAMETER:  1.0207E-01  2.5960E-01  3.5837E-02  3.9391E-02  1.7211E-02  1.0106E-01  2.6177E-01  6.6971E-01  4.9712E-02 -8.6623E-01
             3.2510E-01
 GRADIENT:  -1.1271E+01  1.2244E+01  6.9688E+00  3.0811E+00 -4.0383E+01  5.3166E+00  6.3809E+00  1.3255E+01  5.5873E+00  8.0008E-01
             2.2705E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1651.99531862731        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      526
 NPARAMETR:  1.0077E+00  1.1851E+00  8.7550E-01  9.2713E-01  9.3006E-01  9.8624E-01  1.1377E+00  1.3109E+00  9.4825E-01  4.9236E-01
             1.1934E+00
 PARAMETER:  1.0770E-01  2.6984E-01 -3.2956E-02  2.4341E-02  2.7490E-02  8.6145E-02  2.2897E-01  3.7073E-01  4.6862E-02 -6.0855E-01
             2.7679E-01
 GRADIENT:   1.0288E+00  1.0744E+00 -6.6426E-01  4.5309E+00  6.2596E+00 -7.7641E-01  8.3539E-02 -5.6451E-01 -6.2683E-01 -2.1713E+00
            -3.3957E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.21833298801        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      702
 NPARAMETR:  1.0078E+00  1.3348E+00  8.4583E-01  8.3463E-01  9.8806E-01  9.8797E-01  1.0177E+00  1.3741E+00  1.0445E+00  6.0711E-01
             1.1876E+00
 PARAMETER:  1.0781E-01  3.8880E-01 -6.7437E-02 -8.0765E-02  8.7988E-02  8.7894E-02  1.1752E-01  4.1783E-01  1.4356E-01 -3.9904E-01
             2.7193E-01
 GRADIENT:   6.7397E-01  2.5940E+00  3.2115E-01  2.7143E+00 -1.7829E+00 -2.2573E-01  1.2563E-01  2.6932E-01  1.9620E-01  5.2778E-01
            -1.0691E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.25292187792        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      877
 NPARAMETR:  1.0092E+00  1.4756E+00  6.8148E-01  7.4108E-01  9.7971E-01  9.8873E-01  9.6124E-01  1.3035E+00  1.1109E+00  5.8093E-01
             1.1821E+00
 PARAMETER:  1.0916E-01  4.8903E-01 -2.8348E-01 -1.9965E-01  7.9502E-02  8.8670E-02  6.0472E-02  3.6507E-01  2.0515E-01 -4.4312E-01
             2.6725E-01
 GRADIENT:   1.9043E+00  1.0513E+01  1.3713E+00  6.1155E+00 -5.8435E+00 -3.4328E-01  8.7298E-01  3.0738E-01  1.8188E-01  8.7996E-01
            -1.6320E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1652.83275964917        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1062
 NPARAMETR:  1.0103E+00  1.7242E+00  2.6150E-01  5.3729E-01  8.4283E-01  9.9270E-01  8.9386E-01  7.9097E-01  1.2136E+00  3.6937E-01
             1.1566E+00
 PARAMETER:  1.1029E-01  6.4477E-01 -1.2413E+00 -5.2122E-01 -7.0987E-02  9.2671E-02 -1.2209E-02 -1.3450E-01  2.9358E-01 -8.9597E-01
             2.4552E-01
 GRADIENT:   1.1297E+01  4.4910E+01  1.4337E+01 -8.4015E+00 -4.8130E+01  1.2315E+00  1.0225E+01  7.6812E-02  2.1229E-02  3.0066E+00
            -8.9756E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.51386510213        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1240
 NPARAMETR:  9.9637E-01  1.8833E+00  1.3602E-01  3.8991E-01  8.6552E-01  9.9134E-01  7.8151E-01  5.3720E-01  1.4540E+00  2.5776E-01
             1.1735E+00
 PARAMETER:  9.6362E-02  7.3301E-01 -1.8949E+00 -8.4183E-01 -4.4423E-02  9.1301E-02 -1.4652E-01 -5.2138E-01  4.7435E-01 -1.2557E+00
             2.5995E-01
 GRADIENT:  -5.1296E+00 -4.4651E+00  5.8836E+00 -1.2927E+01 -1.4939E+01  1.4866E+00 -5.4848E+00 -1.2617E+00 -5.1170E+00 -1.5057E-01
             1.5153E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1657.18120643336        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1421
 NPARAMETR:  1.0001E+00  1.9350E+00  1.2844E-01  3.6956E-01  8.8808E-01  9.8988E-01  7.8731E-01  6.1938E-01  1.5976E+00  2.6951E-01
             1.1016E+00
 PARAMETER:  1.0007E-01  7.6011E-01 -1.9523E+00 -8.9545E-01 -1.8698E-02  8.9828E-02 -1.3913E-01 -3.7903E-01  5.6849E-01 -1.2112E+00
             1.9676E-01
 GRADIENT:   4.7680E+00  1.4042E+01  4.5080E+00 -7.2548E-01 -6.9090E+00 -1.4215E-01  1.4978E-01 -2.5870E+00  1.0253E+00 -1.1326E+00
            -6.2517E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1657.42524146737        NO. OF FUNC. EVALS.: 201
 CUMULATIVE NO. OF FUNC. EVALS.:     1622             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0008E+00  1.9276E+00  1.3031E-01  3.7051E-01  8.9050E-01  9.9014E-01  7.8716E-01  6.5249E-01  1.5953E+00  3.4346E-01
             1.1154E+00
 PARAMETER:  1.0083E-01  7.5627E-01 -1.9378E+00 -8.9288E-01 -1.5969E-02  9.0094E-02 -1.3932E-01 -3.2695E-01  5.6705E-01 -9.6870E-01
             2.0925E-01
 GRADIENT:   3.4690E+02  7.6994E+02  1.6274E+01  8.7959E+01 -5.4894E+00  3.3574E+01  1.1170E+01 -2.2734E+00  1.3727E+01  1.8326E+00
             5.1261E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1657.48984464264        NO. OF FUNC. EVALS.: 159
 CUMULATIVE NO. OF FUNC. EVALS.:     1781
 NPARAMETR:  1.0008E+00  1.9276E+00  1.3031E-01  3.7044E-01  8.9303E-01  9.8977E-01  7.8237E-01  6.5249E-01  1.5953E+00  3.2976E-01
             1.1112E+00
 PARAMETER:  1.0083E-01  7.5627E-01 -1.9378E+00 -8.9307E-01 -1.3136E-02  8.9720E-02 -1.4543E-01 -3.2696E-01  5.6706E-01 -1.0094E+00
             2.0544E-01
 GRADIENT:   5.4821E+00 -6.1182E+00  3.9168E+00 -6.7240E-01 -3.5927E+00 -1.1127E-01  7.2549E-01 -2.4803E+00  3.2616E+00  1.1118E-01
             6.3144E-01

0ITERATION NO.:   52    OBJECTIVE VALUE:  -1657.49373662901        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1838
 NPARAMETR:  1.0008E+00  1.9276E+00  1.3031E-01  3.7019E-01  8.9440E-01  9.8998E-01  7.8012E-01  6.5249E-01  1.5953E+00  3.3071E-01
             1.1101E+00
 PARAMETER:  1.0083E-01  7.5627E-01 -1.9378E+00 -8.9373E-01 -1.1600E-02  8.9934E-02 -1.4830E-01 -3.2696E-01  5.6706E-01 -1.0065E+00
             2.0449E-01
 GRADIENT:   5.3954E+00 -8.6652E+00  3.5411E+00  7.2546E-02 -8.1567E-02 -3.8429E-02 -8.6383E-02 -2.5118E+00  2.9694E+00 -1.2384E-02
             7.3499E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1838
 NO. OF SIG. DIGITS IN FINAL EST.:  2.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -1.6357E-03 -1.4614E-02 -4.7819E-03  2.1020E-02 -2.6003E-02
 SE:             2.9833E-02  2.8274E-02  6.4486E-03  2.4252E-02  1.1095E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5627E-01  6.0524E-01  4.5836E-01  3.8609E-01  1.9091E-02

 ETASHRINKSD(%)  5.6176E-02  5.2775E+00  7.8397E+01  1.8753E+01  6.2832E+01
 ETASHRINKVR(%)  1.1232E-01  1.0276E+01  9.5333E+01  3.3989E+01  8.6185E+01
 EBVSHRINKSD(%)  5.1888E-01  6.0429E+00  8.1358E+01  1.6143E+01  6.3940E+01
 EBVSHRINKVR(%)  1.0351E+00  1.1721E+01  9.6525E+01  2.9680E+01  8.6997E+01
 RELATIVEINF(%)  9.8206E+01  2.0422E+01  1.3916E+00  1.1360E+01  2.3766E+00
 EPSSHRINKSD(%)  4.3259E+01
 EPSSHRINKVR(%)  6.7805E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1657.4937366290135     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -922.34291006527530     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    23.68
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.09
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1657.494       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.00E+00  1.93E+00  1.30E-01  3.70E-01  8.94E-01  9.90E-01  7.80E-01  6.52E-01  1.60E+00  3.31E-01  1.11E+00
 


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
+        1.85E+07
 
 TH 2
+       -9.23E+00  8.91E+04
 
 TH 3
+       -5.97E+02  4.23E+02  4.74E+03
 
 TH 4
+        3.75E+02  1.40E+02 -2.23E+06  1.71E+06
 
 TH 5
+        2.09E+07 -4.45E+02 -3.70E+03  2.67E+03  2.36E+07
 
 TH 6
+       -6.08E+02 -1.27E+00 -2.30E+02  1.77E+02 -2.78E+00  1.99E+02
 
 TH 7
+       -5.84E+02  4.44E+00 -3.43E+02  1.81E+02  3.00E+01 -4.02E-01  2.58E+02
 
 TH 8
+       -8.73E+06 -1.89E+00  1.67E+04 -1.29E+04 -9.85E+06  2.82E+02  2.79E+02  4.12E+06
 
 TH 9
+        1.44E+02 -7.46E+00  9.86E+03  1.26E+06 -2.34E+06  6.65E+01  6.88E+01 -9.85E+05  2.33E+05
 
 TH10
+        5.61E+06 -1.01E+01  1.34E+03 -1.02E+03  6.33E+06 -1.82E+02 -1.52E+02 -2.65E+06 -6.30E+05  1.70E+06
 
 TH11
+        8.23E+06 -1.17E+01  2.44E+03 -1.94E+03 -2.31E+03 -2.63E+02 -2.50E+02 -3.88E+06 -9.23E+05  2.50E+06  3.66E+06
 
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
 #CPUT: Total CPU Time in Seconds,       29.831
Stop Time:
Wed Sep 29 16:12:35 CDT 2021
