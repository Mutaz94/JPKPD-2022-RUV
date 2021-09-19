Sat Sep 18 13:20:17 CDT 2021
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
$DATA ../../../../data/spa/S2/dat26.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m26.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1212.49731260598        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.9257E+01 -3.4540E+01  4.4260E+01 -3.5739E+01  8.8347E+00  2.9786E+01 -2.4759E+00 -1.7958E+02 -5.0774E+00 -4.6003E+00
            -7.0833E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1677.26025605280        NO. OF FUNC. EVALS.:  77
 CUMULATIVE NO. OF FUNC. EVALS.:       90
 NPARAMETR:  1.1106E+00  1.0282E+00  1.0295E+00  9.0438E-01  1.0490E+00  8.8449E-01  1.0209E+00  1.1193E+00  9.5333E-01  1.0314E+00
             1.1571E+00
 PARAMETER:  2.0492E-01  1.2779E-01  1.2907E-01 -5.0616E-04  1.4780E-01 -2.2741E-02  1.2067E-01  2.1274E-01  5.2209E-02  1.3090E-01
             2.4594E-01
 GRADIENT:   3.6689E+02 -1.1799E+02 -1.3738E+01 -1.3446E+02 -6.1780E+00 -5.4170E+01 -8.2406E+00  1.1998E+01  4.8867E+00  4.0002E+00
             9.6167E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1689.57633633798        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      168
 NPARAMETR:  1.1077E+00  1.3389E+00  1.0511E+00  7.5582E-01  1.2159E+00  9.1615E-01  1.0907E+00  1.0455E+00  1.0312E+00  1.4291E+00
             1.0952E+00
 PARAMETER:  2.0231E-01  3.9182E-01  1.4982E-01 -1.7995E-01  2.9547E-01  1.2427E-02  1.8680E-01  1.4447E-01  1.3073E-01  4.5703E-01
             1.9093E-01
 GRADIENT:   3.4611E+02 -3.3894E+01 -7.8816E-01 -5.2305E+01 -3.4962E+01 -3.6150E+01  2.5265E+01  5.0729E+00  9.7264E+00  2.9736E+01
             8.4626E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1719.15178758716        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      238
 NPARAMETR:  9.8835E-01  1.0101E+00  1.1474E+00  9.8530E-01  1.1026E+00  8.8644E-01  1.0639E+00  8.7053E-01  8.8108E-01  1.1358E+00
             8.6246E-01
 PARAMETER:  8.8283E-02  1.1007E-01  2.3753E-01  8.5194E-02  1.9770E-01 -2.0538E-02  1.6192E-01 -3.8659E-02 -2.6607E-02  2.2737E-01
            -4.7961E-02
 GRADIENT:   3.7136E+01 -2.8436E+01 -3.0913E+00 -3.1401E+01  5.4255E+00 -2.0718E+01 -6.0582E+00  1.7538E+00 -2.2714E+00 -5.4903E+00
             9.9433E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1720.28300861144        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      313
 NPARAMETR:  9.7051E-01  1.0416E+00  1.1358E+00  9.7936E-01  1.1157E+00  9.4364E-01  1.0944E+00  7.3237E-01  8.9135E-01  1.2058E+00
             8.6093E-01
 PARAMETER:  7.0069E-02  1.4078E-01  2.2737E-01  7.9143E-02  2.0948E-01  4.1991E-02  1.9018E-01 -2.1147E-01 -1.5022E-02  2.8715E-01
            -4.9740E-02
 GRADIENT:  -5.7093E+00 -7.7471E+00 -1.6136E+00 -7.3975E+00  1.7197E+00  4.2710E+00 -7.6566E-01  4.5753E-01  4.6855E-01  1.6163E+00
             1.4783E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1721.33087694869        NO. OF FUNC. EVALS.: 118
 CUMULATIVE NO. OF FUNC. EVALS.:      431
 NPARAMETR:  9.9558E-01  1.1291E+00  1.0444E+00  9.3055E-01  1.1100E+00  9.4710E-01  1.0665E+00  6.0023E-01  9.0791E-01  1.1726E+00
             8.5829E-01
 PARAMETER:  9.5571E-02  2.2146E-01  1.4343E-01  2.8021E-02  2.0435E-01  4.5653E-02  1.6443E-01 -4.1044E-01  3.3914E-03  2.5919E-01
            -5.2813E-02
 GRADIENT:  -2.7663E+00 -1.5333E+00  3.2824E+00 -3.3396E+00 -3.9437E+00 -2.5026E-01 -6.5164E-01 -3.0916E-01 -1.3642E+00 -2.2321E+00
            -2.4554E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1721.72498036934        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      606
 NPARAMETR:  9.9766E-01  1.4277E+00  8.4334E-01  7.4648E-01  1.1813E+00  9.4660E-01  9.1604E-01  4.9498E-01  1.0457E+00  1.2035E+00
             8.5905E-01
 PARAMETER:  9.7654E-02  4.5608E-01 -7.0388E-02 -1.9238E-01  2.6663E-01  4.5118E-02  1.2300E-02 -6.0324E-01  1.4473E-01  2.8520E-01
            -5.1929E-02
 GRADIENT:  -2.3153E+00  1.3142E+01 -4.7476E-01  1.1875E+01 -9.3105E-01 -1.1511E+00 -1.0060E+00  3.1158E-01 -4.9394E-02  1.0156E-02
            -1.2533E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1721.96946243446        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      783
 NPARAMETR:  9.9915E-01  1.5989E+00  7.0560E-01  6.2567E-01  1.2255E+00  9.4964E-01  8.5797E-01  2.4447E-01  1.1531E+00  1.2181E+00
             8.6218E-01
 PARAMETER:  9.9152E-02  5.6933E-01 -2.4871E-01 -3.6894E-01  3.0332E-01  4.8327E-02 -5.3182E-02 -1.3087E+00  2.4245E-01  2.9729E-01
            -4.8292E-02
 GRADIENT:  -4.3118E-01  4.0922E+00 -8.7788E-02  1.7328E+00 -6.9876E-01 -1.1352E-01  2.0332E-01  1.2254E-01  1.7768E-02  5.7120E-02
             2.1344E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1722.02717169809        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      960
 NPARAMETR:  9.9920E-01  1.5752E+00  7.1952E-01  6.4018E-01  1.2184E+00  9.5036E-01  8.6619E-01  1.2792E-01  1.1398E+00  1.2160E+00
             8.6247E-01
 PARAMETER:  9.9198E-02  5.5440E-01 -2.2917E-01 -3.4601E-01  2.9754E-01  4.9084E-02 -4.3651E-02 -1.9564E+00  2.3087E-01  2.9560E-01
            -4.7952E-02
 GRADIENT:  -1.1100E-01  2.4360E+00  5.0287E-02  9.8366E-01 -1.0005E-01  2.0086E-01  3.5688E-01  3.3685E-02 -5.8522E-02 -8.3458E-02
             1.4984E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1722.03150931620        NO. OF FUNC. EVALS.: 168
 CUMULATIVE NO. OF FUNC. EVALS.:     1128
 NPARAMETR:  9.9901E-01  1.5712E+00  7.2163E-01  6.4242E-01  1.2170E+00  9.5028E-01  8.6706E-01  1.1756E-01  1.1367E+00  1.2148E+00
             8.6248E-01
 PARAMETER:  9.9010E-02  5.5182E-01 -2.2624E-01 -3.4252E-01  2.9637E-01  4.9001E-02 -4.2642E-02 -2.0408E+00  2.2813E-01  2.9456E-01
            -4.7948E-02
 GRADIENT:   5.7894E+01  6.4326E+01  3.4461E-01  1.4757E+01  2.2341E+00  7.1909E+00  1.0113E+00  3.0417E-02  7.7223E-01  1.9606E-01
             2.0169E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1722.03198768477        NO. OF FUNC. EVALS.: 123
 CUMULATIVE NO. OF FUNC. EVALS.:     1251
 NPARAMETR:  9.9879E-01  1.5707E+00  7.2171E-01  6.4238E-01  1.2170E+00  9.5028E-01  8.6693E-01  1.1757E-01  1.1373E+00  1.2155E+00
             8.6247E-01
 PARAMETER:  9.8789E-02  5.5152E-01 -2.2613E-01 -3.4257E-01  2.9635E-01  4.8999E-02 -4.2792E-02 -2.0407E+00  2.2869E-01  2.9514E-01
            -4.7950E-02
 GRADIENT:   5.7364E+01  6.3620E+01  2.5418E-01  1.4488E+01  2.2268E+00  7.1941E+00  1.0054E+00  3.0852E-02  8.6505E-01  3.2840E-01
             2.3738E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1722.03275112335        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1430
 NPARAMETR:  9.9920E-01  1.5697E+00  7.2173E-01  6.4245E-01  1.2166E+00  9.5028E-01  8.6602E-01  1.1757E-01  1.1374E+00  1.2160E+00
             8.6247E-01
 PARAMETER:  9.9196E-02  5.5091E-01 -2.2610E-01 -3.4246E-01  2.9607E-01  4.8999E-02 -4.3844E-02 -2.0407E+00  2.2872E-01  2.9560E-01
            -4.7950E-02
 GRADIENT:  -6.5768E-02  7.9319E-02 -3.6650E-02 -2.2821E-01 -1.0370E-01  1.7577E-01  3.7027E-02  2.9072E-02 -8.8892E-02  8.1124E-02
             1.7790E-01

0ITERATION NO.:   58    OBJECTIVE VALUE:  -1722.03290834459        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1567
 NPARAMETR:  9.9920E-01  1.5697E+00  7.2173E-01  6.4245E-01  1.2166E+00  9.5028E-01  8.6602E-01  1.1698E-01  1.1374E+00  1.2160E+00
             8.6247E-01
 PARAMETER:  9.9196E-02  5.5090E-01 -2.2610E-01 -3.4246E-01  2.9607E-01  4.8999E-02 -4.3850E-02 -2.0457E+00  2.2874E-01  2.9559E-01
            -4.7952E-02
 GRADIENT:   1.3327E+05 -2.8677E+05  4.6767E+05  1.0421E+04 -8.4205E+04 -2.6335E+04 -8.3690E+04 -9.2154E+04 -7.6019E+05 -9.8893E+03
            -1.1072E+05

 #TERM:
0MINIMIZATION SUCCESSFUL
 HOWEVER, PROBLEMS OCCURRED WITH THE MINIMIZATION.
 REGARD THE RESULTS OF THE ESTIMATION STEP CAREFULLY, AND ACCEPT THEM ONLY
 AFTER CHECKING THAT THE COVARIANCE STEP PRODUCES REASONABLE OUTPUT.
 NO. OF FUNCTION EVALUATIONS USED:     1567
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -2.5570E-04 -2.3734E-02 -3.9730E-03  1.9596E-02 -3.3142E-02
 SE:             2.9849E-02  2.3323E-02  1.3182E-03  2.1849E-02  2.3659E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9316E-01  3.0886E-01  2.5795E-03  3.6978E-01  1.6127E-01

 ETASHRINKSD(%)  3.0227E-03  2.1865E+01  9.5584E+01  2.6803E+01  2.0740E+01
 ETASHRINKVR(%)  6.0453E-03  3.8950E+01  9.9805E+01  4.6422E+01  3.7178E+01
 EBVSHRINKSD(%)  3.5217E-01  2.0759E+01  9.6288E+01  3.0035E+01  1.7008E+01
 EBVSHRINKVR(%)  7.0311E-01  3.7209E+01  9.9862E+01  5.1049E+01  3.1122E+01
 RELATIVEINF(%)  9.9063E+01  2.6883E+00  1.4704E-02  1.9800E+00  1.5121E+01
 EPSSHRINKSD(%)  4.4356E+01
 EPSSHRINKVR(%)  6.9037E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1722.0329083445927     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -986.88208178085449     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.73
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     8.97
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1722.033       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.99E-01  1.57E+00  7.22E-01  6.42E-01  1.22E+00  9.50E-01  8.66E-01  1.17E-01  1.14E+00  1.22E+00  8.62E-01
 


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
+        7.33E+08
 
 TH 2
+        3.46E+08  1.23E+07
 
 TH 3
+        1.86E+09 -9.54E+07  4.28E+08
 
 TH 4
+        8.22E+08 -2.30E+08  1.41E+09  2.07E+08
 
 TH 5
+        1.57E+09 -2.08E+08  2.92E+08  7.30E+08  1.89E+08
 
 TH 6
+        6.24E+07 -5.03E+08  2.22E+08  3.75E+07 -3.23E+09 -1.74E+03
 
 TH 7
+       -2.83E+04  6.52E+08  3.30E+08  5.93E+07  1.19E+09 -2.98E+08  1.64E+08
 
 TH 8
+        2.01E+09 -2.18E+07 -3.05E+07  9.96E+08 -8.04E+07 -3.73E+08 -3.42E+08  1.57E+09
 
 TH 9
+       -3.77E+07 -2.43E+07  3.90E+07  4.43E+08 -3.98E+08 -1.10E+09  1.25E+08 -2.95E+08  6.03E+08
 
 TH10
+       -1.47E+08 -4.87E+07 -7.10E+08 -2.85E+08 -2.62E+08  1.24E+08 -1.08E+08 -6.35E+08 -1.09E+07  4.44E+07
 
 TH11
+       -1.17E+09 -2.09E+07  2.22E+09 -3.08E+09  7.66E+08  5.40E+08  6.70E+09 -3.23E+08 -1.10E+08 -8.82E+07  6.18E+09
 
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
 #CPUT: Total CPU Time in Seconds,       31.769
Stop Time:
Sat Sep 18 13:20:50 CDT 2021
