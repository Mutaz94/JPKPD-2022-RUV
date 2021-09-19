Sat Sep 18 15:27:52 CDT 2021
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
$DATA ../../../../data/spa/D/dat57.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m57.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   14922.5138391671        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.4143E+01  2.4598E+02  4.1676E+00  1.0274E+02  2.2893E+02 -2.3749E+03 -8.5324E+02 -1.2795E+02 -1.3480E+03 -7.4940E+02
            -2.7743E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -631.691752631155        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.5257E+00  1.1863E+00  9.2920E-01  1.6548E+00  1.1948E+00  1.7300E+00  1.0962E+00  9.6963E-01  1.0314E+00  1.1213E+00
             1.5016E+01
 PARAMETER:  5.2248E-01  2.7081E-01  2.6571E-02  6.0368E-01  2.7796E-01  6.4811E-01  1.9181E-01  6.9161E-02  1.3092E-01  2.1447E-01
             2.8091E+00
 GRADIENT:   3.1796E+01  2.0984E+01 -1.7291E+00  3.5755E+01 -8.9719E+00  2.1578E+01  1.2388E+00  4.3884E+00  1.1263E+01  3.4421E+00
             6.0797E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -642.952130731247        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.4232E+00  1.0117E+00  1.3555E+00  1.6956E+00  2.2375E+00  1.5295E+00  1.6770E+00  2.4938E-01  7.8899E-01  5.0571E+00
             1.4037E+01
 PARAMETER:  4.5292E-01  1.1160E-01  4.0415E-01  6.2803E-01  9.0535E-01  5.2492E-01  6.1701E-01 -1.2888E+00 -1.3700E-01  1.7208E+00
             2.7417E+00
 GRADIENT:   2.7732E+01  2.6334E+01  1.1386E+00  5.5809E+01 -1.6012E+01 -5.5624E+00  1.0173E+00  8.2970E-02  4.7981E+00  1.1989E+01
             1.7669E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -664.525674772582        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.1793E+00  2.1096E-01  1.2062E+00  1.6282E+00  2.2633E+00  1.4344E+00  3.4545E+00  1.0000E-02  1.5687E-01  4.1132E+00
             1.2873E+01
 PARAMETER:  2.6488E-01 -1.4561E+00  2.8745E-01  5.8749E-01  9.1683E-01  4.6076E-01  1.3397E+00 -4.9233E+00 -1.7523E+00  1.5142E+00
             2.6552E+00
 GRADIENT:   2.0010E+00  8.5135E-02  1.3552E+01 -7.2128E+01 -1.3842E+01  1.5869E+01  2.6354E+00  0.0000E+00  8.2630E-01 -9.2960E-02
             4.1278E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -685.215601879211        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.1335E+00  3.6182E-02  6.0331E-01  1.6470E+00  4.6227E+00  1.4153E+00  6.8074E+00  1.0000E-02  1.2909E-02  6.9375E+00
             1.1568E+01
 PARAMETER:  2.2528E-01 -3.2192E+00 -4.0533E-01  5.9897E-01  1.6310E+00  4.4732E-01  2.0180E+00 -1.1301E+01 -4.2499E+00  2.0369E+00
             2.5483E+00
 GRADIENT:   1.2580E+01  1.6670E+00 -1.3324E+01  6.7182E+01 -4.3176E+00 -2.5047E+01  2.1075E-01  0.0000E+00  9.3773E-03  5.8787E+00
            -6.3828E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -723.737007973610        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      369
 NPARAMETR:  8.3076E-01  4.0134E-02  1.1672E-01  9.5806E-01  5.0302E+00  2.0043E+00  3.4372E-01  1.0000E-02  1.0000E-02  9.5245E-01
             1.3117E+01
 PARAMETER: -8.5416E-02 -3.1155E+00 -2.0479E+00  5.7156E-02  1.7155E+00  7.9529E-01 -9.6793E-01 -1.8917E+01 -5.4120E+00  5.1282E-02
             2.6739E+00
 GRADIENT:  -2.2914E+01 -2.4027E+00 -3.3529E+01  1.1382E+02  8.4390E-01  3.3176E+01  1.4587E-02  0.0000E+00  0.0000E+00  2.9123E+00
             6.6870E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -757.376977693523        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      439
 NPARAMETR:  5.0075E-01  2.2154E-02  3.0263E-02  3.3922E-01  9.4569E+00  1.6602E+00  2.6801E-02  1.0000E-02  1.0000E-02  4.4463E-01
             1.0390E+01
 PARAMETER: -5.9165E-01 -3.7097E+00 -3.3978E+00 -9.8110E-01  2.3467E+00  6.0692E-01 -3.5193E+00 -2.8593E+01 -6.1041E+00 -7.1052E-01
             2.4408E+00
 GRADIENT:  -1.2591E+01  6.9962E+00  1.6171E+01 -6.9824E+00 -1.2190E+00 -1.4686E+01  2.0050E-03  0.0000E+00  0.0000E+00  4.1827E-02
            -3.4245E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -758.530285169845        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:      512
 NPARAMETR:  4.2339E-01  1.2297E-02  1.6597E-02  2.1241E-01  1.5407E+01  1.6964E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.9712E-01
             1.0367E+01
 PARAMETER: -7.5946E-01 -4.2984E+00 -3.9985E+00 -1.4493E+00  2.8348E+00  6.2849E-01 -5.5105E+00 -3.4112E+01 -7.4434E+00 -1.1136E+00
             2.4386E+00
 GRADIENT:   4.2693E+01  3.0573E+00 -1.7298E+01  1.1842E+01 -1.4064E-01 -7.9308E+00  0.0000E+00  0.0000E+00  0.0000E+00  4.4315E-04
            -4.0172E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -758.531494664208        NO. OF FUNC. EVALS.: 110
 CUMULATIVE NO. OF FUNC. EVALS.:      622
 NPARAMETR:  4.1813E-01  1.1940E-02  1.6099E-02  2.0719E-01  1.5814E+01  1.6970E+00  1.0000E-02  1.0000E-02  1.0000E-02  2.9217E-01
             1.0364E+01
 PARAMETER: -7.7196E-01 -4.3279E+00 -4.0290E+00 -1.4741E+00  2.8609E+00  6.2885E-01 -5.6161E+00 -3.4399E+01 -7.5089E+00 -1.1304E+00
             2.4383E+00
 GRADIENT:   3.7239E+01  2.9254E+00 -2.5257E+01  8.7230E+00 -1.3199E-01 -8.7459E+00  0.0000E+00  0.0000E+00  0.0000E+00  3.5845E-04
            -4.2061E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -760.716606599875        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      798
 NPARAMETR:  4.2199E-01  1.0527E-02  1.8550E-02  2.3088E-01  1.4461E+01  1.7213E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.0741E-01
             1.0797E+01
 PARAMETER: -7.6278E-01 -4.4538E+00 -3.8873E+00 -1.3658E+00  2.7715E+00  6.4307E-01 -5.5042E+00 -3.3730E+01 -7.6869E+00 -1.0796E+00
             2.4793E+00
 GRADIENT:  -2.1591E+00  7.0885E-01  2.1002E+00 -2.0755E+00 -5.2722E-02 -1.3193E+00  0.0000E+00  0.0000E+00  0.0000E+00  4.2445E-05
            -6.6575E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -760.756287645560        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      976            RESET HESSIAN, TYPE II
 NPARAMETR:  4.2706E-01  1.0000E-02  1.8960E-02  2.3522E-01  1.5554E+01  1.7290E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.0566E-01
             1.0819E+01
 PARAMETER: -7.5083E-01 -4.5493E+00 -3.8654E+00 -1.3472E+00  2.8443E+00  6.4756E-01 -5.4591E+00 -3.3741E+01 -7.7906E+00 -1.0853E+00
             2.4813E+00
 GRADIENT:   5.3145E+00  0.0000E+00  6.4533E+00  3.1619E+00 -7.3677E-04  9.9562E-01  0.0000E+00  0.0000E+00  0.0000E+00  6.0532E-06
             1.8972E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -760.756370808926        NO. OF FUNC. EVALS.: 137
 CUMULATIVE NO. OF FUNC. EVALS.:     1113            RESET HESSIAN, TYPE II
 NPARAMETR:  4.2690E-01  1.0000E-02  1.8961E-02  2.3521E-01  1.5976E+01  1.7289E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.0565E-01
             1.0819E+01
 PARAMETER: -7.5120E-01 -4.5493E+00 -3.8654E+00 -1.3473E+00  2.8711E+00  6.4747E-01 -5.4591E+00 -3.3741E+01 -7.7906E+00 -1.0853E+00
             2.4813E+00
 GRADIENT:   5.0902E+00  0.0000E+00  6.6800E+00  2.9962E+00  1.9543E-04  9.7980E-01  0.0000E+00  0.0000E+00  0.0000E+00  4.3207E-06
             1.9940E+00

0ITERATION NO.:   57    OBJECTIVE VALUE:  -760.756370808926        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1170
 NPARAMETR:  4.2690E-01  1.0000E-02  1.8961E-02  2.3521E-01  1.5976E+01  1.7289E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.0565E-01
             1.0819E+01
 PARAMETER: -7.5120E-01 -4.5493E+00 -3.8654E+00 -1.3473E+00  2.8711E+00  6.4747E-01 -5.4591E+00 -3.3741E+01 -7.7906E+00 -1.0853E+00
             2.4813E+00
 GRADIENT:   9.7253E-03  0.0000E+00 -6.7063E-02  1.2069E-01 -6.2854E-05 -7.6760E-03  0.0000E+00  0.0000E+00  0.0000E+00  3.7975E-06
             8.9498E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1170
 NO. OF SIG. DIGITS IN FINAL EST.:  3.5
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5228E-03  3.2415E-06  9.8565E-05 -1.8874E-04  4.2822E-06
 SE:             2.9114E-02  1.3198E-06  3.1474E-04  3.6919E-04  2.0683E-05
 N:                     100         100         100         100         100

 P VAL.:         9.5829E-01  1.4049E-02  7.5416E-01  6.0918E-01  8.3598E-01

 ETASHRINKSD(%)  2.4655E+00  9.9996E+01  9.8946E+01  9.8763E+01  9.9931E+01
 ETASHRINKVR(%)  4.8702E+00  1.0000E+02  9.9989E+01  9.9985E+01  1.0000E+02
 EBVSHRINKSD(%)  2.7051E+00  9.9992E+01  9.8940E+01  9.8722E+01  9.9931E+01
 EBVSHRINKVR(%)  5.3370E+00  1.0000E+02  9.9989E+01  9.9984E+01  1.0000E+02
 RELATIVEINF(%)  1.4645E+00  9.2018E-08  2.7673E-05  3.0383E-05  9.0645E-07
 EPSSHRINKSD(%)  7.2311E+00
 EPSSHRINKVR(%)  1.3939E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -760.75637080892591     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -25.605544245187730     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    13.22
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.22
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -760.756       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         4.27E-01  1.00E-02  1.90E-02  2.35E-01  1.60E+01  1.73E+00  1.00E-02  1.00E-02  1.00E-02  3.06E-01  1.08E+01
 


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
+        1.91E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.20E+04  0.00E+00  2.16E+06
 
 TH 4
+       -5.01E+02  0.00E+00 -2.07E+05  2.14E+04
 
 TH 5
+        1.46E-01  0.00E+00 -4.10E+00  3.33E-01  1.45E-04
 
 TH 6
+        1.65E+00  0.00E+00  1.04E+03 -1.14E+02  1.98E-03  5.71E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+       -1.71E-02  0.00E+00 -7.31E-02 -4.54E-03  2.60E-04 -9.89E-04  0.00E+00  0.00E+00  0.00E+00  4.55E-04
 
 TH11
+       -1.82E+01  0.00E+00  3.69E+02 -2.46E+01 -2.06E-03  1.02E+00  0.00E+00  0.00E+00  0.00E+00  4.70E-04  3.47E+00
 
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
 #CPUT: Total CPU Time in Seconds,       19.506
Stop Time:
Sat Sep 18 15:28:13 CDT 2021
