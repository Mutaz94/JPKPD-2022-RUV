Sat Sep 18 07:11:23 CDT 2021
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
$DATA ../../../../data/int/D/dat58.csv ignore=@
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
 (2E4.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m58.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   22600.4919836957        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.3161E+01  4.0386E+02  3.2054E+01  2.8050E+02  3.8625E+02 -2.9900E+03 -1.1540E+03 -1.1799E+02 -1.7876E+03 -1.2784E+03
            -4.5663E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1029.88996974997        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  1.7396E+00  1.4034E+00  8.4014E-01  1.4620E+00  7.1619E-01  5.3935E+00  5.1632E+00  9.9105E-01  3.0124E+00  3.3051E+00
             1.0032E+01
 PARAMETER:  6.5366E-01  4.3887E-01 -7.4184E-02  4.7982E-01 -2.3381E-01  1.7852E+00  1.7416E+00  9.1013E-02  1.2027E+00  1.2955E+00
             2.4058E+00
 GRADIENT:   2.2897E+01 -5.1952E+00 -1.9019E+01  4.0310E+01 -6.5154E+01  1.2759E+02  1.2736E+01  4.1305E+00  8.0171E+01  7.4269E+01
             4.3788E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1095.20724300838        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.1985E+00  1.2702E+00  3.3621E+01  1.9716E+00  2.1645E+00  6.1187E+00  1.1484E+01  7.7665E-01  5.5366E+00  2.7968E+00
             9.2817E+00
 PARAMETER:  2.8111E-01  3.3920E-01  3.6151E+00  7.7885E-01  8.7219E-01  1.9113E+00  2.5410E+00 -1.5277E-01  1.8114E+00  1.1285E+00
             2.3280E+00
 GRADIENT:  -3.6321E+00  1.3100E+01 -1.6879E+00 -3.3779E+00 -1.6770E+01  1.5039E+02  8.0173E+01  3.1171E-02  8.5982E+01  8.9393E+01
             4.3067E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1294.21152131677        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      226
 NPARAMETR:  1.0823E+00  4.1424E-01  1.4461E+02  1.5889E+00  2.8021E+00  3.2020E+00  6.9196E+00  2.5175E+00  1.2834E+00  1.3760E+00
             8.2866E+00
 PARAMETER:  1.7912E-01 -7.8130E-01  5.0740E+00  5.6307E-01  1.1304E+00  1.2638E+00  2.0344E+00  1.0233E+00  3.4953E-01  4.1921E-01
             2.2146E+00
 GRADIENT:  -2.5340E+01 -1.4181E+01 -5.1942E-01  2.9532E+01  4.4973E+01  5.4627E+00 -2.8003E+01  1.1347E-02 -9.2997E+00  2.0999E+01
             2.2580E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1328.05426504940        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      296
 NPARAMETR:  1.1587E+00  4.4185E-01  7.3074E+01  1.4169E+00  2.2600E+00  3.3063E+00  6.9722E+00  6.2095E+00  1.6735E+00  1.0196E+00
             7.0396E+00
 PARAMETER:  2.4733E-01 -7.1679E-01  4.3915E+00  4.4846E-01  9.1536E-01  1.2958E+00  2.0419E+00  1.9261E+00  6.1494E-01  1.1944E-01
             2.0515E+00
 GRADIENT:  -6.3852E+00 -1.2365E+01  1.7070E+00 -1.0441E+01 -1.4131E+01 -1.1160E+00  1.3114E-02 -1.2209E+00  6.4479E+00  7.9384E+00
            -4.1618E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1331.64078721279        NO. OF FUNC. EVALS.: 143
 CUMULATIVE NO. OF FUNC. EVALS.:      439             RESET HESSIAN, TYPE I
 NPARAMETR:  1.1689E+00  6.7599E-01  7.1734E+01  1.4146E+00  2.2626E+00  3.3298E+00  7.2444E+00  6.3458E+00  1.6783E+00  1.0173E+00
             7.0712E+00
 PARAMETER:  2.5604E-01 -2.9158E-01  4.3730E+00  4.4688E-01  9.1653E-01  1.3029E+00  2.0802E+00  1.9478E+00  6.1780E-01  1.1719E-01
             2.0560E+00
 GRADIENT:  -4.8014E+00  1.4470E+00  1.5812E+00  2.1716E+00 -1.8073E+01  1.6993E+00  2.5549E+01 -9.1629E-01  6.0887E+00  7.7988E+00
            -1.3602E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1332.20298233605        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      620
 NPARAMETR:  1.2253E+00  6.1165E-01  7.1696E+01  1.4155E+00  2.2890E+00  3.4686E+00  7.5279E+00  6.3469E+00  1.6658E+00  1.0166E+00
             7.0684E+00
 PARAMETER:  3.0315E-01 -3.9160E-01  4.3724E+00  4.4746E-01  9.2813E-01  1.3437E+00  2.1186E+00  1.9480E+00  6.1033E-01  1.1650E-01
             2.0556E+00
 GRADIENT:   2.1744E+00 -4.2591E-01  1.3995E+00 -7.8755E+00 -1.1170E+01  4.2390E+00  2.7811E-01 -8.9391E-01  6.6048E+00  7.9813E+00
            -4.7823E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1332.64381143878        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      795
 NPARAMETR:  1.2187E+00  6.1705E-01  7.1647E+01  1.4247E+00  2.3654E+00  3.4224E+00  7.5000E+00  6.3472E+00  1.5957E+00  1.0121E+00
             7.1052E+00
 PARAMETER:  2.9781E-01 -3.8281E-01  4.3718E+00  4.5395E-01  9.6096E-01  1.3304E+00  2.1149E+00  1.9480E+00  5.6733E-01  1.1204E-01
             2.0608E+00
 GRADIENT:   8.3984E-01 -3.0099E-01  7.6861E-01 -3.2603E+00  4.7211E+00 -6.6709E-01 -1.2029E+00 -5.9878E-01  3.5161E+00  8.1319E+00
            -7.1215E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1333.50579152689        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:      978
 NPARAMETR:  1.2302E+00  5.4638E-01  7.0330E+01  1.4625E+00  2.3477E+00  3.3830E+00  7.8569E+00  6.3632E+00  1.5191E+00  8.9146E-01
             7.1269E+00
 PARAMETER:  3.0718E-01 -5.0445E-01  4.3532E+00  4.8016E-01  9.5342E-01  1.3188E+00  2.1614E+00  1.9505E+00  5.1815E-01 -1.4891E-02
             2.0639E+00
 GRADIENT:   2.5586E+00 -3.5221E-01  1.5378E+00  3.0901E+00 -2.4590E-01 -4.8895E+00 -1.0153E+00 -1.5531E+00 -1.9880E+00  4.9437E+00
            -4.3858E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1333.64950376465        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1162
 NPARAMETR:  1.2294E+00  5.4183E-01  7.0061E+01  1.4674E+00  2.3481E+00  3.3850E+00  7.8817E+00  6.3664E+00  1.5289E+00  8.6844E-01
             7.1293E+00
 PARAMETER:  3.0652E-01 -5.1281E-01  4.3494E+00  4.8348E-01  9.5361E-01  1.3194E+00  2.1645E+00  1.9510E+00  5.2456E-01 -4.1061E-02
             2.0642E+00
 GRADIENT:   2.4189E+00 -2.5043E-01  2.2905E+00  2.0799E+00  7.6399E-01 -4.6642E+00 -8.3613E-01 -2.8420E+00 -2.2068E+00  4.8857E+00
            -3.2704E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1333.67699369921        NO. OF FUNC. EVALS.: 145
 CUMULATIVE NO. OF FUNC. EVALS.:     1307
 NPARAMETR:  1.2296E+00  5.4673E-01  7.0219E+01  1.4670E+00  2.3493E+00  3.4261E+00  7.8906E+00  6.3598E+00  1.5286E+00  8.6827E-01
             7.1373E+00
 PARAMETER:  3.0668E-01 -5.0380E-01  4.3516E+00  4.8325E-01  9.5412E-01  1.3314E+00  2.1657E+00  1.9500E+00  5.2433E-01 -4.1250E-02
             2.0653E+00
 GRADIENT:   2.3881E+00 -5.9336E-02  1.9246E+00  2.6947E+00  2.2439E-01 -1.9312E-02 -2.9706E-01 -2.2089E+00 -1.8653E+00  4.6975E+00
            -1.1375E+00

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1333.72734910358        NO. OF FUNC. EVALS.: 195
 CUMULATIVE NO. OF FUNC. EVALS.:     1502             RESET HESSIAN, TYPE I
 NPARAMETR:  1.2199E+00  5.4725E-01  7.0167E+01  1.4563E+00  2.3466E+00  3.4275E+00  7.8917E+00  6.3694E+00  1.5565E+00  8.6484E-01
             7.1412E+00
 PARAMETER:  2.9876E-01 -5.0285E-01  4.3509E+00  4.7589E-01  9.5297E-01  1.3318E+00  2.1658E+00  1.9515E+00  5.4247E-01 -4.5211E-02
             2.0659E+00
 GRADIENT:   3.4703E+00  4.5387E-01  2.1423E+00 -2.5783E-01  1.1405E+00  1.3376E+01  3.4187E+01 -2.5460E+00  7.0126E-01  4.8623E+00
             7.0769E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1333.76064307898        NO. OF FUNC. EVALS.: 128
 CUMULATIVE NO. OF FUNC. EVALS.:     1630
 NPARAMETR:  1.2085E+00  5.5486E-01  7.0291E+01  1.4566E+00  2.3414E+00  3.4179E+00  7.7556E+00  6.3642E+00  1.5742E+00  8.5850E-01
             7.1472E+00
 PARAMETER:  2.8912E-01 -4.8901E-01  4.3526E+00  4.7612E-01  9.5074E-01  1.3289E+00  2.1484E+00  1.9507E+00  5.5375E-01 -5.2568E-02
             2.0667E+00
 GRADIENT:  -1.0342E+00  4.8285E-01 -1.9209E+02  7.2618E+02 -3.6433E+02 -1.1289E+00 -8.1654E+01  4.3220E+02  6.2704E+02  9.6360E+00
            -1.6165E+02

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1630
 NO. OF SIG. DIGITS UNREPORTABLE

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.4053E-03  5.7094E-02 -5.6775E-03 -8.5977E-02 -7.4900E-03
 SE:             2.9641E-02  2.2170E-02  2.9837E-03  1.6842E-02  1.2885E-02
 N:                     100         100         100         100         100

 P VAL.:         8.2891E-01  1.0016E-02  5.7063E-02  3.3128E-07  5.6103E-01

 ETASHRINKSD(%)  6.9944E-01  2.5728E+01  9.0004E+01  4.3579E+01  5.6835E+01
 ETASHRINKVR(%)  1.3940E+00  4.4837E+01  9.9001E+01  6.8166E+01  8.1368E+01
 EBVSHRINKSD(%)  8.5474E-01  2.2127E+01  9.1292E+01  4.0845E+01  5.5106E+01
 EBVSHRINKVR(%)  1.7022E+00  3.9357E+01  9.9242E+01  6.5007E+01  7.9845E+01
 RELATIVEINF(%)  9.8252E+01  3.2372E+01  2.4750E-01  1.8705E+01  6.5263E+00
 EPSSHRINKSD(%)  9.4322E+00
 EPSSHRINKVR(%)  1.7975E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1333.7606430789788     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       320.32871668943199     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    56.72
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    19.29
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1333.761       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.21E+00  5.55E-01  7.03E+01  1.46E+00  2.34E+00  3.42E+00  7.76E+00  6.36E+00  1.57E+00  8.59E-01  7.15E+00
 


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
+        6.62E+01
 
 TH 2
+        3.79E+01 -5.12E+02
 
 TH 3
+       -1.53E-02  1.08E-01  2.28E+00
 
 TH 4
+        1.95E+01 -2.57E+02 -1.58E+00  1.83E+05
 
 TH 5
+       -4.26E+00  5.93E+01  1.63E-01 -4.28E+02  4.31E+04
 
 TH 6
+       -1.91E-01  1.07E+01  3.00E-03  5.69E+00 -1.53E+00  1.60E+01
 
 TH 7
+       -1.08E+00  2.01E+01  9.13E-02 -4.97E+01  1.19E+01 -3.50E-01  3.17E+02
 
 TH 8
+        5.71E-01 -3.81E+00 -3.34E+01  4.40E+01 -8.06E+00 -5.93E-02 -2.35E+00  1.39E+03
 
 TH 9
+        2.29E+01 -3.29E+02 -4.00E-01  3.68E+02 -8.14E+01  6.47E+00 -6.02E+03  1.16E+01  1.15E+05
 
 TH10
+        1.06E+02 -1.55E+03  8.27E-02 -8.74E+02  2.12E+02  2.95E+01  4.68E+01 -5.35E+00 -9.81E+02 -4.17E+03
 
 TH11
+       -2.97E+00  5.10E+00  2.28E-01 -1.23E+04  3.83E+03  5.92E-01  1.68E+00 -6.67E+00 -4.14E+00  2.59E+01  4.26E+02
 
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
 #CPUT: Total CPU Time in Seconds,       76.132
Stop Time:
Sat Sep 18 07:12:40 CDT 2021
