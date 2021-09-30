Wed Sep 29 13:05:41 CDT 2021
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
$DATA ../../../../data/spa/A2/dat84.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m84.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -591.258971636295        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.7876E+02  3.7049E+01 -2.5661E+01  6.9491E+01  1.0924E+02  8.7546E+01 -5.8354E+01  3.0580E+00 -1.1371E+02 -4.2300E+01
            -1.9924E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1364.50517291967        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.9767E-01  1.0544E+00  1.1984E+00  9.5675E-01  1.0299E+00  7.0685E-01  1.0208E+00  9.4008E-01  1.0150E+00  7.9284E-01
             2.3631E+00
 PARAMETER:  9.7669E-02  1.5296E-01  2.8096E-01  5.5784E-02  1.2948E-01 -2.4694E-01  1.2059E-01  3.8213E-02  1.1488E-01 -1.3213E-01
             9.5997E-01
 GRADIENT:  -6.1164E+01 -2.2375E+01  5.8837E+00 -5.6671E+01 -2.8262E+01 -3.3436E+01 -4.0536E+00  3.4334E+00 -1.2761E+01  1.6469E+01
            -2.3495E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1377.18115512425        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.8251E-01  9.3576E-01  9.9947E-01  1.0430E+00  8.8554E-01  8.6421E-01  8.9239E-01  8.3503E-01  9.9420E-01  2.0855E-01
             2.4683E+00
 PARAMETER:  8.2355E-02  3.3603E-02  9.9473E-02  1.4207E-01 -2.1562E-02 -4.5939E-02 -1.3847E-02 -8.0291E-02  9.4184E-02 -1.4676E+00
             1.0035E+00
 GRADIENT:  -7.2709E+01  4.9085E+00  1.8527E+01 -2.4749E+01 -5.8810E+01  2.8086E+01 -1.6121E+01  2.9869E+00 -1.8868E+01  4.5517E-01
            -2.0043E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1413.05925697051        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  1.0229E+00  9.9820E-01  1.3305E+00  1.0489E+00  1.0752E+00  7.6181E-01  1.0372E+00  3.5194E-01  9.7816E-01  2.1725E-01
             3.1072E+00
 PARAMETER:  1.2260E-01  9.8197E-02  3.8553E-01  1.4771E-01  1.7248E-01 -1.7206E-01  1.3648E-01 -9.4429E-01  7.7913E-02 -1.4267E+00
             1.2337E+00
 GRADIENT:  -3.7893E+00 -3.7950E+00  6.4154E-01 -7.4340E+00 -3.8232E+00  2.0717E+00  8.4411E-01  5.6259E-01  2.5033E-01  1.2125E+00
             6.3171E+00

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1413.98915752636        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  1.0241E+00  1.0056E+00  1.3387E+00  1.0512E+00  1.0846E+00  7.5690E-01  1.0163E+00  7.6116E-02  9.8220E-01  4.4569E-02
             3.0951E+00
 PARAMETER:  1.2381E-01  1.0559E-01  3.9173E-01  1.4989E-01  1.8119E-01 -1.7853E-01  1.1614E-01 -2.4755E+00  8.2044E-02 -3.0107E+00
             1.2298E+00
 GRADIENT:   5.8616E-02  7.0280E-01  5.4639E-01  3.9882E-01 -1.5772E+00 -1.5634E-01 -8.7789E-02  2.6781E-02 -3.2220E-01  5.0314E-02
             8.3589E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1414.92273598157        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      480
 NPARAMETR:  1.0437E+00  8.7526E-01  1.4330E+00  1.1606E+00  1.0714E+00  7.7037E-01  1.1027E+00  2.3230E-02  9.3061E-01  1.0810E-02
             3.1397E+00
 PARAMETER:  1.4275E-01 -3.3232E-02  4.5974E-01  2.4896E-01  1.6900E-01 -1.6089E-01  1.9772E-01 -3.6623E+00  2.8090E-02 -4.4273E+00
             1.2441E+00
 GRADIENT:   6.2113E+00  7.4705E+00  6.8600E-01  1.0146E+01 -2.9952E+00  1.7472E+00 -1.9957E-01  2.5739E-03 -7.1918E-01  2.6197E-03
            -5.2572E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1415.44404138689        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      656
 NPARAMETR:  1.0400E+00  5.4131E-01  1.4824E+00  1.3716E+00  9.7336E-01  7.6823E-01  1.4554E+00  1.0000E-02  8.3544E-01  1.0000E-02
             3.1413E+00
 PARAMETER:  1.3923E-01 -5.1376E-01  4.9363E-01  4.1601E-01  7.2996E-02 -1.6367E-01  4.7528E-01 -1.0913E+01 -7.9800E-02 -1.1541E+01
             1.2446E+00
 GRADIENT:   4.8969E+00  6.2015E+00  1.9979E+00  1.7594E+01 -6.5483E+00  1.1993E+00 -3.9464E-01  0.0000E+00 -1.8977E+00  0.0000E+00
             1.9619E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1415.75054987651        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:      835
 NPARAMETR:  1.0347E+00  3.4540E-01  1.5113E+00  1.4860E+00  9.2802E-01  7.6366E-01  1.8712E+00  1.0000E-02  7.9636E-01  1.0000E-02
             3.1356E+00
 PARAMETER:  1.3408E-01 -9.6305E-01  5.1295E-01  4.9611E-01  2.5297E-02 -1.6964E-01  7.2656E-01 -1.9948E+01 -1.2770E-01 -1.9852E+01
             1.2428E+00
 GRADIENT:  -2.2794E+00  2.3966E+00  9.5072E-01  1.1165E+01 -3.6240E+00 -2.4624E-01 -2.5546E-01  0.0000E+00 -1.5784E+00  0.0000E+00
            -3.1081E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1415.89860232937        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1012
 NPARAMETR:  1.0318E+00  1.9012E-01  1.5517E+00  1.5785E+00  8.9970E-01  7.6224E-01  2.6516E+00  1.0000E-02  7.7161E-01  1.0000E-02
             3.1279E+00
 PARAMETER:  1.3126E-01 -1.5601E+00  5.3934E-01  5.5649E-01 -5.6977E-03 -1.7150E-01  1.0752E+00 -3.3416E+01 -1.5927E-01 -3.1953E+01
             1.2404E+00
 GRADIENT:  -2.8566E+00  9.8501E-01  7.3454E-01  7.2391E+00 -2.3953E+00 -4.7542E-01  8.7023E-02  0.0000E+00 -5.6276E-01  0.0000E+00
            -1.7260E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1415.95747004977        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1191
 NPARAMETR:  1.0312E+00  1.0146E-01  1.5844E+00  1.6330E+00  8.8764E-01  7.6276E-01  3.6824E+00  1.0000E-02  7.5759E-01  1.0000E-02
             3.1360E+00
 PARAMETER:  1.3073E-01 -2.1881E+00  5.6023E-01  5.9043E-01 -1.9186E-02 -1.7082E-01  1.4036E+00 -4.8397E+01 -1.7761E-01 -4.5298E+01
             1.2430E+00
 GRADIENT:  -1.4353E-01  4.6884E-01  5.0859E-01  6.2261E+00 -1.7834E+00  5.5337E-02  7.0978E-02  0.0000E+00 -3.5595E-01  0.0000E+00
             3.7113E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1415.99276120427        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1368
 NPARAMETR:  1.0300E+00  4.5748E-02  1.6088E+00  1.6656E+00  8.8161E-01  7.6211E-01  5.3866E+00  1.0000E-02  7.5008E-01  1.0000E-02
             3.1335E+00
 PARAMETER:  1.2954E-01 -2.9846E+00  5.7549E-01  6.1019E-01 -2.6002E-02 -1.7166E-01  1.7839E+00 -6.7934E+01 -1.8758E-01 -6.2619E+01
             1.2422E+00
 GRADIENT:  -5.3224E-01  1.2621E-01  2.5670E-01  2.7307E+00 -8.1306E-01 -5.0261E-02 -2.7641E-03  0.0000E+00 -1.5637E-01  0.0000E+00
            -1.7600E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1416.00414448375        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1545
 NPARAMETR:  1.0295E+00  1.7476E-02  1.6200E+00  1.6828E+00  8.7805E-01  7.6197E-01  8.3850E+00  1.0000E-02  7.4614E-01  1.0000E-02
             3.1338E+00
 PARAMETER:  1.2910E-01 -3.9469E+00  5.8240E-01  6.2048E-01 -3.0050E-02 -1.7184E-01  2.2264E+00 -9.1890E+01 -1.9285E-01 -8.3811E+01
             1.2423E+00
 GRADIENT:  -3.9794E-01  3.9996E-02  2.3861E-01  2.2767E+00 -7.3607E-01 -3.0577E-02 -1.7389E-02  0.0000E+00 -1.2110E-01  0.0000E+00
            -1.0944E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1416.01062468143        NO. OF FUNC. EVALS.: 164
 CUMULATIVE NO. OF FUNC. EVALS.:     1709
 NPARAMETR:  1.0292E+00  1.0000E-02  1.6208E+00  1.6848E+00  8.7716E-01  7.6182E-01  1.1626E+01  1.0000E-02  7.4519E-01  1.0000E-02
             3.1333E+00
 PARAMETER:  1.2878E-01 -4.5331E+00  5.8294E-01  6.2165E-01 -3.1067E-02 -1.7205E-01  2.5533E+00 -1.0569E+02 -1.9412E-01 -9.6005E+01
             1.2421E+00
 GRADIENT:   5.4521E+01  2.9859E-02  6.9620E-01  1.1313E+02  1.1356E+00  3.7160E+00  8.8150E-02  0.0000E+00  2.7275E+00  0.0000E+00
             9.8865E+00

0ITERATION NO.:   62    OBJECTIVE VALUE:  -1416.01062468143        NO. OF FUNC. EVALS.:  66
 CUMULATIVE NO. OF FUNC. EVALS.:     1775
 NPARAMETR:  1.0295E+00  1.0000E-02  1.6218E+00  1.6850E+00  8.7628E-01  7.6191E-01  1.1581E+01  1.0000E-02  7.4524E-01  1.0000E-02
             3.1333E+00
 PARAMETER:  1.2878E-01 -4.5331E+00  5.8294E-01  6.2165E-01 -3.1067E-02 -1.7205E-01  2.5533E+00 -1.0569E+02 -1.9412E-01 -9.6005E+01
             1.2421E+00
 GRADIENT:  -6.9937E-01  3.0711E-04 -8.4022E-02 -2.6640E-01  2.7724E-01 -2.4355E-02  8.0815E-04  0.0000E+00 -1.0183E-02  0.0000E+00
            -5.4055E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1775
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:        -3.2342E-04 -1.4688E-04  8.4070E-05 -1.0133E-02 -2.3464E-05
 SE:             2.8301E-02  1.7185E-03  8.1765E-05  2.6040E-02  1.5981E-04
 N:                     100         100         100         100         100

 P VAL.:         9.9088E-01  9.3189E-01  3.0386E-01  6.9719E-01  8.8327E-01

 ETASHRINKSD(%)  5.1885E+00  9.4243E+01  9.9726E+01  1.2762E+01  9.9465E+01
 ETASHRINKVR(%)  1.0108E+01  9.9669E+01  9.9999E+01  2.3895E+01  9.9997E+01
 EBVSHRINKSD(%)  5.1369E+00  9.4310E+01  9.9670E+01  1.2602E+01  9.9369E+01
 EBVSHRINKVR(%)  1.0010E+01  9.9676E+01  9.9999E+01  2.3615E+01  9.9996E+01
 RELATIVEINF(%)  7.0384E+01  1.9328E-03  4.0448E-05  5.4979E-01  2.1385E-04
 EPSSHRINKSD(%)  2.1995E+01
 EPSSHRINKVR(%)  3.9152E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1416.0106246814282     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -680.85979811769005     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    22.27
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.33
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1416.011       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.03E+00  1.00E-02  1.62E+00  1.68E+00  8.77E-01  7.62E-01  1.16E+01  1.00E-02  7.45E-01  1.00E-02  3.13E+00
 


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
+        1.63E+03
 
 TH 2
+       -1.89E+01  2.80E+02
 
 TH 3
+        2.04E+00  7.47E+00  4.25E+01
 
 TH 4
+       -1.08E+02  6.94E+01  2.57E+00  5.35E+02
 
 TH 5
+        1.31E+00 -4.27E+01 -1.45E+02 -1.39E+02  5.60E+02
 
 TH 6
+       -4.16E+00 -1.92E+00  4.68E+00 -1.97E+01 -7.09E+00  2.77E+02
 
 TH 7
+        1.15E-04  1.53E-01  2.73E-03 -7.15E-03  5.81E-03  2.79E-03  1.28E-03
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.10E+01 -3.31E+00  5.60E+00 -1.42E+01  1.71E+01 -1.92E+00  3.63E-02  0.00E+00  2.12E+02
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.41E+01 -1.01E+00 -1.71E-01 -8.86E+00  3.62E+00  7.34E+00  3.88E-03  0.00E+00  1.30E+01  0.00E+00  4.35E+01
 
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
 #CPUT: Total CPU Time in Seconds,       28.652
Stop Time:
Wed Sep 29 13:06:11 CDT 2021
