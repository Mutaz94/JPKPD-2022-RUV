Sat Sep 25 10:06:18 CDT 2021
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
$DATA ../../../../data/spa/S1/dat79.csv ignore=@
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
Current Date:       25 SEP 2021
Days until program expires : 204
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
 RAW OUTPUT FILE (FILE): m79.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1665.72270393682        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   8.6313E+00 -4.3630E+01 -1.3200E+00 -4.3290E+01 -2.1723E+01  2.7015E+01 -1.2907E+01 -1.1039E-01  1.5379E+01  9.7545E+00
            -6.4067E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1673.49612855292        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.9674E-01  1.2020E+00  1.0646E+00  9.3660E-01  1.1548E+00  8.8610E-01  1.1904E+00  1.0178E+00  8.4221E-01  8.9856E-01
             1.2586E+00
 PARAMETER:  9.6733E-02  2.8396E-01  1.6256E-01  3.4504E-02  2.4397E-01 -2.0926E-02  2.7433E-01  1.1769E-01 -7.1721E-02 -6.9667E-03
             3.3002E-01
 GRADIENT:  -3.3781E+01  3.1099E+01 -2.8860E+00  3.7030E+01  3.5372E+01 -2.0575E+01  4.3240E+00 -6.7709E+00 -4.3469E+00 -7.8702E+00
             2.9362E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1674.63338965705        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.0042E+00  1.0387E+00  1.2407E+00  1.0306E+00  1.0954E+00  9.0442E-01  1.4576E+00  1.5603E+00  7.3323E-01  7.7806E-01
             1.2343E+00
 PARAMETER:  1.0420E-01  1.3795E-01  3.1566E-01  1.3011E-01  1.9110E-01 -4.5752E-04  4.7677E-01  5.4491E-01 -2.1029E-01 -1.5095E-01
             3.1047E-01
 GRADIENT:  -3.6530E+00  2.3699E+01  3.4552E+00  1.3060E+01 -9.2040E+00 -1.0268E+01  1.6706E+01  3.6269E+00 -1.2742E+00 -7.0160E+00
             2.3998E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1677.79628570209        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      227
 NPARAMETR:  1.0030E+00  1.0762E+00  1.1760E+00  9.9490E-01  1.1159E+00  9.2751E-01  1.2386E+00  1.3318E+00  8.2682E-01  9.1681E-01
             1.1512E+00
 PARAMETER:  1.0303E-01  1.7346E-01  2.6213E-01  9.4888E-02  2.0965E-01  2.4746E-02  3.1402E-01  3.8652E-01 -9.0173E-02  1.3141E-02
             2.4082E-01
 GRADIENT:  -6.1404E-01  4.0886E+00  8.4741E-01  4.1558E+00  6.0226E-02 -5.7162E-01  4.5251E-01 -1.8953E-01  2.4390E-02 -7.2663E-01
             8.8184E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1678.06791476111        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      407
 NPARAMETR:  1.0190E+00  1.1560E+00  1.1196E+00  9.4821E-01  1.1281E+00  9.3913E-01  1.1665E+00  1.3340E+00  8.6161E-01  9.2640E-01
             1.1516E+00
 PARAMETER:  1.1881E-01  2.4496E-01  2.1301E-01  4.6822E-02  2.2052E-01  3.7197E-02  2.5403E-01  3.8818E-01 -4.8955E-02  2.3556E-02
             2.4116E-01
 GRADIENT:   3.4899E+00  3.2989E+00  8.0375E-01  4.4625E+00 -2.5292E+00  1.7327E+00 -8.2213E-01 -1.3465E-02 -1.9156E-01  8.6769E-02
            -8.2528E-03

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1678.26669870349        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      585
 NPARAMETR:  1.0182E+00  1.3577E+00  8.9297E-01  8.1277E-01  1.1321E+00  9.3313E-01  1.0333E+00  1.2013E+00  9.5813E-01  9.0947E-01
             1.1502E+00
 PARAMETER:  1.1800E-01  4.0579E-01 -1.3201E-02 -1.0730E-01  2.2411E-01  3.0791E-02  1.3276E-01  2.8337E-01  5.7225E-02  5.1100E-03
             2.3992E-01
 GRADIENT:  -3.0393E+00  3.3236E+00  9.2523E-01  3.0682E+00 -1.4795E+00 -1.5463E+00  2.1806E-01 -3.3971E-01 -1.7174E-01 -5.4892E-02
            -2.1698E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1678.86521826556        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      763
 NPARAMETR:  1.0216E+00  1.7518E+00  5.0383E-01  5.5579E-01  1.1583E+00  9.3853E-01  8.4530E-01  8.7662E-01  1.2350E+00  8.9306E-01
             1.1498E+00
 PARAMETER:  1.2140E-01  6.6064E-01 -5.8551E-01 -4.8737E-01  2.4698E-01  3.6557E-02 -6.8070E-02 -3.1687E-02  3.1105E-01 -1.3102E-02
             2.3960E-01
 GRADIENT:  -1.3290E-01  2.1311E+01  5.4047E-01  1.1003E+01 -6.2997E+00 -2.4546E-01 -1.4033E+00  2.3605E-01  9.1806E-01 -3.1703E-02
            -6.3432E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1679.40863081739        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      944
 NPARAMETR:  1.0218E+00  1.9664E+00  3.7141E-01  4.0311E-01  1.2602E+00  9.3909E-01  7.7708E-01  7.4091E-01  1.5115E+00  9.6817E-01
             1.1549E+00
 PARAMETER:  1.2155E-01  7.7623E-01 -8.9045E-01 -8.0854E-01  3.3124E-01  3.7153E-02 -1.5222E-01 -1.9987E-01  5.1311E-01  6.7650E-02
             2.4404E-01
 GRADIENT:   3.0310E-01  1.8998E+00 -2.8555E-01  2.0664E+00  6.2860E-01  2.1406E-01 -4.7055E-01  2.2859E-01  6.7303E-01  3.3101E-01
             7.5935E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1679.45901452590        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  1.0216E+00  2.0385E+00  3.1629E-01  3.5248E-01  1.2850E+00  9.3808E-01  7.6086E-01  6.4343E-01  1.6152E+00  9.7789E-01
             1.1570E+00
 PARAMETER:  1.2140E-01  8.1223E-01 -1.0511E+00 -9.4278E-01  3.5079E-01  3.6078E-02 -1.7331E-01 -3.4095E-01  5.7948E-01  7.7645E-02
             2.4585E-01
 GRADIENT:   9.4547E-03 -8.3317E-01 -3.9622E-01 -2.5594E-02  6.0503E-01 -4.0773E-02  3.4029E-01  2.5690E-01  7.8464E-02 -1.7164E-01
             2.8571E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1679.56868652176        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1299
 NPARAMETR:  1.0218E+00  2.0227E+00  3.1259E-01  3.6091E-01  1.2684E+00  9.3842E-01  7.6460E-01  1.1685E-01  1.5949E+00  9.6511E-01
             1.1559E+00
 PARAMETER:  1.2159E-01  8.0443E-01 -1.0629E+00 -9.1913E-01  3.3779E-01  3.6438E-02 -1.6840E-01 -2.0468E+00  5.6680E-01  6.4484E-02
             2.4487E-01
 GRADIENT:   3.7416E-01 -1.8728E+00 -5.0667E-03 -2.9479E-01  1.2541E+00 -3.8656E-02 -1.3188E-01  6.5933E-03 -1.8634E-01 -4.1225E-01
            -1.3574E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1679.57365761707        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1475
 NPARAMETR:  1.0217E+00  2.0316E+00  3.0752E-01  3.5577E-01  1.2711E+00  9.3847E-01  7.6276E-01  2.9819E-02  1.6091E+00  9.6977E-01
             1.1562E+00
 PARAMETER:  1.2145E-01  8.0882E-01 -1.0792E+00 -9.3348E-01  3.3988E-01  3.6491E-02 -1.7082E-01 -3.4126E+00  5.7570E-01  6.9302E-02
             2.4515E-01
 GRADIENT:   1.5960E-02 -7.1653E-02  1.6948E-02 -4.1874E-02 -3.3033E-02  1.9060E-03  5.9804E-03  4.7100E-04  9.4623E-03  2.8071E-03
             3.6528E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1679.57390957236        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     1657
 NPARAMETR:  1.0217E+00  2.0309E+00  3.0763E-01  3.5614E-01  1.2707E+00  9.3847E-01  7.6290E-01  1.0000E-02  1.6080E+00  9.6951E-01
             1.1562E+00
 PARAMETER:  1.2146E-01  8.0850E-01 -1.0789E+00 -9.3244E-01  3.3956E-01  3.6495E-02 -1.7063E-01 -4.8539E+00  5.7502E-01  6.9036E-02
             2.4513E-01
 GRADIENT:   3.2821E-02 -2.2599E-01 -2.1510E-02 -1.3354E-02  4.0955E-02  1.4721E-03  5.0767E-03  0.0000E+00  2.2778E-02  1.8314E-02
             1.3972E-02

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1679.57391349001        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1714
 NPARAMETR:  1.0217E+00  2.0310E+00  3.0764E-01  3.5614E-01  1.2707E+00  9.3847E-01  7.6289E-01  1.0000E-02  1.6080E+00  9.6941E-01
             1.1562E+00
 PARAMETER:  1.2145E-01  8.0852E-01 -1.0788E+00 -9.3243E-01  3.3955E-01  3.6493E-02 -1.7064E-01 -4.8539E+00  5.7497E-01  6.8936E-02
             2.4512E-01
 GRADIENT:   2.2690E-02 -1.2837E-01 -5.9936E-03 -1.5557E-02  1.9765E-02  8.2931E-04  4.9178E-04  0.0000E+00  1.0603E-02  3.9104E-03
             2.6093E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1714
 NO. OF SIG. DIGITS IN FINAL EST.:  3.7
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.2627E-04 -2.7208E-02 -2.0694E-04  3.4539E-02 -4.2008E-02
 SE:             2.9774E-02  2.6331E-02  7.1514E-05  2.0337E-02  2.0594E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9662E-01  3.0145E-01  3.8077E-03  8.9445E-02  4.1371E-02

 ETASHRINKSD(%)  2.5500E-01  1.1789E+01  9.9760E+01  3.1869E+01  3.1006E+01
 ETASHRINKVR(%)  5.0936E-01  2.2188E+01  9.9999E+01  5.3581E+01  5.2398E+01
 EBVSHRINKSD(%)  6.5217E-01  1.2094E+01  9.9794E+01  3.5018E+01  2.9057E+01
 EBVSHRINKVR(%)  1.3001E+00  2.2725E+01  1.0000E+02  5.7774E+01  4.9670E+01
 RELATIVEINF(%)  9.8641E+01  6.9793E+00  3.4796E-05  2.5859E+00  1.3228E+01
 EPSSHRINKSD(%)  4.3134E+01
 EPSSHRINKVR(%)  6.7662E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1679.5739134900064     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -944.42308692626818     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.00
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.14
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1679.574       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  2.03E+00  3.08E-01  3.56E-01  1.27E+00  9.38E-01  7.63E-01  1.00E-02  1.61E+00  9.69E-01  1.16E+00
 


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
+        1.19E+03
 
 TH 2
+       -9.18E+00  3.70E+02
 
 TH 3
+        8.91E+00  1.36E+02  5.80E+02
 
 TH 4
+       -2.58E+01  3.19E+02 -6.31E+02  1.51E+03
 
 TH 5
+       -6.73E+00 -1.15E+02 -3.21E+02  3.47E+02  3.78E+02
 
 TH 6
+        2.75E-01 -1.88E+00  1.61E+00 -8.10E+00 -1.85E+00  2.23E+02
 
 TH 7
+        8.86E-01  8.65E+00 -2.72E+01 -2.25E+01 -1.59E-02 -2.93E+00  2.17E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        1.28E+00 -1.18E+01 -4.53E+01  8.10E+01 -3.36E+00 -4.15E-01  9.40E+00  0.00E+00  2.59E+01
 
 TH10
+        4.57E-01 -1.18E+01 -2.90E+01  5.52E+00 -5.96E+01  1.54E+00  9.97E+00  0.00E+00  5.58E+00  6.62E+01
 
 TH11
+       -9.79E+00 -1.63E+01 -2.65E+01  1.42E+01 -7.86E+00  3.24E+00  1.30E+01  0.00E+00  4.40E+00  1.65E+01  1.62E+02
 
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
 #CPUT: Total CPU Time in Seconds,       27.212
Stop Time:
Sat Sep 25 10:06:47 CDT 2021
