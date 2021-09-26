Sat Sep 25 00:52:14 CDT 2021
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
$DATA ../../../../data/int/SL2/dat9.csv ignore=@
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
 (2E4.0,E19.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m9.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1404.92743898956        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.9414E+02 -6.0165E+01  1.3583E+02  7.9668E+01  1.4750E+02 -4.5362E+00 -8.2822E+01 -2.1359E+02 -9.6011E+01 -1.6979E+01
            -4.4252E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -2842.28960092685        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:       85
 NPARAMETR:  9.3456E-01  1.3819E+00  8.3726E-01  7.8674E-01  1.1382E+00  9.0958E-01  1.0256E+00  1.0196E+00  1.1003E+00  1.0350E+00
             2.1703E+00
 PARAMETER:  3.2321E-02  4.2346E-01 -7.7615E-02 -1.3986E-01  2.2946E-01  5.2278E-03  1.2527E-01  1.1942E-01  1.9559E-01  1.3439E-01
             8.7487E-01
 GRADIENT:  -1.4716E+01  1.3759E+00 -2.2188E+01 -8.5407E+00  8.6200E+00 -2.5597E+01  1.9407E+01 -1.5961E+00 -1.3081E+01 -2.7174E+01
            -1.6855E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -2851.68246547476        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  9.3988E-01  1.6583E+00  8.8464E-01  6.7575E-01  1.3008E+00  9.0965E-01  7.6326E-01  1.0378E+00  1.4044E+00  1.2539E+00
             2.1895E+00
 PARAMETER:  3.7998E-02  6.0582E-01 -2.2571E-02 -2.9193E-01  3.6296E-01  5.3034E-03 -1.7016E-01  1.3707E-01  4.3963E-01  3.2628E-01
             8.8366E-01
 GRADIENT:  -1.7746E+00  1.0165E+02  4.9294E+00  4.6452E+01 -2.6331E+01 -2.5613E+01  2.0534E+00 -5.9451E+00  1.9662E+00 -1.4726E+01
            -1.4075E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -2861.35339273489        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.3983E-01  1.7052E+00  9.6816E-01  5.9900E-01  1.4156E+00  9.7390E-01  6.4053E-01  1.6386E+00  1.4913E+00  1.4184E+00
             2.3187E+00
 PARAMETER:  3.7943E-02  6.3368E-01  6.7642E-02 -4.1250E-01  4.4759E-01  7.3554E-02 -3.4546E-01  5.9385E-01  4.9965E-01  4.4953E-01
             9.4099E-01
 GRADIENT:  -3.8355E+00  1.3620E+01  2.1381E+00  7.1943E+00 -3.6772E-01  2.6948E+00 -1.3370E+00 -2.2335E-01 -2.7606E+00  3.5146E+00
             1.5441E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -2861.86537892340        NO. OF FUNC. EVALS.:  99
 CUMULATIVE NO. OF FUNC. EVALS.:      328
 NPARAMETR:  9.4168E-01  1.7368E+00  9.2741E-01  5.7551E-01  1.4296E+00  9.7290E-01  6.3021E-01  1.6191E+00  1.5547E+00  1.4175E+00
             2.3089E+00
 PARAMETER:  3.9909E-02  6.5203E-01  2.4640E-02 -4.5250E-01  4.5739E-01  7.2525E-02 -3.6170E-01  5.8186E-01  5.4125E-01  4.4892E-01
             9.3679E-01
 GRADIENT:  -6.1818E+00 -6.3492E+00  1.2192E+00  3.2421E+00 -4.5270E+00  1.6225E+00 -1.2165E+00 -3.1111E-01 -1.5482E+00  3.4785E-01
             4.5690E+00

0ITERATION NO.:   25    OBJECTIVE VALUE:  -2863.52246861755        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      504
 NPARAMETR:  9.4341E-01  2.0862E+00  6.4328E-01  3.5144E-01  1.6512E+00  9.7101E-01  6.7432E-01  1.5004E+00  2.0344E+00  1.5206E+00
             2.3103E+00
 PARAMETER:  4.1744E-02  8.3534E-01 -3.4117E-01 -9.4571E-01  6.0151E-01  7.0581E-02 -2.9405E-01  5.0572E-01  8.1022E-01  5.1913E-01
             9.3736E-01
 GRADIENT:  -2.4079E+00  1.2167E+01  1.8820E+00  1.8198E+00 -5.0021E+00  4.8423E-01  2.7265E+00  4.2569E-01  9.7591E-01 -7.9603E-01
             3.3717E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -2864.43784376252        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      680
 NPARAMETR:  9.4545E-01  2.2605E+00  3.4556E-01  2.3529E-01  1.7566E+00  9.6945E-01  6.6886E-01  7.4552E-01  2.4500E+00  1.5885E+00
             2.2994E+00
 PARAMETER:  4.3903E-02  9.1557E-01 -9.6259E-01 -1.3469E+00  6.6341E-01  6.8976E-02 -3.0218E-01 -1.9367E-01  9.9610E-01  5.6276E-01
             9.3266E-01
 GRADIENT:   2.5640E+00  8.0459E+00 -3.7370E+00  4.5946E+00  3.0692E+00 -7.0616E-01 -6.0389E-01  9.0922E-01 -2.7199E+00 -1.6562E+00
            -4.8517E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -2865.09851379439        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      855
 NPARAMETR:  9.4406E-01  2.3314E+00  3.0372E-01  1.8062E-01  1.8059E+00  9.7127E-01  6.6478E-01  3.8785E-01  2.9138E+00  1.6220E+00
             2.3018E+00
 PARAMETER:  4.2431E-02  9.4646E-01 -1.0916E+00 -1.6114E+00  6.9106E-01  7.0845E-02 -3.0830E-01 -8.4713E-01  1.1695E+00  5.8363E-01
             9.3368E-01
 GRADIENT:  -3.8474E-01  2.7741E+00 -3.1440E-01 -5.3156E-01 -2.4181E+00 -3.6365E-02  1.8355E+00  3.3409E-01  4.5730E-01 -3.4002E-01
            -5.6023E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -2865.27566512697        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1032
 NPARAMETR:  9.4371E-01  2.3160E+00  3.1261E-01  1.9111E-01  1.7968E+00  9.7196E-01  6.6178E-01  1.0387E-01  2.8412E+00  1.6172E+00
             2.3021E+00
 PARAMETER:  4.2064E-02  9.3983E-01 -1.0628E+00 -1.5549E+00  6.8602E-01  7.1555E-02 -3.1282E-01 -2.1646E+00  1.1442E+00  5.8072E-01
             9.3382E-01
 GRADIENT:  -1.2635E+00  2.2353E+00  9.5140E-02  3.8379E-01 -6.9628E-01  2.2346E-01  2.8406E-01  2.1047E-02  2.6559E-01  7.6269E-02
             6.5697E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -2865.28684187270        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1216            RESET HESSIAN, TYPE II
 NPARAMETR:  9.4422E-01  2.3176E+00  3.1008E-01  1.8893E-01  1.8005E+00  9.7138E-01  6.6084E-01  4.6104E-02  2.8501E+00  1.6184E+00
             2.3014E+00
 PARAMETER:  4.2607E-02  9.4055E-01 -1.0709E+00 -1.5664E+00  6.8804E-01  7.0962E-02 -3.1425E-01 -2.9769E+00  1.1474E+00  5.8141E-01
             9.3353E-01
 GRADIENT:   7.2996E+00  4.6379E+01  1.0703E-01  1.2419E+00  4.2727E+00  6.3160E-01  5.7769E-01  4.2913E-03  8.9205E-01  9.2333E-01
             1.6835E+00

0ITERATION NO.:   50    OBJECTIVE VALUE:  -2865.28691077506        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:     1394
 NPARAMETR:  9.4422E-01  2.3177E+00  3.0997E-01  1.8901E-01  1.8005E+00  9.7134E-01  6.6086E-01  4.5576E-02  2.8509E+00  1.6183E+00
             2.3014E+00
 PARAMETER:  4.2604E-02  9.4056E-01 -1.0713E+00 -1.5660E+00  6.8806E-01  7.0921E-02 -3.1422E-01 -2.9884E+00  1.1476E+00  5.8138E-01
             9.3353E-01
 GRADIENT:  -6.2623E-03  5.3389E-02 -4.5448E-03  1.3773E-02  6.4627E-03 -1.0340E-03  3.6884E-03  4.1691E-03  1.1173E-02 -2.3504E-03
             8.5826E-04

0ITERATION NO.:   55    OBJECTIVE VALUE:  -2865.28883222280        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1574
 NPARAMETR:  9.4422E-01  2.3170E+00  3.1040E-01  1.8948E-01  1.7998E+00  9.7134E-01  6.6092E-01  1.1041E-02  2.8470E+00  1.6180E+00
             2.3014E+00
 PARAMETER:  4.2605E-02  9.4026E-01 -1.0699E+00 -1.5635E+00  6.8770E-01  7.0924E-02 -3.1412E-01 -4.4061E+00  1.1463E+00  5.8118E-01
             9.3352E-01
 GRADIENT:  -9.2212E-04  1.1251E-02  5.9207E-03  3.0803E-02 -1.9008E-02  6.1031E-04 -1.9950E-03  2.4317E-04 -5.5122E-03  6.6562E-04
            -1.2531E-02

0ITERATION NO.:   60    OBJECTIVE VALUE:  -2865.28889454487        NO. OF FUNC. EVALS.: 169
 CUMULATIVE NO. OF FUNC. EVALS.:     1743
 NPARAMETR:  9.4420E-01  2.3177E+00  3.0982E-01  1.8926E-01  1.8005E+00  9.7136E-01  6.6084E-01  1.0000E-02  2.8538E+00  1.6184E+00
             2.3013E+00
 PARAMETER:  4.2608E-02  9.4053E-01 -1.0714E+00 -1.5660E+00  6.8802E-01  7.0922E-02 -3.1421E-01 -4.5233E+00  1.1475E+00  5.8136E-01
             9.3353E-01
 GRADIENT:   2.5489E-03 -9.5813E-03  1.0750E-03 -7.3803E-03 -2.3714E-03 -3.2142E-04  3.5408E-04  0.0000E+00 -6.7909E-03 -1.0647E-03
             4.4265E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1743
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5485E-03 -3.5505E-02 -1.8292E-05  4.4201E-02 -2.3715E-02
 SE:             2.9486E-02  2.4594E-02  3.3806E-05  1.9366E-02  2.6725E-02
 N:                     100         100         100         100         100

 P VAL.:         9.5812E-01  1.4884E-01  5.8845E-01  2.2469E-02  3.7488E-01

 ETASHRINKSD(%)  1.2177E+00  1.7607E+01  9.9887E+01  3.5120E+01  1.0468E+01
 ETASHRINKVR(%)  2.4207E+00  3.2114E+01  1.0000E+02  5.7906E+01  1.9840E+01
 EBVSHRINKSD(%)  1.4426E+00  1.4799E+01  9.9861E+01  4.2639E+01  8.0635E+00
 EBVSHRINKVR(%)  2.8643E+00  2.7408E+01  1.0000E+02  6.7097E+01  1.5477E+01
 RELATIVEINF(%)  9.7098E+01  2.0812E+01  1.2281E-04  8.9729E+00  5.3376E+01
 EPSSHRINKSD(%)  1.7030E+01
 EPSSHRINKVR(%)  3.1159E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          900
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    1654.0893597684108     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -2865.2888945448672     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -1211.1995347764564     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    45.62
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:    13.65
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -2865.289       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.44E-01  2.32E+00  3.10E-01  1.89E-01  1.80E+00  9.71E-01  6.61E-01  1.00E-02  2.85E+00  1.62E+00  2.30E+00
 


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
+        1.28E+03
 
 TH 2
+       -1.12E+01  3.58E+02
 
 TH 3
+        2.46E+00  5.68E+01  2.59E+02
 
 TH 4
+       -1.38E+01  3.85E+02 -2.15E+02  1.43E+03
 
 TH 5
+       -2.89E+00 -3.34E+01 -1.94E+01  9.92E+01  1.11E+02
 
 TH 6
+        9.24E+00 -2.47E+00  2.89E+00 -9.93E+00 -7.32E-01  1.93E+02
 
 TH 7
+        2.97E+00 -2.17E+01 -1.12E+01  1.99E+01 -7.43E+00 -6.54E-01  2.51E+02
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        8.06E-01 -3.99E+00 -5.67E+00  4.96E+01  5.51E-03 -2.67E-01  1.07E+01  0.00E+00  6.79E+00
 
 TH10
+       -2.80E-01 -3.00E+00  1.10E+01  9.23E+00 -6.23E+00  6.38E-02  9.51E-01  0.00E+00 -2.38E-01  5.27E+01
 
 TH11
+       -1.62E+01 -1.78E+01 -5.62E+00 -2.53E+00 -4.66E-02  4.77E+00  4.46E+00  0.00E+00  2.13E+00  5.60E+00  2.19E+02
 
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
 #CPUT: Total CPU Time in Seconds,       59.362
Stop Time:
Sat Sep 25 00:53:14 CDT 2021
