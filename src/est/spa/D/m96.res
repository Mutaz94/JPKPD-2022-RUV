Sat Sep 18 15:45:21 CDT 2021
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
$DATA ../../../../data/spa/D/dat96.csv ignore=@
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
 (E4.0,E3.0,E22.0,E4.0,2E2.0)

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
 RAW OUTPUT FILE (FILE): m96.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   26676.2521866487        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   5.9367E+02  6.7320E+02  2.1951E+01  7.1168E+02  1.4010E+02 -3.4202E+03 -1.2276E+03 -1.3705E+02 -1.7139E+03 -7.5962E+02
            -4.9076E+04

0ITERATION NO.:    5    OBJECTIVE VALUE:  -395.570690862582        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:       83
 NPARAMETR:  9.8294E-01  7.8344E-01  7.8658E-01  1.1514E+00  1.9695E+00  1.8907E+00  9.8808E-01  9.3768E-01  7.3274E-01  8.6834E-01
             1.4885E+01
 PARAMETER:  8.2791E-02 -1.4407E-01 -1.4006E-01  2.4097E-01  7.7779E-01  7.3694E-01  8.8012E-02  3.5649E-02 -2.1096E-01 -4.1170E-02
             2.8003E+00
 GRADIENT:  -5.9272E+01  5.6196E-01 -2.7520E+00 -2.4725E+01 -8.5530E+00  3.1701E+01  3.2006E-01  2.6870E+00  5.4906E+00 -1.7400E-02
            -1.2103E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -424.728556421279        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      155
 NPARAMETR:  1.1175E+00  4.3875E-01  3.4542E-01  1.5491E+00  8.1011E+00  1.5722E+00  3.9532E-01  4.7091E-02  1.1934E-01  2.2385E+00
             1.7628E+01
 PARAMETER:  2.1107E-01 -7.2383E-01 -9.6298E-01  5.3764E-01  2.1920E+00  5.5247E-01 -8.2805E-01 -2.9557E+00 -2.0258E+00  9.0579E-01
             2.9695E+00
 GRADIENT:  -3.6863E+01  5.3698E+01 -4.2414E+01  1.3818E+02 -1.1060E+01 -9.4366E+01  9.2420E-01  4.4595E-02  2.4943E-01  1.3643E+00
            -2.3996E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -471.236430049034        NO. OF FUNC. EVALS.:  76
 CUMULATIVE NO. OF FUNC. EVALS.:      231
 NPARAMETR:  9.5666E-01  1.4657E-01  1.3425E-01  1.2527E+00  2.7503E+01  2.0962E+00  5.5889E-02  1.0000E-02  1.3922E-02  4.3039E+00
             1.9371E+01
 PARAMETER:  5.5689E-02 -1.8202E+00 -1.9080E+00  3.2529E-01  3.4143E+00  8.4015E-01 -2.7844E+00 -6.3639E+00 -4.1743E+00  1.5595E+00
             3.0638E+00
 GRADIENT:  -1.6099E+01  2.0658E+01 -5.8837E+01  1.2576E+02 -1.6796E+00 -3.1987E+01  7.9738E-02  0.0000E+00  9.7618E-03  9.2852E-01
             5.2802E+01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -521.442180965406        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      302
 NPARAMETR:  4.9183E-01  1.4820E-02  2.4643E-02  3.0023E-01  1.2872E+03  1.8655E+00  1.0000E-02  1.0000E-02  5.3706E-02  7.5387E+01
             1.3825E+01
 PARAMETER: -6.0962E-01 -4.1118E+00 -3.6033E+00 -1.1032E+00  7.2603E+00  7.2353E-01 -8.3120E+00 -5.3459E+00 -2.8242E+00  4.4226E+00
             2.7265E+00
 GRADIENT:   3.1781E+01 -1.3925E+00 -5.7238E+01  7.9632E+01  3.2465E-03 -3.5632E+01  0.0000E+00  0.0000E+00  1.3492E-01 -5.0630E-05
            -5.5946E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -525.909811190205        NO. OF FUNC. EVALS.:  78
 CUMULATIVE NO. OF FUNC. EVALS.:      380
 NPARAMETR:  3.6740E-01  1.0000E-02  1.2458E-02  1.7257E-01  5.0131E+03  1.8879E+00  1.0000E-02  1.0000E-02  1.4132E-02  2.0730E+02
             1.3895E+01
 PARAMETER: -9.0130E-01 -5.0322E+00 -4.2854E+00 -1.6569E+00  8.6198E+00  7.3546E-01 -1.0415E+01 -5.7643E+00 -4.1593E+00  5.4342E+00
             2.7315E+00
 GRADIENT:   1.9379E+01  0.0000E+00 -4.5873E+01  5.7776E+01 -2.9767E-03 -2.5862E+01  0.0000E+00  0.0000E+00  1.0641E-02  5.5822E-03
            -3.3575E+01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -528.525579393916        NO. OF FUNC. EVALS.: 163
 CUMULATIVE NO. OF FUNC. EVALS.:      543
 NPARAMETR:  3.7981E-01  1.0000E-02  1.4001E-02  1.8455E-01  3.7860E+03  2.0224E+00  1.0000E-02  1.0000E-02  1.0000E-02  1.7133E+02
             1.4526E+01
 PARAMETER: -8.6807E-01 -4.9184E+00 -4.1686E+00 -1.5898E+00  8.3391E+00  8.0427E-01 -1.0158E+01 -6.2866E+00 -4.6225E+00  5.2436E+00
             2.7759E+00
 GRADIENT:  -3.5246E-01  0.0000E+00  1.5421E-01  4.8667E-02 -3.1862E-04 -1.6504E-01  0.0000E+00  0.0000E+00  0.0000E+00  2.6559E-04
             1.5775E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -528.526063503493        NO. OF FUNC. EVALS.: 162
 CUMULATIVE NO. OF FUNC. EVALS.:      705
 NPARAMETR:  3.7965E-01  1.0000E-02  1.3959E-02  1.8411E-01  8.2429E+03  2.0233E+00  1.0000E-02  1.0000E-02  1.0000E-02  7.8724E+01
             1.4523E+01
 PARAMETER: -8.6850E-01 -4.9224E+00 -4.1716E+00 -1.5922E+00  9.1171E+00  8.0473E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  4.4659E+00
             2.7757E+00
 GRADIENT:   6.5812E-02  0.0000E+00 -2.7992E-02  1.6086E-02 -4.8324E-05 -1.5123E-02  0.0000E+00  0.0000E+00  0.0000E+00  6.2315E-06
            -6.9382E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -528.526093626631        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:      881
 NPARAMETR:  3.7966E-01  1.0000E-02  1.3964E-02  1.8416E-01  5.6149E+04  2.0234E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3512E+01
             1.4524E+01
 PARAMETER: -8.6847E-01 -4.9224E+00 -4.1713E+00 -1.5920E+00  1.1036E+01  8.0477E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  3.6119E+00
             2.7758E+00
 GRADIENT:  -1.9465E-03  0.0000E+00  9.7578E-03  9.1276E-03 -6.2034E-06  7.6927E-04  0.0000E+00  0.0000E+00  0.0000E+00  2.2505E-08
             2.1365E-03

0ITERATION NO.:   45    OBJECTIVE VALUE:  -528.526100198684        NO. OF FUNC. EVALS.: 190
 CUMULATIVE NO. OF FUNC. EVALS.:     1071
 NPARAMETR:  3.7964E-01  1.0000E-02  1.3963E-02  1.8414E-01  1.5859E+05  2.0233E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3750E+01
             1.4524E+01
 PARAMETER: -8.6854E-01 -4.9224E+00 -4.1713E+00 -1.5920E+00  1.2074E+01  8.0475E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  3.6190E+00
             2.7758E+00
 GRADIENT:  -1.9303E-02  0.0000E+00  4.4239E-02 -2.5954E-02 -2.2038E-06 -2.5850E-03  0.0000E+00  0.0000E+00  0.0000E+00  2.8587E-09
             8.4778E-03

0ITERATION NO.:   50    OBJECTIVE VALUE:  -528.526139038476        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:     1252
 NPARAMETR:  3.7947E-01  1.0000E-02  1.3948E-02  1.8399E-01  2.3231E+10  2.0233E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3676E+01
             1.4524E+01
 PARAMETER: -8.6897E-01 -4.9224E+00 -4.1724E+00 -1.5929E+00  2.3969E+01  8.0475E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  3.6168E+00
             2.7758E+00
 GRADIENT:  -2.1623E-02  0.0000E+00 -4.0985E-02  6.0418E-02 -1.6601E-11 -5.0565E-03  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             1.4674E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -528.526143168821        NO. OF FUNC. EVALS.: 196
 CUMULATIVE NO. OF FUNC. EVALS.:     1448             RESET HESSIAN, TYPE I
 NPARAMETR:  3.7947E-01  1.0000E-02  1.3946E-02  1.8397E-01  1.0256E+11  2.0233E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3444E+01
             1.4524E+01
 PARAMETER: -8.6898E-01 -4.9224E+00 -4.1726E+00 -1.5930E+00  2.5454E+01  8.0475E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  3.6099E+00
             2.7758E+00
 GRADIENT:   4.3146E+00  0.0000E+00  5.6364E+00  2.4083E+00 -4.4664E-12  1.1236E+00  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
             2.2096E+00

0ITERATION NO.:   57    OBJECTIVE VALUE:  -528.526143168821        NO. OF FUNC. EVALS.:  57
 CUMULATIVE NO. OF FUNC. EVALS.:     1505
 NPARAMETR:  3.7947E-01  1.0000E-02  1.3946E-02  1.8397E-01  1.0256E+11  2.0233E+00  1.0000E-02  1.0000E-02  1.0000E-02  3.3444E+01
             1.4524E+01
 PARAMETER: -8.6898E-01 -4.9224E+00 -4.1726E+00 -1.5930E+00  2.5454E+01  8.0475E-01 -1.0166E+01 -6.2933E+00 -4.6267E+00  3.6099E+00
             2.7758E+00
 GRADIENT:   9.3780E-03  0.0000E+00 -6.5174E-02  7.2528E-02 -4.4664E-12 -5.6662E-03  0.0000E+00  0.0000E+00  0.0000E+00  0.0000E+00
            -1.1030E-02

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     1505
 NO. OF SIG. DIGITS IN FINAL EST.:  4.6
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.4548E-03  6.9543E-06  7.9663E-05 -1.6537E-04 -1.5678E-14
 SE:             2.9133E-02  3.6541E-06  2.5720E-04  2.8762E-04  2.9915E-13
 N:                     100         100         100         100         100

 P VAL.:         9.6017E-01  5.7019E-02  7.5677E-01  5.6533E-01  9.5820E-01

 ETASHRINKSD(%)  2.4005E+00  9.9988E+01  9.9138E+01  9.9036E+01  1.0000E+02
 ETASHRINKVR(%)  4.7435E+00  1.0000E+02  9.9993E+01  9.9991E+01  1.0000E+02
 EBVSHRINKSD(%)  2.6751E+00  9.9977E+01  9.9080E+01  9.8921E+01  1.0000E+02
 EBVSHRINKVR(%)  5.2785E+00  1.0000E+02  9.9992E+01  9.9988E+01  1.0000E+02
 RELATIVEINF(%)  8.0830E+00  2.9475E-06  4.0987E-05  5.5714E-05  0.0000E+00
 EPSSHRINKSD(%)  5.3537E+00
 EPSSHRINKVR(%)  1.0421E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -528.52614316882057     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:       206.62468339491761     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    19.59
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.64
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************     -528.526       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         3.79E-01  1.00E-02  1.39E-02  1.84E-01  1.03E+11  2.02E+00  1.00E-02  1.00E-02  1.00E-02  3.34E+01  1.45E+01
 


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
+        1.77E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.40E+04  0.00E+00  2.89E+06
 
 TH 4
+       -4.83E+02  0.00E+00 -2.52E+05  2.40E+04
 
 TH 5
+        1.79E-15  0.00E+00  2.40E-14  7.04E-15  1.76E-28
 
 TH 6
+       -4.07E+00  0.00E+00  1.38E+03 -1.31E+02 -4.86E-17  3.97E+01
 
 TH 7
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 8
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH 9
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH10
+        5.91E-05  0.00E+00 -2.97E-04 -7.37E-05 -1.23E-17 -7.73E-06  0.00E+00  0.00E+00  0.00E+00 -2.29E-07
 
 TH11
+       -1.53E+01  0.00E+00  3.61E+02 -2.02E+01  5.09E-17  2.01E-01  0.00E+00  0.00E+00  0.00E+00 -1.92E-08  1.53E+00
 
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
 #CPUT: Total CPU Time in Seconds,       26.301
Stop Time:
Sat Sep 18 15:45:49 CDT 2021
