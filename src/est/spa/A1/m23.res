Wed Sep 29 12:00:08 CDT 2021
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
$DATA ../../../../data/spa/A1/dat23.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m23.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1475.16674097378        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   3.8970E+02  5.0724E+01  2.7938E+01  8.6204E+01  1.6747E+01  5.6474E+01  4.3619E+00 -2.6187E+01  1.1936E+01  8.6122E+00
            -3.9293E+02

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1561.21528029229        NO. OF FUNC. EVALS.:  95
 CUMULATIVE NO. OF FUNC. EVALS.:      108
 NPARAMETR:  1.0263E+00  1.0296E+00  9.1609E-01  9.6784E-01  1.0151E+00  8.6568E-01  9.3719E-01  1.1171E+00  8.5586E-01  8.0626E-01
             1.9395E+00
 PARAMETER:  1.2594E-01  1.2913E-01  1.2356E-02  6.7313E-02  1.1503E-01 -4.4237E-02  3.5135E-02  2.1074E-01 -5.5650E-02 -1.1534E-01
             7.6242E-01
 GRADIENT:  -1.6331E+01 -4.1918E+01 -2.6263E+01 -1.7714E+01  4.8522E+01 -3.2775E+01 -4.8622E+00 -3.1718E+00 -2.2091E+00  9.1818E+00
             5.0380E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1567.00600795993        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      283
 NPARAMETR:  1.0500E+00  1.1386E+00  9.7737E-01  9.2353E-01  9.8084E-01  9.6130E-01  9.9080E-01  1.6040E+00  9.1052E-01  4.1157E-01
             1.8975E+00
 PARAMETER:  1.4882E-01  2.2981E-01  7.7111E-02  2.0452E-02  8.0657E-02  6.0535E-02  9.0754E-02  5.7250E-01  6.2560E-03 -7.8778E-01
             7.4053E-01
 GRADIENT:   4.1052E+01  3.3657E+01  2.3885E+01  6.6800E+00 -5.6329E+01  6.6253E+00  3.8967E+00 -2.1407E+00  4.8423E+00  7.8616E-01
             3.7967E+01

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1570.57251913233        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      461
 NPARAMETR:  1.0337E+00  1.1892E+00  8.4688E-01  8.7918E-01  9.8302E-01  9.4291E-01  9.2212E-01  1.5284E+00  9.2612E-01  5.2868E-01
             1.7227E+00
 PARAMETER:  1.3313E-01  2.7329E-01 -6.6200E-02 -2.8761E-02  8.2878E-02  4.1218E-02  1.8918E-02  5.2423E-01  2.3243E-02 -5.3737E-01
             6.4389E-01
 GRADIENT:   5.5213E+00  8.1104E+00  3.0313E+00  7.5894E+00 -5.5315E+00 -1.5888E+00 -1.0437E+00 -2.7413E-01 -1.1708E+00  1.5628E-01
            -6.0282E-02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1573.42344884840        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      641
 NPARAMETR:  1.0111E+00  1.4187E+00  2.1895E-01  6.2389E-01  6.7517E-01  9.8262E-01  8.8247E-01  6.1243E-01  9.9558E-01  2.4752E-01
             1.6834E+00
 PARAMETER:  1.1106E-01  4.4973E-01 -1.4189E+00 -3.7178E-01 -2.9279E-01  8.2468E-02 -2.5032E-02 -3.9032E-01  9.5568E-02 -1.2963E+00
             6.2081E-01
 GRADIENT:  -2.6482E+01  6.1136E+01  3.3661E+01 -3.2285E+01 -9.9411E+01  1.6405E+01  1.5860E+01 -1.3989E+00  1.3417E+01  2.5435E+00
             1.0016E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1579.69014729804        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:      821
 NPARAMETR:  1.0260E+00  1.5603E+00  1.6654E-01  5.3624E-01  7.3991E-01  8.9965E-01  7.6045E-01  6.5298E-01  9.3800E-01  1.4086E-01
             1.6771E+00
 PARAMETER:  1.2567E-01  5.4490E-01 -1.6925E+00 -5.2317E-01 -2.0123E-01 -5.7455E-03 -1.7384E-01 -3.2620E-01  3.5993E-02 -1.8600E+00
             6.1704E-01
 GRADIENT:   1.4476E+01  1.7196E+01 -4.1662E-01  2.6470E+01  1.1721E+01 -1.4826E+01 -5.1514E+00 -4.5240E+00 -1.2535E+01  3.1148E-01
            -7.8451E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1599.53415549162        NO. OF FUNC. EVALS.: 184
 CUMULATIVE NO. OF FUNC. EVALS.:     1005             RESET HESSIAN, TYPE I
 NPARAMETR:  1.0176E+00  2.0112E+00  1.0738E-01  2.7854E-01  9.9601E-01  9.1159E-01  6.7927E-01  1.6348E+00  1.4704E+00  1.1223E-02
             1.7378E+00
 PARAMETER:  1.1743E-01  7.9873E-01 -2.1314E+00 -1.1782E+00  9.6001E-02  7.4347E-03 -2.8673E-01  5.9155E-01  4.8555E-01 -4.3898E+00
             6.5263E-01
 GRADIENT:   1.5852E+02  4.2746E+02  2.8487E+01  1.2211E+01 -5.1049E+01  5.5693E+00  4.7662E+00 -1.5893E+01 -6.6845E+00  2.9776E-04
             5.3808E+01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1600.88974436633        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1182
 NPARAMETR:  1.0254E+00  2.0099E+00  1.0733E-01  2.8414E-01  1.0191E+00  9.2482E-01  6.7643E-01  1.6348E+00  1.5845E+00  1.6332E-02
             1.7375E+00
 PARAMETER:  1.2507E-01  7.9811E-01 -2.1318E+00 -1.1583E+00  1.1893E-01  2.1848E-02 -2.9093E-01  5.9153E-01  5.6028E-01 -4.0146E+00
             6.5244E-01
 GRADIENT:   4.0587E+00  2.0714E+01  1.1296E+01  7.1257E-02  6.7929E-01 -1.3500E+00  2.7032E-01 -1.6807E+01 -4.4289E-01 -2.2221E-03
             5.0582E+01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1601.39030366973        NO. OF FUNC. EVALS.: 180
 CUMULATIVE NO. OF FUNC. EVALS.:     1362
 NPARAMETR:  1.0251E+00  2.0009E+00  1.0655E-01  2.8600E-01  1.0161E+00  9.2179E-01  6.7065E-01  1.6356E+00  1.5641E+00  3.0859E-01
             1.7316E+00
 PARAMETER:  1.2477E-01  7.9360E-01 -2.1391E+00 -1.1518E+00  1.1596E-01  1.8559E-02 -2.9951E-01  5.9198E-01  5.4730E-01 -1.0757E+00
             6.4907E-01
 GRADIENT:   3.3129E+00  3.6258E+00  9.9751E+00  2.4384E-02  1.2969E+00 -2.4708E+00 -6.8120E-02 -1.6016E+01 -5.7891E-02  2.1286E-01
             6.0677E+01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1603.26168129729        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1538
 NPARAMETR:  1.0205E+00  1.9220E+00  9.9258E-02  3.0964E-01  9.5253E-01  9.2633E-01  6.8508E-01  1.6434E+00  1.4406E+00  2.2484E-01
             1.6748E+00
 PARAMETER:  1.2029E-01  7.5337E-01 -2.2100E+00 -1.0723E+00  5.1372E-02  2.3477E-02 -2.7822E-01  5.9675E-01  4.6506E-01 -1.3924E+00
             6.1568E-01
 GRADIENT:   2.2524E+00 -3.4659E+01 -8.8907E-01  8.1286E-02  4.2253E-01  3.1016E-01  4.0666E-01 -1.0099E+01 -2.8574E-01 -7.7350E-02
             4.6558E+01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1604.65660824064        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     1717
 NPARAMETR:  1.0165E+00  1.9180E+00  9.2639E-02  3.0947E-01  9.3450E-01  9.2970E-01  6.6748E-01  1.6559E+00  1.4540E+00  4.3667E-01
             1.5918E+00
 PARAMETER:  1.1636E-01  7.5129E-01 -2.2790E+00 -1.0729E+00  3.2255E-02  2.7110E-02 -3.0424E-01  6.0436E-01  4.7430E-01 -7.2858E-01
             5.6487E-01
 GRADIENT:  -3.4492E+00 -3.0961E+01 -5.8600E+00  1.3152E+00 -2.9363E+01  1.7272E+00 -5.5132E+00 -6.2659E+00  1.5625E+00  2.3203E+00
             4.2630E+01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1604.89676803133        NO. OF FUNC. EVALS.: 165
 CUMULATIVE NO. OF FUNC. EVALS.:     1882
 NPARAMETR:  1.0171E+00  1.9125E+00  9.3538E-02  3.1031E-01  9.3497E-01  9.2324E-01  6.8521E-01  1.6515E+00  1.4619E+00  3.3543E-01
             1.5928E+00
 PARAMETER:  1.1691E-01  7.4841E-01 -2.2694E+00 -1.0702E+00  3.2754E-02  2.0133E-02 -2.7804E-01  6.0166E-01  4.7976E-01 -9.9235E-01
             5.6551E-01
 GRADIENT:  -1.4670E+00 -3.4307E+01 -6.2544E+00  1.2473E+00 -1.8453E+01 -9.6850E-01  2.4891E-01 -6.8504E+00  1.4747E+00 -6.0611E-02
             3.5947E+01

0ITERATION NO.:   57    OBJECTIVE VALUE:  -1604.91359163824        NO. OF FUNC. EVALS.:  88
 CUMULATIVE NO. OF FUNC. EVALS.:     1970
 NPARAMETR:  1.0171E+00  1.9131E+00  9.3784E-02  3.1028E-01  9.3534E-01  9.3123E-01  6.8456E-01  1.6615E+00  1.4612E+00  3.4214E-01
             1.5859E+00
 PARAMETER:  1.1691E-01  7.4880E-01 -2.2686E+00 -1.0702E+00  3.2760E-02  2.9019E-02 -2.7898E-01  6.0172E-01  4.7974E-01 -9.7244E-01
             5.6535E-01
 GRADIENT:  -1.6277E+00  1.1490E+03 -5.2688E+00  1.3756E+00 -1.9170E+01  2.2821E+00  6.1774E-02 -6.7882E+00  1.3457E+00  2.5761E-02
             3.5977E+01
 NUMSIGDIG:         1.9         2.3         1.4         2.6         0.8         0.9         3.1         0.3         1.3         2.3
                    0.5

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1970
 NO. OF SIG. DIGITS IN FINAL EST.:  0.3

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.0523E-03 -8.0184E-03 -1.8438E-03  2.1573E-02 -1.6973E-02
 SE:             2.9446E-02  2.7068E-02  1.2469E-02  2.1430E-02  9.8925E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7149E-01  7.6705E-01  8.8244E-01  3.1409E-01  8.6213E-02

 ETASHRINKSD(%)  1.3520E+00  9.3199E+00  5.8226E+01  2.8206E+01  6.6859E+01
 ETASHRINKVR(%)  2.6857E+00  1.7771E+01  8.2550E+01  4.8456E+01  8.9017E+01
 EBVSHRINKSD(%)  1.0940E+00  1.0028E+01  6.1688E+01  2.7900E+01  6.6827E+01
 EBVSHRINKVR(%)  2.1761E+00  1.9050E+01  8.5322E+01  4.8016E+01  8.8996E+01
 RELATIVEINF(%)  9.6421E+01  1.6094E+01  6.1993E+00  8.0296E+00  1.9298E+00
 EPSSHRINKSD(%)  4.3266E+01
 EPSSHRINKVR(%)  6.7813E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1604.9135916382388     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -869.76276507450063     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    26.66
0R MATRIX ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
 BUT NONSINGULAR
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.71
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1604.914       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  1.91E+00  9.36E-02  3.10E-01  9.35E-01  9.31E-01  6.85E-01  1.65E+00  1.46E+00  3.42E-01  1.59E+00
 


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
+        1.21E+03
 
 TH 2
+       -1.29E+01  1.12E+04
 
 TH 3
+       -1.11E+02  3.12E+03  7.33E+03
 
 TH 4
+       -4.06E+01  1.07E+02 -3.14E+03  2.82E+03
 
 TH 5
+       -4.76E+01 -4.76E+02 -1.31E+03  9.52E+02  1.24E+03
 
 TH 6
+        2.05E+00  9.86E+00 -9.99E+00 -1.38E+01 -7.43E+00  2.17E+02
 
 TH 7
+        1.38E+00  5.22E+01 -1.42E+01  1.44E+01  2.62E+00  4.90E-02  2.85E+02
 
 TH 8
+       -1.84E+05  1.53E+04 -1.01E+05 -6.61E+04 -2.00E+01  1.02E+00  1.21E+00  2.18E+04
 
 TH 9
+        1.13E+00 -2.83E+01 -1.88E+01  7.03E+01 -2.11E+01 -9.55E-01  1.12E+01 -3.13E+04  2.79E+01
 
 TH10
+       -2.63E+00 -8.41E+00  3.65E+01 -6.88E+00 -3.08E+01  1.82E-01  8.07E+00  6.12E-01  5.36E+00  2.29E+01
 
 TH11
+       -1.09E+01  1.10E+01  1.21E+02 -7.25E+00 -1.81E+01  2.63E+00  1.31E+01 -2.43E+04  1.23E+01  2.66E+01  7.50E+01
 
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
 #CPUT: Total CPU Time in Seconds,       33.434
Stop Time:
Wed Sep 29 12:00:43 CDT 2021
