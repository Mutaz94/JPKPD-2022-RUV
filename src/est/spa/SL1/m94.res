Sat Sep 18 11:58:47 CDT 2021
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
$DATA ../../../../data/spa/SL1/dat94.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m94.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:  -1660.02200253561        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   4.7922E+01 -6.2607E+01 -3.4328E+01 -7.5293E+01  4.9385E+01 -6.4455E+00 -8.0379E+00  1.0124E+01 -3.9997E+01  2.5801E+01
             1.1950E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1669.63635785122        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  1.0080E+00  9.4644E-01  1.0113E+00  1.0863E+00  9.1122E-01  9.7555E-01  9.7874E-01  9.0104E-01  1.2053E+00  6.9449E-01
             1.0370E+00
 PARAMETER:  1.0801E-01  4.4952E-02  1.1121E-01  1.8275E-01  7.0331E-03  7.5250E-02  7.8511E-02 -4.2021E-03  2.8672E-01 -2.6457E-01
             1.3633E-01
 GRADIENT:   6.6956E+01  1.0369E+01  1.6702E+01  1.3895E+01 -7.9913E+00 -1.7882E+01  6.4598E-01 -1.4297E+00  1.6709E+01 -5.1687E+00
             1.9601E+01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1671.82262314597        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      157
 NPARAMETR:  1.0056E+00  9.5666E-01  7.3708E-01  1.0687E+00  7.9217E-01  9.8549E-01  1.2028E+00  6.1430E-01  1.1346E+00  5.8578E-01
             9.9404E-01
 PARAMETER:  1.0554E-01  5.5693E-02 -2.0506E-01  1.6649E-01 -1.3298E-01  8.5388E-02  2.8463E-01 -3.8728E-01  2.2625E-01 -4.3481E-01
             9.4025E-02
 GRADIENT:   5.7911E+01  6.9120E-01 -2.2252E+01  3.6304E+01  2.7493E+01 -1.4595E+01  5.5126E+00  3.1500E+00  1.4949E+01  6.8294E-01
             9.3061E+00

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1672.85267820686        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      228
 NPARAMETR:  9.8339E-01  9.8880E-01  7.1225E-01  1.0266E+00  7.8919E-01  1.0137E+00  1.1562E+00  5.0690E-01  1.0938E+00  6.0729E-01
             9.7286E-01
 PARAMETER:  8.3250E-02  8.8739E-02 -2.3933E-01  1.2629E-01 -1.3675E-01  1.1363E-01  2.4517E-01 -5.7944E-01  1.8962E-01 -3.9875E-01
             7.2485E-02
 GRADIENT:   8.6663E+00 -6.4454E-01 -4.9626E+00  2.4481E+00  4.9605E+00 -2.2022E+00  1.3956E+00  1.5984E+00  2.1087E+00  1.4387E+00
             9.5066E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1672.86749340462        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  9.8056E-01  9.7394E-01  6.9972E-01  1.0323E+00  7.7582E-01  1.0171E+00  1.1711E+00  4.4702E-01  1.0774E+00  6.0171E-01
             9.7151E-01
 PARAMETER:  8.0372E-02  7.3598E-02 -2.5707E-01  1.3175E-01 -1.5383E-01  1.1699E-01  2.5795E-01 -7.0516E-01  1.7459E-01 -4.0798E-01
             7.1093E-02
 GRADIENT:   2.6169E+00 -8.7343E-01 -2.9417E+00  6.7518E-01  2.8023E+00 -9.0555E-01  7.0112E-01  1.0541E+00  6.0908E-01  1.0448E+00
             4.1138E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1672.87106094781        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      368
 NPARAMETR:  9.7970E-01  9.7170E-01  6.8444E-01  1.0310E+00  7.6635E-01  1.0187E+00  1.1769E+00  3.8773E-01  1.0711E+00  5.9596E-01
             9.7095E-01
 PARAMETER:  7.9495E-02  7.1292E-02 -2.7915E-01  1.3055E-01 -1.6612E-01  1.1848E-01  2.6287E-01 -8.4745E-01  1.6871E-01 -4.1758E-01
             7.0521E-02
 GRADIENT:   5.8500E-01 -8.3054E-01 -1.8725E+00  4.6502E-02  1.6507E+00 -3.6876E-01  3.7468E-01  6.9195E-01  1.1142E-01  7.2717E-01
             2.0302E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1672.87148098048        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      438
 NPARAMETR:  9.7948E-01  9.7450E-01  6.7057E-01  1.0274E+00  7.5979E-01  1.0194E+00  1.1780E+00  3.3006E-01  1.0689E+00  5.9176E-01
             9.7068E-01
 PARAMETER:  7.9264E-02  7.4170E-02 -2.9963E-01  1.2701E-01 -1.7471E-01  1.1920E-01  2.6385E-01 -1.0085E+00  1.6663E-01 -4.2466E-01
             7.0244E-02
 GRADIENT:  -1.2341E-01 -6.3175E-01 -1.1986E+00 -1.4296E-01  9.2181E-01 -1.3053E-01  2.0065E-01  4.4524E-01 -6.0019E-02  4.9911E-01
             1.1339E-01

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1672.87253642592        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      508
 NPARAMETR:  9.7953E-01  9.7886E-01  6.5956E-01  1.0234E+00  7.5549E-01  1.0197E+00  1.1771E+00  2.7704E-01  1.0687E+00  5.8901E-01
             9.7056E-01
 PARAMETER:  7.9313E-02  7.8634E-02 -3.1618E-01  1.2309E-01 -1.8039E-01  1.1950E-01  2.6309E-01 -1.1836E+00  1.6641E-01 -4.2932E-01
             7.0114E-02
 GRADIENT:  -2.1470E-01 -4.1956E-01 -7.9048E-01 -1.3099E-01  5.1481E-01 -4.9020E-02  1.1495E-01  2.8444E-01 -7.8088E-02  3.4011E-01
             7.9133E-02

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1672.87329733139        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      578
 NPARAMETR:  9.7964E-01  9.8299E-01  6.5136E-01  1.0199E+00  7.5267E-01  1.0198E+00  1.1757E+00  2.2966E-01  1.0691E+00  5.8729E-01
             9.7049E-01
 PARAMETER:  7.9433E-02  8.2841E-02 -3.2869E-01  1.1968E-01 -1.8412E-01  1.1961E-01  2.6190E-01 -1.3712E+00  1.6683E-01 -4.3224E-01
             7.0042E-02
 GRADIENT:  -1.1168E-01 -2.5862E-01 -5.4289E-01 -6.6824E-02  3.0671E-01 -3.0833E-02  7.2968E-02  1.8114E-01 -4.6790E-02  2.3025E-01
             6.4848E-02

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1672.87340515365        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      649
 NPARAMETR:  9.7970E-01  9.8478E-01  6.4811E-01  1.0184E+00  7.5163E-01  1.0198E+00  1.1751E+00  2.0734E-01  1.0694E+00  5.8670E-01
             9.7046E-01
 PARAMETER:  7.9487E-02  8.4659E-02 -3.3370E-01  1.1823E-01 -1.8552E-01  1.1965E-01  2.6132E-01 -1.4734E+00  1.6709E-01 -4.3324E-01
             7.0012E-02
 GRADIENT:  -5.5082E-02 -1.9805E-01 -4.4388E-01 -3.8851E-02  2.4080E-01 -2.8947E-02  5.8368E-02  1.4279E-01 -2.9894E-02  1.8462E-01
             5.7332E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1673.35173219742        NO. OF FUNC. EVALS.: 119
 CUMULATIVE NO. OF FUNC. EVALS.:      768
 NPARAMETR:  9.8937E-01  9.5766E-01  6.9003E-01  1.0443E+00  7.6418E-01  1.0371E+00  1.1872E+00  2.2324E-01  1.0611E+00  6.1418E-01
             9.7189E-01
 PARAMETER:  8.9312E-02  5.6741E-02 -2.7102E-01  1.4332E-01 -1.6895E-01  1.3639E-01  2.7161E-01 -1.3995E+00  1.5931E-01 -3.8747E-01
             7.1491E-02
 GRADIENT:  -2.1240E+01 -8.1134E-02  4.5495E+00 -4.6328E+00 -4.1137E+00  1.3678E+00 -2.5520E+00 -1.3373E-01 -2.0984E+00 -1.8550E+00
            -8.4651E-01

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1673.57948423191        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      943
 NPARAMETR:  1.0007E+00  8.7139E-01  7.5266E-01  1.1083E+00  7.6282E-01  1.0293E+00  1.2731E+00  4.4528E-01  1.0406E+00  6.3057E-01
             9.7325E-01
 PARAMETER:  1.0069E-01 -3.7668E-02 -1.8415E-01  2.0282E-01 -1.7073E-01  1.2890E-01  3.4144E-01 -7.0905E-01  1.3983E-01 -3.6113E-01
             7.2881E-02
 GRADIENT:   3.7541E+00  1.5663E+00 -1.1004E+00  1.7001E+00 -7.7805E-01 -1.0815E+00  5.3387E-01  2.7509E-01  1.1865E+00  5.5454E-01
             6.5429E-01

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1673.65190086722        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1120
 NPARAMETR:  9.9688E-01  7.5817E-01  8.2948E-01  1.1844E+00  7.5968E-01  1.0317E+00  1.3650E+00  5.5904E-01  9.9836E-01  6.4571E-01
             9.7184E-01
 PARAMETER:  9.6871E-02 -1.7685E-01 -8.6957E-02  2.6927E-01 -1.7485E-01  1.3117E-01  4.1118E-01 -4.8154E-01  9.8361E-02 -3.3741E-01
             7.1438E-02
 GRADIENT:  -1.2779E+00  2.1755E+00  3.1575E+00  8.0297E-01 -3.4043E+00  3.1562E-01 -4.7072E-01 -2.3520E-01 -9.9240E-01 -6.1714E-01
            -5.4368E-01

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1673.68932782657        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1296
 NPARAMETR:  9.9585E-01  6.5687E-01  8.7016E-01  1.2484E+00  7.4663E-01  1.0298E+00  1.4894E+00  6.1062E-01  9.6950E-01  6.5926E-01
             9.7204E-01
 PARAMETER:  9.5841E-02 -3.2027E-01 -3.9073E-02  3.2184E-01 -1.9218E-01  1.2932E-01  4.9840E-01 -3.9328E-01  6.9020E-02 -3.1663E-01
             7.1645E-02
 GRADIENT:  -7.5451E-01  6.9135E-01  2.0997E-01  7.4102E-01 -7.5752E-01  2.1108E-02  2.8818E-01  5.2791E-02 -1.2379E-01  1.7514E-01
             6.7974E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1673.69256918694        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1471
 NPARAMETR:  9.9549E-01  6.1082E-01  8.9678E-01  1.2781E+00  7.4512E-01  1.0290E+00  1.5292E+00  6.4395E-01  9.5972E-01  6.6660E-01
             9.7200E-01
 PARAMETER:  9.5476E-02 -3.9296E-01 -8.9428E-03  3.4540E-01 -1.9421E-01  1.2858E-01  5.2474E-01 -3.4013E-01  5.8885E-02 -3.0557E-01
             7.1603E-02
 GRADIENT:   7.9641E-02  1.2635E-01  1.6246E-01  2.0145E-01 -2.2912E-01 -5.2835E-04 -1.2160E-02 -1.2877E-02  1.7776E-02 -2.9946E-02
            -8.6882E-03

0ITERATION NO.:   75    OBJECTIVE VALUE:  -1673.69269589153        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1646
 NPARAMETR:  9.9476E-01  5.7471E-01  9.1903E-01  1.3030E+00  7.4406E-01  1.0284E+00  1.5640E+00  6.7152E-01  9.5144E-01  6.7445E-01
             9.7201E-01
 PARAMETER:  9.4744E-02 -4.5388E-01  1.5568E-02  3.6465E-01 -1.9563E-01  1.2801E-01  5.4728E-01 -2.9821E-01  5.0219E-02 -2.9386E-01
             7.1613E-02
 GRADIENT:  -8.5072E-02  1.1563E+00  1.2871E+00  1.4359E+00 -2.3122E+00  1.4290E-03  2.0319E-01  8.3957E-03 -1.5407E-01  1.5401E-01
             2.3963E-02

0ITERATION NO.:   80    OBJECTIVE VALUE:  -1673.69611529445        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     1829
 NPARAMETR:  9.9328E-01  5.2560E-01  9.5126E-01  1.3375E+00  7.4334E-01  1.0275E+00  1.6139E+00  7.0929E-01  9.3959E-01  6.8940E-01
             9.7202E-01
 PARAMETER:  9.3253E-02 -5.4321E-01  5.0029E-02  3.9077E-01 -1.9660E-01  1.2713E-01  5.7866E-01 -2.4349E-01  3.7690E-02 -2.7193E-01
             7.1625E-02
 GRADIENT:  -1.1981E+00  3.4453E+00  3.8536E+00  3.8796E+00 -7.4098E+00 -1.7020E-02  9.1423E-01  2.0020E-01 -8.4332E-01  1.0485E+00
             2.0187E-01

0ITERATION NO.:   85    OBJECTIVE VALUE:  -1673.93179302194        NO. OF FUNC. EVALS.: 179
 CUMULATIVE NO. OF FUNC. EVALS.:     2008
 NPARAMETR:  9.8774E-01  3.3600E-01  9.9010E-01  1.4532E+00  7.1109E-01  1.0260E+00  1.9798E+00  7.5598E-01  8.9262E-01  7.1085E-01
             9.7131E-01
 PARAMETER:  8.7661E-02 -9.9065E-01  9.0050E-02  4.7376E-01 -2.4095E-01  1.2563E-01  7.8300E-01 -1.7974E-01 -1.3592E-02 -2.4129E-01
             7.0893E-02
 GRADIENT:  -5.5510E+00  3.6515E+00  3.2423E+00  8.9757E+00 -1.0284E+01  3.8086E-01  1.1352E+00  9.3838E-01 -1.4860E+00  2.8418E+00
             7.1103E-01

0ITERATION NO.:   90    OBJECTIVE VALUE:  -1674.25791355747        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2183
 NPARAMETR:  9.8759E-01  1.8836E-01  1.0467E+00  1.5401E+00  7.0185E-01  1.0219E+00  2.3705E+00  8.3652E-01  8.6535E-01  6.9241E-01
             9.7082E-01
 PARAMETER:  8.7510E-02 -1.5694E+00  1.4561E-01  5.3187E-01 -2.5403E-01  1.2163E-01  9.6309E-01 -7.8509E-02 -4.4626E-02 -2.6758E-01
             7.0391E-02
 GRADIENT:   1.2505E+00  1.1150E+00  5.2255E+00  6.3792E+00 -5.9564E+00 -3.3803E-01  1.6014E-01 -5.6601E-01  3.8455E-01 -4.3109E-01
            -4.3580E-01

0ITERATION NO.:   95    OBJECTIVE VALUE:  -1674.36736161099        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2358
 NPARAMETR:  9.8517E-01  1.0675E-01  1.0322E+00  1.5807E+00  6.7993E-01  1.0218E+00  2.5971E+00  8.2814E-01  8.4809E-01  6.9845E-01
             9.7078E-01
 PARAMETER:  8.5060E-02 -2.1373E+00  1.3172E-01  5.5784E-01 -2.8577E-01  1.2152E-01  1.0544E+00 -8.8568E-02 -6.4764E-02 -2.5889E-01
             7.0344E-02
 GRADIENT:  -8.0970E-02  1.4317E-01  4.1912E-01  9.9174E-01 -1.4263E+00  2.8317E-02 -8.3273E-02 -1.4989E-01 -7.3415E-02  1.0607E-01
             5.0632E-02

0ITERATION NO.:  100    OBJECTIVE VALUE:  -1674.39811243790        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2533
 NPARAMETR:  9.8387E-01  5.7073E-02  1.0558E+00  1.6116E+00  6.7942E-01  1.0208E+00  2.8382E+00  8.6467E-01  8.3561E-01  6.9442E-01
             9.7070E-01
 PARAMETER:  8.3742E-02 -2.7634E+00  1.5430E-01  5.7723E-01 -2.8651E-01  1.2063E-01  1.1432E+00 -4.5404E-02 -7.9588E-02 -2.6468E-01
             7.0257E-02
 GRADIENT:  -2.3439E-01  7.6835E-02  3.2367E-01  1.9760E+00 -1.7074E-01 -3.9672E-02 -3.8937E-02 -7.0824E-02 -2.5031E-01 -9.1563E-02
            -4.1286E-02

0ITERATION NO.:  105    OBJECTIVE VALUE:  -1674.41627093481        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     2708
 NPARAMETR:  9.8316E-01  2.3051E-02  1.0537E+00  1.6299E+00  6.7117E-01  1.0204E+00  3.1459E+00  8.6772E-01  8.2759E-01  6.9257E-01
             9.7048E-01
 PARAMETER:  8.3021E-02 -3.6701E+00  1.5227E-01  5.8855E-01 -2.9874E-01  1.2023E-01  1.2461E+00 -4.1882E-02 -8.9233E-02 -2.6735E-01
             7.0038E-02
 GRADIENT:   3.2056E-02  3.9359E-02  7.3720E-01  1.8691E+00 -1.4522E+00 -1.2254E-02 -1.0258E-02  3.7497E-02 -5.8657E-02  7.4467E-02
             6.3907E-03

0ITERATION NO.:  110    OBJECTIVE VALUE:  -1674.42496782983        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     2884
 NPARAMETR:  9.8284E-01  1.0370E-02  1.0540E+00  1.6360E+00  6.6917E-01  1.0203E+00  3.4167E+00  8.6899E-01  8.2432E-01  6.9170E-01
             9.7041E-01
 PARAMETER:  8.2691E-02 -4.4689E+00  1.5260E-01  5.9223E-01 -3.0171E-01  1.2008E-01  1.3287E+00 -4.0424E-02 -9.3191E-02 -2.6861E-01
             6.9963E-02
 GRADIENT:   4.6612E-02  5.5198E-03 -8.8482E-02 -2.7381E-01  2.7100E-01  4.5798E-03 -2.4592E-03 -9.3173E-03  3.5178E-02 -3.0727E-02
            -7.8253E-03

0ITERATION NO.:  115    OBJECTIVE VALUE:  -1674.42550990915        NO. OF FUNC. EVALS.: 182
 CUMULATIVE NO. OF FUNC. EVALS.:     3066
 NPARAMETR:  9.8280E-01  1.0000E-02  1.0528E+00  1.6361E+00  6.6865E-01  1.0202E+00  3.8077E+00  8.6784E-01  8.2415E-01  6.9183E-01
             9.7041E-01
 PARAMETER:  8.2654E-02 -4.5177E+00  1.5150E-01  5.9234E-01 -3.0250E-01  1.2003E-01  1.4370E+00 -4.1750E-02 -9.3405E-02 -2.6842E-01
             6.9965E-02
 GRADIENT:  -3.9028E-02  0.0000E+00 -2.0371E-01 -1.8683E-02  2.7995E-01 -1.0090E-02 -2.8514E-03  1.2303E-02 -7.6106E-03 -1.3551E-03
             7.6500E-03

0ITERATION NO.:  120    OBJECTIVE VALUE:  -1674.42866200015        NO. OF FUNC. EVALS.: 183
 CUMULATIVE NO. OF FUNC. EVALS.:     3249            RESET HESSIAN, TYPE II
 NPARAMETR:  9.8282E-01  1.0000E-02  1.0523E+00  1.6361E+00  6.6890E-01  1.0203E+00  8.7333E+00  8.6876E-01  8.2368E-01  6.9115E-01
             9.7046E-01
 PARAMETER:  8.2672E-02 -4.5177E+00  1.5101E-01  5.9231E-01 -3.0213E-01  1.2006E-01  2.2671E+00 -4.0690E-02 -9.3971E-02 -2.6939E-01
             7.0013E-02
 GRADIENT:   4.2300E+01  0.0000E+00 -7.4638E-01  1.1495E+02  5.1443E+00  5.3261E+00  4.4655E-02  1.3147E-01  1.4932E+00  1.1156E-01
             1.4047E-01

0ITERATION NO.:  125    OBJECTIVE VALUE:  -1674.42972596247        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     3426
 NPARAMETR:  9.8282E-01  1.0000E-02  1.0529E+00  1.6362E+00  6.6854E-01  1.0202E+00  9.3568E+00  8.6800E-01  8.2355E-01  6.9103E-01
             9.7038E-01
 PARAMETER:  8.2675E-02 -4.5177E+00  1.5154E-01  5.9237E-01 -3.0266E-01  1.2005E-01  2.3361E+00 -4.1560E-02 -9.4126E-02 -2.6958E-01
             6.9937E-02
 GRADIENT:   3.9359E-02  0.0000E+00  4.8469E-02 -2.8630E-02  4.2436E-02  1.1863E-02 -1.6361E-04 -3.9895E-03 -4.0375E-03 -6.2240E-03
            -5.2004E-03

0ITERATION NO.:  126    OBJECTIVE VALUE:  -1674.42972596247        NO. OF FUNC. EVALS.:  22
 CUMULATIVE NO. OF FUNC. EVALS.:     3448
 NPARAMETR:  9.8282E-01  1.0000E-02  1.0529E+00  1.6362E+00  6.6854E-01  1.0202E+00  9.3568E+00  8.6800E-01  8.2355E-01  6.9103E-01
             9.7038E-01
 PARAMETER:  8.2675E-02 -4.5177E+00  1.5154E-01  5.9237E-01 -3.0266E-01  1.2005E-01  2.3361E+00 -4.1560E-02 -9.4126E-02 -2.6958E-01
             6.9937E-02
 GRADIENT:   3.9359E-02  0.0000E+00  4.8469E-02 -2.8630E-02  4.2436E-02  1.1863E-02 -1.6361E-04 -3.9895E-03 -4.0375E-03 -6.2240E-03
            -5.2004E-03

 #TERM:
0MINIMIZATION SUCCESSFUL
 NO. OF FUNCTION EVALUATIONS USED:     3448
 NO. OF SIG. DIGITS IN FINAL EST.:  3.3
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         1.5629E-04 -3.1450E-04 -1.7641E-02 -3.7673E-03 -2.2767E-02
 SE:             2.9845E-02  1.4960E-03  1.8174E-02  2.9406E-02  2.0483E-02
 N:                     100         100         100         100         100

 P VAL.:         9.9582E-01  8.3349E-01  3.3173E-01  8.9806E-01  2.6635E-01

 ETASHRINKSD(%)  1.6451E-02  9.4988E+01  3.9113E+01  1.4858E+00  3.1380E+01
 ETASHRINKVR(%)  3.2899E-02  9.9749E+01  6.2928E+01  2.9494E+00  5.2913E+01
 EBVSHRINKSD(%)  3.9198E-01  9.5078E+01  4.0424E+01  1.8788E+00  3.0335E+01
 EBVSHRINKVR(%)  7.8243E-01  9.9758E+01  6.4507E+01  3.7223E+00  5.1468E+01
 RELATIVEINF(%)  9.2747E+01  7.9792E-03  4.0650E+00  5.3778E+00  2.2626E+00
 EPSSHRINKSD(%)  4.4796E+01
 EPSSHRINKVR(%)  6.9525E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1674.4297259624657     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -939.27889939872750     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    38.41
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     5.55
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1674.430       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.83E-01  1.00E-02  1.05E+00  1.64E+00  6.69E-01  1.02E+00  9.36E+00  8.68E-01  8.24E-01  6.91E-01  9.70E-01
 


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
+        1.12E+03
 
 TH 2
+        0.00E+00  0.00E+00
 
 TH 3
+       -1.10E+00  0.00E+00  4.88E+02
 
 TH 4
+       -7.85E+00  0.00E+00 -3.95E+01  5.90E+02
 
 TH 5
+        6.77E-01  0.00E+00 -1.00E+03 -9.47E+01  2.56E+03
 
 TH 6
+        1.03E+01  0.00E+00 -7.10E-01 -1.77E+00 -1.50E+00  1.95E+02
 
 TH 7
+        2.26E-02  0.00E+00  5.89E-03 -3.40E-03 -1.21E-02 -6.33E-03  6.31E-04
 
 TH 8
+       -2.93E+00  0.00E+00 -4.23E+01 -2.55E+00 -1.01E+01  9.26E+00 -1.89E-02  6.52E+01
 
 TH 9
+        1.54E+01  0.00E+00  1.14E+01 -1.12E+00  5.79E+00  1.63E+00 -1.43E-04 -6.84E-02  2.91E+02
 
 TH10
+       -3.09E+00  0.00E+00 -1.54E+00  1.40E-01 -1.30E+02 -1.13E-01  1.92E-02  5.06E+01  4.96E+00  1.10E+02
 
 TH11
+       -2.74E+00  0.00E+00 -1.30E+01 -5.89E+00 -1.20E+01  4.56E+00  3.51E-02  3.27E+00  6.59E+00  2.09E+01  2.29E+02
 
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
 #CPUT: Total CPU Time in Seconds,       44.019
Stop Time:
Sat Sep 18 11:59:32 CDT 2021
