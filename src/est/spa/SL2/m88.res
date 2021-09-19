Sat Sep 18 12:31:24 CDT 2021
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
 GRADIENT:   3.1567E+01 -8.5124E+01 -4.8228E+01 -5.1057E+01  7.9535E+01  8.3117E+00 -1.1385E+01 -4.5724E+00  8.3677E-01  2.1383E+01
            -7.6821E+01

0ITERATION NO.:    5    OBJECTIVE VALUE:  -1649.54478597676        NO. OF FUNC. EVALS.:  73
 CUMULATIVE NO. OF FUNC. EVALS.:       86
 NPARAMETR:  9.9518E-01  1.1485E+00  1.1187E+00  9.7239E-01  1.0073E+00  9.7220E-01  1.1103E+00  1.1246E+00  9.6843E-01  6.9568E-01
             1.2022E+00
 PARAMETER:  9.5166E-02  2.3847E-01  2.1216E-01  7.2001E-02  1.0727E-01  7.1806E-02  2.0461E-01  2.1746E-01  6.7922E-02 -2.6286E-01
             2.8416E-01
 GRADIENT:   3.1823E+00  3.0238E+01  2.3636E+01  9.5101E+00 -1.4598E+01 -3.0967E+00 -1.1234E+00 -1.5550E+01 -2.9635E+00 -2.7873E+00
             9.9864E-01

0ITERATION NO.:   10    OBJECTIVE VALUE:  -1651.84339024655        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      160
 NPARAMETR:  9.9948E-01  1.0435E+00  1.0890E+00  1.0376E+00  9.5049E-01  9.8170E-01  1.2338E+00  1.3896E+00  9.1096E-01  5.8592E-01
             1.1892E+00
 PARAMETER:  9.9479E-02  1.4256E-01  1.8525E-01  1.3687E-01  4.9221E-02  8.1529E-02  3.1013E-01  4.2902E-01  6.7386E-03 -4.3458E-01
             2.7326E-01
 GRADIENT:   1.4641E+01  1.7275E+01  3.4668E+00  2.2305E+01 -9.7019E+00  8.9267E-01  2.6016E+00 -2.9346E-01  6.6086E-01 -3.9171E-01
             7.8348E-02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1652.10558191673        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:      298
 NPARAMETR:  1.0057E+00  1.0602E+00  1.1480E+00  1.0175E+00  9.8785E-01  9.8637E-01  1.1714E+00  1.4449E+00  9.3117E-01  6.4390E-01
             1.1902E+00
 PARAMETER:  1.0571E-01  1.5847E-01  2.3804E-01  1.1732E-01  8.7774E-02  8.6279E-02  2.5820E-01  4.6807E-01  2.8682E-02 -3.4021E-01
             2.7411E-01
 GRADIENT:  -5.5379E-01 -1.6496E+00  4.4827E-01 -2.2033E+00  8.9375E-02 -8.7849E-03 -3.3214E-01 -6.4542E-01 -5.2314E-01 -1.1817E-02
             2.4579E-01

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1652.15878859170        NO. OF FUNC. EVALS.: 178
 CUMULATIVE NO. OF FUNC. EVALS.:      476
 NPARAMETR:  1.0067E+00  1.2007E+00  9.8058E-01  9.2083E-01  9.8303E-01  9.8731E-01  1.0890E+00  1.3848E+00  9.8060E-01  6.1684E-01
             1.1907E+00
 PARAMETER:  1.0671E-01  2.8290E-01  8.0389E-02  1.7517E-02  8.2883E-02  8.7226E-02  1.8528E-01  4.2552E-01  8.0413E-02 -3.8315E-01
             2.7452E-01
 GRADIENT:  -1.0720E+00 -1.5625E+00  1.5736E+00 -3.4593E+00 -1.8295E+00 -1.3266E-01 -3.6482E-01 -6.8009E-01 -6.3970E-01  3.6555E-02
             5.5176E-01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1652.24739225828        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:      653
 NPARAMETR:  1.0095E+00  1.3752E+00  7.1280E-01  8.0345E-01  9.4213E-01  9.9011E-01  1.0193E+00  1.2508E+00  1.0423E+00  5.2947E-01
             1.1858E+00
 PARAMETER:  1.0949E-01  4.1860E-01 -2.3856E-01 -1.1884E-01  4.0383E-02  9.0058E-02  1.1908E-01  3.2382E-01  1.4145E-01 -5.3588E-01
             2.7046E-01
 GRADIENT:   1.4395E+00  9.4429E+00  1.3025E+00  6.4620E+00 -4.9551E+00  1.6610E-01  5.8615E-01  1.6923E-01  2.2456E-01  1.1744E-01
            -4.5551E-01

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1653.67541892582        NO. OF FUNC. EVALS.: 181
 CUMULATIVE NO. OF FUNC. EVALS.:      834
 NPARAMETR:  1.0022E+00  1.6459E+00  2.4262E-01  5.6667E-01  7.9414E-01  9.8818E-01  9.3371E-01  6.8335E-01  1.1496E+00  2.7611E-01
             1.1692E+00
 PARAMETER:  1.0219E-01  5.9830E-01 -1.3163E+00 -4.6797E-01 -1.3049E-01  8.8112E-02  3.1414E-02 -2.8075E-01  2.3945E-01 -1.1869E+00
             2.5632E-01
 GRADIENT:  -5.0894E+00  1.9230E+01  3.1172E+00  2.6072E+00 -1.0648E+01 -3.7189E-01  1.7500E+01  4.8885E-01 -1.0080E-01  2.1586E+00
             7.0135E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1656.76694991014        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1009
 NPARAMETR:  9.9990E-01  1.7873E+00  1.4952E-01  4.3879E-01  8.1940E-01  9.9082E-01  8.0660E-01  5.3352E-01  1.3970E+00  1.9786E-01
             1.1085E+00
 PARAMETER:  9.9903E-02  6.8068E-01 -1.8003E+00 -7.2372E-01 -9.9188E-02  9.0775E-02 -1.1493E-01 -5.2827E-01  4.3433E-01 -1.5202E+00
             2.0298E-01
 GRADIENT:   4.4388E+00 -1.6103E+01  4.3395E+00 -5.9280E+00 -1.3660E+00  8.4727E-01 -5.2606E+00 -1.8593E+00  2.2146E+00  1.4124E-01
            -4.6366E+00

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1657.39275554376        NO. OF FUNC. EVALS.: 186
 CUMULATIVE NO. OF FUNC. EVALS.:     1195
 NPARAMETR:  9.9653E-01  1.8632E+00  1.3680E-01  4.0186E-01  8.4965E-01  9.8788E-01  8.1751E-01  6.9965E-01  1.4697E+00  1.5880E-01
             1.1324E+00
 PARAMETER:  9.6526E-02  7.2229E-01 -1.8893E+00 -8.1164E-01 -6.2928E-02  8.7805E-02 -1.0149E-01 -2.5718E-01  4.8504E-01 -1.7401E+00
             2.2438E-01
 GRADIENT:  -3.2356E+00  9.7935E+00  2.2613E+00 -6.0837E+00 -6.5759E+00 -7.1971E-01  4.9796E+00 -2.1632E+00 -2.5359E+00 -2.9007E-01
             5.8735E+00

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1657.66035416388        NO. OF FUNC. EVALS.: 111
 CUMULATIVE NO. OF FUNC. EVALS.:     1306
 NPARAMETR:  9.9653E-01  1.8632E+00  1.3680E-01  4.0187E-01  8.5689E-01  9.8424E-01  7.9472E-01  6.9965E-01  1.4697E+00  2.7328E-01
             1.1042E+00
 PARAMETER:  9.6527E-02  7.2229E-01 -1.8893E+00 -8.1164E-01 -5.4447E-02  8.4111E-02 -1.2977E-01 -2.5718E-01  4.8504E-01 -1.1973E+00
             1.9911E-01
 GRADIENT:   3.0850E+01  6.4715E+01  2.4130E+00  8.5686E+00  4.7173E+00  1.0979E+00  5.6089E-02 -2.5922E+00 -6.9942E-01  1.0621E-01
            -5.7643E-01

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1657.67284857227        NO. OF FUNC. EVALS.: 135
 CUMULATIVE NO. OF FUNC. EVALS.:     1441
 NPARAMETR:  9.9653E-01  1.8632E+00  1.3680E-01  4.0187E-01  8.5536E-01  9.8997E-01  7.9787E-01  6.9965E-01  1.4697E+00  2.6352E-01
             1.1064E+00
 PARAMETER:  9.6527E-02  7.2229E-01 -1.8893E+00 -8.1164E-01 -5.6232E-02  8.9921E-02 -1.2581E-01 -2.5718E-01  4.8504E-01 -1.2336E+00
             2.0108E-01
 GRADIENT:  -3.2940E+00 -3.0654E+00  1.6207E+00 -2.5404E+00  2.0452E-01 -2.3399E-02  1.9307E-01 -2.5634E+00 -1.5310E+00 -6.7980E-03
            -3.1728E-02

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1658.33594524090        NO. OF FUNC. EVALS.: 185
 CUMULATIVE NO. OF FUNC. EVALS.:     1626
 NPARAMETR:  9.9688E-01  1.8638E+00  1.3680E-01  4.0218E-01  8.5547E-01  9.9015E-01  7.9718E-01  1.0264E+00  1.4721E+00  2.6309E-01
             1.1050E+00
 PARAMETER:  9.6870E-02  7.2264E-01 -1.8892E+00 -8.1086E-01 -5.6099E-02  9.0103E-02 -1.2668E-01  1.2609E-01  4.8667E-01 -1.2353E+00
             1.9984E-01
 GRADIENT:  -1.6331E+00 -4.5113E+00 -4.1780E+00 -2.2861E+00 -7.5196E+00 -5.6992E-01  1.0049E+00  1.1140E-01 -1.6734E+00  3.1344E-01
             7.8241E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1658.42779535500        NO. OF FUNC. EVALS.: 176
 CUMULATIVE NO. OF FUNC. EVALS.:     1802
 NPARAMETR:  9.9706E-01  1.8644E+00  1.3779E-01  4.0255E-01  8.5678E-01  9.9079E-01  7.9648E-01  1.0480E+00  1.4763E+00  2.5860E-01
             1.0982E+00
 PARAMETER:  9.7056E-02  7.2295E-01 -1.8820E+00 -8.0993E-01 -5.4577E-02  9.0750E-02 -1.2755E-01  1.4685E-01  4.8957E-01 -1.2525E+00
             1.9368E-01
 GRADIENT:  -1.3016E+00 -4.0750E+00 -4.0194E+00 -1.6081E+00 -5.2247E+00 -4.3032E-01  4.6656E-01  7.5000E-02 -1.6549E+00  1.4484E-01
             5.5963E+00

0ITERATION NO.:   63    OBJECTIVE VALUE:  -1658.42918420934        NO. OF FUNC. EVALS.: 115
 CUMULATIVE NO. OF FUNC. EVALS.:     1917
 NPARAMETR:  9.9706E-01  1.8644E+00  1.3779E-01  4.0255E-01  8.5678E-01  9.9091E-01  7.9636E-01  1.0480E+00  1.4763E+00  2.5848E-01
             1.0980E+00
 PARAMETER:  9.7056E-02  7.2296E-01 -1.8820E+00 -8.0993E-01 -5.4578E-02  9.0766E-02 -1.2757E-01  1.4685E-01  4.8957E-01 -1.2529E+00
             1.9345E-01
 GRADIENT:   5.0027E+05  6.9192E+04  5.0805E+01  6.1749E+04  5.0027E+05 -4.1669E-01  4.4178E-01 -3.4072E+05  1.0216E+05 -3.9929E+04
            -2.5861E+05
 NUMSIGDIG:         3.3         3.3         6.3         3.3         3.3         1.7         1.7         3.3         3.3         3.3
                    3.3

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1917
 NO. OF SIG. DIGITS IN FINAL EST.:  1.7

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         7.8761E-04 -1.0048E-02 -4.5688E-03  1.8923E-02 -2.1040E-02
 SE:             2.9855E-02  2.8299E-02  9.2989E-03  2.4806E-02  8.9382E-03
 N:                     100         100         100         100         100

 P VAL.:         9.7895E-01  7.2253E-01  6.2319E-01  4.4555E-01  1.8578E-02

 ETASHRINKSD(%)  1.0000E-10  5.1952E+00  6.8848E+01  1.6897E+01  7.0056E+01
 ETASHRINKVR(%)  1.0000E-10  1.0121E+01  9.0295E+01  3.0939E+01  9.1033E+01
 EBVSHRINKSD(%)  5.2053E-01  5.5723E+00  6.8683E+01  1.6047E+01  7.1007E+01
 EBVSHRINKVR(%)  1.0383E+00  1.0834E+01  9.0192E+01  2.9519E+01  9.1594E+01
 RELATIVEINF(%)  9.8235E+01  1.9678E+01  3.9466E+00  1.1120E+01  1.4432E+00
 EPSSHRINKSD(%)  4.4374E+01
 EPSSHRINKVR(%)  6.9058E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1658.4291842093410     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -923.27835764560280     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    25.73
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0S MATRIX UNOBTAINABLE
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.81
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1658.429       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         9.97E-01  1.86E+00  1.38E-01  4.03E-01  8.57E-01  9.91E-01  7.96E-01  1.05E+00  1.48E+00  2.58E-01  1.10E+00
 


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
+        1.26E+09
 
 TH 2
+        3.39E+02  6.88E+06
 
 TH 3
+       -4.83E+08 -3.57E+07  3.72E+08
 
 TH 4
+        1.42E+03 -2.85E+07 -1.48E+08  1.18E+08
 
 TH 5
+        4.79E+04  3.05E+03 -1.77E+03  1.59E+04  1.70E+09
 
 TH 6
+        2.96E+03  2.17E+02  3.95E+00  8.96E+02  3.43E+03  2.01E+02
 
 TH 7
+        6.94E+03  5.23E+02 -9.56E+01  2.12E+03  8.10E+03 -6.24E-01  1.21E+09
 
 TH 8
+       -3.07E+05 -2.27E+04 -2.96E+02 -9.39E+04 -3.11E+04 -1.91E+03 -4.49E+03  5.28E+08
 
 TH 9
+        6.49E+02 -6.52E+02 -6.67E+07 -2.86E+04  6.57E+03  4.07E+02  9.61E+02 -4.24E+04  2.39E+07
 
 TH10
+       -2.22E+04 -1.65E+03  4.46E+01 -6.84E+03 -1.49E+04 -9.09E+02 -3.80E+08 -2.51E+08 -3.06E+03  1.19E+08
 
 TH11
+       -3.85E+04 -2.85E+03  2.27E+08 -1.18E+04  6.87E+08 -1.38E+03 -3.25E+03 -3.83E+08 -5.28E+03  1.05E+04  2.77E+08
 
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
 #CPUT: Total CPU Time in Seconds,       32.611
Stop Time:
Sat Sep 18 12:31:58 CDT 2021
