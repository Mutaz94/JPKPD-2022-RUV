Sat Sep 25 09:23:48 CDT 2021
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
$DATA ../../../../data/spa/A3/dat62.csv ignore=@
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
 RAW OUTPUT FILE (FILE): m62.ext
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


0ITERATION NO.:    0    OBJECTIVE VALUE:   606.112687478183        NO. OF FUNC. EVALS.:  13
 CUMULATIVE NO. OF FUNC. EVALS.:       13
 NPARAMETR:  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00  1.0000E+00
             1.0000E+00
 PARAMETER:  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01  1.0000E-01
             1.0000E-01
 GRADIENT:   1.7323E+01  1.0781E+02  1.6062E+02 -3.6004E+01  1.9686E+02  2.9663E+01 -8.2707E+01 -8.7129E+01 -1.9188E+02 -2.2796E+02
            -3.9573E+03

0ITERATION NO.:    5    OBJECTIVE VALUE:  -869.653556233740        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:       84
 NPARAMETR:  1.3010E+00  8.2864E-01  6.2240E-01  1.3951E+00  5.9503E-01  5.2486E-01  1.1632E+00  1.1885E+00  1.2483E+00  1.5751E+00
             1.3972E+01
 PARAMETER:  3.6314E-01 -8.7972E-02 -3.7418E-01  4.3299E-01 -4.1915E-01 -5.4463E-01  2.5118E-01  2.7272E-01  3.2182E-01  5.5431E-01
             2.7370E+00
 GRADIENT:  -9.7026E-01 -3.4421E+01 -2.4792E+01 -4.4357E+01  1.0725E+01 -9.5256E+00  8.0060E+00  8.7416E+00  4.2142E+01  3.8518E+01
             4.2158E+02

0ITERATION NO.:   10    OBJECTIVE VALUE:  -952.041287853057        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      154
 NPARAMETR:  1.0329E+00  3.6072E-01  5.6566E-02  1.0708E+00  8.8939E-02  7.1544E-01  1.8233E+00  1.9261E+00  1.8843E+00  1.3814E+00
             7.8187E+00
 PARAMETER:  1.3238E-01 -9.1965E-01 -2.7723E+00  1.6840E-01 -2.3198E+00 -2.3486E-01  7.0064E-01  7.5551E-01  7.3356E-01  4.2311E-01
             2.1565E+00
 GRADIENT:  -2.1955E+02  1.3011E+02  9.5045E+01  8.8863E+00 -3.3003E+02 -1.7513E+01  1.6545E+01  1.7603E+01  1.4110E+01 -1.0028E+01
             1.8879E+02

0ITERATION NO.:   15    OBJECTIVE VALUE:  -1195.36331533229        NO. OF FUNC. EVALS.:  75
 CUMULATIVE NO. OF FUNC. EVALS.:      229
 NPARAMETR:  9.4124E-01  2.8584E-01  4.4817E-02  1.0271E+00  1.3803E-01  8.9656E-01  6.0458E-01  4.7554E-01  3.3437E+00  8.2761E-01
             4.5715E+00
 PARAMETER:  3.9440E-02 -1.1523E+00 -3.0052E+00  1.2671E-01 -1.8803E+00 -9.1922E-03 -4.0322E-01 -6.4331E-01  1.3071E+00 -8.9216E-02
             1.6198E+00
 GRADIENT:  -1.4383E+02 -1.6205E+01 -2.2773E+01  2.5998E+01  7.4633E+01 -2.2335E+00  9.6688E+00  1.4277E+00 -7.2248E-01  3.3294E+01
             1.2043E+02

0ITERATION NO.:   20    OBJECTIVE VALUE:  -1234.08657607153        NO. OF FUNC. EVALS.:  72
 CUMULATIVE NO. OF FUNC. EVALS.:      301
 NPARAMETR:  1.0239E+00  3.6134E-01  6.4621E-02  8.1806E-01  1.6005E-01  9.2152E-01  2.9309E-01  1.6998E+00  2.1053E+00  5.4771E-01
             3.9200E+00
 PARAMETER:  1.2357E-01 -9.1795E-01 -2.6392E+00 -1.0082E-01 -1.7323E+00  1.8275E-02 -1.1273E+00  6.3053E-01  8.4444E-01 -5.0202E-01
             1.4661E+00
 GRADIENT:   7.9116E+00  5.7312E+01 -7.1313E+00  3.9301E+00 -7.3349E+01  1.2026E+01  1.9681E+00  1.3898E+01  1.5639E+01 -1.6349E-01
             4.7120E+01

0ITERATION NO.:   25    OBJECTIVE VALUE:  -1249.19675118388        NO. OF FUNC. EVALS.:  71
 CUMULATIVE NO. OF FUNC. EVALS.:      372
 NPARAMETR:  1.0363E+00  3.7038E-01  8.7197E-02  8.8625E-01  1.8181E-01  9.1602E-01  1.0962E-01  1.5754E+00  1.6083E+00  5.5121E-01
             3.3513E+00
 PARAMETER:  1.3566E-01 -8.9323E-01 -2.3396E+00 -2.0761E-02 -1.6048E+00  1.2282E-02 -2.1107E+00  5.5451E-01  5.7518E-01 -4.9564E-01
             1.3093E+00
 GRADIENT:   3.4312E+01  2.1660E+00 -1.2428E+01  1.0542E+01 -1.3412E+01  9.4289E+00  2.7065E-01 -4.7196E+00 -9.5314E+00  3.2759E+00
            -9.8197E+00

0ITERATION NO.:   30    OBJECTIVE VALUE:  -1257.71741784068        NO. OF FUNC. EVALS.:  70
 CUMULATIVE NO. OF FUNC. EVALS.:      442
 NPARAMETR:  1.0174E+00  4.8861E-01  1.2408E-01  9.6032E-01  2.3717E-01  8.6781E-01  4.9332E-02  1.8514E+00  1.6045E+00  3.0978E-01
             3.2978E+00
 PARAMETER:  1.1726E-01 -6.1618E-01 -1.9868E+00  5.9511E-02 -1.3390E+00 -4.1786E-02 -2.9092E+00  7.1595E-01  5.7284E-01 -1.0719E+00
             1.2933E+00
 GRADIENT:  -1.3489E+00  6.6744E-01 -1.2273E+00  3.1131E+00  9.5666E+00 -2.4193E+00  3.5345E-02 -6.9109E-01  5.7498E+00  2.6439E+00
             2.3094E+00

0ITERATION NO.:   35    OBJECTIVE VALUE:  -1258.33820592065        NO. OF FUNC. EVALS.:  74
 CUMULATIVE NO. OF FUNC. EVALS.:      516
 NPARAMETR:  1.0209E+00  4.7819E-01  1.2879E-01  9.6518E-01  2.3810E-01  8.7508E-01  4.1653E-02  1.8912E+00  1.5416E+00  2.0818E-01
             3.3169E+00
 PARAMETER:  1.2067E-01 -6.3776E-01 -1.9495E+00  6.4559E-02 -1.3350E+00 -3.3443E-02 -3.0784E+00  7.3720E-01  5.3285E-01 -1.4693E+00
             1.2990E+00
 GRADIENT:   1.1754E+00 -2.3084E+00 -3.5639E-01 -8.9282E-01  2.7292E+00  6.9969E-01  1.3021E-02  6.1357E-01 -6.5281E-01  2.3555E-01
             9.6855E-01

0ITERATION NO.:   40    OBJECTIVE VALUE:  -1258.72919196494        NO. OF FUNC. EVALS.: 158
 CUMULATIVE NO. OF FUNC. EVALS.:      674
 NPARAMETR:  1.0219E+00  5.0476E-01  1.3740E-01  9.7558E-01  2.5088E-01  8.7104E-01  3.2932E-02  1.9223E+00  1.5363E+00  1.3055E-01
             3.3187E+00
 PARAMETER:  1.2165E-01 -5.8368E-01 -1.8848E+00  7.5277E-02 -1.2828E+00 -3.8063E-02 -3.3133E+00  7.5351E-01  5.2935E-01 -1.9360E+00
             1.2996E+00
 GRADIENT:   1.0416E+00 -9.2854E-01 -1.6323E-01 -4.3676E-01  1.7272E+00 -1.2395E-01  2.8565E-03 -6.0524E-01  2.1397E-01  1.5627E-01
             7.1874E-01

0ITERATION NO.:   45    OBJECTIVE VALUE:  -1258.76666427304        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:      849
 NPARAMETR:  1.0218E+00  5.0668E-01  1.3975E-01  9.7994E-01  2.5315E-01  8.7082E-01  2.5634E-02  1.9403E+00  1.5327E+00  5.4949E-02
             3.3180E+00
 PARAMETER:  1.2156E-01 -5.7988E-01 -1.8679E+00  7.9739E-02 -1.2738E+00 -3.8324E-02 -3.5638E+00  7.6284E-01  5.2705E-01 -2.8013E+00
             1.2994E+00
 GRADIENT:   9.6583E-02 -5.9667E-01 -6.7544E-02 -2.5687E-03  1.0322E+00  9.0698E-03 -2.7600E-04 -8.6859E-02  3.4304E-02  6.1181E-03
             6.7171E-02

0ITERATION NO.:   50    OBJECTIVE VALUE:  -1258.76749157857        NO. OF FUNC. EVALS.: 175
 CUMULATIVE NO. OF FUNC. EVALS.:     1024
 NPARAMETR:  1.0218E+00  5.0694E-01  1.3991E-01  9.8016E-01  2.5332E-01  8.7074E-01  2.2812E-02  1.9423E+00  1.5327E+00  3.4597E-02
             3.3192E+00
 PARAMETER:  1.2156E-01 -5.7936E-01 -1.8667E+00  7.9963E-02 -1.2731E+00 -3.8417E-02 -3.6805E+00  7.6387E-01  5.2700E-01 -3.2640E+00
             1.2997E+00
 GRADIENT:  -2.6755E-03 -2.3525E-03  2.4832E-03  2.4831E-03  1.9091E-03  6.8925E-04 -5.5929E-04  1.6465E-03  2.2208E-03  1.9800E-04
             4.6832E-03

0ITERATION NO.:   55    OBJECTIVE VALUE:  -1258.76932034376        NO. OF FUNC. EVALS.: 142
 CUMULATIVE NO. OF FUNC. EVALS.:     1166
 NPARAMETR:  1.0218E+00  5.0697E-01  1.3990E-01  9.8015E-01  2.5332E-01  8.7070E-01  7.6177E-02  1.9422E+00  1.5325E+00  3.1130E-02
             3.3188E+00
 PARAMETER:  1.2154E-01 -5.7930E-01 -1.8668E+00  7.9948E-02 -1.2731E+00 -3.8455E-02 -2.4747E+00  7.6384E-01  5.2688E-01 -3.3696E+00
             1.2996E+00
 GRADIENT:   3.3419E+00  2.0603E-01  1.8355E+00  4.0416E-01  7.6457E+00  1.4222E-01  1.1071E-03  4.5320E-01  6.9196E-01  3.6086E-03
             1.2686E+00

0ITERATION NO.:   60    OBJECTIVE VALUE:  -1258.77027728003        NO. OF FUNC. EVALS.: 138
 CUMULATIVE NO. OF FUNC. EVALS.:     1304
 NPARAMETR:  1.0219E+00  5.0734E-01  1.3986E-01  9.8002E-01  2.5326E-01  8.7083E-01  1.0331E-01  1.9419E+00  1.5316E+00  2.4261E-02
             3.3177E+00
 PARAMETER:  1.2163E-01 -5.7858E-01 -1.8671E+00  7.9816E-02 -1.2733E+00 -3.8309E-02 -2.1701E+00  7.6367E-01  5.2631E-01 -3.6189E+00
             1.2993E+00
 GRADIENT:   2.4418E-01 -5.2408E-02  2.2526E-01  1.0595E-01 -2.9184E-01  5.9946E-04  1.3650E-03  2.3315E-02  1.5576E-02  1.3813E-03
            -2.2881E-02

0ITERATION NO.:   65    OBJECTIVE VALUE:  -1258.77105251328        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1481
 NPARAMETR:  1.0218E+00  5.0795E-01  1.3979E-01  9.7971E-01  2.5343E-01  8.7082E-01  1.0684E-01  1.9418E+00  1.5314E+00  1.0000E-02
             3.3189E+00
 PARAMETER:  1.2154E-01 -5.7737E-01 -1.8676E+00  7.9500E-02 -1.2727E+00 -3.8315E-02 -2.1365E+00  7.6361E-01  5.2616E-01 -4.7318E+00
             1.2996E+00
 GRADIENT:   1.2394E-01 -4.1610E-02  2.8621E-02  1.8192E-02 -2.2952E-02  1.6113E-02  1.3937E-04  5.7536E-03 -3.5980E-04  0.0000E+00
             4.6905E-02

0ITERATION NO.:   70    OBJECTIVE VALUE:  -1258.77107100250        NO. OF FUNC. EVALS.: 177
 CUMULATIVE NO. OF FUNC. EVALS.:     1658
 NPARAMETR:  1.0217E+00  5.0813E-01  1.3977E-01  9.7961E-01  2.5347E-01  8.7077E-01  1.0776E-01  1.9418E+00  1.5315E+00  1.0000E-02
             3.3187E+00
 PARAMETER:  1.2147E-01 -5.7701E-01 -1.8678E+00  7.9396E-02 -1.2725E+00 -3.8378E-02 -2.1278E+00  7.6362E-01  5.2624E-01 -4.9435E+00
             1.2996E+00
 GRADIENT:  -1.1584E-02 -7.3492E-03 -2.1543E-03 -3.0497E-03  9.7505E-03 -2.6585E-03  4.0920E-05 -1.6532E-03  2.7007E-03  0.0000E+00
            -2.6683E-03

0ITERATION NO.:   71    OBJECTIVE VALUE:  -1258.77107100250        NO. OF FUNC. EVALS.:  28
 CUMULATIVE NO. OF FUNC. EVALS.:     1686
 NPARAMETR:  1.0218E+00  5.0817E-01  1.3977E-01  9.7966E-01  2.5346E-01  8.7084E-01  1.0753E-01  1.9422E+00  1.5313E+00  1.0000E-02
             3.3188E+00
 PARAMETER:  1.2147E-01 -5.7701E-01 -1.8678E+00  7.9396E-02 -1.2725E+00 -3.8378E-02 -2.1278E+00  7.6362E-01  5.2624E-01 -4.9435E+00
             1.2996E+00
 GRADIENT:  -1.2673E-02 -7.4313E-03 -1.4147E-04 -2.6276E-03  9.0617E-03 -1.6856E-03  9.5671E-06 -2.2687E-03  2.0537E-03  0.0000E+00
            -1.9388E-03

 #TERM:
0MINIMIZATION TERMINATED
 DUE TO ROUNDING ERRORS (ERROR=134)
 NO. OF FUNCTION EVALUATIONS USED:     1686
 NO. OF SIG. DIGITS UNREPORTABLE
0PARAMETER ESTIMATE IS NEAR ITS BOUNDARY

 ETABAR IS THE ARITHMETIC MEAN OF THE ETA-ESTIMATES,
 AND THE P-VALUE IS GIVEN FOR THE NULL HYPOTHESIS THAT THE TRUE MEAN IS 0.

 ETABAR:         6.1751E-03 -6.5757E-03  1.7216E-02 -3.5197E-03  8.0028E-04
 SE:             2.8321E-02  2.2966E-03  2.3505E-02  2.5858E-02  3.7642E-04
 N:                     100         100         100         100         100

 P VAL.:         8.2740E-01  4.1936E-03  4.6391E-01  8.9173E-01  3.3502E-02

 ETASHRINKSD(%)  5.1225E+00  9.2306E+01  2.1255E+01  1.3371E+01  9.8739E+01
 ETASHRINKVR(%)  9.9826E+00  9.9408E+01  3.7992E+01  2.4954E+01  9.9984E+01
 EBVSHRINKSD(%)  5.2774E+00  9.2399E+01  2.1389E+01  1.1991E+01  9.8797E+01
 EBVSHRINKVR(%)  1.0276E+01  9.9422E+01  3.8204E+01  2.2544E+01  9.9986E+01
 RELATIVEINF(%)  8.6245E+01  4.3190E-02  2.5072E+01  6.2421E+01  8.7905E-04
 EPSSHRINKSD(%)  3.1352E+01
 EPSSHRINKVR(%)  5.2875E+01

  
 TOTAL DATA POINTS NORMALLY DISTRIBUTED (N):          400
 N*LOG(2PI) CONSTANT TO OBJECTIVE FUNCTION:    735.15082656373818     
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -1258.7710710024985     
 OBJECTIVE FUNCTION VALUE WITH CONSTANT:      -523.62024443876032     
 REPORTED OBJECTIVE FUNCTION DOES NOT CONTAIN CONSTANT
  
 TOTAL EFFECTIVE ETAS (NIND*NETA):                           500
  
 #TERE:
 Elapsed estimation  time in seconds:    21.02
0R MATRIX ALGORITHMICALLY SINGULAR
 AND ALGORITHMICALLY NON-POSITIVE-SEMIDEFINITE
0R MATRIX IS OUTPUT
0COVARIANCE STEP ABORTED
 Elapsed covariance  time in seconds:     6.90
 Elapsed postprocess time in seconds:     0.00
1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 #OBJT:**************                       MINIMUM VALUE OF OBJECTIVE FUNCTION                      ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 





 #OBJV:********************************************    -1258.771       **************************************************
1
 ************************************************************************************************************************
 ********************                                                                                ********************
 ********************               FIRST ORDER CONDITIONAL ESTIMATION WITH INTERACTION              ********************
 ********************                             FINAL PARAMETER ESTIMATE                           ********************
 ********************                                                                                ********************
 ************************************************************************************************************************
 


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   *********


         TH 1      TH 2      TH 3      TH 4      TH 5      TH 6      TH 7      TH 8      TH 9      TH10      TH11     
 
         1.02E+00  5.08E-01  1.40E-01  9.80E-01  2.53E-01  8.71E-01  1.08E-01  1.94E+00  1.53E+00  1.00E-02  3.32E+00
 


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
+        1.27E+03
 
 TH 2
+       -1.15E+01  2.60E+03
 
 TH 3
+       -3.81E+02  2.69E+03  9.09E+03
 
 TH 4
+       -4.51E+01  2.35E+02 -2.54E+02  3.77E+02
 
 TH 5
+        3.33E+02 -7.98E+03 -1.08E+04 -5.31E+02  2.83E+04
 
 TH 6
+       -3.66E+00 -2.27E+01 -2.03E+01 -7.87E+00  1.01E+02  2.10E+02
 
 TH 7
+       -6.55E-01 -3.85E+01 -9.88E+00 -6.93E-01  1.08E+02  1.08E-02  2.11E+00
 
 TH 8
+        3.36E+00 -3.39E+00 -3.82E+01  7.42E-01 -2.62E+01  3.56E+00  4.88E-01  2.30E+01
 
 TH 9
+        1.54E+01 -6.63E+01  2.93E+01 -1.31E+01  2.09E+02  1.20E+00  2.40E+00  1.18E-01  4.21E+01
 
 TH10
+        0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00  0.00E+00
 
 TH11
+       -2.19E+01 -6.63E+01  2.12E+01 -5.94E-01  1.78E+02  4.09E+00  2.21E+00  4.13E+00  8.84E+00  0.00E+00  2.65E+01
 
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
 #CPUT: Total CPU Time in Seconds,       27.984
Stop Time:
Sat Sep 25 09:24:18 CDT 2021
